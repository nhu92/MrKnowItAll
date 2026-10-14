"""Run multiple real mixtures to genus on Zaratan from a login node.

The login node only fetches Kew files and manages Slurm jobs. FASTQ processing,
reference construction, and phylogenetic analysis run on compute nodes.
Restarting this command resumes submitted and completed stages.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import shutil
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from contextlib import contextmanager
from datetime import UTC, datetime
from pathlib import Path

TERMINAL_FAILURES = {
    "BOOT_FAIL",
    "CANCELLED",
    "DEADLINE",
    "FAILED",
    "NODE_FAIL",
    "OUT_OF_MEMORY",
    "PREEMPTED",
    "TIMEOUT",
}
PREFETCH_LOCK = threading.Lock()
PRINT_LOCK = threading.Lock()
NUMERIC_THREAD_VARIABLES = (
    "OPENBLAS_NUM_THREADS",
    "OMP_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
)


def say(message: str) -> None:
    with PRINT_LOCK:
        print(f"{datetime.now(UTC).isoformat()} {message}", flush=True)


def run_command(args: list[str], *, cwd: Path | None = None, env: dict | None = None) -> str:
    command_env = (env if env is not None else os.environ).copy()
    # Login-node helper commands only need one BLAS thread. Unbounded NumPy
    # imports can otherwise exhaust Zaratan's per-user process/thread limit.
    for name in NUMERIC_THREAD_VARIABLES:
        command_env[name] = "1"
    result = subprocess.run(
        args, cwd=cwd, env=command_env, capture_output=True, text=True, check=False
    )
    if result.returncode:
        raise RuntimeError(
            f"Command failed ({result.returncode}): {' '.join(args)}\n"
            f"stdout: {result.stdout[-2000:]}\nstderr: {result.stderr[-2000:]}"
        )
    return result.stdout.strip()


def top_aggregate_taxon(predictions: Path) -> str:
    with predictions.open(encoding="utf-8-sig", newline="") as handle:
        return next(csv.DictReader(handle), {}).get("row_name", "")


def sample_pairs(sample_root: Path, selected: set[str], excluded: set[str]) -> dict[str, tuple[Path, Path]]:
    pairs = {}
    for directory in sorted(sample_root.glob("EMseq*")):
        if not directory.is_dir() or directory.name in excluded:
            continue
        if selected and directory.name not in selected:
            continue
        reads1 = sorted(directory.glob("*_R1_*.fastq.gz"))
        reads2 = sorted(directory.glob("*_R2_*.fastq.gz"))
        if len(reads1) != 1 or len(reads2) != 1:
            raise ValueError(f"Expected one FASTQ pair in {directory}, found {len(reads1)}/{len(reads2)}")
        if reads1[0].name.replace("_R1_", "_R2_") != reads2[0].name:
            raise ValueError(f"FASTQ pair names do not match in {directory}")
        pairs[directory.name] = (reads1[0].resolve(), reads2[0].resolve())
    if selected - pairs.keys():
        raise ValueError(f"Samples not found: {', '.join(sorted(selected - pairs.keys()))}")
    if not pairs:
        raise ValueError(f"No FASTQ pairs found under {sample_root}")
    return pairs


def slurm_state(job_id: str) -> str | None:
    result = subprocess.run(
        ["sacct", "-n", "-P", "-j", job_id, "--format=JobIDRaw,State"],
        capture_output=True, text=True, check=False,
    )
    if result.returncode:
        raise RuntimeError(f"sacct failed for {job_id}: {result.stderr.strip()}")
    for line in result.stdout.splitlines():
        fields = line.split("|")
        if len(fields) >= 2 and fields[0] == job_id:
            return fields[1].split()[0].rstrip("+")
    return None


def wait_for_stage(stage: str, job_id: str, result_dir: Path, poll_seconds: int) -> None:
    last_state = None
    while True:
        if (result_dir / "SUCCESS").exists():
            say(f"{stage}: complete ({job_id})")
            return
        state = slurm_state(job_id)
        if state != last_state:
            say(f"{stage}: job {job_id} {state or 'not yet in accounting'}")
            last_state = state
        if state == "COMPLETED":
            raise RuntimeError(f"{stage}: Slurm completed without {result_dir / 'SUCCESS'}")
        if state in TERMINAL_FAILURES:
            raise RuntimeError(f"{stage}: Slurm job {job_id} ended {state}; inspect its logs")
        time.sleep(poll_seconds)


@contextmanager
def submission_lock(stage_dir: Path):
    """Prevent simultaneous controllers from submitting the same Slurm stage."""
    with (stage_dir / ".submit.lock").open("a+", encoding="utf-8") as handle:
        if os.name == "posix":
            import fcntl

            fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
        try:
            yield
        finally:
            if os.name == "posix":
                fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


@contextmanager
def controller_lock(run_root: Path):
    """Allow only one full pipeline controller per shared run directory."""
    with (run_root / ".controller.lock").open("a+", encoding="utf-8") as handle:
        if os.name == "posix":
            import fcntl

            try:
                fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError as exc:
                raise RuntimeError(
                    f"Another controller is already running for {run_root}"
                ) from exc
        try:
            yield
        finally:
            if os.name == "posix":
                fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


def queue_stage(
    stage: str,
    source_script: Path,
    stage_dir: Path,
    environment: dict[str, str],
    memory: str,
) -> tuple[Path, str | None]:
    result_dir = stage_dir / "results"
    stage_dir.mkdir(parents=True, exist_ok=True)
    with submission_lock(stage_dir):
        if (result_dir / "SUCCESS").exists():
            say(f"{stage}: reusing completed result")
            return result_dir, None
        (stage_dir / "logs").mkdir(exist_ok=True)
        script = stage_dir / source_script.name
        shutil.copy2(source_script, script)
        job_file = stage_dir / "JOB_ID"
        if job_file.exists():
            job_id = job_file.read_text(encoding="utf-8").strip()
        else:
            if result_dir.exists():
                raise RuntimeError(f"{stage}: partial results exist without a recorded job: {result_dir}")
            env = {
                key: value for key, value in os.environ.items()
                if not key.startswith(("SPROUT_", "MRKNOWITALL_"))
            }
            env.update(environment)
            job_id = run_command(
                ["sbatch", "--parsable", "--export=ALL", f"--mem={memory}", script.name],
                cwd=stage_dir,
                env=env,
            ).split(";")[0]
            if not re.fullmatch(r"\d+", job_id):
                raise RuntimeError(f"Unexpected sbatch response: {job_id}")
            job_file.write_text(job_id + "\n", encoding="utf-8")
            say(f"{stage}: submitted job {job_id}")
    return result_dir, job_id


def submit_stage(
    stage: str,
    source_script: Path,
    stage_dir: Path,
    environment: dict[str, str],
    memory: str,
    poll_seconds: int,
) -> Path:
    result_dir, job_id = queue_stage(
        stage, source_script, stage_dir, environment, memory
    )
    if job_id is None:
        return result_dir
    wait_for_stage(stage, job_id, result_dir, poll_seconds)
    return result_dir


def prefetch(
    stage: str,
    level: str,
    parent_file: Path,
    report: Path,
    cache: Path,
    backbone: Path,
    allowed_taxa_file: Path | None = None,
    representatives: int | None = None,
    pool_multiplier: int = 2,
) -> None:
    if report.exists() and report.stat().st_size > 0:
        say(f"{stage}: reusing Kew prefetch report")
        return
    report.parent.mkdir(parents=True, exist_ok=True)
    representatives = representatives or (4 if level == "family" else 2)
    command = [
        sys.executable, "-m", "sprout_refbuilder.cli",
        "--cache", str(cache), "prefetch-kew-panel",
        "--level", level, "--parent-file", str(parent_file),
        "--backbone", str(backbone), "--report", str(report),
        "--representatives", str(representatives),
        "--pool-multiplier", str(pool_multiplier),
        "--minimum-recovered-loci", "50", "--download-workers", "8",
    ]
    if allowed_taxa_file is not None:
        command.extend(["--allowed-taxa-file", str(allowed_taxa_file)])
    # A shared cache avoids duplicate downloads; one writer at a time prevents
    # two samples from writing the same recovery file concurrently.
    with PREFETCH_LOCK:
        if report.exists() and report.stat().st_size > 0:
            return
        say(f"{stage}: downloading Kew recoveries")
        try:
            output = run_command(command)
        except Exception:
            report.unlink(missing_ok=True)
            raise
        summary = json.loads(output)
        say(
            f"{stage}: cached {summary.get('usable_recoveries', '?')} usable "
            f"recoveries from release {summary.get('kew_release', '?')}"
        )


def common_environment(args: argparse.Namespace, project: str) -> dict[str, str]:
    return {
        "SPROUT_PROJECT": project,
        "SPROUT_SOURCE": str(args.sprout),
        "MRKNOWITALL_SOURCE": str(args.repo),
        "MRKNOWITALL_BACKBONE": str(args.backbone),
        "SPROUT_REF_CACHE": str(args.cache),
        "SPROUT_CONDA_ENV": args.conda_env,
    }


def process_sample(
    sample: str,
    reads: tuple[Path, Path],
    args: argparse.Namespace,
) -> str:
    sample_dir = args.run_root / sample
    sample_dir.mkdir(parents=True, exist_ok=True)
    if (sample_dir / "branch" / "plan.json").exists():
        raise RuntimeError(
            f"{sample}: this run root uses the earlier single-pass branch format; "
            "choose a fresh run root for progressive screening"
        )
    project = sample.lower()
    common = common_environment(args, project)
    scripts = args.repo / "examples" / "hpc"
    order_dir = sample_dir / "order"
    old_sample = args.reuse_upstream_from / sample if args.reuse_upstream_from else None
    old_order = old_sample / "order" / "results" if old_sample else None
    if old_order is not None and (old_order / "SUCCESS").is_file():
        order = old_order
        say(f"{sample}/order: reusing completed upstream result at {order}")
    else:
        order = submit_stage(
            f"{sample}/order", scripts / "zaratan_query_order.sbatch", order_dir,
            common | {
                "SPROUT_READ1": str(reads[0]),
                "SPROUT_READ2": str(reads[1]),
                "SPROUT_STAGE_OUTPUT": str(order_dir / "results"),
            },
            args.order_memory, args.poll_seconds,
        )
    family_dir = sample_dir / "family"
    old_family_dir = old_sample / "family" if old_sample else None
    old_family = old_family_dir / "results" if old_family_dir else None
    if old_family is not None and (old_family / "SUCCESS").is_file():
        family = old_family
        say(f"{sample}/family: reusing completed upstream result at {family}")
    elif old_family_dir is not None and (old_family_dir / "JOB_ID").is_file():
        old_job_id = (old_family_dir / "JOB_ID").read_text(encoding="utf-8").strip()
        if not re.fullmatch(r"\d+", old_job_id):
            raise RuntimeError(f"{sample}/family: invalid upstream job ID")
        say(f"{sample}/family: waiting for existing upstream job {old_job_id}")
        wait_for_stage(f"{sample}/family", old_job_id, old_family, args.poll_seconds)
        family = old_family
    else:
        prefetch(
            f"{sample}/family", "family", order / "candidates.txt",
            family_dir / "prefetch.json", args.cache, args.backbone,
        )
        family = submit_stage(
            f"{sample}/family", scripts / "zaratan_refine_stage.sbatch", family_dir,
            common | {
                "SPROUT_LEVEL": "family",
                "SPROUT_PARENT_FILE": str(order / "candidates.txt"),
                "SPROUT_QUERY_EXON": str(order / "exon"),
                "SPROUT_QUERY_GENE_LIST": str(order / "gene_list.txt"),
                "SPROUT_STAGE_OUTPUT": str(family_dir / "results"),
            },
            args.refine_memory, args.poll_seconds,
        )
    (sample_dir / "UPSTREAM_SOURCES.json").write_text(
        json.dumps({"order_results": str(order), "family_results": str(family)}, indent=2),
        encoding="utf-8",
    )
    branch_root = sample_dir / "branch"
    branch_root.mkdir(exist_ok=True)
    genus_allowed: Path | None = None
    branch_plan: Path | None = None
    plan: dict = {}
    for round_number in range(1, args.max_branch_rounds + 1):
        branch_dir = branch_root / f"round_{round_number:02d}"
        branch_dir.mkdir(exist_ok=True)
        branch_plan = branch_dir / "plan.json"
        if not branch_plan.exists():
            command = [
                sys.executable, "-m", "sprout_refbuilder.cli",
                "--cache", str(args.cache), "plan-genus-branches",
                "--parent-file", str(family / "candidates.txt"),
                "--backbone", str(args.backbone), "--output", str(branch_plan),
                "--threshold", str(min(args.branch_threshold, args.max_final_genera)),
                "--max-group-genera", str(args.genus_shard_size),
            ]
            if genus_allowed is not None:
                command.extend(["--allowed-genera-file", str(genus_allowed)])
            # Planning reads/initializes the shared Kew manifest and species tree.
            # Serialize it with prefetch so concurrent samples never see partial files.
            with PREFETCH_LOCK:
                if not branch_plan.exists():
                    run_command(command)
        plan = json.loads(branch_plan.read_text(encoding="utf-8"))
        source_count = len(plan["all_genera"])
        say(
            f"{sample}/branch round {round_number}: {source_count} genera; "
            f"screening {'enabled' if plan['triggered'] else 'not needed'}"
        )
        if source_count <= args.max_final_genera or not plan["triggered"]:
            break
        anchor_file = branch_dir / "anchors.txt"
        anchor_file.write_text(
            "".join(f"{name}\n" for name in plan["anchor_genera"]),
            encoding="utf-8",
        )
        prefetch(
            f"{sample}/branch round {round_number}", "genus",
            family / "candidates.txt", branch_dir / "prefetch.json",
            args.cache, args.backbone, allowed_taxa_file=anchor_file,
            representatives=1, pool_multiplier=1,
        )
        branch_results = submit_stage(
            f"{sample}/branch round {round_number}",
            scripts / "zaratan_refine_stage.sbatch", branch_dir,
            common | {
                "SPROUT_LEVEL": "genus",
                "SPROUT_PARENT_FILE": str(family / "candidates.txt"),
                "SPROUT_QUERY_EXON": str(order / "exon"),
                "SPROUT_QUERY_GENE_LIST": str(order / "gene_list.txt"),
                "SPROUT_STAGE_OUTPUT": str(branch_dir / "results"),
                "SPROUT_ALLOWED_TAXA_FILE": str(anchor_file),
                "SPROUT_REF_REPRESENTATIVES": "1",
                "SPROUT_REF_POOL_MULTIPLIER": "1",
                "SPROUT_STAGE_ID": f"branch_round_{round_number:02d}",
            },
            args.refine_memory, args.poll_seconds,
        )
        next_allowed = branch_dir / "allowed_genera.txt"
        if not next_allowed.exists():
            run_command([
                sys.executable, "-m", "sprout_refbuilder.cli", "expand-genus-branches",
                "--plan", str(branch_plan),
                "--evidence", str(branch_results / "candidate_evidence.csv"),
                "--output", str(next_allowed),
                "--report", str(branch_dir / "expansion.json"),
            ])
        expanded = json.loads((branch_dir / "expansion.json").read_text(encoding="utf-8"))
        remaining = expanded["selected_genera"]
        say(f"{sample}/branch round {round_number}: {remaining}/{source_count} genera retained")
        genus_allowed = next_allowed
        if remaining >= source_count or remaining <= args.max_final_genera:
            break

    remaining_count = (
        len(genus_allowed.read_text(encoding="utf-8").splitlines())
        if genus_allowed is not None else len(plan["all_genera"])
    )
    if remaining_count > args.max_postscreen_genera:
        reason = (
            f"Tree screening stopped at {remaining_count} genera, above the "
            f"{args.max_postscreen_genera}-genus post-screen limit; no full-panel fallback"
        )
        unresolved = {
            "status": "genus_unresolved_overflow",
            "remaining_genera": remaining_count,
            "limit": args.max_postscreen_genera,
            "reason": reason,
        }
        (branch_root / "UNRESOLVED.json").write_text(
            json.dumps(unresolved, indent=2), encoding="utf-8"
        )
        summary = sample_dir / "RESULT.txt"
        summary.write_text(
            f"sample={sample}\n"
            f"order={','.join((order / 'candidates.txt').read_text().splitlines())}\n"
            f"family={','.join((family / 'candidates.txt').read_text().splitlines())}\n"
            f"genus_status=unresolved_overflow\nremaining_genera={remaining_count}\n"
            f"reason={reason}\n",
            encoding="utf-8",
        )
        run_command([
            sys.executable, str(args.repo / "examples/hpc/export_sample_report.py"),
            "--sample-dir", str(sample_dir), "--sprout", str(args.sprout),
        ])
        say(f"{sample}: {reason}")
        return summary.read_text(encoding="utf-8").strip()
    if genus_allowed is None and remaining_count > args.max_final_genera:
        genus_allowed = branch_root / "all_small_genera.txt"
        genus_allowed.write_text(
            "".join(f"{name}\n" for name in plan["all_genera"]), encoding="utf-8"
        )
    if genus_allowed is not None:
        allowed_count = len(genus_allowed.read_text(encoding="utf-8").splitlines())
        if allowed_count > args.max_final_genera:
            shard_root = sample_dir / "genus_shards"
            shard_plan = shard_root / "shards.json"
            if not shard_plan.exists():
                if branch_plan is None:
                    raise RuntimeError("Genus shard planning requires a tree plan")
                run_command([
                    sys.executable, "-m", "sprout_refbuilder.cli", "plan-genus-shards",
                    "--plan", str(branch_plan),
                    "--allowed-taxa-file", str(genus_allowed),
                    "--output-dir", str(shard_root),
                    "--max-genera-per-shard", str(args.genus_shard_size),
                ])
            shard_manifest = json.loads(shard_plan.read_text(encoding="utf-8"))
            if set(shard_manifest["allowed_genera"]) != set(genus_allowed.read_text().split()):
                raise RuntimeError(f"{sample}: existing shard plan does not match the genus set")
            shards = shard_manifest["shards"]
            say(
                f"{sample}/genus_shards: {allowed_count} genera in {len(shards)} "
                "independent Slurm panels"
            )

            def run_shard(shard: dict) -> None:
                shard_id = shard["id"]
                shard_dir = shard_root / shard_id
                allowed_file = shard_root / f"{shard_id}.txt"
                prefetch(
                    f"{sample}/{shard_id}", "genus", family / "candidates.txt",
                    shard_dir / "prefetch.json", args.cache, args.backbone,
                    allowed_taxa_file=allowed_file, representatives=1,
                    pool_multiplier=1,
                )
                submit_stage(
                    f"{sample}/{shard_id}", scripts / "zaratan_refine_stage.sbatch",
                    shard_dir,
                    common | {
                        "SPROUT_LEVEL": "genus",
                        "SPROUT_PARENT_FILE": str(family / "candidates.txt"),
                        "SPROUT_QUERY_EXON": str(order / "exon"),
                        "SPROUT_QUERY_GENE_LIST": str(order / "gene_list.txt"),
                        "SPROUT_STAGE_OUTPUT": str(shard_dir / "results"),
                        "SPROUT_ALLOWED_TAXA_FILE": str(allowed_file),
                        "SPROUT_REF_REPRESENTATIVES": "1",
                        "SPROUT_REF_POOL_MULTIPLIER": "1",
                    },
                    args.shard_memory, args.poll_seconds,
                )

            with ThreadPoolExecutor(max_workers=args.max_parallel_shards) as executor:
                futures = [executor.submit(run_shard, shard) for shard in shards]
                for future in as_completed(futures):
                    future.result()
            shortlist = shard_root / "shortlist.txt"
            if not shortlist.exists():
                run_command([
                    sys.executable, "-m", "sprout_refbuilder.cli",
                    "shortlist-genus-shards", "--plan", str(shard_plan),
                    "--shard-runs", str(shard_root), "--output", str(shortlist),
                    "--report", str(shard_root / "shortlist.json"),
                    "--max-final-genera", str(args.max_final_genera),
                ])
            shortlisted = len(shortlist.read_text(encoding="utf-8").splitlines())
            say(f"{sample}/genus_shards: {shortlisted} genera retained for joint validation")
            genus_allowed = shortlist
    genus_dir = sample_dir / "genus"
    prefetch(
        f"{sample}/genus", "genus", family / "candidates.txt",
        genus_dir / "prefetch.json", args.cache, args.backbone,
        allowed_taxa_file=genus_allowed,
    )
    genus = submit_stage(
        f"{sample}/genus", scripts / "zaratan_refine_stage.sbatch", genus_dir,
        common | {
            "SPROUT_LEVEL": "genus",
            "SPROUT_PARENT_FILE": str(family / "candidates.txt"),
            "SPROUT_QUERY_EXON": str(order / "exon"),
            "SPROUT_QUERY_GENE_LIST": str(order / "gene_list.txt"),
            "SPROUT_STAGE_OUTPUT": str(genus_dir / "results"),
            **({"SPROUT_ALLOWED_TAXA_FILE": str(genus_allowed)} if genus_allowed else {}),
        },
        args.refine_memory, args.poll_seconds,
    )
    summary = sample_dir / "RESULT.txt"
    evidence_candidates = (genus / "candidates.txt").read_text(encoding="utf-8").splitlines()
    summary.write_text(
        f"sample={sample}\norder={','.join((order / 'candidates.txt').read_text().splitlines())}"
        f"\nfamily={','.join((family / 'candidates.txt').read_text().splitlines())}"
        f"\ngenus={','.join(evidence_candidates)}"
        f"\ngenus_top_by_evidence={evidence_candidates[0] if evidence_candidates else ''}"
        f"\ngenus_top_by_sprout_score={top_aggregate_taxon(genus / 'predictions.csv')}\n",
        encoding="utf-8",
    )
    run_command([
        sys.executable, str(args.repo / "examples/hpc/export_sample_report.py"),
        "--sample-dir", str(sample_dir), "--sprout", str(args.sprout),
    ])
    say(f"{sample}: finished order → family → genus")
    return summary.read_text(encoding="utf-8").strip()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-root", required=True, type=Path)
    parser.add_argument("--run-root", required=True, type=Path)
    parser.add_argument("--repo", required=True, type=Path)
    parser.add_argument("--sprout", required=True, type=Path)
    parser.add_argument("--backbone", required=True, type=Path)
    parser.add_argument("--cache", required=True, type=Path)
    parser.add_argument(
        "--reuse-upstream-from", type=Path,
        help="Reuse completed order/family results (and wait for a recorded family job) from an earlier run",
    )
    parser.add_argument("--sample", action="append", default=[])
    parser.add_argument("--exclude", action="append", default=[])
    parser.add_argument(
        "--validate-only", action="store_true",
        help="Check selected FASTQ pairs and upstream reuse without creating directories or submitting jobs",
    )
    parser.add_argument(
        "--submit-orders-only", action="store_true",
        help="Queue every selected sample's independent order job, then exit",
    )
    parser.add_argument("--max-parallel-samples", type=int, default=3)
    parser.add_argument("--poll-seconds", type=int, default=30)
    parser.add_argument("--order-memory", default="96G")
    parser.add_argument("--refine-memory", default="192G")
    parser.add_argument("--shard-memory", default="96G")
    parser.add_argument("--genus-shard-size", type=int, default=40)
    parser.add_argument("--max-parallel-shards", type=int, default=2)
    parser.add_argument("--max-final-genera", type=int, default=80)
    parser.add_argument(
        "--max-branch-rounds", type=int, default=5,
        help="Maximum coarse-to-fine tree screens before reporting unresolved overflow",
    )
    parser.add_argument(
        "--max-postscreen-genera", type=int, default=160,
        help="Maximum remaining genera eligible for small final shards",
    )
    parser.add_argument(
        "--branch-threshold", type=int, default=120,
        help="Screen Kew tree branches when the genus count exceeds this limit",
    )
    parser.add_argument("--conda-env", default="sprout")
    args = parser.parse_args()
    if min(
        args.max_parallel_samples, args.max_parallel_shards, args.poll_seconds,
        args.branch_threshold, args.genus_shard_size, args.max_final_genera,
        args.max_branch_rounds, args.max_postscreen_genera,
    ) < 1:
        parser.error("Parallel counts, limits, and poll interval must be positive")
    if args.max_postscreen_genera < args.max_final_genera:
        parser.error("--max-postscreen-genera must be at least --max-final-genera")
    for path in (args.sample_root, args.repo, args.sprout, args.backbone):
        if not path.exists():
            parser.error(f"Required path does not exist: {path}")
    if args.reuse_upstream_from is not None and not args.reuse_upstream_from.is_dir():
        parser.error(f"Upstream run root does not exist: {args.reuse_upstream_from}")
    if args.reuse_upstream_from is not None and args.reuse_upstream_from.resolve() == args.run_root.resolve():
        parser.error("--reuse-upstream-from must differ from --run-root")
    os.environ["PYTHONPATH"] = str(args.repo / "src") + os.pathsep + os.environ.get(
        "PYTHONPATH", ""
    )
    pairs = sample_pairs(args.sample_root, set(args.sample), set(args.exclude))
    if args.validate_only:
        for sample, reads in pairs.items():
            old = args.reuse_upstream_from / sample if args.reuse_upstream_from else None
            order_done = bool(old and (old / "order" / "results" / "SUCCESS").is_file())
            family_done = bool(old and (old / "family" / "results" / "SUCCESS").is_file())
            family_job = bool(old and (old / "family" / "JOB_ID").is_file())
            say(
                f"{sample}: R1={reads[0].name} R2={reads[1].name}; "
                f"order={'reuse' if order_done else 'new'}; "
                f"family={'reuse' if family_done else 'wait' if family_job else 'new'}"
            )
        say(f"Validated {len(pairs)} sample workflows; no jobs submitted")
        return 0
    args.run_root.mkdir(parents=True, exist_ok=True)
    args.cache.mkdir(parents=True, exist_ok=True)
    if args.submit_orders_only:
        for sample, reads in pairs.items():
            order_dir = args.run_root / sample / "order"
            queue_stage(
                f"{sample}/order",
                args.repo / "examples/hpc/zaratan_query_order.sbatch",
                order_dir,
                common_environment(args, sample.lower()) | {
                    "SPROUT_READ1": str(reads[0]),
                    "SPROUT_READ2": str(reads[1]),
                    "SPROUT_STAGE_OUTPUT": str(order_dir / "results"),
                },
                args.order_memory,
            )
        say(f"Queued or reused order stages for {len(pairs)} samples")
        return 0
    with controller_lock(args.run_root):
        say(f"Starting {len(pairs)} samples with {args.max_parallel_samples} concurrent pipelines")
        failed = []
        with ThreadPoolExecutor(max_workers=args.max_parallel_samples) as executor:
            futures = {
                executor.submit(process_sample, sample, reads, args): sample
                for sample, reads in pairs.items()
            }
            for future in as_completed(futures):
                sample = futures[future]
                try:
                    say(future.result())
                except (OSError, RuntimeError, ValueError) as exc:
                    failed.append(sample)
                    say(f"{sample}: FAILED: {exc}")
        if failed:
            say(f"Failed samples: {', '.join(sorted(failed))}")
            return 1
        say("All selected samples reached genus")
        return 0


if __name__ == "__main__":
    sys.exit(main())
