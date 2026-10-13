"""Run multiple real mixtures to genus on Zaratan from a login node.

The login node only fetches Kew files and manages Slurm jobs. FASTQ processing,
reference construction, and phylogenetic analysis run on compute nodes.
Restarting this command resumes submitted and completed stages.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
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


def say(message: str) -> None:
    with PRINT_LOCK:
        print(f"{datetime.now(UTC).isoformat()} {message}", flush=True)


def run_command(args: list[str], *, cwd: Path | None = None, env: dict | None = None) -> str:
    result = subprocess.run(
        args, cwd=cwd, env=env, capture_output=True, text=True, check=False
    )
    if result.returncode:
        raise RuntimeError(
            f"Command failed ({result.returncode}): {' '.join(args)}\n"
            f"stdout: {result.stdout[-2000:]}\nstderr: {result.stderr[-2000:]}"
        )
    return result.stdout.strip()


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


def submit_stage(
    stage: str,
    source_script: Path,
    stage_dir: Path,
    environment: dict[str, str],
    memory: str,
    poll_seconds: int,
) -> Path:
    result_dir = stage_dir / "results"
    if (result_dir / "SUCCESS").exists():
        say(f"{stage}: reusing completed result")
        return result_dir
    stage_dir.mkdir(parents=True, exist_ok=True)
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


def process_sample(
    sample: str,
    reads: tuple[Path, Path],
    args: argparse.Namespace,
) -> str:
    sample_dir = args.run_root / sample
    sample_dir.mkdir(parents=True, exist_ok=True)
    project = sample.lower()
    common = {
        "SPROUT_PROJECT": project,
        "SPROUT_SOURCE": str(args.sprout),
        "MRKNOWITALL_SOURCE": str(args.repo),
        "MRKNOWITALL_BACKBONE": str(args.backbone),
        "SPROUT_REF_CACHE": str(args.cache),
        "SPROUT_CONDA_ENV": args.conda_env,
    }
    scripts = args.repo / "examples" / "hpc"
    order_dir = sample_dir / "order"
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
    branch_dir = sample_dir / "branch"
    branch_dir.mkdir(exist_ok=True)
    branch_plan = branch_dir / "plan.json"
    if not branch_plan.exists():
        say(f"{sample}/branch: counting Kew genus candidates")
        run_command([
            sys.executable, "-m", "sprout_refbuilder.cli",
            "--cache", str(args.cache), "plan-genus-branches",
            "--parent-file", str(family / "candidates.txt"),
            "--backbone", str(args.backbone), "--output", str(branch_plan),
            "--threshold", str(args.branch_threshold),
        ])
    plan = json.loads(branch_plan.read_text(encoding="utf-8"))
    say(
        f"{sample}/branch: {plan['pool_records']} planned Kew recoveries; "
        f"screening {'enabled' if plan['triggered'] else 'skipped'}"
    )
    genus_allowed: Path | None = None
    if plan["triggered"]:
        anchor_file = branch_dir / "anchors.txt"
        anchor_file.write_text(
            "".join(f"{name}\n" for name in plan["anchor_genera"]),
            encoding="utf-8",
        )
        prefetch(
            f"{sample}/branch", "genus", family / "candidates.txt",
            branch_dir / "prefetch.json", args.cache, args.backbone,
            allowed_taxa_file=anchor_file, representatives=1, pool_multiplier=1,
        )
        branch_results = submit_stage(
            f"{sample}/branch", scripts / "zaratan_refine_stage.sbatch", branch_dir,
            common | {
                "SPROUT_LEVEL": "genus",
                "SPROUT_PARENT_FILE": str(family / "candidates.txt"),
                "SPROUT_QUERY_EXON": str(order / "exon"),
                "SPROUT_QUERY_GENE_LIST": str(order / "gene_list.txt"),
                "SPROUT_STAGE_OUTPUT": str(branch_dir / "results"),
                "SPROUT_ALLOWED_TAXA_FILE": str(anchor_file),
                "SPROUT_REF_REPRESENTATIVES": "1",
                "SPROUT_REF_POOL_MULTIPLIER": "1",
            },
            args.refine_memory, args.poll_seconds,
        )
        genus_allowed = branch_dir / "allowed_genera.txt"
        if not genus_allowed.exists():
            run_command([
                sys.executable, "-m", "sprout_refbuilder.cli", "expand-genus-branches",
                "--plan", str(branch_plan),
                "--evidence", str(branch_results / "candidate_evidence.csv"),
                "--output", str(genus_allowed),
                "--report", str(branch_dir / "expansion.json"),
            ])
        expanded = json.loads((branch_dir / "expansion.json").read_text(encoding="utf-8"))
        say(
            f"{sample}/branch: {expanded['mode']} selected "
            f"{expanded['selected_genera']}/{expanded['total_genera']} genera"
        )
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
    summary.write_text(
        f"sample={sample}\norder={','.join((order / 'candidates.txt').read_text().splitlines())}"
        f"\nfamily={','.join((family / 'candidates.txt').read_text().splitlines())}"
        f"\ngenus={','.join((genus / 'candidates.txt').read_text().splitlines())}\n",
        encoding="utf-8",
    )
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
    parser.add_argument("--sample", action="append", default=[])
    parser.add_argument("--exclude", action="append", default=[])
    parser.add_argument("--max-parallel-samples", type=int, default=3)
    parser.add_argument("--poll-seconds", type=int, default=30)
    parser.add_argument("--order-memory", default="96G")
    parser.add_argument("--refine-memory", default="192G")
    parser.add_argument(
        "--branch-threshold", type=int, default=120,
        help="Screen Kew tree branches when the genus pool exceeds this many recovery FASTAs",
    )
    parser.add_argument("--conda-env", default="sprout")
    args = parser.parse_args()
    if args.max_parallel_samples < 1 or args.poll_seconds < 1 or args.branch_threshold < 1:
        parser.error("Parallel sample count and poll interval must be positive")
    for path in (args.sample_root, args.repo, args.sprout, args.backbone):
        if not path.exists():
            parser.error(f"Required path does not exist: {path}")
    os.environ["PYTHONPATH"] = str(args.repo / "src") + os.pathsep + os.environ.get(
        "PYTHONPATH", ""
    )
    pairs = sample_pairs(args.sample_root, set(args.sample), set(args.exclude))
    args.run_root.mkdir(parents=True, exist_ok=True)
    args.cache.mkdir(parents=True, exist_ok=True)
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
