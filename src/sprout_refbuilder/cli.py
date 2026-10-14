from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

from .assembly import assemble_with_hybpiper
from .backbone import Backbone, install_backbone_archive
from .branches import expand_genus_branches, plan_genus_branches
from .builder import ReferenceBuilder
from .candidates import select_mixture_candidates
from .clade import BalancedKewPanelBuilder, read_taxon_candidates
from .kew import DEFAULT_BASE_URL, KewClient
from .ml import train_anomaly_model
from .models import Taxonomy
from .panel import HierarchicalPanelBuilder, run_mix7_demo
from .shards import plan_genus_shards, shortlist_genus_shards


def _cache(value: str | None) -> Path:
    return Path(value or os.environ.get("SPROUT_REF_CACHE", Path.home() / ".cache/sprout-ref"))


def parser() -> argparse.ArgumentParser:
    root = argparse.ArgumentParser(prog="sprout-ref", description=__doc__)
    root.add_argument("--cache", help="Kew metadata/recovery cache directory")
    root.add_argument("--kew-base-url", default=DEFAULT_BASE_URL)
    commands = root.add_subparsers(dest="command", required=True)

    commands.add_parser("sync", help="Cache the current Kew release manifest and recovery index")
    search = commands.add_parser("search", help="Search taxa in the cached/current Kew release")
    search.add_argument("taxon")
    search.add_argument("--limit", type=int, default=20)

    inspect = commands.add_parser("inspect-backbone", help="Validate and summarize alignments")
    inspect.add_argument("backbone", type=Path)

    install = commands.add_parser(
        "install-backbone", help="Safely install the curated Angiosperms353 alignment ZIP"
    )
    install.add_argument("--archive", required=True, type=Path)
    install.add_argument("--destination", required=True, type=Path)
    install.add_argument("--overwrite", action="store_true")
    install.add_argument("--expected-loci", type=int, default=353)
    install.add_argument(
        "--taxonomy-corrections",
        type=Path,
        help="Override TSV; defaults to the corrections bundled with MrKnowItAll",
    )
    install.add_argument(
        "--no-taxonomy-corrections",
        action="store_true",
        help="Install original headers without bundled replacements or exclusions",
    )

    build = commands.add_parser("build", help="Build a SPrOUT reference from Kew assemblies")
    _build_arguments(build)

    local = commands.add_parser("build-local", help="Build from a gene-labelled assembled FASTA")
    _build_arguments(local, include_taxonomy=True)
    local.add_argument("--recovery-fasta", required=True, type=Path)
    local.add_argument("--sequence-id", default="local")

    reads = commands.add_parser("build-reads", help="Recover loci with HybPiper, then build")
    _build_arguments(reads, include_taxonomy=True)
    reads.add_argument("--read1", required=True, type=Path)
    reads.add_argument("--read2", required=True, type=Path)
    reads.add_argument("--target-file", required=True, type=Path)
    reads.add_argument("--threads", type=int, default=8)

    panel = commands.add_parser(
        "build-panel", help="Build an order/family/genus SPrOUT reference panel"
    )
    panel.add_argument("--level", required=True, choices=["order", "family", "genus"])
    panel.add_argument("--backbone", required=True, type=Path)
    panel.add_argument("--summary", required=True, type=Path)
    panel.add_argument("--output", required=True, type=Path)
    panel.add_argument(
        "--parent-taxon", action="append", default=[], help="Detected parent order/family"
    )
    panel.add_argument(
        "--include-species", action="append", default=[], help="Exact Kew species to integrate"
    )
    panel.add_argument("--representatives", type=int, default=2)
    panel.add_argument("--python-aligner", action="store_true")

    kew_panel = commands.add_parser(
        "build-kew-panel",
        help="Build a balanced family/genus/species panel from Kew assembled recoveries",
    )
    kew_panel.add_argument("--level", required=True, choices=["family", "genus", "species"])
    kew_panel.add_argument("--parent-taxon", action="append", default=[])
    kew_panel.add_argument(
        "--parent-file", type=Path, help="SPrOUT candidate TXT or ranked prediction CSV"
    )
    kew_panel.add_argument(
        "--minimum-z", type=float, help="Optional z-score filter when --parent-file is CSV"
    )
    kew_panel.add_argument(
        "--top", type=int, default=0, help="Use only the top N parents from --parent-file"
    )
    kew_panel.add_argument("--backbone", required=True, type=Path)
    kew_panel.add_argument("--output", required=True, type=Path)
    kew_panel.add_argument("--allowed-taxa-file", type=Path)
    kew_panel.add_argument("--representatives", type=int, default=4)
    kew_panel.add_argument("--pool-multiplier", type=int, default=2)
    kew_panel.add_argument("--minimum-recovered-loci", type=int, default=50)
    kew_panel.add_argument("--minimum-score", type=float, default=0.50)
    kew_panel.add_argument("--download-workers", type=int, default=4)
    kew_panel.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Parallel locus workers for QC and alignment (use allocated HPC CPUs)",
    )
    kew_panel.add_argument("--python-aligner", action="store_true")
    kew_panel.add_argument(
        "--full-backbone-alignment",
        action="store_true",
        help="Use all curated backbone rows during MAFFT addition (comparison mode)",
    )

    prefetch = commands.add_parser(
        "prefetch-kew-panel",
        help="Download a balanced Kew candidate pool for a later offline HPC job",
    )
    prefetch.add_argument("--level", required=True, choices=["family", "genus", "species"])
    prefetch.add_argument("--parent-taxon", action="append", default=[])
    prefetch.add_argument("--parent-file", type=Path)
    prefetch.add_argument("--minimum-z", type=float)
    prefetch.add_argument("--top", type=int, default=0)
    prefetch.add_argument("--backbone", required=True, type=Path)
    prefetch.add_argument("--report", required=True, type=Path)
    prefetch.add_argument("--allowed-taxa-file", type=Path)
    prefetch.add_argument("--representatives", type=int, default=4)
    prefetch.add_argument("--pool-multiplier", type=int, default=2)
    prefetch.add_argument("--minimum-recovered-loci", type=int, default=50)
    prefetch.add_argument("--download-workers", type=int, default=4)

    branch_plan = commands.add_parser(
        "plan-genus-branches", help="Use Kew's species tree to plan large-family screening"
    )
    branch_plan.add_argument("--parent-file", required=True, type=Path)
    branch_plan.add_argument("--backbone", required=True, type=Path)
    branch_plan.add_argument("--output", required=True, type=Path)
    branch_plan.add_argument("--threshold", type=int, default=120)
    branch_plan.add_argument("--max-group-genera", type=int, default=40)
    branch_plan.add_argument("--max-anchor-genera", type=int, default=120)
    branch_plan.add_argument("--max-branches-per-family", type=int, default=6)
    branch_plan.add_argument("--allowed-genera-file", type=Path)

    branch_expand = commands.add_parser(
        "expand-genus-branches", help="Expand supported anchor branches to a genus list"
    )
    branch_expand.add_argument("--plan", required=True, type=Path)
    branch_expand.add_argument("--evidence", required=True, type=Path)
    branch_expand.add_argument("--output", required=True, type=Path)
    branch_expand.add_argument("--report", required=True, type=Path)

    shard_plan = commands.add_parser(
        "plan-genus-shards", help="Split a large genus set along Kew tree branches"
    )
    shard_plan.add_argument("--plan", required=True, type=Path)
    shard_plan.add_argument("--allowed-taxa-file", required=True, type=Path)
    shard_plan.add_argument("--output-dir", required=True, type=Path)
    shard_plan.add_argument("--max-genera-per-shard", type=int, default=40)

    shard_shortlist = commands.add_parser(
        "shortlist-genus-shards", help="Combine within-shard candidates for one small joint panel"
    )
    shard_shortlist.add_argument("--plan", required=True, type=Path)
    shard_shortlist.add_argument("--shard-runs", required=True, type=Path)
    shard_shortlist.add_argument("--output", required=True, type=Path)
    shard_shortlist.add_argument("--report", required=True, type=Path)
    shard_shortlist.add_argument("--max-final-genera", type=int, default=80)

    candidates = commands.add_parser(
        "select-candidates",
        help="Combine aggregate scores with per-query nearest-reference locus votes",
    )
    candidates.add_argument("--matrix-dir", required=True, type=Path)
    candidates.add_argument("--predictions", required=True, type=Path)
    candidates.add_argument("--level", required=True, choices=["order", "family", "genus"])
    candidates.add_argument("--project", required=True)
    candidates.add_argument("--output", required=True, type=Path)
    candidates.add_argument("--report", required=True, type=Path)
    candidates.add_argument("--minimum-locus-votes", type=int, default=2)
    candidates.add_argument("--minimum-z", type=float, default=0.0)
    candidates.add_argument("--always-top", type=int, default=1)
    candidates.add_argument("--max-candidates", type=int, default=12)

    demo = commands.add_parser("demo-mix7", help="Build the controlled mix7 hierarchy demo")
    demo.add_argument("--backbone", required=True, type=Path)
    demo.add_argument("--summary", required=True, type=Path)
    demo.add_argument("--output", required=True, type=Path)
    demo.add_argument("--python-aligner", action="store_true")

    ml = commands.add_parser("train-qc", help="Fit an optional IsolationForest on the backbone")
    ml.add_argument("--backbone", required=True, type=Path)
    ml.add_argument("--output", required=True, type=Path)

    serve = commands.add_parser("web", help="Start the local Web app")
    serve.add_argument("--host", default="127.0.0.1")
    serve.add_argument("--port", type=int, default=8000)
    return root


def _build_arguments(command: argparse.ArgumentParser, include_taxonomy: bool = False) -> None:
    command.add_argument("--taxon", required=True, help="Accepted binomial, e.g. Abatia rugosa")
    command.add_argument("--backbone", required=True, type=Path)
    command.add_argument("--output", required=True, type=Path)
    command.add_argument("--minimum-score", type=float, default=0.50)
    command.add_argument("--python-aligner", action="store_true", help="Do not invoke MAFFT")
    command.add_argument("--ml-model", type=Path, help="Optional trained IsolationForest joblib")
    command.add_argument("--ml-min-score", type=float, default=-0.05)
    command.add_argument("--no-archive", action="store_true", help="Skip the output ZIP")
    if include_taxonomy:
        command.add_argument("--order", required=True)
        command.add_argument("--family", required=True)


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    cache = _cache(args.cache)
    if args.command == "inspect-backbone":
        print(json.dumps(Backbone(args.backbone).metadata(), indent=2))
        return 0
    if args.command == "install-backbone":
        report = install_backbone_archive(
            args.archive,
            args.destination,
            overwrite=args.overwrite,
            expected_loci=args.expected_loci,
            taxonomy_corrections=args.taxonomy_corrections,
            apply_builtin_corrections=not args.no_taxonomy_corrections,
        )
        print(json.dumps(report, indent=2))
        return 0
    if args.command == "train-qc":
        print(json.dumps(train_anomaly_model(args.backbone, args.output), indent=2))
        return 0
    if args.command == "select-candidates":
        report = select_mixture_candidates(
            args.matrix_dir,
            args.predictions,
            args.output,
            args.report,
            level=args.level,
            project=args.project,
            minimum_locus_votes=args.minimum_locus_votes,
            minimum_z=args.minimum_z,
            always_top=args.always_top,
            max_candidates=args.max_candidates,
        )
        print(json.dumps(report, indent=2))
        return 0
    if args.command == "expand-genus-branches":
        print(
            json.dumps(
                expand_genus_branches(args.plan, args.evidence, args.output, args.report),
                indent=2,
            )
        )
        return 0
    if args.command == "plan-genus-shards":
        print(json.dumps(plan_genus_shards(
            args.plan, args.allowed_taxa_file, args.output_dir,
            max_genera_per_shard=args.max_genera_per_shard,
        ), indent=2))
        return 0
    if args.command == "shortlist-genus-shards":
        print(json.dumps(shortlist_genus_shards(
            args.plan, args.shard_runs, args.output, args.report,
            max_final_genera=args.max_final_genera,
        ), indent=2))
        return 0
    if args.command == "web":
        import uvicorn

        uvicorn.run("sprout_refbuilder.web:app", host=args.host, port=args.port, reload=False)
        return 0

    offline = os.environ.get("SPROUT_REF_OFFLINE", "").casefold() in {"1", "true", "yes"}
    with KewClient(cache, args.kew_base_url, offline=offline) as client:
        if args.command == "sync":
            print(json.dumps(client.sync(force=True), indent=2))
            return 0
        if args.command == "search":
            print(json.dumps([record.to_dict() for record in client.search(args.taxon, args.limit)], indent=2))
            return 0
        if args.command == "build-panel":
            report = HierarchicalPanelBuilder(
                args.backbone,
                args.summary,
                args.output,
                representatives_per_taxon=args.representatives,
                prefer_mafft=not args.python_aligner,
            ).build(
                args.level,
                client,
                parent_taxa=args.parent_taxon,
                include_species=args.include_species,
            )
            print(json.dumps(report, indent=2))
            return 0
        if args.command == "build-kew-panel":
            parents = list(args.parent_taxon)
            if args.parent_file:
                parents.extend(
                    read_taxon_candidates(args.parent_file, args.minimum_z, args.top)
                )
            parents = list(dict.fromkeys(parents))
            report = BalancedKewPanelBuilder(
                args.backbone,
                args.output,
                representatives_per_taxon=args.representatives,
                pool_multiplier=args.pool_multiplier,
                minimum_recovered_loci=args.minimum_recovered_loci,
                minimum_score=args.minimum_score,
                prefer_mafft=not args.python_aligner,
                download_workers=args.download_workers,
                workers=args.threads,
                full_backbone_alignment=args.full_backbone_alignment,
            ).build(
                args.level,
                parents,
                client,
                allowed_groups=set(read_taxon_candidates(args.allowed_taxa_file))
                if args.allowed_taxa_file
                else None,
            )
            print(json.dumps(report, indent=2))
            return 0
        if args.command == "prefetch-kew-panel":
            parents = list(args.parent_taxon)
            if args.parent_file:
                parents.extend(
                    read_taxon_candidates(args.parent_file, args.minimum_z, args.top)
                )
            parents = list(dict.fromkeys(parents))
            report = BalancedKewPanelBuilder(
                args.backbone,
                args.report.parent / ".prefetch-only",
                representatives_per_taxon=args.representatives,
                pool_multiplier=args.pool_multiplier,
                minimum_recovered_loci=args.minimum_recovered_loci,
                download_workers=args.download_workers,
            ).prefetch(
                args.level,
                parents,
                client,
                allowed_groups=set(read_taxon_candidates(args.allowed_taxa_file))
                if args.allowed_taxa_file
                else None,
            )
            args.report.parent.mkdir(parents=True, exist_ok=True)
            args.report.write_text(json.dumps(report, indent=2), encoding="utf-8")
            print(json.dumps(report, indent=2))
            return 0
        if args.command == "plan-genus-branches":
            report = plan_genus_branches(
                client,
                args.backbone,
                args.parent_file,
                args.output,
                threshold=args.threshold,
                max_group_genera=args.max_group_genera,
                max_anchor_genera=args.max_anchor_genera,
                max_branches_per_family=args.max_branches_per_family,
                allowed_genera_file=args.allowed_genera_file,
            )
            print(json.dumps({key: value for key, value in report.items() if key != "branches"}, indent=2))
            return 0
        if args.command == "demo-mix7":
            report = run_mix7_demo(
                args.backbone,
                args.summary,
                args.output,
                client,
                prefer_mafft=not args.python_aligner,
            )
            print(json.dumps(report, indent=2))
            return 0
        builder = ReferenceBuilder(
            args.backbone,
            args.output,
            minimum_score=args.minimum_score,
            prefer_mafft=not args.python_aligner,
            ml_model=args.ml_model,
            ml_min_score=args.ml_min_score,
            create_archive=not args.no_archive,
        )
        if args.command == "build":
            report = builder.build_from_kew(args.taxon, client)
        else:
            genus, epithet = args.taxon.replace("_", " ").split()[:2]
            taxonomy = Taxonomy(args.order, args.family, genus, epithet)
            if args.command == "build-local":
                recovery = args.recovery_fasta
            else:
                recovery = assemble_with_hybpiper(
                    args.read1,
                    args.read2,
                    args.target_file,
                    f"{genus}_{epithet}",
                    args.output.parent / f"{args.output.name}_hybpiper",
                    args.threads,
                )
            sequence_id = args.sequence_id if hasattr(args, "sequence_id") else "local"
            report = builder.build_from_fasta(args.taxon, taxonomy, recovery, sequence_id)
        print(json.dumps(report.to_dict(), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
