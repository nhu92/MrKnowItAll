from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

from .assembly import assemble_with_hybpiper
from .backbone import Backbone
from .builder import ReferenceBuilder
from .kew import DEFAULT_BASE_URL, KewClient
from .ml import train_anomaly_model
from .models import Taxonomy


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

    ml = commands.add_parser("train-qc", help="Fit an optional IsolationForest on the 871 backbone")
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
    if include_taxonomy:
        command.add_argument("--order", required=True)
        command.add_argument("--family", required=True)


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    cache = _cache(args.cache)
    if args.command == "inspect-backbone":
        print(json.dumps(Backbone(args.backbone).metadata(), indent=2))
        return 0
    if args.command == "train-qc":
        print(json.dumps(train_anomaly_model(args.backbone, args.output), indent=2))
        return 0
    if args.command == "web":
        import uvicorn

        uvicorn.run("sprout_refbuilder.web:app", host=args.host, port=args.port, reload=False)
        return 0

    with KewClient(cache, args.kew_base_url) as client:
        if args.command == "sync":
            print(json.dumps(client.sync(force=True), indent=2))
            return 0
        if args.command == "search":
            print(json.dumps([record.to_dict() for record in client.search(args.taxon, args.limit)], indent=2))
            return 0
        builder = ReferenceBuilder(
            args.backbone,
            args.output,
            minimum_score=args.minimum_score,
            prefer_mafft=not args.python_aligner,
            ml_model=args.ml_model,
            ml_min_score=args.ml_min_score,
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
