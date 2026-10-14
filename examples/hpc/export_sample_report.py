#!/usr/bin/env python3
"""Combine one MrKnowItAll sample into a readable and web-ready SPrOUT report.

New HPC stages already contain per-gene evidence. Completed older stages can be
backfilled from their rankings, with gene evidence explicitly marked unavailable.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path


def export_sample(sample_dir: Path, sprout: Path) -> tuple[Path, Path]:
    exporter = sprout / "export_pea.py"
    if not exporter.is_file():
        raise FileNotFoundError(f"SPrOUT exporter not found: {exporter}")
    project = sample_dir.name.lower()
    stage_files = []
    upstream_file = sample_dir / "UPSTREAM_SOURCES.json"
    upstream = (
        json.loads(upstream_file.read_text(encoding="utf-8"))
        if upstream_file.is_file() else {}
    )
    stages = [
        ("order", "order", 0.9, Path(upstream.get("order_results", sample_dir / "order" / "results"))),
        ("family", "family", 0.0, Path(upstream.get("family_results", sample_dir / "family" / "results"))),
    ]
    stages.extend(
        (f"branch_{round_dir.name}", "genus", 0.0, round_dir / "results")
        for round_dir in sorted((sample_dir / "branch").glob("round_*"))
        if round_dir.is_dir()
    )
    stages.extend([
        ("branch", "genus", 0.0, sample_dir / "branch" / "results"),
        ("genus", "genus", 0.0, sample_dir / "genus" / "results"),
    ])
    for stage_id, level, z, results in stages:
        if not (results / "SUCCESS").is_file():
            continue
        stage_file = results / "stage.pea.json"
        if not stage_file.is_file():
            if results.parent.parent != sample_dir:
                export_dir = sample_dir / "upstream_exports"
                export_dir.mkdir(exist_ok=True)
                stage_file = export_dir / f"{stage_id}.pea.json"
            predictions = results / "predictions.csv"
            if not predictions.is_file():
                raise FileNotFoundError(f"Completed stage lacks predictions: {results}")
            command = [
                sys.executable, str(exporter), "stage",
                "--project", project, "--level", level, "--stage-id", stage_id,
                "--predictions", str(predictions),
                "--z-threshold", str(z), "--output", str(stage_file),
                "--markdown", str(stage_file.parent / f"{stage_id}.report.md"),
            ]
            for flag, name in (
                ("--evidence", "candidate_evidence.csv"),
                ("--candidates", "candidates.txt"),
                ("--gene-list", "gene_list.txt"),
            ):
                if (results / name).is_file():
                    command.extend([flag, str(results / name)])
            subprocess.run(command, check=True)
        stage_files.append(stage_file)
    if not stage_files:
        raise RuntimeError(f"No completed stages in {sample_dir}")
    package = sample_dir / "RESULT.pea.json.gz"
    readable = sample_dir / "RESULT.md"
    subprocess.run(
        [sys.executable, str(exporter), "merge", "--project", project,
         "--input", *map(str, stage_files), "--output", str(package),
         "--markdown", str(readable)],
        check=True,
    )
    return package, readable


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-dir", required=True, type=Path)
    parser.add_argument("--sprout", required=True, type=Path)
    args = parser.parse_args()
    package, readable = export_sample(args.sample_dir, args.sprout)
    print(f"Web report: {package}")
    print(f"Readable report: {readable}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
