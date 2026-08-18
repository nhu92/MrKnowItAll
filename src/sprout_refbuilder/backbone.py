from __future__ import annotations

import csv
import hashlib
import json
import re
import shutil
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from statistics import median
from zipfile import ZipFile

from .fasta import FastaRecord, iter_fasta_files, read_fasta, ungap

FASTA_SUFFIXES = {".fasta", ".fa", ".fas", ".fna"}
SUMMARY_FILES = {"gene_summary.csv", "species_summary.csv"}


@dataclass(frozen=True)
class LocusBackbone:
    gene_id: str
    path: Path
    records: tuple[FastaRecord, ...]
    alignment_length: int
    median_ungapped_length: float


class Backbone:
    def __init__(self, directory: Path) -> None:
        if not directory.is_dir():
            raise FileNotFoundError(f"Backbone directory does not exist: {directory}")
        loci: dict[str, LocusBackbone] = {}
        for path in iter_fasta_files(directory):
            records = tuple(read_fasta(path))
            if not records:
                continue
            widths = {len(record.sequence) for record in records}
            if len(widths) != 1:
                raise ValueError(f"Backbone file is not aligned (unequal widths): {path}")
            gene_id = path.stem
            lengths = [len(ungap(record.sequence)) for record in records]
            loci[gene_id] = LocusBackbone(
                gene_id=gene_id,
                path=path,
                records=records,
                alignment_length=widths.pop(),
                median_ungapped_length=median(lengths),
            )
        if not loci:
            raise ValueError(f"No FASTA alignments were found in {directory}")
        self.directory = directory
        self.loci = loci

    def metadata(self) -> dict:
        all_records = [record for locus in self.loci.values() for record in locus.records]
        specimens = {re.sub(r"-\d+$", "", record.id) for record in all_records}
        species = {"_".join(record.id.split("_")[2:4]) for record in all_records}
        metadata = {
            "directory": str(self.directory),
            "locus_count": len(self.loci),
            "sequence_count": len(all_records),
            "specimen_count": len(specimens),
            "species_count": len(species),
            "loci": sorted(self.loci),
        }
        summary = self.directory / "species_summary.csv"
        if summary.is_file():
            with summary.open(encoding="utf-8-sig", newline="") as source:
                rows = list(csv.DictReader(source))
            if rows and all("Species" in row for row in rows):
                metadata.update(
                    specimen_count=len(rows),
                    species_count=len({row["Species"] for row in rows}),
                    order_count=len({row["Order"] for row in rows if row.get("Order")}),
                    family_count=len({row["Family"] for row in rows if row.get("Family")}),
                )
        return metadata

    def write_metadata(self, path: Path) -> None:
        path.write_text(json.dumps(self.metadata(), indent=2), encoding="utf-8")


def _archive_members(names: Sequence[str]) -> tuple[list[str], list[str]]:
    fasta_members: list[str] = []
    summary_members: list[str] = []
    basenames: set[str] = set()
    for name in names:
        archive_path = Path(name)
        if archive_path.name in {"", ".", ".."}:
            continue
        is_fasta = archive_path.suffix.lower() in FASTA_SUFFIXES
        is_summary = archive_path.name in SUMMARY_FILES
        if not is_fasta and not is_summary:
            continue
        if archive_path.name in basenames:
            raise ValueError(f"Archive has duplicate output filename: {archive_path.name}")
        basenames.add(archive_path.name)
        (fasta_members if is_fasta else summary_members).append(name)
    return sorted(fasta_members), sorted(summary_members)


def install_backbone_archive(
    archive: Path,
    destination: Path,
    *,
    overwrite: bool = False,
    expected_loci: int = 353,
) -> dict:
    """Install a curated alignment ZIP without preserving untrusted archive paths."""
    archive = archive.resolve()
    if not archive.is_file():
        raise FileNotFoundError(f"Backbone archive does not exist: {archive}")

    digest = hashlib.sha256()
    with archive.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)

    with ZipFile(archive) as bundle:
        fasta_members, summary_members = _archive_members(bundle.namelist())
        if len(fasta_members) != expected_loci:
            raise ValueError(
                f"Expected {expected_loci} FASTA loci, found {len(fasta_members)} in {archive.name}"
            )
        outputs = [destination / Path(member).name for member in fasta_members + summary_members]
        conflicts = [path for path in outputs if path.exists()]
        if conflicts and not overwrite:
            preview = ", ".join(path.name for path in conflicts[:5])
            raise FileExistsError(f"Destination files already exist ({preview}); use --overwrite")

        destination.mkdir(parents=True, exist_ok=True)
        for member, output in zip(fasta_members + summary_members, outputs, strict=True):
            with bundle.open(member) as source, output.open("wb") as target:
                shutil.copyfileobj(source, target)

    metadata = Backbone(destination).metadata()
    metadata.pop("directory", None)
    metadata.pop("loci", None)
    manifest = {
        "schema_version": 1,
        "source_archive": archive.name,
        "source_sha256": digest.hexdigest(),
        "summary_files": [Path(member).name for member in summary_members],
        **metadata,
    }
    (destination / "backbone_install.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    return manifest
