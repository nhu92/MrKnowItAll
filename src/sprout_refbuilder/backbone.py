from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from statistics import median

from .fasta import FastaRecord, iter_fasta_files, read_fasta, ungap


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
        species = {record.id for locus in self.loci.values() for record in locus.records}
        return {
            "directory": str(self.directory),
            "locus_count": len(self.loci),
            "header_count": len(species),
            "loci": sorted(self.loci),
        }

    def write_metadata(self, path: Path) -> None:
        path.write_text(json.dumps(self.metadata(), indent=2), encoding="utf-8")

