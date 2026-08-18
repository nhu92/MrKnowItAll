from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path

from Bio import Align

from .backbone import LocusBackbone
from .fasta import FastaRecord, read_fasta, ungap, write_fasta


def add_to_alignment(
    locus: LocusBackbone,
    candidate: FastaRecord,
    prefer_mafft: bool = True,
) -> tuple[list[FastaRecord], str]:
    if prefer_mafft and shutil.which("mafft"):
        return _mafft_add(locus, candidate), "mafft-addfragments-keeplength"
    projected = _project_pairwise(locus, candidate)
    return [*locus.records, projected], "biopython-pairwise-projection"


def _mafft_add(locus: LocusBackbone, candidate: FastaRecord) -> list[FastaRecord]:
    with tempfile.TemporaryDirectory(prefix="sprout-ref-") as temporary:
        directory = Path(temporary)
        backbone_file = directory / "backbone.fasta"
        candidate_file = directory / "candidate.fasta"
        output_file = directory / "output.fasta"
        write_fasta(backbone_file, locus.records)
        write_fasta(candidate_file, [candidate])
        command = [
            "mafft",
            "--quiet",
            "--preservecase",
            "--addfragments",
            str(candidate_file),
            "--keeplength",
            str(backbone_file),
        ]
        result = subprocess.run(command, check=True, capture_output=True, text=True)
        output_file.write_text(result.stdout, encoding="utf-8")
        records = read_fasta(output_file)
    if len(records) != len(locus.records) + 1:
        raise RuntimeError(f"MAFFT did not add exactly one sequence for locus {locus.gene_id}")
    return records


def _project_pairwise(locus: LocusBackbone, candidate: FastaRecord) -> FastaRecord:
    candidate_sequence = ungap(candidate.sequence)
    template = max(locus.records, key=lambda record: _quick_similarity(candidate_sequence, record.sequence))
    template_sequence = ungap(template.sequence)
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -4.0
    aligner.extend_gap_score = -0.5
    alignment = aligner.align(template_sequence, candidate_sequence)[0]

    template_columns = [i for i, base in enumerate(template.sequence) if base != "-"]
    projected = ["-"] * locus.alignment_length
    coordinates = alignment.coordinates
    for segment in range(coordinates.shape[1] - 1):
        t0, t1 = int(coordinates[0, segment]), int(coordinates[0, segment + 1])
        c0, c1 = int(coordinates[1, segment]), int(coordinates[1, segment + 1])
        if t1 > t0 and c1 > c0:
            width = min(t1 - t0, c1 - c0)
            for offset in range(width):
                projected[template_columns[t0 + offset]] = candidate_sequence[c0 + offset]
    return FastaRecord(candidate.id, candidate.description, "".join(projected))


def _quick_similarity(query: str, aligned_reference: str) -> float:
    reference = ungap(aligned_reference)
    k = 7
    left = {query[i : i + k] for i in range(max(0, len(query) - k + 1))}
    right = {reference[i : i + k] for i in range(max(0, len(reference) - k + 1))}
    union = left | right
    return len(left & right) / len(union) if union else 0.0

