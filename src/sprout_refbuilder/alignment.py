from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path

from Bio import Align

from .backbone import LocusBackbone
from .fasta import FastaRecord, parse_fasta_text, read_fasta, ungap, write_fasta


def add_to_alignment(
    locus: LocusBackbone,
    candidate: FastaRecord,
    prefer_mafft: bool = True,
) -> tuple[list[FastaRecord], str]:
    if prefer_mafft and shutil.which("mafft"):
        return _mafft_add(locus, candidate), "mafft-addfragments-keeplength"
    projected = _project_pairwise(locus, candidate)
    return [*locus.records, projected], "biopython-pairwise-projection"


def add_many_to_alignment(
    locus: LocusBackbone,
    candidates: list[FastaRecord],
    prefer_mafft: bool = True,
) -> tuple[list[FastaRecord], str]:
    """Project several references into one curated locus alignment in one operation.

    Batch addition is important for clade panels: invoking MAFFT once per species would make a
    family panel need tens of thousands of processes.  The curated alignment is used only as a
    coordinate scaffold; callers may retain or discard its original records afterwards.
    """
    if not candidates:
        return list(locus.records), "none"
    if prefer_mafft and shutil.which("mafft"):
        return _mafft_add_many(locus, candidates), "mafft-addfragments-keeplength"
    projected = _project_many_pairwise(locus, candidates)
    return [*locus.records, *projected], "biopython-pairwise-projection"


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
            "--thread",
            "1",
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


def _mafft_add_many(locus: LocusBackbone, candidates: list[FastaRecord]) -> list[FastaRecord]:
    with tempfile.TemporaryDirectory(prefix="sprout-ref-batch-") as temporary:
        directory = Path(temporary)
        backbone_file = directory / "backbone.fasta"
        candidate_file = directory / "candidates.fasta"
        write_fasta(backbone_file, locus.records)
        write_fasta(candidate_file, candidates)
        command = [
            "mafft",
            "--quiet",
            "--thread",
            "1",
            "--preservecase",
            "--addfragments",
            str(candidate_file),
            "--keeplength",
            str(backbone_file),
        ]
        result = subprocess.run(command, check=True, capture_output=True, text=True)
        records = parse_fasta_text(result.stdout)
    if len(records) != len(locus.records) + len(candidates):
        raise RuntimeError(
            f"MAFFT added {len(records) - len(locus.records)} sequences instead of "
            f"{len(candidates)} for locus {locus.gene_id}"
        )
    return records


def _project_pairwise(locus: LocusBackbone, candidate: FastaRecord) -> FastaRecord:
    candidate_sequence = ungap(candidate.sequence)
    template = max(locus.records, key=lambda record: _quick_similarity(candidate_sequence, record.sequence))
    return _project_to_template(locus, candidate, template)


def _project_many_pairwise(
    locus: LocusBackbone, candidates: list[FastaRecord]
) -> list[FastaRecord]:
    """Pairwise fallback with cached reference k-mers for batch-panel construction."""
    k = 7
    reference_kmers = [(_kmer_set(record.sequence, k), record) for record in locus.records]
    projected: list[FastaRecord] = []
    for candidate in candidates:
        query = _kmer_set(candidate.sequence, k)
        template = max(
            reference_kmers,
            key=lambda item: _set_similarity(query, item[0]),
        )[1]
        projected.append(_project_to_template(locus, candidate, template))
    return projected


def _project_to_template(
    locus: LocusBackbone, candidate: FastaRecord, template: FastaRecord
) -> FastaRecord:
    candidate_sequence = ungap(candidate.sequence)
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
    k = 7
    return _set_similarity(_kmer_set(query, k), _kmer_set(aligned_reference, k))


def _kmer_set(sequence: str, k: int) -> set[str]:
    clean = ungap(sequence)
    return {clean[i : i + k] for i in range(max(0, len(clean) - k + 1))}


def _set_similarity(left: set[str], right: set[str]) -> float:
    intersection = len(left & right)
    union_size = len(left) + len(right) - intersection
    return intersection / union_size if union_size else 0.0
