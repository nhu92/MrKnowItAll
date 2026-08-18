from __future__ import annotations

import re
from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from pathlib import Path

DNA_RE = re.compile(r"^[ACGTRYSWKMBDHVN.\-]+$", re.IGNORECASE)


@dataclass(frozen=True)
class FastaRecord:
    id: str
    description: str
    sequence: str


def parse_fasta_text(text: str) -> list[FastaRecord]:
    records: list[FastaRecord] = []
    header: str | None = None
    chunks: list[str] = []
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append(_record(header, chunks))
            header, chunks = line[1:].strip(), []
        elif header is None:
            raise ValueError("FASTA sequence encountered before a header")
        else:
            chunks.append(line)
    if header is not None:
        records.append(_record(header, chunks))
    return records


def _record(header: str, chunks: list[str]) -> FastaRecord:
    sequence = "".join(chunks).upper().replace(".", "-")
    if not sequence or not DNA_RE.fullmatch(sequence):
        raise ValueError(f"Invalid or empty DNA sequence for {header!r}")
    return FastaRecord(header.split()[0], header, sequence)


def read_fasta(path: Path) -> list[FastaRecord]:
    return parse_fasta_text(path.read_text(encoding="utf-8"))


def format_fasta(records: Iterable[FastaRecord], width: int = 80) -> str:
    lines: list[str] = []
    for record in records:
        lines.append(f">{record.description}")
        lines.extend(record.sequence[i : i + width] for i in range(0, len(record.sequence), width))
    return "\n".join(lines) + ("\n" if lines else "")


def write_fasta(path: Path, records: Iterable[FastaRecord]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(format_fasta(records), encoding="utf-8")


def iter_fasta_files(directory: Path) -> Iterator[Path]:
    for pattern in ("*.fasta", "*.fa", "*.fas", "*.fna"):
        yield from sorted(directory.glob(pattern))


def ungap(sequence: str) -> str:
    return sequence.replace("-", "").replace(".", "")


def parse_kew_gene_id(record: FastaRecord) -> str:
    # Kew: >gene_id Gene_Name:... Species:... Repository:... Sequence_ID:...
    return record.id

