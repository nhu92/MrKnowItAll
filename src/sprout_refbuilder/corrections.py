from __future__ import annotations

import csv
import hashlib
from collections import Counter
from dataclasses import dataclass
from importlib.resources import files
from pathlib import Path

from .fasta import FastaRecord, iter_fasta_files, read_fasta, write_fasta


@dataclass(frozen=True)
class TaxonomyCorrection:
    action: str
    original: str
    corrected: str
    notes: str


def builtin_corrections_resource():
    return files("sprout_refbuilder").joinpath("resources/taxonomy_corrections.tsv")


def load_taxonomy_corrections(path: Path | None = None) -> tuple[list[TaxonomyCorrection], str]:
    resource = path or builtin_corrections_resource()
    digest = hashlib.sha256()
    rows: list[TaxonomyCorrection] = []
    with resource.open("rb") as raw:
        content = raw.read()
    digest.update(content)
    text = content.decode("utf-8-sig")
    reader = csv.DictReader(text.splitlines(), delimiter="\t")
    required = {"action", "original", "corrected", "notes"}
    if not reader.fieldnames or not required.issubset(reader.fieldnames):
        raise ValueError(f"Taxonomy correction file must contain {sorted(required)}")
    for row in reader:
        correction = TaxonomyCorrection(
            row["action"].strip().casefold(),
            row["original"].strip(),
            row["corrected"].strip(),
            row["notes"].strip(),
        )
        if correction.action not in {"replace", "exclude"}:
            raise ValueError(f"Unknown taxonomy correction action: {correction.action}")
        if not correction.original or (correction.action == "replace" and not correction.corrected):
            raise ValueError(f"Incomplete taxonomy correction: {row}")
        rows.append(correction)
    return rows, digest.hexdigest()


def apply_taxonomy_corrections(
    directory: Path, corrections: list[TaxonomyCorrection]
) -> dict:
    replacements = [row for row in corrections if row.action == "replace"]
    exclusions = [row for row in corrections if row.action == "exclude"]
    replacement_hits: Counter[str] = Counter()
    exclusion_hits: Counter[str] = Counter()

    for path in iter_fasta_files(directory):
        output: list[FastaRecord] = []
        changed = False
        for record in read_fasta(path):
            excluded = next(
                (row for row in exclusions if _prefix_suffix(record.id, row.original) is not None),
                None,
            )
            if excluded:
                exclusion_hits[excluded.original] += 1
                changed = True
                continue
            updated = record
            for row in replacements:
                suffix = _prefix_suffix(updated.id, row.original)
                if suffix is None:
                    continue
                new_id = row.corrected + suffix
                _, separator, remainder = updated.description.partition(" ")
                new_description = new_id + (separator + remainder if separator else "")
                updated = FastaRecord(new_id, new_description, updated.sequence)
                replacement_hits[row.original] += 1
                changed = True
                break
            output.append(updated)
        if changed:
            write_fasta(path, output)

    _write_summaries(directory)
    shadowed_replacements = sorted(
        row.original
        for row in replacements
        if not replacement_hits[row.original]
        and any(
            exclusion_hits[excluded.original]
            and _prefix_suffix(excluded.original, row.original) is not None
            for excluded in exclusions
        )
    )
    return {
        "rules": len(corrections),
        "replacement_rules": len(replacements),
        "exclusion_rules": len(exclusions),
        "replaced_locus_sequences": sum(replacement_hits.values()),
        "excluded_locus_sequences": sum(exclusion_hits.values()),
        "shadowed_by_exclusion": shadowed_replacements,
        "unmatched_replacements": sorted(
            row.original
            for row in replacements
            if not replacement_hits[row.original] and row.original not in shadowed_replacements
        ),
        "unmatched_exclusions": sorted(
            row.original for row in exclusions if not exclusion_hits[row.original]
        ),
    }


def _prefix_suffix(identifier: str, prefix: str) -> str | None:
    if identifier == prefix:
        return ""
    for separator in ("_", "-"):
        marker = prefix + separator
        if identifier.startswith(marker):
            return identifier[len(prefix) :]
    return None


def _write_summaries(directory: Path) -> None:
    specimen_loci: Counter[tuple[str, str, str, str]] = Counter()
    gene_rows: list[dict[str, str | int]] = []
    for path in iter_fasta_files(directory):
        records = read_fasta(path)
        if records:
            gene_rows.append(
                {
                    "Gene": path.stem,
                    "Aligned Sequences": len(records),
                    "Alignment Length (1st Seq)": len(records[0].sequence),
                }
            )
        seen: set[tuple[str, str, str, str]] = set()
        for record in records:
            parsed = _summary_key(record.id)
            if parsed:
                seen.add(parsed)
        specimen_loci.update(seen)

    with (directory / "gene_summary.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["Gene", "Aligned Sequences", "Alignment Length (1st Seq)"]
        )
        writer.writeheader()
        writer.writerows(gene_rows)
    with (directory / "species_summary.csv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["Species", "Order", "Family", "Data Source", "Gene Count"]
        )
        writer.writeheader()
        for (species, order, family, source), count in sorted(specimen_loci.items()):
            writer.writerow(
                {
                    "Species": species,
                    "Order": order,
                    "Family": family,
                    "Data Source": source,
                    "Gene Count": count,
                }
            )


def _summary_key(identifier: str) -> tuple[str, str, str, str] | None:
    fields = identifier.split("_", 4)
    if len(fields) != 5:
        return None
    order, family, genus, species, _tail = fields
    source = data_source_from_identifier(identifier)
    return f"{genus}_{species}", order, family, source


def data_source_from_identifier(identifier: str) -> str:
    fields = identifier.split("_", 4)
    if len(fields) != 5:
        return ""
    raw_source = fields[4].rsplit("-", 1)[0]
    return "Timilsena" if raw_source.startswith("Timilsena") else raw_source
