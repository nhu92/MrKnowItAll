from __future__ import annotations

import csv
import math
from collections import defaultdict
from pathlib import Path
from statistics import median


LEVEL_INDEX = {"order": 0, "family": 1, "genus": 2}


def taxon_from_sprout_header(header: str, level: str) -> str | None:
    """Extract a taxon from a SPrOUT reference leaf name."""
    value = header.strip().strip("'")
    if value.startswith("_R_"):
        value = value[3:]
    parts = value.split("_")
    index = LEVEL_INDEX[level]
    return parts[index] if len(parts) > index and parts[index] else None


def select_mixture_candidates(
    matrix_dir: Path,
    predictions_csv: Path,
    output_path: Path,
    report_path: Path,
    *,
    level: str,
    project: str,
    minimum_locus_votes: int = 2,
    minimum_z: float = 0.0,
    always_top: int = 1,
    max_candidates: int = 12,
) -> dict:
    """Select high-recall hierarchical candidates for a mixed sample.

    SPrOUT's aggregate score is useful for dominant components but a global z-score can suppress
    low-abundance mixture members. This selector unions that evidence with nearest-reference
    votes from individual query exons, counting a taxon at most once per tree/locus matrix.
    """
    if level not in LEVEL_INDEX:
        raise ValueError(f"level must be one of {sorted(LEVEL_INDEX)}")
    if minimum_locus_votes < 1 or always_top < 0 or max_candidates < 1:
        raise ValueError("candidate selection counts must be positive")
    matrices = sorted(matrix_dir.glob("*.cleaned.csv"))
    if not matrices:
        raise FileNotFoundError(f"No cleaned distance matrices found in {matrix_dir}")

    query_votes: dict[str, int] = defaultdict(int)
    locus_votes: dict[str, int] = defaultdict(int)
    winning_distances: dict[str, list[float]] = defaultdict(list)
    query_columns_seen = 0
    for matrix in matrices:
        with matrix.open(encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle)
            if not reader.fieldnames:
                continue
            label_column = reader.fieldnames[0]
            query_columns = [name for name in reader.fieldnames[1:] if project in name]
            rows = list(reader)
        query_columns_seen += len(query_columns)
        locus_winners: set[str] = set()
        for query in query_columns:
            scored: list[tuple[float, str]] = []
            for row in rows:
                label = (row.get(label_column) or "").strip()
                if not label or project in label or "NODE" in label:
                    continue
                taxon = taxon_from_sprout_header(label, level)
                try:
                    distance = float(row.get(query, ""))
                except (TypeError, ValueError):
                    continue
                if taxon and math.isfinite(distance):
                    scored.append((distance, taxon))
            if not scored:
                continue
            best_distance = min(distance for distance, _ in scored)
            winners = {
                taxon
                for distance, taxon in scored
                if math.isclose(distance, best_distance, rel_tol=1e-9, abs_tol=1e-12)
            }
            for taxon in winners:
                query_votes[taxon] += 1
                locus_winners.add(taxon)
                winning_distances[taxon].append(best_distance)
        for taxon in locus_winners:
            locus_votes[taxon] += 1

    if query_columns_seen == 0:
        raise RuntimeError(
            f"No distance-matrix query column contains project token {project!r}"
        )

    predictions: dict[str, dict[str, float | int]] = {}
    with predictions_csv.open(encoding="utf-8-sig", newline="") as handle:
        for rank, row in enumerate(csv.DictReader(handle), start=1):
            taxon = (row.get("row_name") or "").strip()
            if not taxon:
                continue
            score = _optional_float(row.get("sum_of_total_value")) or 0.0
            z_score = _optional_float(row.get("z_score"))
            predictions[taxon] = {"rank": rank, "score": score, "z_score": z_score}

    all_taxa = set(predictions) | set(locus_votes)
    top_taxa = {
        taxon
        for taxon, _ in sorted(
            predictions.items(), key=lambda item: int(item[1]["rank"])
        )[:always_top]
    }
    selected_by: dict[str, list[str]] = defaultdict(list)
    for taxon in all_taxa:
        if taxon in top_taxa:
            selected_by[taxon].append("aggregate_top")
        prediction = predictions.get(taxon, {})
        z_score = prediction.get("z_score")
        if isinstance(z_score, float) and z_score >= minimum_z:
            selected_by[taxon].append("aggregate_z")
        if locus_votes[taxon] >= minimum_locus_votes:
            selected_by[taxon].append("locus_votes")

    def priority(taxon: str) -> tuple:
        prediction = predictions.get(taxon, {})
        return (
            -locus_votes[taxon],
            -query_votes[taxon],
            -float(prediction.get("score", 0.0)),
            int(prediction.get("rank", 10**9)),
            taxon.casefold(),
        )

    selected = sorted((taxon for taxon in all_taxa if selected_by[taxon]), key=priority)
    selected = selected[:max_candidates]
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("".join(f"{taxon}\n" for taxon in selected), encoding="utf-8")

    rows = []
    for taxon in sorted(all_taxa, key=priority):
        prediction = predictions.get(taxon, {})
        rows.append(
            {
                "taxon": taxon,
                "selected": taxon in selected,
                "selected_by": ";".join(selected_by[taxon]),
                "locus_votes": locus_votes[taxon],
                "query_votes": query_votes[taxon],
                "median_winning_distance": (
                    median(winning_distances[taxon]) if winning_distances[taxon] else ""
                ),
                "aggregate_rank": prediction.get("rank", ""),
                "aggregate_score": prediction.get("score", ""),
                "aggregate_z": prediction.get("z_score", ""),
            }
        )
    with report_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return {
        "level": level,
        "matrices": len(matrices),
        "query_columns": query_columns_seen,
        "selected": selected,
        "report": str(report_path),
    }


def _optional_float(value: str | None) -> float | None:
    if value in (None, ""):
        return None
    try:
        result = float(value)
    except ValueError:
        return None
    return result if math.isfinite(result) else None
