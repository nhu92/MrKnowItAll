from __future__ import annotations

import math
from collections.abc import Iterable
from functools import lru_cache

from .backbone import LocusBackbone
from .fasta import FastaRecord, ungap
from .models import QCResult, SequenceFeatures, Taxonomy
from .taxonomy import taxonomy_from_header

STOP_CODONS = {"TAA", "TAG", "TGA"}
COMPLEMENT = str.maketrans("ACGTRYMKSWBDHVN", "TGCAYRKMSWVHDBN")


@lru_cache(maxsize=4096)
def kmers(sequence: str, k: int = 9) -> frozenset[str]:
    clean = ungap(sequence).upper()
    return frozenset(
        clean[index : index + k]
        for index in range(max(0, len(clean) - k + 1))
        if set(clean[index : index + k]) <= {"A", "C", "G", "T"}
    )


def jaccard(left: frozenset[str], right: frozenset[str]) -> float:
    intersection = len(left & right)
    union_size = len(left) + len(right) - intersection
    return intersection / union_size if union_size else 0.0


def nearest_similarity(sequence: str, references: Iterable[FastaRecord], k: int = 9) -> float:
    query = kmers(sequence, k)
    return max((jaccard(query, kmers(record.sequence, k)) for record in references), default=0.0)


def clade_support(
    sequence: str,
    references: Iterable[FastaRecord],
    expected: Taxonomy | None,
    k: int = 9,
    neighbors: int = 25,
) -> tuple[float, float, float, float]:
    """Return nearest similarity and weighted order/family/genus neighbor support.

    A value of -1 means the expected clade is absent from this backbone locus, so it must not
    be used to reject the candidate. This is a fast alignment-neighborhood proxy; production
    deployments can add EPA-ng placement as a stricter second stage.
    """
    query = kmers(sequence, k)
    scored: list[tuple[float, Taxonomy]] = []
    all_taxa: list[Taxonomy] = []
    for record in references:
        taxonomy = taxonomy_from_header(record.description)
        if taxonomy is None:
            continue
        all_taxa.append(taxonomy)
        scored.append((jaccard(query, kmers(record.sequence, k)), taxonomy))
    scored.sort(key=lambda item: item[0], reverse=True)
    nearest = scored[0][0] if scored else 0.0
    if expected is None:
        return nearest, -1.0, -1.0, -1.0
    top = scored[:neighbors]
    total = sum(score for score, _ in top)
    if total <= 0:
        return nearest, 0.0, 0.0, 0.0
    has_order = any(t.order.casefold() == expected.order.casefold() for t in all_taxa)
    has_family = any(t.family.casefold() == expected.family.casefold() for t in all_taxa)
    has_genus = any(t.genus.casefold() == expected.genus.casefold() for t in all_taxa)
    order = sum(
        score for score, taxon in top if taxon.order.casefold() == expected.order.casefold()
    ) / total
    family = sum(
        score for score, taxon in top if taxon.family.casefold() == expected.family.casefold()
    ) / total
    genus = sum(
        score for score, taxon in top if taxon.genus.casefold() == expected.genus.casefold()
    ) / total
    return (
        nearest,
        order if has_order else -1.0,
        family if has_family else -1.0,
        genus if has_genus else -1.0,
    )


def best_stop_fraction(sequence: str) -> float:
    clean = "".join(base for base in ungap(sequence).upper() if base in "ACGT")
    if len(clean) < 3:
        return 1.0
    reverse = clean.translate(COMPLEMENT)[::-1]
    fractions = []
    for strand in (clean, reverse):
        for frame in range(3):
            codons = [strand[i : i + 3] for i in range(frame, len(strand) - 2, 3)]
            fractions.append(sum(c in STOP_CODONS for c in codons) / max(1, len(codons)))
    return min(fractions)


def evaluate_candidate(
    candidate: FastaRecord,
    locus: LocusBackbone,
    duplicate_count: int = 1,
    minimum_score: float = 0.50,
    expected_taxonomy: Taxonomy | None = None,
) -> QCResult:
    sequence = ungap(candidate.sequence)
    length_ratio = len(sequence) / max(1.0, locus.median_ungapped_length)
    ambiguity = sum(base not in "ACGT" for base in sequence.upper()) / max(1, len(sequence))
    similarity, order_support, family_support, genus_support = clade_support(
        sequence, locus.records, expected_taxonomy
    )
    stops = best_stop_fraction(sequence)
    features = SequenceFeatures(
        length_ratio,
        ambiguity,
        similarity,
        stops,
        order_support,
        family_support,
        genus_support,
        duplicate_count,
    )

    reasons: list[str] = []
    if length_ratio < 0.35:
        reasons.append("sequence shorter than 35% of the backbone locus median")
    if length_ratio > 2.0:
        reasons.append("sequence longer than 200% of the backbone locus median")
    if ambiguity > 0.15:
        reasons.append("more than 15% ambiguous bases")
    if similarity < 0.015:
        reasons.append("no credible 9-mer neighborhood in the curated backbone")
    if stops > 0.08:
        reasons.append("excess stop codons in every translated frame")
    if order_support >= 0 and order_support < 0.15:
        reasons.append("top backbone neighbors give less than 15% support to the expected order")
    # Historical backbones can use a former family circumscription (for example Allium in
    # Alliaceae versus current APG/Kew Asparagaceae). Strong same-genus support supersedes a
    # family-name mismatch; otherwise the family gate remains active.
    if family_support >= 0 and family_support < 0.05 and genus_support < 0.05:
        reasons.append("top backbone neighbors give less than 5% support to the expected family")

    length_score = math.exp(-abs(math.log(max(length_ratio, 1e-6))))
    ambiguity_score = max(0.0, 1.0 - ambiguity / 0.15)
    stop_score = max(0.0, 1.0 - stops / 0.08)
    duplicate_penalty = 1.0 / math.sqrt(max(1, duplicate_count))
    clade_score = max(
        value for value in (order_support, family_support, genus_support, 0.5) if value >= 0
    )
    score = (
        0.25 * length_score
        + 0.30 * min(1.0, similarity / 0.20)
        + 0.18 * ambiguity_score
        + 0.10 * stop_score
        + 0.12 * clade_score
        + 0.05 * duplicate_penalty
    )
    if score < minimum_score:
        reasons.append(f"composite QC score {score:.3f} is below {minimum_score:.3f}")
    return QCResult(not reasons, round(score, 6), tuple(reasons), features)


def choose_candidate(
    candidates: list[FastaRecord],
    locus: LocusBackbone,
    minimum_score: float = 0.50,
    expected_taxonomy: Taxonomy | None = None,
) -> tuple[FastaRecord | None, QCResult | None, list[tuple[FastaRecord, QCResult]]]:
    evaluated = [
        (
            candidate,
            evaluate_candidate(
                candidate, locus, len(candidates), minimum_score, expected_taxonomy
            ),
        )
        for candidate in candidates
    ]
    evaluated.sort(key=lambda item: item[1].score, reverse=True)
    if not evaluated:
        return None, None, []
    best_candidate, best_qc = evaluated[0]
    return (best_candidate if best_qc.accepted else None), best_qc, evaluated
