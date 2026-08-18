from __future__ import annotations

from dataclasses import asdict, dataclass, field
from pathlib import Path


@dataclass(frozen=True)
class KewRecord:
    repository: str
    sequence_id: str
    sequence_type: str
    species: str
    project: str
    recovery_file: str | None = None

    def to_dict(self) -> dict[str, str | None]:
        return asdict(self)


@dataclass(frozen=True)
class Taxonomy:
    order: str
    family: str
    genus: str
    species: str

    @property
    def sprout_prefix(self) -> str:
        return f"{self.order}_{self.family}_{self.genus}_{self.species}"


@dataclass(frozen=True)
class SequenceFeatures:
    length_ratio: float
    ambiguity_fraction: float
    nearest_kmer_similarity: float
    stop_codon_fraction: float
    order_support: float
    family_support: float
    genus_support: float
    duplicate_count: int

    def vector(self) -> list[float]:
        return [
            self.length_ratio,
            self.ambiguity_fraction,
            self.nearest_kmer_similarity,
            self.stop_codon_fraction,
            self.order_support,
            self.family_support,
            self.genus_support,
            float(self.duplicate_count),
        ]


@dataclass(frozen=True)
class QCResult:
    accepted: bool
    score: float
    reasons: tuple[str, ...]
    features: SequenceFeatures

    def to_dict(self) -> dict:
        result = asdict(self)
        result["reasons"] = list(self.reasons)
        return result


@dataclass
class BuildReport:
    requested_taxon: str
    resolved_taxon: str
    source: str
    output_dir: Path
    release: str = "unknown"
    accepted_loci: list[str] = field(default_factory=list)
    rejected_loci: dict[str, list[str]] = field(default_factory=dict)
    missing_loci: list[str] = field(default_factory=list)
    provenance: dict[str, str] = field(default_factory=dict)

    def to_dict(self) -> dict:
        return {
            **asdict(self),
            "output_dir": str(self.output_dir),
        }
