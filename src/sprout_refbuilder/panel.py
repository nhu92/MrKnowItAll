from __future__ import annotations

import csv
import json
import shutil
import tempfile
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import ClassVar

from .builder import ReferenceBuilder
from .fasta import FastaRecord, iter_fasta_files, read_fasta, write_fasta
from .kew import KewClient
from .models import Taxonomy
from .taxonomy import kew_taxonomy_index, taxonomy_from_header


@dataclass(frozen=True)
class PanelCandidate:
    species_name: str
    taxonomy: Taxonomy
    source: str
    gene_count: int

    @property
    def group_species(self) -> str:
        return f"{self.taxonomy.genus}_{self.taxonomy.species}"


class HierarchicalPanelBuilder:
    """Subsample a curated alignment skeleton and optionally add exact Kew recoveries."""

    LEVELS: ClassVar[frozenset[str]] = frozenset({"order", "family", "genus"})

    def __init__(
        self,
        backbone_dir: Path,
        summary_csv: Path,
        output_dir: Path,
        representatives_per_taxon: int = 2,
        prefer_mafft: bool = True,
    ) -> None:
        self.backbone_dir = backbone_dir
        self.summary_csv = summary_csv
        self.output_dir = output_dir
        self.representatives_per_taxon = representatives_per_taxon
        self.prefer_mafft = prefer_mafft

    def build(
        self,
        level: str,
        client: KewClient,
        parent_taxa: list[str] | None = None,
        include_species: list[str] | None = None,
    ) -> dict:
        level = level.casefold()
        if level not in self.LEVELS:
            raise ValueError(f"level must be one of {sorted(self.LEVELS)}")
        if self.output_dir.exists() and any(self.output_dir.iterdir()):
            raise FileExistsError(f"Output directory is not empty: {self.output_dir}")
        self.output_dir.parent.mkdir(parents=True, exist_ok=True)
        metadata = client.sync()
        taxonomy_index = kew_taxonomy_index(client.fetch_species_tree())
        candidates = self._load_candidates(taxonomy_index)
        selected = self._select(candidates, level, parent_taxa or [])
        augmentations: list[dict] = []

        with tempfile.TemporaryDirectory(
            prefix=f"sprout-{level}-panel-", dir=self.output_dir.parent
        ) as temporary:
            temporary_dir = Path(temporary)
            base_ref = temporary_dir / "base_ref"
            self._write_selected_alignments(selected, base_ref)
            current_ref = base_ref
            for index, species in enumerate(include_species or [], start=1):
                stage = temporary_dir / f"augment_{index}"
                report = ReferenceBuilder(
                    current_ref,
                    stage,
                    prefer_mafft=self.prefer_mafft,
                    create_archive=False,
                ).build_from_kew(species, client)
                augmentations.append(report.to_dict())
                current_ref = stage / "ref"
            shutil.copytree(current_ref, self.output_dir / "ref")

        gene_ids = sorted(path.stem for path in iter_fasta_files(self.output_dir / "ref"))
        (self.output_dir / "gene.list.txt").write_text(
            "".join(f"{gene}\n" for gene in gene_ids), encoding="utf-8"
        )
        selection = [
            {
                "species_name": candidate.species_name,
                "source": candidate.source,
                "gene_count": candidate.gene_count,
                **asdict(candidate.taxonomy),
            }
            for candidate in selected
        ]
        report = {
            "level": level,
            "parent_taxa": parent_taxa or [],
            "representatives_per_taxon": self.representatives_per_taxon,
            "selected_specimens": len(selected),
            "selected_groups": len({self._group_key(item, level) for item in selected}),
            "loci": len(gene_ids),
            "kew_release": metadata.get("release", "unknown"),
            "kew_augmentations": augmentations,
            "selection": selection,
            "sprout_reference_argument": str(self.output_dir / "ref"),
        }
        (self.output_dir / "panel_report.json").write_text(
            json.dumps(report, indent=2), encoding="utf-8"
        )
        shutil.make_archive(str(self.output_dir) + "_bundle", "zip", self.output_dir)
        return report

    def _load_candidates(
        self, taxonomy_index: dict[str, list]
    ) -> list[PanelCandidate]:
        candidates: list[PanelCandidate] = []
        with self.summary_csv.open(encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle):
                species_name = row["Species"].replace("_", " ")
                indexed = taxonomy_index.get(species_name.casefold(), [])
                if indexed:
                    taxonomy = indexed[0].taxonomy
                else:
                    genus, species = row["Species"].split("_", 1)
                    taxonomy = Taxonomy(row["Order"], row["Family"], genus, species)
                candidates.append(
                    PanelCandidate(
                        species_name=species_name,
                        taxonomy=taxonomy,
                        source=row["Data Source"],
                        gene_count=int(row["Gene Count"]),
                    )
                )
        return candidates

    def _select(
        self, candidates: list[PanelCandidate], level: str, parent_taxa: list[str]
    ) -> list[PanelCandidate]:
        parent = {value.casefold() for value in parent_taxa}
        if level == "family" and parent:
            candidates = [c for c in candidates if c.taxonomy.order.casefold() in parent]
        if level == "genus" and parent:
            candidates = [c for c in candidates if c.taxonomy.family.casefold() in parent]
        grouped: dict[str, list[PanelCandidate]] = {}
        for candidate in candidates:
            grouped.setdefault(self._group_key(candidate, level), []).append(candidate)

        selected: list[PanelCandidate] = []
        for group in sorted(grouped):
            rows = sorted(grouped[group], key=lambda item: item.gene_count, reverse=True)
            child_seen: set[str] = set()
            first_pass: list[PanelCandidate] = []
            for row in rows:
                child = row.taxonomy.family if level == "order" else row.taxonomy.genus
                if child.casefold() not in child_seen:
                    child_seen.add(child.casefold())
                    first_pass.append(row)
            ordered = first_pass + [row for row in rows if row not in first_pass]
            selected.extend(ordered[: self.representatives_per_taxon])
        return selected

    @staticmethod
    def _group_key(candidate: PanelCandidate, level: str) -> str:
        return getattr(candidate.taxonomy, level).casefold()

    def _write_selected_alignments(
        self, selected: list[PanelCandidate], output_dir: Path
    ) -> None:
        wanted = {
            (candidate.group_species.casefold(), candidate.source.casefold()): candidate
            for candidate in selected
        }
        output_dir.mkdir(parents=True, exist_ok=True)
        for path in iter_fasta_files(self.backbone_dir):
            records: list[FastaRecord] = []
            for record in read_fasta(path):
                old_taxonomy = taxonomy_from_header(record.description)
                if old_taxonomy is None:
                    continue
                species_key = f"{old_taxonomy.genus}_{old_taxonomy.species}".casefold()
                tail = record.id.split("_", 4)[4] if len(record.id.split("_", 4)) == 5 else ""
                source = tail.split("-", 1)[0].casefold()
                candidate = wanted.get((species_key, source))
                if candidate is None:
                    continue
                first, separator, remainder = record.description.partition(" ")
                first_tail = first.split("_", 4)[4] if len(first.split("_", 4)) == 5 else tail
                new_first = f"{candidate.taxonomy.sprout_prefix}_{first_tail}"
                description = new_first + (separator + remainder if separator else "")
                records.append(FastaRecord(new_first, description, record.sequence))
            if records:
                write_fasta(output_dir / path.name, records)


def run_mix7_demo(
    backbone_dir: Path,
    summary_csv: Path,
    output_dir: Path,
    client: KewClient,
    prefer_mafft: bool = True,
) -> dict:
    """Build a controlled hierarchical demo from the supplied mix7 truth table."""
    truth = ["Allium sativum", "Asparagus officinalis", "Brassica oleracea", "Artocarpus heterophyllus"]
    orders = ["Asparagales", "Brassicales", "Rosales"]
    families = ["Asparagaceae", "Brassicaceae", "Moraceae"]
    reports = {}
    reports["order"] = HierarchicalPanelBuilder(
        backbone_dir, summary_csv, output_dir / "01_order", 2, prefer_mafft
    ).build("order", client)
    reports["family"] = HierarchicalPanelBuilder(
        backbone_dir, summary_csv, output_dir / "02_family", 2, prefer_mafft
    ).build("family", client, parent_taxa=orders)
    reports["genus"] = HierarchicalPanelBuilder(
        backbone_dir, summary_csv, output_dir / "03_genus", 2, prefer_mafft
    ).build("genus", client, parent_taxa=families, include_species=truth)
    demo = {
        "sample": "mix7 / EM-10",
        "ground_truth": truth,
        "simulated_order_calls": orders,
        "simulated_family_calls": families,
        "truth_usage": (
            "Ground truth is not used to select order/family skeletons. Correct upstream calls are "
            "simulated because no mix7 FASTQ was supplied. Exact truth species are added only to "
            "the controlled genus-panel validation to test Kew integration and output coverage."
        ),
        "reports": reports,
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "demo_report.json").write_text(json.dumps(demo, indent=2), encoding="utf-8")
    return demo
