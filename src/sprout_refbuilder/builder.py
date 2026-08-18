from __future__ import annotations

import json
import shutil
from collections import defaultdict
from pathlib import Path

from .alignment import add_to_alignment
from .backbone import Backbone
from .fasta import FastaRecord, format_fasta, parse_kew_gene_id, read_fasta, write_fasta
from .kew import KewClient
from .models import BuildReport, KewRecord, Taxonomy
from .qc import choose_candidate
from .taxonomy import taxonomy_from_kew_tree


class ReferenceBuilder:
    def __init__(
        self,
        backbone_dir: Path,
        output_dir: Path,
        minimum_score: float = 0.50,
        prefer_mafft: bool = True,
        ml_model: Path | None = None,
        ml_min_score: float = -0.05,
        create_archive: bool = True,
    ) -> None:
        self.backbone = Backbone(backbone_dir)
        self.output_dir = output_dir
        self.minimum_score = minimum_score
        self.prefer_mafft = prefer_mafft
        self.ml_model = ml_model
        self.ml_min_score = ml_min_score
        self.create_archive = create_archive

    def build_from_kew(self, taxon: str, client: KewClient) -> BuildReport:
        metadata = client.sync()
        matches = client.resolve_exact(taxon)
        if not matches:
            suggestions = ", ".join(r.species for r in client.search(taxon, limit=5))
            raise LookupError(
                f"No exact Kew recovery for {taxon!r}. Closest indexed names: {suggestions}. "
                "A species-labelled reference cannot be synthesized from relatives alone."
            )
        recovery_paths = [(record, client.fetch_recovery(record)) for record in matches]
        tree_path = client.fetch_species_tree()
        taxonomy = taxonomy_from_kew_tree(tree_path, matches[0].species, matches[0].sequence_id)
        report = self._build(
            requested_taxon=taxon,
            resolved_taxon=matches[0].species,
            source="kew-paftol-assembled-recovery",
            taxonomy=taxonomy,
            recoveries=recovery_paths,
            release=str(metadata.get("release", "unknown")),
        )
        report.provenance["kew_base_url"] = client.base_url
        report.provenance["sequence_ids"] = ",".join(record.sequence_id for record in matches)
        self._write_report(report)
        return report

    def build_from_fasta(
        self,
        taxon: str,
        taxonomy: Taxonomy,
        recovery_fasta: Path,
        sequence_id: str = "local",
    ) -> BuildReport:
        record = KewRecord("local", sequence_id, "assembled", taxon, "local", recovery_fasta.name)
        report = self._build(
            requested_taxon=taxon,
            resolved_taxon=taxon,
            source="local-assembled-recovery",
            taxonomy=taxonomy,
            recoveries=[(record, recovery_fasta)],
            release="local",
        )
        self._write_report(report)
        return report

    def _build(
        self,
        requested_taxon: str,
        resolved_taxon: str,
        source: str,
        taxonomy: Taxonomy,
        recoveries: list[tuple[KewRecord, Path]],
        release: str,
    ) -> BuildReport:
        if self.output_dir.exists() and any(self.output_dir.iterdir()):
            raise FileExistsError(f"Output directory is not empty: {self.output_dir}")
        reference_dir = self.output_dir / "ref"
        reference_dir.mkdir(parents=True, exist_ok=True)
        candidates: dict[str, list[tuple[FastaRecord, KewRecord]]] = defaultdict(list)
        for source_record, path in recoveries:
            for record in read_fasta(path):
                candidates[parse_kew_gene_id(record)].append((record, source_record))

        report = BuildReport(requested_taxon, resolved_taxon, source, self.output_dir, release)
        combined: list[FastaRecord] = []
        qc_rows: list[dict] = []
        alignment_methods: set[str] = set()
        for gene_id, locus in sorted(self.backbone.loci.items()):
            gene_candidates = candidates.get(gene_id, [])
            if not gene_candidates:
                shutil.copyfile(locus.path, reference_dir / f"{gene_id}.fasta")
                report.missing_loci.append(gene_id)
                continue
            candidate_only = [item[0] for item in gene_candidates]
            selected, best_qc, evaluated = choose_candidate(
                candidate_only,
                locus,
                minimum_score=self.minimum_score,
                expected_taxonomy=taxonomy,
            )
            selected_qc_row: dict | None = None
            for candidate, qc in evaluated:
                row = {"gene_id": gene_id, "candidate": candidate.id, **qc.to_dict()}
                qc_rows.append(row)
                if candidate is selected:
                    selected_qc_row = row
            if selected is None or best_qc is None:
                shutil.copyfile(locus.path, reference_dir / f"{gene_id}.fasta")
                report.rejected_loci[gene_id] = list(best_qc.reasons if best_qc else ("no candidate",))
                continue
            if self.ml_model:
                from .ml import anomaly_score

                ml_score = anomaly_score(self.ml_model, selected, locus.median_ungapped_length)
                if selected_qc_row is not None:
                    selected_qc_row["ml_anomaly_score"] = ml_score
                if ml_score < self.ml_min_score:
                    shutil.copyfile(locus.path, reference_dir / f"{gene_id}.fasta")
                    report.rejected_loci[gene_id] = [
                        f"ML anomaly {ml_score:.4f} is below {self.ml_min_score:.4f}"
                    ]
                    continue
            selected_index = candidate_only.index(selected)
            source_record = gene_candidates[selected_index][1]
            header = f"{taxonomy.sprout_prefix}_{source_record.sequence_id}"
            normalized = FastaRecord(header, header, selected.sequence)
            aligned, method = add_to_alignment(locus, normalized, prefer_mafft=self.prefer_mafft)
            alignment_methods.add(method)
            write_fasta(reference_dir / f"{gene_id}.fasta", aligned)
            combined.append(
                FastaRecord(
                    f"{header}|gene={gene_id}",
                    f"{header}|gene={gene_id}|source={source_record.repository}:{source_record.sequence_id}",
                    selected.sequence.replace("-", ""),
                )
            )
            report.accepted_loci.append(gene_id)

        (self.output_dir / "gene.list.txt").write_text(
            "".join(f"{gene}\n" for gene in sorted(self.backbone.loci)), encoding="utf-8"
        )
        (self.output_dir / "requested_species_sequences.fasta").write_text(
            format_fasta(combined), encoding="utf-8"
        )
        (self.output_dir / "qc.json").write_text(json.dumps(qc_rows, indent=2), encoding="utf-8")
        report.provenance["alignment_method"] = ",".join(sorted(alignment_methods)) or "none"
        report.provenance["backbone"] = str(self.backbone.directory)
        report.provenance["ml_model"] = str(self.ml_model) if self.ml_model else "disabled"
        report.provenance["ml_min_score"] = str(self.ml_min_score)
        report.provenance["sprout_reference_argument"] = str(reference_dir)
        return report

    def _write_report(self, report: BuildReport) -> None:
        (self.output_dir / "report.json").write_text(
            json.dumps(report.to_dict(), indent=2), encoding="utf-8"
        )
        # Write the report first and keep the archive outside its source directory.
        if self.create_archive:
            shutil.make_archive(str(self.output_dir) + "_bundle", "zip", self.output_dir)
