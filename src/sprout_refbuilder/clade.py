from __future__ import annotations

import csv
import json
import shutil
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from statistics import median
from typing import ClassVar

from .alignment import add_many_to_alignment
from .backbone import Backbone, LocusBackbone
from .fasta import FastaRecord, parse_kew_gene_id, read_fasta, ungap, write_fasta
from .kew import KewClient
from .models import KewRecord, Taxonomy
from .qc import choose_candidate
from .taxonomy import kew_taxonomy_index


@dataclass(frozen=True)
class CladeRecovery:
    record: KewRecord
    taxonomy: Taxonomy
    path: Path
    recovered_loci: int


def _load_locus(path: Path) -> LocusBackbone:
    """Load one locus inside a worker instead of copying the full backbone to every process."""
    records = tuple(read_fasta(path))
    if not records:
        raise ValueError(f"Backbone locus is empty: {path}")
    widths = {len(record.sequence) for record in records}
    if len(widths) != 1:
        raise ValueError(f"Backbone file is not aligned (unequal widths): {path}")
    lengths = [len(ungap(record.sequence)) for record in records]
    return LocusBackbone(path.stem, path, records, widths.pop(), median(lengths))


def _evaluate_locus(
    job: tuple[
        str,
        Path,
        str,
        float,
        list[tuple[CladeRecovery, list[FastaRecord]]],
    ],
) -> tuple[str, list[tuple[FastaRecord, str]], list[dict]]:
    """Evaluate all Kew candidates for one locus in one process.

    Keeping every candidate for a locus together lets the process-local k-mer cache reuse the
    curated backbone features across all recoveries.
    """
    gene_id, locus_path, level, minimum_score, entries = job
    locus = _load_locus(locus_path)
    accepted: list[tuple[FastaRecord, str]] = []
    qc_rows: list[dict] = []
    for recovery, raw in entries:
        chosen, best_qc, evaluated = choose_candidate(
            raw,
            locus,
            minimum_score=minimum_score,
            expected_taxonomy=recovery.taxonomy,
        )
        for sequence, qc in evaluated:
            qc_rows.append(
                {
                    "gene_id": gene_id,
                    "sequence_id": recovery.record.sequence_id,
                    "scientific_name": recovery.record.species,
                    "target_group": getattr(recovery.taxonomy, level),
                    "candidate": sequence.id,
                    **qc.to_dict(),
                }
            )
        if chosen is None or best_qc is None:
            continue
        header = f"{recovery.taxonomy.sprout_prefix}_KEW_{recovery.record.sequence_id}"
        accepted.append((FastaRecord(header, header, chosen.sequence), recovery.record.sequence_id))
    return gene_id, accepted, qc_rows


def _align_locus(
    job: tuple[str, Path, list[FastaRecord], bool, Path],
) -> tuple[str, str, int]:
    """Align and write one final panel locus in an isolated process."""
    gene_id, locus_path, candidates, prefer_mafft, output_path = job
    locus = _load_locus(locus_path)
    aligned, method = add_many_to_alignment(locus, candidates, prefer_mafft)
    wanted_ids = {sequence.id for sequence in candidates}
    panel_records = [sequence for sequence in aligned if sequence.id in wanted_ids]
    if len(panel_records) != len(candidates):
        raise RuntimeError(f"Lost a Kew reference while aligning locus {gene_id}")
    write_fasta(output_path, panel_records)
    return gene_id, method, len(panel_records)


def _progress(stage: str, completed: int, total: int) -> None:
    if completed == 1 or completed == total or completed % 10 == 0:
        timestamp = datetime.now().astimezone().isoformat(timespec="seconds")
        print(f"[{timestamp}] {stage}: {completed}/{total} loci completed", flush=True)


def read_taxon_candidates(path: Path, minimum_z: float | None = None, top: int = 0) -> list[str]:
    """Read SPrOUT's candidate TXT or prediction CSV without losing ranking order."""
    if path.suffix.casefold() == ".csv":
        with path.open(encoding="utf-8-sig", newline="") as handle:
            rows = list(csv.DictReader(handle))
        names = []
        for row in rows:
            name = (row.get("row_name") or row.get("taxon_level") or "").strip()
            if not name:
                continue
            if (
                minimum_z is not None
                and row.get("z_score") not in (None, "")
                and float(row["z_score"]) < minimum_z
            ):
                continue
            names.append(name)
    else:
        names = [line.strip() for line in path.read_text(encoding="utf-8-sig").splitlines()]
        names = [name for name in names if name and not name.startswith("#")]
    deduplicated = list(dict.fromkeys(names))
    return deduplicated[:top] if top > 0 else deduplicated


class BalancedKewPanelBuilder:
    """Build a count-balanced SPrOUT panel directly from Kew assembled recoveries.

    Each target taxon receives the same requested number of recovery records.  This avoids the
    reference-count bias caused by SPrOUT's current sum-over-reference prediction statistic.
    Candidate recoveries are aligned against the curated backbone, but only the selected Kew
    records are written to the final panel.
    """

    PARENT_RANK: ClassVar[dict[str, str]] = {
        "family": "order",
        "genus": "family",
        "species": "genus",
    }

    def __init__(
        self,
        backbone_dir: Path,
        output_dir: Path,
        representatives_per_taxon: int = 4,
        pool_multiplier: int = 2,
        minimum_recovered_loci: int = 50,
        minimum_score: float = 0.50,
        prefer_mafft: bool = True,
        download_workers: int = 4,
        workers: int = 1,
    ) -> None:
        if representatives_per_taxon < 1:
            raise ValueError("representatives_per_taxon must be positive")
        if pool_multiplier < 1:
            raise ValueError("pool_multiplier must be positive")
        self.backbone = Backbone(backbone_dir)
        self.output_dir = output_dir
        self.representatives_per_taxon = representatives_per_taxon
        self.pool_multiplier = pool_multiplier
        self.minimum_recovered_loci = minimum_recovered_loci
        self.minimum_score = minimum_score
        self.prefer_mafft = prefer_mafft
        self.download_workers = max(1, download_workers)
        self.workers = max(1, workers)

    def build(self, level: str, parent_taxa: list[str], client: KewClient) -> dict:
        level = level.casefold()
        if level not in self.PARENT_RANK:
            raise ValueError(f"level must be one of {sorted(self.PARENT_RANK)}")
        parents = {value.strip().casefold() for value in parent_taxa if value.strip()}
        if not parents:
            raise ValueError(f"At least one parent {self.PARENT_RANK[level]} is required")
        if self.output_dir.exists() and any(self.output_dir.iterdir()):
            raise FileExistsError(f"Output directory is not empty: {self.output_dir}")

        metadata = client.sync()
        taxonomy_index = kew_taxonomy_index(client.fetch_species_tree())
        candidates = self._records_in_parents(level, parents, client.records(), taxonomy_index)
        pools = self._balanced_pools(level, candidates)
        recoveries = self._download_pool(pools, client)
        if not recoveries:
            raise LookupError(
                f"Kew has no usable recoveries for {self.PARENT_RANK[level]} candidates: "
                f"{', '.join(parent_taxa)}"
            )

        self.output_dir.mkdir(parents=True, exist_ok=True)
        recovery_sequences: dict[str, dict[str, list[FastaRecord]]] = {}
        for recovery in recoveries:
            per_gene: dict[str, list[FastaRecord]] = defaultdict(list)
            for sequence in read_fasta(recovery.path):
                per_gene[parse_kew_gene_id(sequence)].append(sequence)
            recovery_sequences[recovery.record.sequence_id] = per_gene

        recovery_by_id = {item.record.sequence_id: item for item in recoveries}
        qc_rows: list[dict] = []
        accepted_per_record: dict[str, int] = defaultdict(int)
        accepted_by_gene: dict[str, list[tuple[FastaRecord, CladeRecovery]]] = defaultdict(list)
        evaluation_jobs = []
        for gene_id, locus in sorted(self.backbone.loci.items()):
            entries = [
                (
                    recovery,
                    recovery_sequences[recovery.record.sequence_id].get(gene_id, []),
                )
                for recovery in recoveries
            ]
            evaluation_jobs.append((gene_id, locus.path, level, self.minimum_score, entries))

        print(
            f"QC stage: {len(evaluation_jobs)} loci using {self.workers} worker(s)", flush=True
        )
        evaluation_results = self._run_locus_jobs(_evaluate_locus, evaluation_jobs, "QC")
        for gene_id, accepted, locus_qc_rows in evaluation_results:
            qc_rows.extend(locus_qc_rows)
            for sequence, sequence_id in accepted:
                accepted_by_gene[gene_id].append((sequence, recovery_by_id[sequence_id]))
                accepted_per_record[sequence_id] += 1

        final_recoveries = self._select_after_qc(level, recoveries, accepted_per_record)
        if not final_recoveries:
            raise RuntimeError("No Kew recovery passed the specimen-level accepted-locus gate")
        final_ids = {item.record.sequence_id for item in final_recoveries}
        reference_dir = self.output_dir / "ref"
        reference_dir.mkdir()
        alignment_methods: set[str] = set()
        written_genes: list[str] = []
        alignment_jobs = []
        for gene_id, locus in sorted(self.backbone.loci.items()):
            accepted = [
                item
                for item in accepted_by_gene.get(gene_id, [])
                if item[1].record.sequence_id in final_ids
            ]
            groups = {getattr(recovery.taxonomy, level).casefold() for _, recovery in accepted}
            if len(accepted) < 2 or len(groups) < 2:
                continue
            alignment_jobs.append(
                (
                    gene_id,
                    locus.path,
                    [sequence for sequence, _ in accepted],
                    self.prefer_mafft,
                    reference_dir / f"{gene_id}.fasta",
                )
            )

        print(
            f"Alignment stage: {len(alignment_jobs)} loci using {self.workers} worker(s)",
            flush=True,
        )
        alignment_results = self._run_locus_jobs(_align_locus, alignment_jobs, "alignment")
        for gene_id, method, _ in alignment_results:
            alignment_methods.add(method)
            written_genes.append(gene_id)

        if not written_genes:
            raise RuntimeError("No locus retained at least two target groups after automatic QC")
        written_genes.sort(
            key=lambda value: (not value.isdigit(), int(value) if value.isdigit() else value)
        )
        qc_rows.sort(key=lambda row: (row["gene_id"], row["sequence_id"], row["candidate"]))
        (self.output_dir / "gene.list.txt").write_text(
            "".join(f"{gene}\n" for gene in written_genes), encoding="utf-8"
        )
        (self.output_dir / "qc.json").write_text(
            json.dumps(qc_rows, indent=2), encoding="utf-8"
        )

        selection = []
        for recovery in final_recoveries:
            selection.append(
                {
                    **recovery.record.to_dict(),
                    "scientific_name": recovery.record.species,
                    "order": recovery.taxonomy.order,
                    "family": recovery.taxonomy.family,
                    "genus": recovery.taxonomy.genus,
                    "species_epithet": recovery.taxonomy.species,
                    "recovered_loci": recovery.recovered_loci,
                    "accepted_loci": accepted_per_record[recovery.record.sequence_id],
                }
            )
        report = {
            "schema_version": 1,
            "level": level,
            "parent_rank": self.PARENT_RANK[level],
            "parent_taxa": parent_taxa,
            "selection_policy": "equal recovery count per target taxon; child-diversity first",
            "representatives_per_taxon": self.representatives_per_taxon,
            "minimum_recovered_loci": self.minimum_recovered_loci,
            "workers": self.workers,
            "candidate_pool_recoveries": len(recoveries),
            "selected_recoveries": len(final_recoveries),
            "selected_groups": len(
                {getattr(item.taxonomy, level) for item in final_recoveries}
            ),
            "omitted_groups": sorted(
                set(pools)
                - {getattr(item.taxonomy, level) for item in final_recoveries}
            ),
            "loci": len(written_genes),
            "kew_release": metadata.get("release", "unknown"),
            "alignment_method": ",".join(sorted(alignment_methods)),
            "selection": selection,
            "sprout_reference_argument": str(reference_dir),
            "gene_list": str(self.output_dir / "gene.list.txt"),
        }
        (self.output_dir / "panel_report.json").write_text(
            json.dumps(report, indent=2), encoding="utf-8"
        )
        self._write_selection_tsv(selection)
        shutil.make_archive(str(self.output_dir) + "_bundle", "zip", self.output_dir)
        return report

    def _run_locus_jobs(self, function, jobs: list[tuple], stage: str) -> list[tuple]:
        """Run independent locus jobs, preserving a low-overhead serial path for one worker."""
        if not jobs:
            return []
        results = []
        if self.workers == 1:
            for completed, job in enumerate(jobs, start=1):
                results.append(function(job))
                _progress(stage, completed, len(jobs))
            return results
        with ProcessPoolExecutor(max_workers=self.workers) as executor:
            futures = [executor.submit(function, job) for job in jobs]
            for completed, future in enumerate(as_completed(futures), start=1):
                results.append(future.result())
                _progress(stage, completed, len(futures))
        return results

    def _records_in_parents(
        self,
        level: str,
        parents: set[str],
        records: list[KewRecord],
        taxonomy_index: dict[str, list],
    ) -> list[tuple[KewRecord, Taxonomy]]:
        parent_rank = self.PARENT_RANK[level]
        result: list[tuple[KewRecord, Taxonomy]] = []
        for record in records:
            if not record.recovery_file:
                continue
            binomial = " ".join(record.species.replace("_", " ").split()[:2]).casefold()
            indexed = taxonomy_index.get(binomial, [])
            if not indexed:
                continue
            exact = next(
                (
                    item
                    for item in indexed
                    if item.sequence_id.casefold() == record.sequence_id.casefold()
                ),
                indexed[0],
            )
            if getattr(exact.taxonomy, parent_rank).casefold() in parents:
                result.append((record, exact.taxonomy))
        return result

    def _balanced_pools(
        self, level: str, candidates: list[tuple[KewRecord, Taxonomy]]
    ) -> dict[str, list[tuple[KewRecord, Taxonomy]]]:
        grouped: dict[str, list[tuple[KewRecord, Taxonomy]]] = defaultdict(list)
        for item in candidates:
            grouped[getattr(item[1], level)].append(item)
        pool_size = self.representatives_per_taxon * self.pool_multiplier
        pools: dict[str, list[tuple[KewRecord, Taxonomy]]] = {}
        for group, rows in sorted(grouped.items()):
            rows.sort(key=self._record_priority)
            child_seen: set[str] = set()
            diverse: list[tuple[KewRecord, Taxonomy]] = []
            for row in rows:
                child = row[1].genus if level == "family" else row[0].species
                if child.casefold() not in child_seen:
                    child_seen.add(child.casefold())
                    diverse.append(row)
            ordered = diverse + [row for row in rows if row not in diverse]
            pools[group] = ordered[:pool_size]
        return pools

    def _download_pool(
        self,
        pools: dict[str, list[tuple[KewRecord, Taxonomy]]],
        client: KewClient,
    ) -> list[CladeRecovery]:
        flat = [(group, record, taxonomy) for group, rows in pools.items() for record, taxonomy in rows]

        def fetch(item: tuple[str, KewRecord, Taxonomy]) -> tuple[str, CladeRecovery]:
            group, record, taxonomy = item
            path = client.fetch_recovery(record)
            gene_count = len({parse_kew_gene_id(seq) for seq in read_fasta(path)})
            return group, CladeRecovery(record, taxonomy, path, gene_count)

        downloaded: dict[str, list[CladeRecovery]] = defaultdict(list)
        with ThreadPoolExecutor(max_workers=self.download_workers) as executor:
            for group, recovery in executor.map(fetch, flat):
                if recovery.recovered_loci >= self.minimum_recovered_loci:
                    downloaded[group].append(recovery)

        return [row for group in sorted(downloaded) for row in downloaded[group]]

    def _select_after_qc(
        self,
        level: str,
        recoveries: list[CladeRecovery],
        accepted_per_record: dict[str, int],
    ) -> list[CladeRecovery]:
        grouped: dict[str, list[CladeRecovery]] = defaultdict(list)
        for recovery in recoveries:
            if accepted_per_record[recovery.record.sequence_id] >= self.minimum_recovered_loci:
                grouped[getattr(recovery.taxonomy, level)].append(recovery)
        selected: list[CladeRecovery] = []
        for group in sorted(grouped):
            rows = sorted(
                grouped[group],
                key=lambda item: (
                    -accepted_per_record[item.record.sequence_id],
                    -item.recovered_loci,
                    self._record_priority((item.record, item.taxonomy)),
                ),
            )
            child_seen: set[str] = set()
            diverse: list[CladeRecovery] = []
            for row in rows:
                child = row.taxonomy.genus if level == "family" else row.record.species
                if child.casefold() not in child_seen:
                    child_seen.add(child.casefold())
                    diverse.append(row)
            ordered = diverse + [row for row in rows if row not in diverse]
            selected.extend(ordered[: self.representatives_per_taxon])
        return selected

    @staticmethod
    def _record_priority(item: tuple[KewRecord, Taxonomy]) -> tuple:
        record, _ = item
        sequence_type = record.sequence_type.casefold()
        type_rank = 0 if sequence_type == "annotated_genome" else 1 if "genome" in sequence_type else 2
        project_rank = 0 if "paftol" in record.project.casefold() else 1
        return type_rank, project_rank, record.species.casefold(), record.sequence_id.casefold()

    def _write_selection_tsv(self, rows: list[dict]) -> None:
        path = self.output_dir / "selection.tsv"
        if not rows:
            path.write_text("", encoding="utf-8")
            return
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
