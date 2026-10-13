"""Conservative, tree-guided screening for unusually large genus panels."""

from __future__ import annotations

import csv
import json
from collections import defaultdict
from pathlib import Path

from Bio import Phylo

from .clade import BalancedKewPanelBuilder, read_taxon_candidates
from .kew import KewClient
from .taxonomy import kew_taxonomy_index


def _tip_fields(name: str | None) -> tuple[str, str, str, str] | None:
    fields = (name or "").split("_", 4)
    if len(fields) != 5:
        return None
    return fields[0], fields[1], fields[2], fields[4]


def plan_genus_branches(
    client: KewClient,
    backbone: Path,
    parent_file: Path,
    output: Path,
    *,
    threshold: int = 120,
    max_group_genera: int = 20,
    max_anchor_genera: int = 120,
) -> dict:
    """Count the exact genus pool and choose divergent anchors per tree branch.

    This reads only Kew metadata and its species tree. It does not download recovery
    FASTAs, so a small family can proceed directly to the ordinary genus panel.
    """
    parents = read_taxon_candidates(parent_file)
    if not parents:
        raise ValueError("The family candidate file is empty")
    builder = BalancedKewPanelBuilder(
        backbone, output.parent / ".unused", representatives_per_taxon=2, pool_multiplier=2
    )
    tree_path = client.fetch_species_tree()
    indexed = kew_taxonomy_index(tree_path)
    candidates = builder._records_in_parents(
        "genus", {name.casefold() for name in parents}, client.records(), indexed
    )
    pools = builder._balanced_pools("genus", candidates)
    pool_records = sum(len(rows) for rows in pools.values())
    all_genera = sorted(pools)
    report: dict = {
        "schema_version": 1,
        "kew_release": client.sync().get("release", "unknown"),
        "parent_families": parents,
        "trigger_threshold": threshold,
        "pool_records": pool_records,
        "candidate_genera": len(all_genera),
        "all_genera": all_genera,
        "triggered": False,
        "branches": [],
        "anchor_genera": [],
    }
    if pool_records <= threshold or len(all_genera) < 2:
        report["reason"] = "candidate pool is below threshold"
        _write_json(output, report)
        return report

    tree = Phylo.read(tree_path, "newick")
    parent_node: dict[int, object] = {}
    root_distance = {id(tree.root): 0.0}

    stack = [tree.root]
    while stack:
        node = stack.pop()
        for child in node.clades:
            parent_node[id(child)] = node
            root_distance[id(child)] = root_distance[id(node)] + (
                child.branch_length if child.branch_length is not None else 1.0
            )
            stack.append(child)

    def tip_distance(left, right) -> float:
        ancestors = set()
        node = left
        while node is not None:
            ancestors.add(id(node))
            node = parent_node.get(id(node))
        node = right
        while id(node) not in ancestors:
            node = parent_node[id(node)]
        return root_distance[id(left)] + root_distance[id(right)] - 2 * root_distance[id(node)]

    candidate_index = {
        (
            taxonomy.family.casefold(),
            taxonomy.genus.casefold(),
            record.sequence_id.casefold(),
        ): record
        for record, taxonomy in candidates
    }
    family_tips: dict[str, list] = defaultdict(list)
    tip_record: dict[int, object] = {}
    for tip in tree.get_terminals():
        fields = _tip_fields(tip.name)
        if fields is None:
            continue
        _, family, genus, sequence_id = fields
        record = candidate_index.get(
            (family.casefold(), genus.casefold(), sequence_id.casefold())
        )
        if record is not None:
            family_tips[family].append(tip)
            tip_record[id(tip)] = record

    branch_nodes = []
    for family in parents:
        tips = family_tips.get(family, [])
        if not tips:
            continue
        family_root = tree.common_ancestor(tips)
        descendant_genera: dict[int, set[str]] = {}
        descendant_tips: dict[int, list] = {}

        def collect(node) -> set[str]:
            if node.is_terminal():
                record = tip_record.get(id(node))
                names = {record.species.split()[0]} if record is not None else set()
                descendants = [node] if record is not None else []
            else:
                names = set()
                descendants = []
                for child in node.clades:
                    names.update(collect(child))
                    descendants.extend(descendant_tips[id(child)])
            descendant_genera[id(node)] = names
            descendant_tips[id(node)] = descendants
            return names

        collect(family_root)

        def cut(node, limit: int) -> list:
            names = descendant_genera[id(node)]
            if not names:
                return []
            if len(names) <= limit or node.is_terminal():
                return [node]
            return [part for child in node.clades for part in cut(child, limit)]

        limit = max_group_genera
        groups = cut(family_root, limit)
        while len(groups) * 2 > max_anchor_genera and limit < len(all_genera):
            limit += 5
            groups = cut(family_root, limit)
        branch_nodes.extend((family, group, descendant_tips[id(group)]) for group in groups)

    if not branch_nodes or len(branch_nodes) > max_anchor_genera:
        report["reason"] = "tree branches could not be screened within anchor budget"
        _write_json(output, report)
        return report

    branch_nodes.sort(
        key=lambda row: (row[0].casefold(), sorted(tip_record[id(tip)].species for tip in row[2])[0])
    )
    anchors_by_branch: dict[str, list] = {}
    for index, (family, _, tips) in enumerate(branch_nodes, start=1):
        branch_id = f"B{index:03d}"
        by_genus: dict[str, object] = {}
        for tip in tips:
            genus = tip_record[id(tip)].species.split()[0]
            by_genus.setdefault(genus, tip)
        choices = list(by_genus.values())
        first = sorted(choices, key=lambda tip: tip_record[id(tip)].species)[0]
        anchors = [first]
        if len(choices) > 1:
            second = max(
                (tip for tip in choices if tip is not first),
                key=lambda tip: tip_distance(first, tip),
            )
            anchors.append(second)
        anchors_by_branch[branch_id] = anchors
        report["branches"].append(
            {
                "id": branch_id,
                "family": family,
                "genera": sorted(by_genus),
                "anchor_genera": [tip_record[id(tip)].species.split()[0] for tip in anchors],
                "neighbor": None,
            }
        )

    for branch in report["branches"]:
        peers = [
            other for other in report["branches"]
            if other["id"] != branch["id"] and other["family"] == branch["family"]
        ]
        if peers:
            branch["neighbor"] = min(
                peers,
                key=lambda other: min(
                    tip_distance(left, right)
                    for left in anchors_by_branch[branch["id"]]
                    for right in anchors_by_branch[other["id"]]
                ),
            )["id"]
    report["anchor_genera"] = sorted(
        {genus for branch in report["branches"] for genus in branch["anchor_genera"]}
    )
    covered_genera = {genus for branch in report["branches"] for genus in branch["genera"]}
    report["unplaced_genera"] = sorted(set(all_genera) - covered_genera)
    report["triggered"] = True
    report["reason"] = "candidate pool exceeds threshold"
    _write_json(output, report)
    return report


def expand_genus_branches(
    plan_path: Path,
    evidence_csv: Path,
    output: Path,
    report_path: Path,
    *,
    max_selected_branches: int = 8,
) -> dict:
    """Keep supported anchor branches and neighbors, or fail open to the full family."""
    plan = json.loads(plan_path.read_text(encoding="utf-8"))
    branches = plan["branches"]
    all_genera = set(plan["all_genera"])
    unplaced_genera = set(plan.get("unplaced_genera", []))
    selected_branches = set()
    with evidence_csv.open(encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            if row.get("selected", "").casefold() != "true":
                continue
            genus = row.get("taxon", "")
            selected_branches.update(
                branch["id"] for branch in branches if genus in branch["anchor_genera"]
            )
    by_family: dict[str, list] = defaultdict(list)
    for branch in branches:
        by_family[branch["family"]].append(branch)
    missing = [
        family for family, rows in by_family.items()
        if not selected_branches.intersection({row["id"] for row in rows})
    ]
    if not selected_branches or len(selected_branches) > max_selected_branches or missing:
        chosen = all_genera
        mode = "full_fallback"
    else:
        retained = selected_branches | {
            branch["neighbor"] for branch in branches
            if branch["id"] in selected_branches and branch["neighbor"]
        }
        chosen = unplaced_genera | {
            genus for branch in branches if branch["id"] in retained
            for genus in branch["genera"]
        }
        if len(chosen) > 0.8 * len(all_genera):
            chosen = all_genera
            mode = "full_fallback"
        else:
            mode = "branch_screened"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("".join(f"{genus}\n" for genus in sorted(chosen)), encoding="utf-8")
    report = {
        "mode": mode,
        "selected_anchor_branches": sorted(selected_branches),
        "selected_genera": len(chosen),
        "total_genera": len(all_genera),
        "allowed_taxa_file": str(output),
    }
    _write_json(report_path, report)
    return report


def _write_json(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2), encoding="utf-8")
