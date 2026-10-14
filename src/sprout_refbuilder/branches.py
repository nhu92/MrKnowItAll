"""Progressive, tree-guided screening for large genus panels."""

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
    max_group_genera: int = 40,
    max_anchor_genera: int = 120,
    max_branches_per_family: int = 6,
    allowed_genera_file: Path | None = None,
) -> dict:
    """Cut each large family into a few coarse tree groups for one screen round.

    A subsequent round can use ``allowed_genera_file`` to subdivide only the
    previously retained groups. This reads metadata and the tree, not FASTAs.
    """
    if min(threshold, max_group_genera, max_anchor_genera, max_branches_per_family) < 1:
        raise ValueError("Tree-screen limits must be positive")
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
    allowed = None
    if allowed_genera_file is not None:
        allowed = set(read_taxon_candidates(allowed_genera_file))
        if not allowed:
            raise ValueError("The allowed genus file is empty")
        candidates = [
            (record, taxonomy) for record, taxonomy in candidates
            if taxonomy.genus in allowed
        ]
    pools = builder._balanced_pools("genus", candidates)
    pool_records = sum(len(rows) for rows in pools.values())
    all_genera = sorted(pools)
    if allowed is not None and allowed - set(all_genera):
        raise ValueError("Allowed genera are missing from the selected Kew families")
    report: dict = {
        "schema_version": 2,
        "kew_release": client.sync().get("release", "unknown"),
        "parent_families": parents,
        "trigger_threshold": threshold,
        "pool_records": pool_records,
        "candidate_genera": len(all_genera),
        "all_genera": all_genera,
        "max_group_genera": max_group_genera,
        "max_branches_per_family": max_branches_per_family,
        "triggered": False,
        "branches": [],
        "anchor_genera": [],
        "passthrough_genera": [],
    }
    if len(all_genera) <= threshold or len(all_genera) < 2:
        report["reason"] = "candidate genus count is below threshold"
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
    family_tips: dict[str, dict[str, object]] = defaultdict(dict)
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
            family_tips[family].setdefault(genus, tip)
            tip_record[id(tip)] = record

    large_families = [
        family for family in parents
        if len(family_tips.get(family, {})) > max_group_genera
    ]
    if large_families and max_anchor_genera < 2 * len(large_families):
        raise ValueError("Anchor budget cannot cover two anchors per large family")
    family_group_limit = min(
        max_branches_per_family,
        max(1, max_anchor_genera // (2 * max(1, len(large_families)))),
    )
    branch_nodes = []
    passthrough = set()
    for family in parents:
        tips = list(family_tips.get(family, {}).values())
        if not tips:
            continue
        if len(tips) <= max_group_genera:
            passthrough.update(family_tips[family])
            continue
        family_root = tree.common_ancestor(tips)
        descendant_tips: dict[int, list] = {}
        selected_ids = {id(tip) for tip in tips}

        def collect(node) -> list:
            if node.is_terminal():
                descendants = [node] if id(node) in selected_ids else []
            else:
                descendants = []
                for child in node.clades:
                    descendants.extend(collect(child))
            descendant_tips[id(node)] = descendants
            return descendants

        collect(family_root)

        def group_tips(nodes: list) -> list:
            return [tip for node in nodes for tip in descendant_tips[id(node)]]

        def split(nodes: list) -> tuple[list, list] | None:
            parts = nodes
            while len(parts) == 1:
                children = [
                    child for child in parts[0].clades
                    if descendant_tips.get(id(child))
                ]
                if len(children) == 1:
                    parts = children
                elif len(children) > 1:
                    parts = children
                    break
                else:
                    return None
            if len(parts) < 2:
                return None
            bins: tuple[list, list] = ([], [])
            weights = [0, 0]
            for part in sorted(
                parts,
                key=lambda node: (-len(descendant_tips[id(node)]), str(node.name or "")),
            ):
                side = 0 if weights[0] <= weights[1] else 1
                bins[side].append(part)
                weights[side] += len(descendant_tips[id(part)])
            return bins

        groups = [[family_root]]
        while len(groups) < family_group_limit:
            choices = [
                (len(group_tips(group)), index, halves)
                for index, group in enumerate(groups)
                if len(group_tips(group)) > max_group_genera
                if (halves := split(group)) is not None
            ]
            if not choices:
                break
            _, index, halves = max(choices, key=lambda item: (item[0], -item[1]))
            groups[index:index + 1] = list(halves)
        branch_nodes.extend(
            (family, group_tips(group)) for group in groups
        )

    if not branch_nodes:
        report["passthrough_genera"] = sorted(passthrough)
        report["unplaced_genera"] = sorted(set(all_genera) - passthrough)
        report["reason"] = "no large family has screenable tree tips"
        _write_json(output, report)
        return report

    branch_nodes.sort(
        key=lambda row: (row[0].casefold(), min(tip_record[id(tip)].species for tip in row[1]))
    )
    anchors_by_branch: dict[str, list] = {}
    for index, (family, tips) in enumerate(branch_nodes, start=1):
        branch_id = f"B{index:03d}"
        by_genus: dict[str, object] = {}
        for tip in tips:
            genus = tip_record[id(tip)].species.split()[0]
            by_genus.setdefault(genus, tip)
        choices = list(by_genus.values())
        first = min(choices, key=lambda tip: tip_record[id(tip)].species)
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
    report["passthrough_genera"] = sorted(passthrough)
    report["unplaced_genera"] = sorted(set(all_genera) - covered_genera - passthrough)
    report["family_group_limit"] = family_group_limit
    report["triggered"] = True
    report["reason"] = "candidate genus count exceeds threshold"
    _write_json(output, report)
    return report


def expand_genus_branches(
    plan_path: Path,
    evidence_csv: Path,
    output: Path,
    report_path: Path,
    *,
    max_selected_branches_per_family: int = 3,
) -> dict:
    """Keep the strongest supported branches per family, never all genera by default.

    Anchor panels are coarse screens, not final genus calls. A weak family can
    retain two exploratory branches, but the report marks that evidence weak.
    """
    if max_selected_branches_per_family < 1:
        raise ValueError("max_selected_branches_per_family must be positive")
    plan = json.loads(plan_path.read_text(encoding="utf-8"))
    branches = plan["branches"]
    all_genera = set(plan["all_genera"])
    unplaced_genera = set(plan.get("unplaced_genera", []))
    passthrough_genera = set(plan.get("passthrough_genera", []))
    evidence_by_genus: dict[str, dict] = {}
    with evidence_csv.open(encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            if row.get("taxon"):
                evidence_by_genus[row["taxon"]] = row

    def metric(row: dict, key: str) -> float:
        try:
            return float(row.get(key) or 0)
        except (TypeError, ValueError):
            return 0.0

    by_family: dict[str, list] = defaultdict(list)
    for branch in branches:
        by_family[branch["family"]].append(branch)
    selected_branches: set[str] = set()
    weak_families = []
    family_details = []
    for family, family_branches in by_family.items():
        ranked = []
        for branch in family_branches:
            rows = [evidence_by_genus.get(name, {}) for name in branch["anchor_genera"]]
            votes = sum(metric(row, "locus_votes") for row in rows)
            score = max((metric(row, "aggregate_score") for row in rows), default=0.0)
            selected = any(row.get("selected", "").casefold() == "true" for row in rows)
            ranked.append((votes, score, selected, branch))
        ranked.sort(key=lambda item: (-item[0], -item[1], -int(item[2]), item[3]["id"]))
        top_votes = ranked[0][0]
        if top_votes > 0:
            vote_floor = max(2, top_votes * 0.1)
            supported = [row for row in ranked if row[0] >= vote_floor]
            picked = supported[:max_selected_branches_per_family]
            if not picked:
                picked = ranked[:1]
        else:
            weak_families.append(family)
            marked = [row for row in ranked if row[2]]
            picked = (
                marked[:max_selected_branches_per_family]
                if marked else ranked[:min(2, max_selected_branches_per_family)]
            )
        selected_branches.update(row[3]["id"] for row in picked)
        family_details.append({
            "family": family,
            "top_anchor_locus_votes": top_votes,
            "selected_branches": [row[3]["id"] for row in picked],
            "weak_evidence": top_votes == 0,
        })

    retained = selected_branches | {
        branch["neighbor"] for branch in branches
        if branch["id"] in selected_branches and branch["neighbor"]
    }
    chosen = passthrough_genera | unplaced_genera | {
        genus for branch in branches if branch["id"] in retained
        for genus in branch["genera"]
    }
    if not chosen or chosen - all_genera:
        raise RuntimeError("Branch expansion produced an invalid genus subset")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("".join(f"{genus}\n" for genus in sorted(chosen)), encoding="utf-8")
    report = {
        "schema_version": 2,
        "mode": "branch_screened" if len(chosen) < len(all_genera) else "unresolved_no_reduction",
        "selected_anchor_branches": sorted(selected_branches),
        "retained_branches": sorted(retained),
        "weak_families": sorted(weak_families),
        "families": family_details,
        "selected_genera": len(chosen),
        "total_genera": len(all_genera),
        "allowed_taxa_file": str(output),
    }
    _write_json(report_path, report)
    return report


def _write_json(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2), encoding="utf-8")
