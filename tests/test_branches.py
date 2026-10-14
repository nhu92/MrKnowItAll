import json
from pathlib import Path

from sprout_refbuilder.branches import expand_genus_branches, plan_genus_branches
from sprout_refbuilder.clade import BalancedKewPanelBuilder
from sprout_refbuilder.models import KewRecord


class TinyKew:
    def __init__(self, root: Path) -> None:
        root.mkdir()
        self.root = root
        self.tree = root / "tree.nwk"
        leaves = [f"Rosales_FamA_Genus{letter}_species_{letter}1" for letter in "ABCDEF"]
        self.tree.write_text(
            f"((({leaves[0]},{leaves[1]}),({leaves[2]},{leaves[3]})),"
            f"({leaves[4]},{leaves[5]}));\n",
            encoding="utf-8",
        )
        self._records = [
            KewRecord("INSDC", f"{letter}1", "read", f"Genus{letter} species", "PAFTOL", f"{letter}.fasta")
            for letter in "ABCDEF"
        ]

    def fetch_species_tree(self) -> Path:
        return self.tree

    def records(self) -> list[KewRecord]:
        return self._records

    def sync(self) -> dict:
        return {"release": "test"}

    def fetch_recovery(self, record: KewRecord) -> Path:
        assert record.recovery_file
        return self.root / record.recovery_file


def test_branch_threshold_and_expansion(tmp_path: Path) -> None:
    backbone = tmp_path / "backbone"
    backbone.mkdir()
    (backbone / "4691.fasta").write_text(
        ">Rosales_FamA_GenusA_species_X\nATGACG\n", encoding="utf-8"
    )
    parents = tmp_path / "parents.txt"
    parents.write_text("FamA\n", encoding="utf-8")
    client = TinyKew(tmp_path / "kew")
    plan_path = tmp_path / "plan.json"

    small = plan_genus_branches(client, backbone, parents, plan_path, threshold=6)
    assert small["pool_records"] == 6
    assert not small["triggered"]

    plan = plan_genus_branches(
        client, backbone, parents, plan_path,
        threshold=5, max_group_genera=2, max_anchor_genera=6,
    )
    assert plan["triggered"]
    assert len(plan["branches"]) == 3
    assert set(plan["all_genera"]) == {f"Genus{letter}" for letter in "ABCDEF"}
    assert not plan["unplaced_genera"]

    selected_anchor = plan["branches"][0]["anchor_genera"][0]
    evidence = tmp_path / "evidence.csv"
    evidence.write_text(
        f"taxon,selected\n{selected_anchor},True\n", encoding="utf-8"
    )
    allowed = tmp_path / "allowed.txt"
    result = expand_genus_branches(
        plan_path, evidence, allowed, tmp_path / "expansion.json"
    )
    assert result["mode"] == "branch_screened"
    assert result["selected_genera"] == 4
    assert selected_anchor in allowed.read_text(encoding="utf-8")


def test_prefetch_allowed_groups(tmp_path: Path) -> None:
    backbone = tmp_path / "backbone"
    backbone.mkdir()
    (backbone / "4691.fasta").write_text(
        ">Rosales_FamA_GenusA_species_X\nATGACG\n", encoding="utf-8"
    )
    client = TinyKew(tmp_path / "kew")
    for letter in "ABCDEF":
        (tmp_path / "kew" / f"{letter}.fasta").write_text(
            ">4691 Gene_Name:test\nATGACG\n", encoding="utf-8"
        )
    report = BalancedKewPanelBuilder(
        backbone, tmp_path / "unused", representatives_per_taxon=1,
        minimum_recovered_loci=1, download_workers=1,
    ).prefetch("genus", ["FamA"], client, allowed_groups={"GenusA"})  # type: ignore[arg-type]
    assert report["candidate_groups"] == 1
    assert report["usable_recoveries"] == 1


def test_dominant_anchor_narrows_weak_selected_branches(tmp_path: Path) -> None:
    branches = [
        {
            "id": f"B{index:03d}", "family": "FamA", "genera": [f"Genus{index}"],
            "anchor_genera": [f"Genus{index}"], "neighbor": None,
        }
        for index in range(10)
    ]
    plan = tmp_path / "plan.json"
    plan.write_text(json.dumps({
        "branches": branches, "all_genera": [f"Genus{index}" for index in range(10)],
        "unplaced_genera": [],
    }), encoding="utf-8")
    evidence = tmp_path / "evidence.csv"
    header = "taxon,selected,locus_votes,aggregate_score,aggregate_z\n"
    rows = ["Genus0,True,47,152,8"]
    rows.extend(f"Genus{index},True,{10 if index == 1 else 2},{40 if index == 1 else 10},1"
                for index in range(1, 10))
    evidence.write_text(header + "\n".join(rows) + "\n", encoding="utf-8")
    allowed = tmp_path / "allowed.txt"
    result = expand_genus_branches(plan, evidence, allowed, tmp_path / "report.json")
    assert result["mode"] == "branch_screened"
    assert allowed.read_text(encoding="utf-8").splitlines() == ["Genus0", "Genus1"]


def test_diffuse_anchors_never_fall_back_to_full_family(tmp_path: Path) -> None:
    branches = [
        {
            "id": f"B{index:03d}", "family": "FamA", "genera": [f"Genus{index}"],
            "anchor_genera": [f"Genus{index}"], "neighbor": None,
        }
        for index in range(10)
    ]
    plan = tmp_path / "plan.json"
    plan.write_text(json.dumps({
        "branches": branches, "all_genera": [f"Genus{index}" for index in range(10)],
        "unplaced_genera": [],
    }), encoding="utf-8")
    evidence = tmp_path / "evidence.csv"
    evidence.write_text(
        "taxon,selected,locus_votes,aggregate_score,aggregate_z\n"
        + "".join(f"Genus{index},True,10,20,1\n" for index in range(10)),
        encoding="utf-8",
    )
    result = expand_genus_branches(
        plan, evidence, tmp_path / "allowed.txt", tmp_path / "report.json"
    )
    assert result["mode"] == "branch_screened"
    assert result["selected_genera"] == 3
    assert (tmp_path / "allowed.txt").read_text().splitlines() == [
        "Genus0", "Genus1", "Genus2",
    ]


def test_large_family_is_cut_progressively_without_global_anchor_fallback(
    tmp_path: Path,
) -> None:
    backbone = tmp_path / "backbone"
    backbone.mkdir()
    (backbone / "4691.fasta").write_text(
        ">Rosales_FamA_GenusA_species_X\nATGACG\n", encoding="utf-8"
    )
    parents = tmp_path / "parents.txt"
    parents.write_text("FamA\n", encoding="utf-8")
    root = tmp_path / "kew"
    root.mkdir()
    tree = root / "tree.nwk"
    genera = [f"Genus{i:02d}" for i in range(24)]
    tips = [f"Rosales_FamA_{genus}_species_{i:02d}" for i, genus in enumerate(genera)]
    tree.write_text("(" + ",".join(tips) + ");\n", encoding="utf-8")

    class ManyKew:
        def fetch_species_tree(self) -> Path:
            return tree

        def records(self) -> list[KewRecord]:
            return [
                KewRecord("INSDC", f"{i:02d}", "read", f"{genus} species", "PAFTOL", f"{i}.fasta")
                for i, genus in enumerate(genera)
            ]

        def sync(self) -> dict:
            return {"release": "test"}

    client = ManyKew()
    plan_path = tmp_path / "round1.json"
    plan = plan_genus_branches(
        client, backbone, parents, plan_path,
        threshold=8, max_group_genera=4, max_anchor_genera=8,
        max_branches_per_family=4,
    )
    assert plan["triggered"]
    assert len(plan["branches"]) == 4
    assert len(plan["anchor_genera"]) <= 8
    assert {genus for branch in plan["branches"] for genus in branch["genera"]} == set(genera)

    winner = plan["branches"][0]["anchor_genera"][0]
    evidence = tmp_path / "evidence.csv"
    evidence.write_text(
        "taxon,selected,locus_votes,aggregate_score\n"
        f"{winner},True,20,100\n", encoding="utf-8"
    )
    allowed = tmp_path / "allowed.txt"
    expanded = expand_genus_branches(
        plan_path, evidence, allowed, tmp_path / "expanded.json"
    )
    assert 0 < expanded["selected_genera"] < len(genera)
    second = plan_genus_branches(
        client, backbone, parents, tmp_path / "round2.json",
        threshold=8, max_group_genera=4, max_anchor_genera=8,
        max_branches_per_family=4, allowed_genera_file=allowed,
    )
    assert set(second["all_genera"]) == set(allowed.read_text().splitlines())
