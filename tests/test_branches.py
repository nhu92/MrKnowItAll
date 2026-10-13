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
