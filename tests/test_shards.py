import json
from pathlib import Path

import pytest

from sprout_refbuilder.shards import plan_genus_shards, shortlist_genus_shards


def test_phylogenetic_shards_cover_each_genus_once(tmp_path: Path) -> None:
    plan = tmp_path / "branches.json"
    plan.write_text(json.dumps({
        "all_genera": ["A", "B", "C", "D", "E", "F", "G"],
        "branches": [
            {"id": "B001", "genera": ["A", "B", "C"]},
            {"id": "B002", "genera": ["C", "D"]},
            {"id": "B003", "genera": ["E", "F"]},
        ],
    }), encoding="utf-8")
    allowed = tmp_path / "allowed.txt"
    allowed.write_text("A\nB\nC\nD\nE\nF\nG\n", encoding="utf-8")
    report = plan_genus_shards(plan, allowed, tmp_path / "shards", max_genera_per_shard=4)
    assert [row["genera"] for row in report["shards"]] == [
        ["A", "B", "C", "D"], ["E", "F", "G"]
    ]
    assert sorted(g for row in report["shards"] for g in row["genera"]) == report["allowed_genera"]
    assert (tmp_path / "shards" / "S001.txt").read_text(encoding="utf-8") == "A\nB\nC\nD\n"


def test_shortlist_uses_equal_within_shard_quota(tmp_path: Path) -> None:
    shard_root = tmp_path / "shards"
    shard_root.mkdir()
    plan = shard_root / "shards.json"
    plan.write_text(json.dumps({
        "allowed_genera": ["A", "B", "C", "D", "E", "F"],
        "shards": [
            {"id": "S001", "genera": ["A", "B", "C"]},
            {"id": "S002", "genera": ["D", "E", "F"]},
        ],
    }), encoding="utf-8")
    for shard, candidates in (("S001", "B\nA\nC\n"), ("S002", "F\nD\nE\n")):
        result = shard_root / shard / "results"
        result.mkdir(parents=True)
        (result / "candidates.txt").write_text(candidates, encoding="utf-8")
    output = shard_root / "shortlist.txt"
    report = shortlist_genus_shards(
        plan, shard_root, output, shard_root / "shortlist.json", max_final_genera=4
    )
    assert output.read_text(encoding="utf-8").splitlines() == ["B", "A", "F", "D"]
    assert report["quota_per_shard"] == 2
    assert "scores-not-compared" in report["selection_policy"]


def test_shortlist_rejects_too_many_shards_for_final_panel(tmp_path: Path) -> None:
    plan = tmp_path / "shards.json"
    plan.write_text(json.dumps({
        "allowed_genera": ["A", "B"],
        "shards": [{"id": "S001", "genera": ["A"]}, {"id": "S002", "genera": ["B"]}],
    }), encoding="utf-8")
    with pytest.raises(ValueError, match="at least one genus per shard"):
        shortlist_genus_shards(plan, tmp_path, tmp_path / "out.txt", tmp_path / "report.json",
                                max_final_genera=1)
