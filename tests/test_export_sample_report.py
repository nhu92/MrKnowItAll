"""Regression tests for combining multiple genus-screening rounds."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest


SCRIPT = Path(__file__).resolve().parents[1] / "examples/hpc/export_sample_report.py"
SPEC = importlib.util.spec_from_file_location("export_sample_report", SCRIPT)
assert SPEC and SPEC.loader
report = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(report)


def write_stage(directory: Path, sample: str, stage_id: str) -> Path:
    directory.mkdir(parents=True)
    (directory / "SUCCESS").touch()
    path = directory / "stage.pea.json"
    path.write_text(
        json.dumps({"sample_id": sample, "stages": [{"id": stage_id}]}),
        encoding="utf-8",
    )
    return path


def test_multi_round_report_uses_unique_ids_without_changing_stage_outputs(tmp_path, monkeypatch):
    sample = tmp_path / "EMseq001"
    source_order = write_stage(sample / "order/results", "emseq001", "order")
    source_family = write_stage(sample / "family/results", "emseq001", "family")
    source_round_1 = write_stage(sample / "branch/round_01/results", "emseq001", "branch")
    source_round_2 = write_stage(sample / "branch/round_02/results", "emseq001", "branch")
    sprout = tmp_path / "SPrOUT"
    sprout.mkdir()
    (sprout / "export_pea.py").touch()

    def fake_merge(command, check):
        assert check
        assert command[2] == "merge"
        inputs = command[command.index("--input") + 1 : command.index("--output")]
        ids = [json.loads(Path(path).read_text(encoding="utf-8"))["stages"][0]["id"]
               for path in inputs]
        assert ids == ["order", "family", "branch_round_01", "branch_round_02"]

    monkeypatch.setattr(report.subprocess, "run", fake_merge)
    report.export_sample(sample, sprout)
    report.export_sample(sample, sprout)
    assert json.loads(source_round_1.read_text(encoding="utf-8"))["stages"][0]["id"] == "branch"
    assert json.loads(source_round_2.read_text(encoding="utf-8"))["stages"][0]["id"] == "branch"
    assert json.loads(source_order.read_text(encoding="utf-8"))["stages"][0]["id"] == "order"
    assert json.loads(source_family.read_text(encoding="utf-8"))["stages"][0]["id"] == "family"


def test_unexpected_stage_id_is_rejected(tmp_path):
    source = write_stage(tmp_path / "results", "emseq001", "genus")
    with pytest.raises(ValueError, match="Unexpected stage ID"):
        report.stage_file_for_merge(source, "branch_round_01", tmp_path)
