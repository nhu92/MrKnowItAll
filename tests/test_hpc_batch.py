"""Tests for safe login-node orchestration on Zaratan."""

from __future__ import annotations

import importlib.util
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pytest


SCRIPT = Path(__file__).resolve().parents[1] / "examples/hpc/zaratan_batch_to_genus.py"
SPEC = importlib.util.spec_from_file_location("zaratan_batch_to_genus", SCRIPT)
assert SPEC and SPEC.loader
batch = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(batch)


def test_helper_commands_limit_blas_threads(monkeypatch):
    seen = {}

    def fake_run(args, **kwargs):
        seen.update(kwargs["env"])
        return subprocess.CompletedProcess(args, 0, stdout="ok\n", stderr="")

    monkeypatch.setattr(batch.subprocess, "run", fake_run)
    original = {"OPENBLAS_NUM_THREADS": "128"}
    assert batch.run_command(["true"], env=original) == "ok"
    assert original["OPENBLAS_NUM_THREADS"] == "128"
    assert all(seen[name] == "1" for name in batch.NUMERIC_THREAD_VARIABLES)


@pytest.mark.skipif(os.name != "posix", reason="fcntl locking is used on Zaratan")
def test_concurrent_controllers_submit_stage_only_once(tmp_path, monkeypatch):
    source = tmp_path / "stage.sbatch"
    source.write_text("#!/bin/bash\n", encoding="utf-8")
    stage_dir = tmp_path / "stage"
    submitted = []

    def fake_run(args, **kwargs):
        submitted.append(args)
        return "12345"

    monkeypatch.setattr(batch, "run_command", fake_run)
    monkeypatch.setattr(batch, "wait_for_stage", lambda *args: None)

    def submit():
        return batch.submit_stage("sample/order", source, stage_dir, {}, "1G", 1)

    with ThreadPoolExecutor(max_workers=2) as executor:
        assert list(executor.map(lambda _: submit(), range(2))) == [
            stage_dir / "results", stage_dir / "results"
        ]
    assert len(submitted) == 1
    assert (stage_dir / "JOB_ID").read_text(encoding="utf-8") == "12345\n"
