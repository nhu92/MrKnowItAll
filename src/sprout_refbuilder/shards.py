"""Phylogenetic genus shards and a balanced shortlist for a final joint panel."""

from __future__ import annotations

import json
from pathlib import Path


def plan_genus_shards(
    branch_plan: Path,
    allowed_file: Path,
    output_dir: Path,
    *,
    max_genera_per_shard: int = 40,
) -> dict:
    """Partition selected genera without splitting tree branches when possible."""
    if max_genera_per_shard < 1:
        raise ValueError("max_genera_per_shard must be positive")
    plan = json.loads(branch_plan.read_text(encoding="utf-8"))
    allowed = set(allowed_file.read_text(encoding="utf-8-sig").split())
    possible = set(plan["all_genera"])
    if not allowed or allowed - possible:
        raise ValueError("Allowed genera must be a nonempty subset of the branch plan")

    units: list[tuple[str, list[str]]] = []
    assigned: set[str] = set()
    for branch in plan["branches"]:
        genera = sorted((set(branch["genera"]) & allowed) - assigned)
        if genera:
            units.append((branch["id"], genera))
            assigned.update(genera)
    remaining = sorted(allowed - assigned)
    for start in range(0, len(remaining), max_genera_per_shard):
        units.append(("unplaced", remaining[start : start + max_genera_per_shard]))

    shards = []
    current: list[str] = []
    sources: list[str] = []

    def flush() -> None:
        if current:
            shards.append({
                "id": f"S{len(shards) + 1:03d}",
                "genera": list(current),
                "source_branches": list(dict.fromkeys(sources)),
            })
            current.clear()
            sources.clear()

    for source, genera in units:
        for start in range(0, len(genera), max_genera_per_shard):
            part = genera[start : start + max_genera_per_shard]
            if current and len(current) + len(part) > max_genera_per_shard:
                flush()
            current.extend(part)
            sources.append(source)
    flush()
    if set().union(*(set(shard["genera"]) for shard in shards)) != allowed:
        raise RuntimeError("Genus shard plan lost candidates")

    output_dir.mkdir(parents=True, exist_ok=True)
    for shard in shards:
        (output_dir / f"{shard['id']}.txt").write_text(
            "".join(f"{name}\n" for name in shard["genera"]), encoding="utf-8"
        )
    report = {
        "schema_version": 1,
        "allowed_genera": sorted(allowed),
        "max_genera_per_shard": max_genera_per_shard,
        "shards": shards,
    }
    (output_dir / "shards.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    return report


def shortlist_genus_shards(
    shard_plan: Path,
    shard_runs: Path,
    output: Path,
    report_path: Path,
    *,
    max_final_genera: int = 80,
) -> dict:
    """Keep a fair within-shard shortlist; do not compare scores across panels."""
    plan = json.loads(shard_plan.read_text(encoding="utf-8"))
    shards = plan["shards"]
    if max_final_genera < len(shards):
        raise ValueError("max_final_genera must cover at least one genus per shard")
    quota = max_final_genera // len(shards)
    chosen: list[str] = []
    details = []
    for shard in shards:
        candidate_path = shard_runs / shard["id"] / "results" / "candidates.txt"
        allowed = set(shard["genera"])
        candidates = list(dict.fromkeys(candidate_path.read_text(encoding="utf-8").split()))
        candidates = [name for name in candidates if name in allowed]
        if not candidates:
            raise ValueError(f"No genus candidates in shard {shard['id']}")
        selected = candidates[:quota]
        chosen.extend(selected)
        details.append({"id": shard["id"], "selected": selected, "available": len(candidates)})
    if len(chosen) != len(set(chosen)):
        raise RuntimeError("A genus appeared in multiple shard shortlists")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("".join(f"{name}\n" for name in chosen), encoding="utf-8")
    report = {
        "schema_version": 1,
        "selection_policy": "equal-quota-within-shard; scores-not-compared-across-shards",
        "quota_per_shard": quota,
        "shortlisted_genera": len(chosen),
        "source_genera": len(plan["allowed_genera"]),
        "shards": details,
    }
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    return report
