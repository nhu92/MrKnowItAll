#!/bin/bash
#SBATCH --job-name=prepare_871ref
#SBATCH --partition=nocona
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -Eeuo pipefail
trap 'code=$?; echo "FAILED exit=${code} line=${LINENO} $(date -Is)"; exit "$code"' ERR

ARCHIVE_DIR=${MRKNOWITALL_ARCHIVE_DIR:-/lustre/scratch/nhu/202609/871ref}
REPO=${MRKNOWITALL_SOURCE:-$ARCHIVE_DIR/MrKnowItAll}
OUT=${MRKNOWITALL_BACKBONE:-$ARCHIVE_DIR/backbone_corrected}
CONDA_SH=${SPROUT_CONDA_SH:-$HOME/miniconda3/etc/profile.d/conda.sh}
CONDA_ENV=${SPROUT_CONDA_ENV:-sprout}

ARCHIVE=$(find "$ARCHIVE_DIR" -maxdepth 1 -type f \
  -name 'angiosperms_353_v2_interim_targetfile_gene_alignments.zip' -print -quit)
[[ -n "$ARCHIVE" ]] || {
  echo "Backbone ZIP was not found directly inside $ARCHIVE_DIR" >&2
  exit 2
}
[[ -f "$REPO/pyproject.toml" ]] || { echo "MrKnowItAll checkout not found: $REPO" >&2; exit 2; }
if [[ -d "$OUT" ]] && find "$OUT" -mindepth 1 -print -quit | grep -q .; then
  echo "Refusing to overwrite non-empty backbone directory: $OUT" >&2
  exit 4
fi

source "$CONDA_SH"
conda activate "$CONDA_ENV"
python -m pip install -e "${REPO}[dev]"

python -m pytest -q "$REPO/tests"
python -m ruff check "$REPO/src" "$REPO/tests"
sprout-ref install-backbone --archive "$ARCHIVE" --destination "$OUT"
sprout-ref inspect-backbone "$OUT" | tee backbone_inspection.json

python - "$OUT/backbone_install.json" <<'PY'
import json
import sys
from pathlib import Path

manifest = json.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
corrections = manifest.get("taxonomy_corrections") or {}
summary = {
    "schema_version": manifest.get("schema_version"),
    "locus_count": manifest.get("locus_count"),
    "sequence_count": manifest.get("sequence_count"),
    "specimen_count": manifest.get("specimen_count"),
    "species_count": manifest.get("species_count"),
    "correction_rules": corrections.get("rules"),
    "replaced_locus_sequences": corrections.get("replaced_locus_sequences"),
    "excluded_locus_sequences": corrections.get("excluded_locus_sequences"),
    "shadowed_by_exclusion": corrections.get("shadowed_by_exclusion"),
    "unmatched_replacements": corrections.get("unmatched_replacements"),
    "unmatched_exclusions": corrections.get("unmatched_exclusions"),
}
print(json.dumps(summary, indent=2))
if manifest.get("locus_count") != 353 or corrections.get("rules") != 57:
    raise SystemExit("Backbone invariant or correction inventory failed")
if not corrections.get("replaced_locus_sequences") or not corrections.get("excluded_locus_sequences"):
    raise SystemExit("No taxonomy corrections or exclusions matched the backbone")
PY

{
  echo "finished=$(date -Is)"
  echo "archive=$ARCHIVE"
  echo "repository_commit=$(git -C "$REPO" rev-parse HEAD)"
  echo "backbone=$OUT"
} | tee RUN_SUMMARY.txt
touch SUCCESS
