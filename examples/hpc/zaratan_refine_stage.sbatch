#!/bin/bash
#SBATCH --job-name=mkia_refine
#SBATCH --partition=standard
#SBATCH --account=jmeng-prj-aac
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G
#SBATCH --tmp=102400
#SBATCH --time=08:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -Eeuo pipefail
trap 'code=$?; echo "FAILED exit=${code} line=${LINENO} $(date -Is)"; exit "$code"' ERR

LEVEL=${SPROUT_LEVEL:?Set SPROUT_LEVEL to family or genus}
PROJECT=${SPROUT_PROJECT:?Set SPROUT_PROJECT to the query token}
PARENT_FILE=${SPROUT_PARENT_FILE:?Set SPROUT_PARENT_FILE to the previous candidate list}
QUERY_EXON=${SPROUT_QUERY_EXON:?Set SPROUT_QUERY_EXON to the query exon directory}
QUERY_GENE_LIST=${SPROUT_QUERY_GENE_LIST:?Set SPROUT_QUERY_GENE_LIST}
OUTPUT=${SPROUT_STAGE_OUTPUT:?Set SPROUT_STAGE_OUTPUT to a new persistent directory}
SPROUT_SOURCE=${SPROUT_SOURCE:?Set SPROUT_SOURCE}
MRKNOWITALL_SOURCE=${MRKNOWITALL_SOURCE:?Set MRKNOWITALL_SOURCE}
MRKNOWITALL_BACKBONE=${MRKNOWITALL_BACKBONE:?Set MRKNOWITALL_BACKBONE}
KEW_CACHE=${SPROUT_REF_CACHE:?Set SPROUT_REF_CACHE to the prefetched Kew cache}
CONDA_SH=${SPROUT_CONDA_SH:-$HOME/miniconda3/etc/profile.d/conda.sh}
CONDA_ENV=${SPROUT_CONDA_ENV:-sprout}
TREE_WORKERS=${SPROUT_TREE_WORKERS:-4}
BUILD_WORKERS=${SPROUT_REF_BUILD_WORKERS:-${SLURM_CPUS_PER_TASK:-32}}
MIN_LOCUS_VOTES=${SPROUT_MIN_LOCUS_VOTES:-2}

case "$LEVEL" in
  family) REPRESENTATIVES=4; MAX_CANDIDATES=12 ;;
  genus) REPRESENTATIVES=2; MAX_CANDIDATES=20 ;;
  *) echo "SPROUT_LEVEL must be family or genus" >&2; exit 2 ;;
esac
REPRESENTATIVES=${SPROUT_REF_REPRESENTATIVES:-$REPRESENTATIVES}
POOL_MULTIPLIER=${SPROUT_REF_POOL_MULTIPLIER:-2}
ALLOWED_TAXA_FILE=${SPROUT_ALLOWED_TAXA_FILE:-}
if [[ -n "$ALLOWED_TAXA_FILE" ]]; then
  [[ -s "$ALLOWED_TAXA_FILE" ]] || { echo "Missing allowed taxa: $ALLOWED_TAXA_FILE" >&2; exit 2; }
fi

for path in "$PARENT_FILE" "$QUERY_EXON" "$QUERY_GENE_LIST" "$SPROUT_SOURCE" \
  "$SPROUT_SOURCE/export_pea.py" "$MRKNOWITALL_SOURCE/pyproject.toml" \
  "$MRKNOWITALL_BACKBONE" "$KEW_CACHE"; do
  [[ -e "$path" ]] || { echo "Missing required input: $path" >&2; exit 2; }
done
[[ ! -e "$OUTPUT" ]] || { echo "Refusing to overwrite: $OUTPUT" >&2; exit 4; }

work_root=${SLURM_TMPDIR:-/tmp/mrknowitall-${SLURM_JOB_ID:-manual}}
work="$work_root/$PROJECT-$LEVEL"
mkdir -p "$work"
cp -a "$QUERY_EXON" "$work/exon"
cp "$QUERY_GENE_LIST" "$work/query_gene_list.txt"
cp "$PARENT_FILE" "$work/parent_candidates.txt"
allowed_args=()
if [[ -n "$ALLOWED_TAXA_FILE" ]]; then
  cp "$ALLOWED_TAXA_FILE" "$work/allowed_taxa.txt"
  allowed_args=(--allowed-taxa-file allowed_taxa.txt)
fi
cd "$work"

source "$CONDA_SH"
conda activate "$CONDA_ENV"
export PYTHONPATH="$MRKNOWITALL_SOURCE/src${PYTHONPATH:+:$PYTHONPATH}"
export SPROUT_REF_OFFLINE=1

python -m sprout_refbuilder.cli --cache "$KEW_CACHE" build-kew-panel \
  --level "$LEVEL" --parent-file parent_candidates.txt \
  --backbone "$MRKNOWITALL_BACKBONE" --output reference \
  "${allowed_args[@]}" \
  --representatives "$REPRESENTATIVES" --pool-multiplier "$POOL_MULTIPLIER" \
  --minimum-recovered-loci 50 --download-workers 8 --threads "$BUILD_WORKERS"

grep -F -x -f query_gene_list.txt reference/gene.list.txt > gene_list.txt
[[ -s gene_list.txt ]] || { echo "No shared genes for $LEVEL" >&2; exit 7; }
python "$SPROUT_SOURCE/02_exon_trees.py" \
  -t "$TREE_WORKERS" -p "$PROJECT" -g gene_list.txt \
  -e exon -r reference/ref --output_phylo phylo
python "$SPROUT_SOURCE/03_distance_matrices.py" \
  -t "$TREE_WORKERS" -p "$PROJECT" -g gene_list.txt --threshold 1 \
  --input_phylo phylo --output_tree matrix
cp "$PROJECT.cumulative_dist.csv" "$PROJECT.$LEVEL.cumulative_dist.csv"
python "$SPROUT_SOURCE/04_prediction.py" \
  -i "$PROJECT.$LEVEL.cumulative_dist.csv" -o predictions.csv \
  -tl "${LEVEL:0:1}" -z 0 -to z_candidates.txt
python -m sprout_refbuilder.cli select-candidates \
  --matrix-dir matrix --predictions predictions.csv \
  --level "$LEVEL" --project "$PROJECT" \
  --minimum-locus-votes "$MIN_LOCUS_VOTES" --minimum-z 0 \
  --always-top 1 --max-candidates "$MAX_CANDIDATES" \
  --output candidates.txt --report candidate_evidence.csv
stage_id=$LEVEL
if [[ -n "$ALLOWED_TAXA_FILE" && "$LEVEL" == genus && "$REPRESENTATIVES" == 1 ]]; then
  stage_id=branch
fi
python "$SPROUT_SOURCE/export_pea.py" stage \
  --project "$PROJECT" --level "$LEVEL" --stage-id "$stage_id" \
  --predictions predictions.csv --evidence candidate_evidence.csv \
  --candidates candidates.txt --gene-list gene_list.txt --matrix-dir matrix --phylo-dir phylo \
  --context-tree "$KEW_CACHE/treeoflife.current.tree" \
  --z-threshold 0 --output stage.pea.json --markdown stage.report.md

mkdir -p "$OUTPUT"
cp candidates.txt predictions.csv candidate_evidence.csv z_candidates.txt \
  stage.pea.json stage.report.md \
  "$PROJECT.$LEVEL.cumulative_dist.csv" gene_list.txt "$OUTPUT/"
cp reference/panel_report.json reference/selection.tsv reference/qc.json "$OUTPUT/"
cp -a reference/ref "$OUTPUT/reference_ref"
{
  echo "finished=$(date -Is)"
  echo "job_id=${SLURM_JOB_ID:-none}"
  echo "level=$LEVEL"
  echo "project=$PROJECT"
  echo "parents:"
  cat parent_candidates.txt
  echo "candidates:"
  cat candidates.txt
} | tee "$OUTPUT/RUN_SUMMARY.txt"
touch "$OUTPUT/SUCCESS"
