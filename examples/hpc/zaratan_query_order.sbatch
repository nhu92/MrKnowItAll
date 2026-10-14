#!/bin/bash
#SBATCH --job-name=mkia_order
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

PROJECT=${SPROUT_PROJECT:?Set SPROUT_PROJECT}
READ1=${SPROUT_READ1:?Set SPROUT_READ1}
READ2=${SPROUT_READ2:?Set SPROUT_READ2}
OUTPUT=${SPROUT_STAGE_OUTPUT:?Set SPROUT_STAGE_OUTPUT}
SPROUT_SOURCE=${SPROUT_SOURCE:?Set SPROUT_SOURCE}
MRKNOWITALL_SOURCE=${MRKNOWITALL_SOURCE:?Set MRKNOWITALL_SOURCE}
TARGET_FILE=${SPROUT_TARGET_FILE:-$SPROUT_SOURCE/sample_data/50targetfiles.fasta}
SOURCE_GENE_LIST=${SPROUT_GENE_LIST:-$SPROUT_SOURCE/sample_data/gene.list.txt}
ORDER_REFERENCE=${SPROUT_ORDER_REFERENCE:-$SPROUT_SOURCE/sample_data/order106_refs_50genes}
CONDA_SH=${SPROUT_CONDA_SH:-$HOME/miniconda3/etc/profile.d/conda.sh}
CONDA_ENV=${SPROUT_CONDA_ENV:-sprout}
TREE_WORKERS=${SPROUT_TREE_WORKERS:-4}
MIN_LOCUS_VOTES=${SPROUT_MIN_LOCUS_VOTES:-2}
MAX_ORDER_CANDIDATES=${SPROUT_MAX_ORDER_CANDIDATES:-8}

[[ "$PROJECT" =~ ^[A-Za-z0-9_.-]+$ ]] || { echo "Invalid project token" >&2; exit 2; }
for path in "$READ1" "$READ2" "$TARGET_FILE" "$SOURCE_GENE_LIST" \
  "$ORDER_REFERENCE" "$SPROUT_SOURCE/01_exons_assembly.py" \
  "$SPROUT_SOURCE/02_exon_trees.py" "$SPROUT_SOURCE/03_distance_matrices.py" \
  "$SPROUT_SOURCE/04_prediction.py" "$SPROUT_SOURCE/export_pea.py" \
  "$MRKNOWITALL_SOURCE/pyproject.toml"; do
  [[ -e "$path" ]] || { echo "Missing input: $path" >&2; exit 2; }
done
[[ ! -e "$OUTPUT" ]] || { echo "Refusing to overwrite: $OUTPUT" >&2; exit 4; }

work_root=${SLURM_TMPDIR:-/tmp/mrknowitall-${SLURM_JOB_ID:-manual}}
work="$work_root/$PROJECT-order"
mkdir -p "$work/query" "$work/order"

source "$CONDA_SH"
conda activate "$CONDA_ENV"
export PYTHONPATH="$MRKNOWITALL_SOURCE/src${PYTHONPATH:+:$PYTHONPATH}"

fastq_suffix() {
  case "$1" in
    *.fastq.gz) printf '.fastq.gz' ;;
    *.fq.gz) printf '.fq.gz' ;;
    *.fastq) printf '.fastq' ;;
    *.fq) printf '.fq' ;;
    *) echo "Unsupported FASTQ suffix: $1" >&2; return 2 ;;
  esac
}
input_r1="$work/query/$PROJECT.R1$(fastq_suffix "$READ1")"
input_r2="$work/query/$PROJECT.R2$(fastq_suffix "$READ2")"
ln -s "$READ1" "$input_r1"
ln -s "$READ2" "$input_r2"
cp "$SOURCE_GENE_LIST" "$work/query/gene_list.txt"
seqkit stats -a "$input_r1" "$input_r2" > "$work/query/input_reads.stats.txt"

cd "$work/query"
python "$SPROUT_SOURCE/01_exons_assembly.py" \
  -t "${SLURM_CPUS_PER_TASK:-32}" -r1 "$input_r1" -r2 "$input_r2" \
  -p "$PROJECT" -g gene_list.txt -m "$TARGET_FILE" \
  --output_hyb hyb --output_exon exon

cd "$work/order"
python "$SPROUT_SOURCE/02_exon_trees.py" \
  -t "$TREE_WORKERS" -p "$PROJECT" -g "$work/query/gene_list.txt" \
  -e "$work/query/exon" -r "$ORDER_REFERENCE" --output_phylo phylo
python "$SPROUT_SOURCE/03_distance_matrices.py" \
  -t "$TREE_WORKERS" -p "$PROJECT" -g "$work/query/gene_list.txt" \
  --threshold 1 --input_phylo phylo --output_tree matrix
cp "$PROJECT.cumulative_dist.csv" "$PROJECT.order.cumulative_dist.csv"
python "$SPROUT_SOURCE/04_prediction.py" \
  -i "$PROJECT.order.cumulative_dist.csv" -o predictions.csv \
  -tl o -z 0.9 -to z_candidates.txt
python -m sprout_refbuilder.cli select-candidates \
  --matrix-dir matrix --predictions predictions.csv \
  --level order --project "$PROJECT" \
  --minimum-locus-votes "$MIN_LOCUS_VOTES" --minimum-z 0.9 \
  --always-top 1 --max-candidates "$MAX_ORDER_CANDIDATES" \
  --output candidates.txt --report candidate_evidence.csv
[[ -s candidates.txt ]] || { echo "No order candidate" >&2; exit 6; }
context_args=()
if [[ -f "${SPROUT_REF_CACHE:-}/treeoflife.current.tree" ]]; then
  context_args=(--context-tree "$SPROUT_REF_CACHE/treeoflife.current.tree")
fi
python "$SPROUT_SOURCE/export_pea.py" stage \
  --project "$PROJECT" --level order --predictions predictions.csv \
  --evidence candidate_evidence.csv --candidates candidates.txt \
  --gene-list "$work/query/gene_list.txt" --matrix-dir matrix --phylo-dir phylo \
  "${context_args[@]}" \
  --z-threshold 0.9 --output stage.pea.json --markdown stage.report.md

mkdir -p "$OUTPUT"
cp -a "$work/query/exon" "$OUTPUT/exon"
cp "$work/query/gene_list.txt" "$work/query/input_reads.stats.txt" "$OUTPUT/"
cp candidates.txt predictions.csv candidate_evidence.csv z_candidates.txt \
  stage.pea.json stage.report.md \
  "$PROJECT.order.cumulative_dist.csv" "$OUTPUT/"
{
  echo "finished=$(date -Is)"
  echo "job_id=${SLURM_JOB_ID:-none}"
  echo "project=$PROJECT"
  echo "read1=$READ1"
  echo "read2=$READ2"
  echo "order_candidates:"
  cat candidates.txt
} | tee "$OUTPUT/RUN_SUMMARY.txt"
touch "$OUTPUT/SUCCESS"
