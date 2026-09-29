#!/bin/bash
#SBATCH --job-name=sprout_order
#SBATCH --partition=nocona
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -Eeuo pipefail
trap 'code=$?; echo "FAILED exit=${code} line=${LINENO} $(date -Is)"; exit "$code"' ERR

RUN=${SPROUT_RUN_DIR:-$PWD}
SRC=${SPROUT_SOURCE:-/lustre/scratch/nhu/202606/SPrOUT}
CONDA_SH=${SPROUT_CONDA_SH:-$HOME/miniconda3/etc/profile.d/conda.sh}
CONDA_ENV=${SPROUT_CONDA_ENV:-sprout}
THREADS=${SLURM_CPUS_PER_TASK:-32}
TREE_WORKERS=${SPROUT_TREE_WORKERS:-4}
SAMPLE_FRACTION=${SPROUT_SAMPLE_FRACTION:-0.04}
SEED=${SPROUT_SAMPLE_SEED:-20260929}
PROJECT=${SPROUT_PROJECT:-sprout_test}

cd "$RUN"
if [[ -e SUCCESS || -e hyb || -e exon || -e order_phylo || -e order_matrix ]]; then
  echo "Refusing to reuse an order directory containing outputs: $RUN" >&2
  exit 4
fi
for path in \
  "$SRC/01_exons_assembly.py" \
  "$SRC/02_exon_trees.py" \
  "$SRC/03_distance_matrices.py" \
  "$SRC/04_prediction.py" \
  "$SRC/01x02x03.R1.fastq" \
  "$SRC/01x02x03.R2.fastq" \
  "$SRC/sample_data/gene.list.txt" \
  "$SRC/sample_data/50targetfiles.fasta" \
  "$SRC/sample_data/order106_refs_50genes"; do
  [[ -e "$path" ]] || { echo "Missing required input: $path" >&2; exit 2; }
done

source "$CONDA_SH"
conda activate "$CONDA_ENV"
for tool in python fastp hybpiper mafft fasttree trimal seqkit bwa; do
  command -v "$tool" || { echo "Missing required tool: $tool" >&2; exit 127; }
done

{
  echo "started=$(date -Is)"
  echo "host=$(hostname -f)"
  echo "job_id=${SLURM_JOB_ID:-none}"
  echo "source=$SRC"
  echo "sample_fraction=$SAMPLE_FRACTION"
  echo "sample_seed=$SEED"
  echo "assembly_threads=$THREADS"
  echo "tree_workers=$TREE_WORKERS"
  echo "project=$PROJECT"
  git -C "$SRC" rev-parse HEAD 2>/dev/null | sed 's/^/sprout_commit=/' || true
} | tee RUN_PROVENANCE.txt

seqkit sample -p "$SAMPLE_FRACTION" -s "$SEED" \
  "$SRC/01x02x03.R1.fastq" -o "$PROJECT.R1.fastq"
seqkit sample -p "$SAMPLE_FRACTION" -s "$SEED" \
  "$SRC/01x02x03.R2.fastq" -o "$PROJECT.R2.fastq"
seqkit stats -a "$PROJECT.R1.fastq" "$PROJECT.R2.fastq" | tee sampled_reads.stats.txt
r1_count=$(awk 'END {print NR/4}' "$PROJECT.R1.fastq")
r2_count=$(awk 'END {print NR/4}' "$PROJECT.R2.fastq")
[[ "$r1_count" -eq "$r2_count" ]] || {
  echo "Paired read counts differ: R1=$r1_count R2=$r2_count" >&2
  exit 3
}

cp "$SRC/sample_data/gene.list.txt" gene_list.txt
python "$SRC/01_exons_assembly.py" \
  -t "$THREADS" \
  -r1 "$PROJECT.R1.fastq" \
  -r2 "$PROJECT.R2.fastq" \
  -p "$PROJECT" \
  -g gene_list.txt \
  -m "$SRC/sample_data/50targetfiles.fasta" \
  --output_hyb hyb \
  --output_exon exon

python "$SRC/02_exon_trees.py" \
  -t "$TREE_WORKERS" -p "$PROJECT" -g gene_list.txt \
  -e exon -r "$SRC/sample_data/order106_refs_50genes" \
  --output_phylo order_phylo
if ! grep -R -F -q --include='*.tre' "$PROJECT" order_phylo; then
  echo "No query leaf contains project token '$PROJECT'" >&2
  exit 5
fi
python "$SRC/03_distance_matrices.py" \
  -t "$TREE_WORKERS" -p "$PROJECT" -g gene_list.txt --threshold 1 \
  --input_phylo order_phylo --output_tree order_matrix
cp "$PROJECT.cumulative_dist.csv" "$PROJECT.order.cumulative_dist.csv"
python "$SRC/04_prediction.py" \
  -i "$PROJECT.order.cumulative_dist.csv" \
  -o "$PROJECT.predictions_order.csv" \
  -tl o -z 0.9 \
  -to "$PROJECT.order_candidates.txt"
[[ -s "$PROJECT.order_candidates.txt" ]] || {
  echo "No order passed z > 0.9; stopping before Kew refinement" >&2
  exit 6
}

{
  echo "finished=$(date -Is)"
  echo "r1_reads=$r1_count"
  echo "r2_reads=$r2_count"
  echo "exon_files=$(find exon -type f -name '*.fasta' | wc -l)"
  echo "order_trees=$(find order_phylo -type f -name '*.tre' | wc -l)"
  echo "order_ranking:"
  cat "$PROJECT.predictions_order.csv"
  echo "order_candidates:"
  cat "$PROJECT.order_candidates.txt"
} | tee RUN_SUMMARY.txt
touch SUCCESS
