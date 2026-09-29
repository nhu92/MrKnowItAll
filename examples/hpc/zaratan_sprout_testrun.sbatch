#!/bin/bash
#SBATCH --job-name=sprout_test
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -Eeuo pipefail
trap 'code=$?; echo "FAILED exit=${code} line=${LINENO} $(date -Is)"; exit "$code"' ERR

RUN=${SPROUT_RUN_DIR:-$PWD}
SRC=${SPROUT_SOURCE:-$(cd "$RUN/.." && pwd)}
CONDA_SH=${SPROUT_CONDA_SH:-$HOME/miniconda3/etc/profile.d/conda.sh}
CONDA_ENV=${SPROUT_CONDA_ENV:-sprout}
THREADS=${SLURM_CPUS_PER_TASK:-16}
TREE_WORKERS=${SPROUT_TREE_WORKERS:-4}
SAMPLE_FRACTION=${SPROUT_SAMPLE_FRACTION:-0.04}
SEED=${SPROUT_SAMPLE_SEED:-20260929}
PROJECT=${SPROUT_PROJECT:-sprout_test}

cd "$RUN"
if [[ -e SUCCESS || -e hyb || -e exon || -e order_phylo || -e family_phylo ]]; then
  echo "Refusing to reuse a run directory containing pipeline outputs: $RUN" >&2
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
  "$SRC/sample_data/50targetfiles.fasta"; do
  [[ -e "$path" ]] || { echo "Missing required input: $path" >&2; exit 2; }
done

source "$CONDA_SH"
conda activate "$CONDA_ENV"

{
  echo "started=$(date -Is)"
  echo "host=$(hostname -f)"
  echo "run=$RUN"
  echo "source=$SRC"
  echo "job_id=${SLURM_JOB_ID:-none}"
  echo "threads=$THREADS"
  echo "tree_workers=$TREE_WORKERS"
  echo "sample_fraction=$SAMPLE_FRACTION"
  echo "sample_seed=$SEED"
  echo "conda_prefix=$CONDA_PREFIX"
  git -C "$SRC" rev-parse HEAD 2>/dev/null | sed 's/^/sprout_commit=/' || true
} | tee RUN_PROVENANCE.txt

python --version
for tool in python fastp hybpiper mafft fasttree trimal seqkit bwa; do
  command -v "$tool" || { echo "Missing required tool: $tool" >&2; exit 127; }
done
conda list | grep -Ei \
  'hybpiper|mafft|fasttree|fastp|seqkit|trimal|biopython|scikit-learn|pandas|numpy' \
  > software_versions.txt || true

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

python "$SRC/01_exons_assembly.py" \
  -t "$THREADS" \
  -r1 "$PROJECT.R1.fastq" \
  -r2 "$PROJECT.R2.fastq" \
  -p "$PROJECT" \
  -g "$SRC/sample_data/gene.list.txt" \
  -m "$SRC/sample_data/50targetfiles.fasta" \
  --output_hyb hyb \
  --output_exon exon

python "$SRC/02_exon_trees.py" \
  -t "$TREE_WORKERS" \
  -p "$PROJECT" \
  -g gene_list.txt \
  -e exon \
  -r "$SRC/sample_data/order106_refs_50genes" \
  --output_phylo order_phylo
python "$SRC/03_distance_matrices.py" \
  -t "$TREE_WORKERS" \
  -p "$PROJECT" \
  -g gene_list.txt \
  --threshold 1 \
  --input_phylo order_phylo \
  --output_tree order_matrix
cp "$PROJECT.cumulative_dist.csv" "$PROJECT.order.cumulative_dist.csv"
python "$SRC/04_prediction.py" \
  -i "$PROJECT.order.cumulative_dist.csv" \
  -o "$PROJECT.predictions_order.csv" \
  -tl o \
  -z 0.9 \
  -to "$PROJECT.order_candidates.txt"

if [[ -s "$PROJECT.order_candidates.txt" ]]; then
  mkdir ref_family
  while IFS= read -r gene; do
    python "$SRC/misc/pick_match_list.py" \
      "$SRC/sample_data/family298_refs_50genes/${gene}.fasta" \
      "ref_family/${gene}.fasta" \
      "$PROJECT.order_candidates.txt"
  done < gene_list.txt
  python "$SRC/02_exon_trees.py" \
    -t "$TREE_WORKERS" \
    -p "$PROJECT" \
    -g gene_list.txt \
    -e exon \
    -r ref_family \
    --output_phylo family_phylo
  python "$SRC/03_distance_matrices.py" \
    -t "$TREE_WORKERS" \
    -p "$PROJECT" \
    -g gene_list.txt \
    --threshold 1 \
    --input_phylo family_phylo \
    --output_tree family_matrix
  cp "$PROJECT.cumulative_dist.csv" "$PROJECT.family.cumulative_dist.csv"
  python "$SRC/04_prediction.py" \
    -i "$PROJECT.family.cumulative_dist.csv" \
    -o "$PROJECT.predictions_family.csv" \
    -tl f \
    -z 0 \
    -to "$PROJECT.family_candidates.txt"
else
  echo "No order passed z > 0.9; family stage skipped." | tee family_stage_skipped.txt
fi

{
  echo "finished=$(date -Is)"
  echo "r1_reads=$r1_count"
  echo "r2_reads=$r2_count"
  echo "exon_files=$(find exon -type f -name '*.fasta' | wc -l)"
  echo "order_trees=$(find order_phylo -type f -name '*.tre' | wc -l)"
  echo "order_candidates:"
  cat "$PROJECT.order_candidates.txt"
  if [[ -f "$PROJECT.family_candidates.txt" ]]; then
    echo "family_candidates:"
    cat "$PROJECT.family_candidates.txt"
  fi
} | tee RUN_SUMMARY.txt
touch SUCCESS
