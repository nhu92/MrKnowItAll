#!/bin/bash
#SBATCH --job-name=kew_family
#SBATCH --partition=nocona
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -Eeuo pipefail
trap 'code=$?; echo "FAILED exit=${code} line=${LINENO} $(date -Is)"; exit "$code"' ERR

SPROUT_SOURCE=${SPROUT_SOURCE:-/lustre/scratch/nhu/202606/SPrOUT}
SPROUT_TEST_RUN=${SPROUT_TEST_RUN:-$SPROUT_SOURCE/testrun_20260929}
MRKNOWITALL_SOURCE=${MRKNOWITALL_SOURCE:-$SPROUT_SOURCE/MrKnowItAll}
MRKNOWITALL_BACKBONE=${MRKNOWITALL_BACKBONE:-$MRKNOWITALL_SOURCE/data/backbone}
CONDA_SH=${SPROUT_CONDA_SH:-$HOME/miniconda3/etc/profile.d/conda.sh}
CONDA_ENV=${SPROUT_CONDA_ENV:-sprout}
PROJECT=${SPROUT_PROJECT:-sprout_test_kew_family}
THREADS=${SPROUT_TREE_WORKERS:-4}
RUN=${SPROUT_REFINEMENT_RUN:-$PWD}

cd "$RUN"
if [[ -e SUCCESS || -e family_reference || -e family_phylo || -e family_matrix ]]; then
  echo "Refusing to reuse a refinement directory containing outputs: $RUN" >&2
  exit 4
fi
for path in \
  "$SPROUT_SOURCE/02_exon_trees.py" \
  "$SPROUT_SOURCE/03_distance_matrices.py" \
  "$SPROUT_SOURCE/04_prediction.py" \
  "$SPROUT_TEST_RUN/exon" \
  "$SPROUT_TEST_RUN/gene_list.txt" \
  "$SPROUT_TEST_RUN/sprout_test.order_candidates.txt" \
  "$MRKNOWITALL_SOURCE/pyproject.toml"; do
  [[ -e "$path" ]] || { echo "Missing required input: $path" >&2; exit 2; }
done
backbone_count=$(find "$MRKNOWITALL_BACKBONE" -maxdepth 1 -type f -name '*.fasta' | wc -l)
[[ "$backbone_count" -eq 353 ]] || {
  echo "Expected 353 installed backbone alignments at $MRKNOWITALL_BACKBONE; found $backbone_count" >&2
  exit 2
}

source "$CONDA_SH"
conda activate "$CONDA_ENV"
python -m pip install -e "$MRKNOWITALL_SOURCE"
for tool in python mafft fasttree trimal; do
  command -v "$tool" || { echo "Missing required tool: $tool" >&2; exit 127; }
done

{
  echo "started=$(date -Is)"
  echo "host=$(hostname -f)"
  echo "job_id=${SLURM_JOB_ID:-none}"
  echo "sprout_source=$SPROUT_SOURCE"
  echo "sprout_test_run=$SPROUT_TEST_RUN"
  echo "mrknowitall_source=$MRKNOWITALL_SOURCE"
  echo "backbone=$MRKNOWITALL_BACKBONE"
  git -C "$SPROUT_SOURCE" rev-parse HEAD | sed 's/^/sprout_commit=/'
  git -C "$MRKNOWITALL_SOURCE" rev-parse HEAD | sed 's/^/mrknowitall_commit=/'
} | tee RUN_PROVENANCE.txt

sprout-ref build-kew-panel \
  --level family \
  --parent-file "$SPROUT_TEST_RUN/sprout_test.order_candidates.txt" \
  --backbone "$MRKNOWITALL_BACKBONE" \
  --output family_reference \
  --representatives 4 \
  --pool-multiplier 2 \
  --minimum-recovered-loci 50 \
  --download-workers 4

cp "$SPROUT_TEST_RUN/gene_list.txt" gene_list.txt
python "$SPROUT_SOURCE/02_exon_trees.py" \
  -t "$THREADS" -p "$PROJECT" -g gene_list.txt \
  -e "$SPROUT_TEST_RUN/exon" \
  -r family_reference/ref \
  --output_phylo family_phylo
python "$SPROUT_SOURCE/03_distance_matrices.py" \
  -t "$THREADS" -p "$PROJECT" -g gene_list.txt --threshold 1 \
  --input_phylo family_phylo --output_tree family_matrix
cp "$PROJECT.cumulative_dist.csv" "$PROJECT.family.cumulative_dist.csv"
python "$SPROUT_SOURCE/04_prediction.py" \
  -i "$PROJECT.family.cumulative_dist.csv" \
  -o "$PROJECT.predictions_family.csv" \
  -tl f -z 0 \
  -to "$PROJECT.family_candidates.txt"

{
  echo "finished=$(date -Is)"
  echo "reference_loci=$(find family_reference/ref -type f -name '*.fasta' | wc -l)"
  echo "family_trees=$(find family_phylo -type f -name '*.tre' | wc -l)"
  echo "family_ranking:"
  cat "$PROJECT.predictions_family.csv"
} | tee RUN_SUMMARY.txt
touch SUCCESS
