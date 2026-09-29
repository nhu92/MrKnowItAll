#!/bin/bash
#SBATCH --job-name=kew_species
#SBATCH --partition=nocona
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -Eeuo pipefail
trap 'code=$?; echo "FAILED exit=${code} line=${LINENO} $(date -Is)"; exit "$code"' ERR

SPROUT_SOURCE=${SPROUT_SOURCE:-/lustre/scratch/nhu/202606/SPrOUT}
SPROUT_TEST_RUN=${SPROUT_TEST_RUN:-$SPROUT_SOURCE/testrun_20260929}
GENUS_RUN=${SPROUT_GENUS_RUN:-$SPROUT_SOURCE/kew_genus_refinement_20260929}
GENUS_CANDIDATES=${SPROUT_GENUS_CANDIDATES:-$GENUS_RUN/sprout_test.genus_candidates.txt}
MRKNOWITALL_SOURCE=${MRKNOWITALL_SOURCE:-/lustre/scratch/nhu/202609/871ref/MrKnowItAll}
MRKNOWITALL_BACKBONE=${MRKNOWITALL_BACKBONE:-/lustre/scratch/nhu/202609/871ref/backbone_corrected}
CONDA_SH=${SPROUT_CONDA_SH:-$HOME/miniconda3/etc/profile.d/conda.sh}
CONDA_ENV=${SPROUT_CONDA_ENV:-sprout}
PROJECT=${SPROUT_PROJECT:-sprout_test}
TREE_WORKERS=${SPROUT_TREE_WORKERS:-4}
BUILD_WORKERS=${SPROUT_REF_BUILD_WORKERS:-${SLURM_CPUS_PER_TASK:-32}}
RUN=${SPROUT_REFINEMENT_RUN:-$PWD}

cd "$RUN"
if [[ -e SUCCESS || -e species_reference || -e species_phylo || -e species_matrix ]]; then
  echo "Refusing to reuse a species directory containing outputs: $RUN" >&2
  exit 4
fi
for path in \
  "$SPROUT_SOURCE/02_exon_trees.py" \
  "$SPROUT_SOURCE/03_distance_matrices.py" \
  "$SPROUT_SOURCE/04_prediction.py" \
  "$SPROUT_TEST_RUN/exon" \
  "$SPROUT_TEST_RUN/gene_list.txt" \
  "$GENUS_CANDIDATES" \
  "$MRKNOWITALL_SOURCE/pyproject.toml"; do
  [[ -e "$path" ]] || { echo "Missing required input: $path" >&2; exit 2; }
done
backbone_count=$(find "$MRKNOWITALL_BACKBONE" -maxdepth 1 -type f -name '*.fasta' | wc -l)
[[ "$backbone_count" -eq 353 ]] || {
  echo "Expected 353 backbone alignments; found $backbone_count" >&2
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
  echo "genus_candidates=$GENUS_CANDIDATES"
  echo "genera=$(paste -sd, "$GENUS_CANDIDATES")"
  echo "reference_build_workers=$BUILD_WORKERS"
  echo "tree_workers=$TREE_WORKERS"
  git -C "$SPROUT_SOURCE" rev-parse HEAD | sed 's/^/sprout_commit=/'
  git -C "$MRKNOWITALL_SOURCE" rev-parse HEAD | sed 's/^/mrknowitall_commit=/'
} | tee RUN_PROVENANCE.txt

sprout-ref build-kew-panel \
  --level species \
  --parent-file "$GENUS_CANDIDATES" \
  --backbone "$MRKNOWITALL_BACKBONE" \
  --output species_reference \
  --representatives 1 \
  --pool-multiplier 2 \
  --minimum-recovered-loci 50 \
  --download-workers 8 \
  --threads "$BUILD_WORKERS"

grep -F -x -f "$SPROUT_TEST_RUN/gene_list.txt" \
  species_reference/gene.list.txt > gene_list.txt
gene_count=$(wc -l < gene_list.txt)
[[ "$gene_count" -gt 0 ]] || { echo "No shared query/reference genes" >&2; exit 3; }

python "$SPROUT_SOURCE/02_exon_trees.py" \
  -t "$TREE_WORKERS" -p "$PROJECT" -g gene_list.txt \
  -e "$SPROUT_TEST_RUN/exon" \
  -r species_reference/ref \
  --output_phylo species_phylo
if ! grep -R -F -q --include='*.tre' "$PROJECT" species_phylo; then
  echo "No query leaf contains project token '$PROJECT'; refusing all-zero scores" >&2
  exit 5
fi

python "$SPROUT_SOURCE/03_distance_matrices.py" \
  -t "$TREE_WORKERS" -p "$PROJECT" -g gene_list.txt --threshold 1 \
  --input_phylo species_phylo --output_tree species_matrix
cp "$PROJECT.cumulative_dist.csv" "$PROJECT.species.cumulative_dist.csv"
python "$SPROUT_SOURCE/04_prediction.py" \
  -i "$PROJECT.species.cumulative_dist.csv" \
  -o "$PROJECT.predictions_species.csv" \
  -tl s -z 0 \
  -to "$PROJECT.species_candidates.txt"

{
  echo "finished=$(date -Is)"
  echo "shared_genes=$gene_count"
  echo "reference_loci=$(find species_reference/ref -type f -name '*.fasta' | wc -l)"
  echo "species_trees=$(find species_phylo -type f -name '*.tre' | wc -l)"
  echo "species_ranking:"
  cat "$PROJECT.predictions_species.csv"
  echo "species_candidates:"
  cat "$PROJECT.species_candidates.txt"
} | tee RUN_SUMMARY.txt
touch SUCCESS
