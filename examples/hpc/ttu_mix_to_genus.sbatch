#!/bin/bash
#SBATCH --job-name=mix_to_genus
#SBATCH --partition=nocona
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -Eeuo pipefail
trap 'code=$?; echo "FAILED exit=${code} line=${LINENO} $(date -Is)"; exit "$code"' ERR

SPROUT_SOURCE=${SPROUT_SOURCE:-/lustre/scratch/nhu/202606/SPrOUT}
MRKNOWITALL_SOURCE=${MRKNOWITALL_SOURCE:-/lustre/scratch/nhu/202609/871ref/MrKnowItAll}
MRKNOWITALL_BACKBONE=${MRKNOWITALL_BACKBONE:-/lustre/scratch/nhu/202609/871ref/backbone_corrected}
CONDA_SH=${SPROUT_CONDA_SH:-$HOME/miniconda3/etc/profile.d/conda.sh}
CONDA_ENV=${SPROUT_CONDA_ENV:-sprout}
READ1=${SPROUT_READ1:?Set SPROUT_READ1 to the R1 FASTQ/FASTQ.GZ}
READ2=${SPROUT_READ2:?Set SPROUT_READ2 to the R2 FASTQ/FASTQ.GZ}
PROJECT=${SPROUT_PROJECT:-mix_sample}
TARGET_FILE=${SPROUT_TARGET_FILE:-$SPROUT_SOURCE/sample_data/50targetfiles.fasta}
SOURCE_GENE_LIST=${SPROUT_GENE_LIST:-$SPROUT_SOURCE/sample_data/gene.list.txt}
ORDER_REFERENCE=${SPROUT_ORDER_REFERENCE:-$SPROUT_SOURCE/sample_data/order106_refs_50genes}
SAMPLE_FRACTION=${SPROUT_SAMPLE_FRACTION:-1.0}
SAMPLE_SEED=${SPROUT_SAMPLE_SEED:-20260929}
ASSEMBLY_THREADS=${SLURM_CPUS_PER_TASK:-32}
TREE_WORKERS=${SPROUT_TREE_WORKERS:-4}
BUILD_WORKERS=${SPROUT_REF_BUILD_WORKERS:-${SLURM_CPUS_PER_TASK:-32}}
MIN_LOCUS_VOTES=${SPROUT_MIN_LOCUS_VOTES:-2}
MAX_ORDER_CANDIDATES=${SPROUT_MAX_ORDER_CANDIDATES:-8}
MAX_FAMILY_CANDIDATES=${SPROUT_MAX_FAMILY_CANDIDATES:-12}
MAX_GENUS_CANDIDATES=${SPROUT_MAX_GENUS_CANDIDATES:-20}
RUN=${SPROUT_RUN_DIR:-$PWD}

[[ "$PROJECT" =~ ^[A-Za-z0-9_.-]+$ ]] || {
  echo "SPROUT_PROJECT may contain only letters, digits, dot, underscore and hyphen" >&2
  exit 2
}
for path in \
  "$READ1" "$READ2" "$TARGET_FILE" "$SOURCE_GENE_LIST" "$ORDER_REFERENCE" \
  "$SPROUT_SOURCE/01_exons_assembly.py" \
  "$SPROUT_SOURCE/02_exon_trees.py" \
  "$SPROUT_SOURCE/03_distance_matrices.py" \
  "$SPROUT_SOURCE/04_prediction.py" \
  "$MRKNOWITALL_SOURCE/pyproject.toml"; do
  [[ -e "$path" ]] || { echo "Missing required input: $path" >&2; exit 2; }
done
backbone_count=$(find "$MRKNOWITALL_BACKBONE" -maxdepth 1 -type f -name '*.fasta' | wc -l)
[[ "$backbone_count" -eq 353 ]] || {
  echo "Expected 353 backbone alignments; found $backbone_count" >&2
  exit 2
}

cd "$RUN"
if [[ -e SUCCESS || -e 00_query || -e 01_order || -e 02_family || -e 03_genus ]]; then
  echo "Refusing to reuse a directory containing pipeline outputs: $RUN" >&2
  exit 4
fi
mkdir -p 00_query 01_order 02_family 03_genus

source "$CONDA_SH"
conda activate "$CONDA_ENV"
python -m pip install -e "$MRKNOWITALL_SOURCE"
for tool in python fastp hybpiper mafft fasttree trimal seqkit bwa sprout-ref; do
  command -v "$tool" || { echo "Missing required tool: $tool" >&2; exit 127; }
done

{
  echo "started=$(date -Is)"
  echo "host=$(hostname -f)"
  echo "job_id=${SLURM_JOB_ID:-none}"
  echo "project=$PROJECT"
  echo "read1=$READ1"
  echo "read2=$READ2"
  echo "sample_fraction=$SAMPLE_FRACTION"
  echo "target_file=$TARGET_FILE"
  echo "source_gene_list=$SOURCE_GENE_LIST"
  echo "order_reference=$ORDER_REFERENCE"
  echo "backbone=$MRKNOWITALL_BACKBONE"
  echo "minimum_locus_votes=$MIN_LOCUS_VOTES"
  git -C "$SPROUT_SOURCE" rev-parse HEAD 2>/dev/null | sed 's/^/sprout_commit=/' || true
  git -C "$MRKNOWITALL_SOURCE" rev-parse HEAD | sed 's/^/mrknowitall_commit=/'
} | tee RUN_PROVENANCE.txt

query_dir="$RUN/00_query"
if [[ "$SAMPLE_FRACTION" == "1" || "$SAMPLE_FRACTION" == "1.0" ]]; then
  input_r1="$READ1"
  input_r2="$READ2"
else
  input_r1="$query_dir/$PROJECT.R1.fastq.gz"
  input_r2="$query_dir/$PROJECT.R2.fastq.gz"
  seqkit sample -p "$SAMPLE_FRACTION" -s "$SAMPLE_SEED" "$READ1" -o "$input_r1"
  seqkit sample -p "$SAMPLE_FRACTION" -s "$SAMPLE_SEED" "$READ2" -o "$input_r2"
fi
seqkit stats -a "$input_r1" "$input_r2" | tee "$query_dir/input_reads.stats.txt"
cp "$SOURCE_GENE_LIST" "$query_dir/gene_list.txt"

cd "$query_dir"
python "$SPROUT_SOURCE/01_exons_assembly.py" \
  -t "$ASSEMBLY_THREADS" \
  -r1 "$input_r1" -r2 "$input_r2" \
  -p "$PROJECT" \
  -g gene_list.txt -m "$TARGET_FILE" \
  --output_hyb hyb --output_exon exon

order_dir="$RUN/01_order"
cd "$order_dir"
python "$SPROUT_SOURCE/02_exon_trees.py" \
  -t "$TREE_WORKERS" -p "$PROJECT" -g "$query_dir/gene_list.txt" \
  -e "$query_dir/exon" -r "$ORDER_REFERENCE" --output_phylo order_phylo
python "$SPROUT_SOURCE/03_distance_matrices.py" \
  -t "$TREE_WORKERS" -p "$PROJECT" -g "$query_dir/gene_list.txt" --threshold 1 \
  --input_phylo order_phylo --output_tree order_matrix
cp "$PROJECT.cumulative_dist.csv" "$PROJECT.order.cumulative_dist.csv"
python "$SPROUT_SOURCE/04_prediction.py" \
  -i "$PROJECT.order.cumulative_dist.csv" -o order_predictions.csv \
  -tl o -z 0.9 -to order_z_candidates.txt
sprout-ref select-candidates \
  --matrix-dir order_matrix --predictions order_predictions.csv \
  --level order --project "$PROJECT" \
  --minimum-locus-votes "$MIN_LOCUS_VOTES" --minimum-z 0.9 \
  --always-top 1 --max-candidates "$MAX_ORDER_CANDIDATES" \
  --output order_candidates.txt --report order_candidate_evidence.csv
[[ -s order_candidates.txt ]] || { echo "No order candidate survived" >&2; exit 6; }

family_dir="$RUN/02_family"
cd "$family_dir"
sprout-ref build-kew-panel \
  --level family --parent-file "$order_dir/order_candidates.txt" \
  --backbone "$MRKNOWITALL_BACKBONE" --output family_reference \
  --representatives 4 --pool-multiplier 2 --minimum-recovered-loci 50 \
  --download-workers 8 --threads "$BUILD_WORKERS"
grep -F -x -f "$query_dir/gene_list.txt" \
  family_reference/gene.list.txt > gene_list.txt
[[ -s gene_list.txt ]] || { echo "No shared family genes" >&2; exit 7; }
python "$SPROUT_SOURCE/02_exon_trees.py" \
  -t "$TREE_WORKERS" -p "$PROJECT" -g gene_list.txt \
  -e "$query_dir/exon" -r family_reference/ref --output_phylo family_phylo
python "$SPROUT_SOURCE/03_distance_matrices.py" \
  -t "$TREE_WORKERS" -p "$PROJECT" -g gene_list.txt --threshold 1 \
  --input_phylo family_phylo --output_tree family_matrix
cp "$PROJECT.cumulative_dist.csv" "$PROJECT.family.cumulative_dist.csv"
python "$SPROUT_SOURCE/04_prediction.py" \
  -i "$PROJECT.family.cumulative_dist.csv" -o family_predictions.csv \
  -tl f -z 0 -to family_z_candidates.txt
sprout-ref select-candidates \
  --matrix-dir family_matrix --predictions family_predictions.csv \
  --level family --project "$PROJECT" \
  --minimum-locus-votes "$MIN_LOCUS_VOTES" --minimum-z 0 \
  --always-top 1 --max-candidates "$MAX_FAMILY_CANDIDATES" \
  --output family_candidates.txt --report family_candidate_evidence.csv
[[ -s family_candidates.txt ]] || { echo "No family candidate survived" >&2; exit 8; }

genus_dir="$RUN/03_genus"
cd "$genus_dir"
sprout-ref build-kew-panel \
  --level genus --parent-file "$family_dir/family_candidates.txt" \
  --backbone "$MRKNOWITALL_BACKBONE" --output genus_reference \
  --representatives 2 --pool-multiplier 2 --minimum-recovered-loci 50 \
  --download-workers 8 --threads "$BUILD_WORKERS"
grep -F -x -f "$query_dir/gene_list.txt" \
  genus_reference/gene.list.txt > gene_list.txt
[[ -s gene_list.txt ]] || { echo "No shared genus genes" >&2; exit 9; }
python "$SPROUT_SOURCE/02_exon_trees.py" \
  -t "$TREE_WORKERS" -p "$PROJECT" -g gene_list.txt \
  -e "$query_dir/exon" -r genus_reference/ref --output_phylo genus_phylo
python "$SPROUT_SOURCE/03_distance_matrices.py" \
  -t "$TREE_WORKERS" -p "$PROJECT" -g gene_list.txt --threshold 1 \
  --input_phylo genus_phylo --output_tree genus_matrix
cp "$PROJECT.cumulative_dist.csv" "$PROJECT.genus.cumulative_dist.csv"
python "$SPROUT_SOURCE/04_prediction.py" \
  -i "$PROJECT.genus.cumulative_dist.csv" -o genus_predictions.csv \
  -tl g -z 0 -to genus_z_candidates.txt
sprout-ref select-candidates \
  --matrix-dir genus_matrix --predictions genus_predictions.csv \
  --level genus --project "$PROJECT" \
  --minimum-locus-votes "$MIN_LOCUS_VOTES" --minimum-z 0 \
  --always-top 1 --max-candidates "$MAX_GENUS_CANDIDATES" \
  --output genus_candidates.txt --report genus_candidate_evidence.csv

cp genus_candidates.txt "$RUN/FINAL_GENUS_CANDIDATES.txt"
cp genus_candidate_evidence.csv "$RUN/FINAL_GENUS_EVIDENCE.csv"
cp genus_predictions.csv "$RUN/FINAL_GENUS_RANKING.csv"
{
  echo "finished=$(date -Is)"
  echo "project=$PROJECT"
  echo "order_candidates:"
  cat "$order_dir/order_candidates.txt"
  echo "family_candidates:"
  cat "$family_dir/family_candidates.txt"
  echo "genus_candidates:"
  cat genus_candidates.txt
} | tee "$RUN/RUN_SUMMARY.txt"
touch "$RUN/SUCCESS"
