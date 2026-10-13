#!/bin/bash
#SBATCH --job-name=pea_check
#SBATCH --partition=standard
#SBATCH --account=jmeng-prj-aac
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --tmp=20480
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -Eeuo pipefail
SPROUT_SOURCE=${SPROUT_SOURCE:?Set SPROUT_SOURCE}
QUERY_EXON=${QUERY_EXON:?Set QUERY_EXON}
QUERY_GENE_LIST=${QUERY_GENE_LIST:?Set QUERY_GENE_LIST}
REPORT_DIR=${REPORT_DIR:?Set REPORT_DIR to a new directory}
PROJECT=${SPROUT_PROJECT:-emseq003}
[[ ! -e "$REPORT_DIR" ]] || { echo "Refusing to overwrite $REPORT_DIR" >&2; exit 4; }
for path in "$SPROUT_SOURCE/export_pea.py" "$QUERY_EXON" "$QUERY_GENE_LIST"; do
  [[ -e "$path" ]] || { echo "Missing $path" >&2; exit 2; }
done
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate sprout
work=${SLURM_TMPDIR:-/tmp}/pea-check-${SLURM_JOB_ID:-manual}
mkdir -p "$work"
cd "$work"
python "$SPROUT_SOURCE/02_exon_trees.py" -t 2 -p "$PROJECT" \
  -g "$QUERY_GENE_LIST" -e "$QUERY_EXON" \
  -r "$SPROUT_SOURCE/sample_data/order106_refs_50genes" --output_phylo phylo
python "$SPROUT_SOURCE/03_distance_matrices.py" -t 4 -p "$PROJECT" \
  -g "$QUERY_GENE_LIST" --threshold 1 --input_phylo phylo --output_tree matrix
python "$SPROUT_SOURCE/04_prediction.py" -i "$PROJECT.cumulative_dist.csv" \
  -o predictions.csv -tl o -z 0.9 -to candidates.txt
python "$SPROUT_SOURCE/export_pea.py" stage --project "$PROJECT" --level order \
  --predictions predictions.csv --candidates candidates.txt \
  --gene-list "$QUERY_GENE_LIST" --matrix-dir matrix --z-threshold 0.9 \
  --output "$PROJECT.pea.json.gz" --markdown "$PROJECT.report.md"
python - "$PROJECT.pea.json.gz" <<'PY'
import gzip, json, sys
with gzip.open(sys.argv[1], 'rt', encoding='utf-8') as handle:
    report = json.load(handle)
stage = report['stages'][0]
assert stage['evidence_available'] and stage['genes'], 'Missing real per-gene evidence'
assert stage['counts']['query_exons'] > 0, 'No matched query exons'
print('INTEGRATION_OK', report['sample_id'], len(stage['genes']), stage['counts'])
PY
mkdir "$REPORT_DIR"
cp "$PROJECT.pea.json.gz" "$PROJECT.report.md" "$REPORT_DIR/"
touch "$REPORT_DIR/SUCCESS"
