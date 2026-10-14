#!/usr/bin/env bash
set -eo pipefail

# Run this controller inside tmux on a Zaratan login node. Only lightweight
# metadata/download orchestration happens here; analysis runs via sbatch.
source /home/nhu/miniconda3/etc/profile.d/conda.sh
conda activate sprout

work=/scratch/zt1/project/jmeng-prj/user/nhu/mrknowitall_hierarchy_20261013
repo="$work/repo"
run_root="$work/realmix_runs"
old_root=/scratch/zt1/project/jmeng-prj/user/nhu/mrknowitall_test_20261013/realmix_runs

test -d "$repo/.git" || { echo "Repository not found: $repo" >&2; exit 2; }
mkdir -p "$run_root"
export PYTHONPATH="$repo/src${PYTHONPATH:+:$PYTHONPATH}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$repo"

python examples/hpc/zaratan_batch_to_genus.py \
  --sample-root /scratch/zt1/project/jmeng-prj/user/nhu/202610/realmixes \
  --run-root "$run_root" \
  --reuse-upstream-from "$old_root" \
  --repo "$repo" \
  --sprout /home/nhu/mrknowitall_test_20261013/SPrOUT \
  --backbone /home/nhu/mrknowitall_test_20261013/backbone_corrected \
  --cache "$work/kew_cache" \
  --sample EMseq001 --sample EMseq002 --sample EMseq003 \
  --sample EMseq004 --sample EMseq005 --sample EMseq006 \
  --sample EMseq007 --sample EMseq008 --sample EMseq009 \
  --sample EMseq010 --sample EMseq011 \
  --max-parallel-samples 11 \
  --max-parallel-shards 2 \
  --max-branch-rounds 5 \
  --max-postscreen-genera 160 \
  --poll-seconds 45 \
  2>&1 | tee -a "$work/controller.log"
