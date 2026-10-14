# HPC examples

These Slurm scripts run MrKnowItAll with SPrOUT on a computing cluster.

## Mixed reads to genus

`ttu_mix_to_genus.sbatch` is the main script. It accepts paired reads and runs the complete
order → family → genus workflow.

```bash
run=/path/to/new/run
mkdir -p "$run/logs"
cp ttu_mix_to_genus.sbatch "$run/"
cd "$run"

sbatch --export=ALL,\
SPROUT_READ1=/path/to/sample_R1.fastq.gz,\
SPROUT_READ2=/path/to/sample_R2.fastq.gz,\
SPROUT_PROJECT=my_sample \
ttu_mix_to_genus.sbatch
```

The main results are:

- `FINAL_GENUS_CANDIDATES.txt`
- `FINAL_GENUS_RANKING.csv`
- `FINAL_GENUS_EVIDENCE.csv`
- `RUN_SUMMARY.txt`
- `RESULT.md` (readable ranking, Z values, and evidence counts)
- `RESULT.pea.json.gz` (small file for [SPrOUT-PEA](https://nhu92.github.io/SPrOUT-PEA/))

Set `SPROUT_SAMPLE_FRACTION=0.04` for a small test or `1.0` for all reads.

## Prepare the 353-gene backbone

The alignment ZIP is too large for GitHub. On TTU, copy `ttu_prepare_backbone.sbatch` into a
new run directory and submit it:

```bash
mkdir -p prepare_backbone_run/logs
cp ttu_prepare_backbone.sbatch prepare_backbone_run/
cd prepare_backbone_run
sbatch ttu_prepare_backbone.sbatch
```

## Other scripts

- `zaratan_sprout_testrun.sbatch`: small SPrOUT test on UMD Zaratan.
- `zaratan_batch_to_genus.py`: one login-node command for multiple real samples on Zaratan.
- `zaratan_query_order.sbatch`: query assembly and order prediction on a compute node.
- `zaratan_refine_stage.sbatch`: offline family, branch-anchor, or genus refinement.
- `ttu_kew_family_refinement.sbatch`: build family references and repeat SPrOUT.
- `ttu_kew_genus_refinement.sbatch`: build genus references and repeat SPrOUT.
- `ttu_kew_species_refinement.sbatch`: optional species-level experiment.

The separate refinement scripts are useful for debugging. For a normal sample, use
`ttu_mix_to_genus.sbatch`.

## Zaratan: several real samples

Zaratan compute nodes cannot reach Kew. The batch controller downloads Kew recoveries on
the login node and submits each analysis stage to Slurm. It can run several samples at once
and resumes stages that already succeeded or have a recorded Slurm job ID.

To put every independent order analysis in the Slurm queue immediately, add
`--submit-orders-only` and select samples with repeated `--sample` options (or
omit them to select every sample). Run the normal controller afterward; it will
reuse those order job IDs. Family and genus jobs cannot be queued until the
preceding stage identifies candidate taxa and the Kew sequences have been cached.
`--max-parallel-samples` limits active sample pipelines, not pre-submitted order jobs.
Only one full controller should use a run directory at a time; newer controllers
enforce this with a shared lock across Zaratan login nodes.

```bash
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate sprout
python examples/hpc/zaratan_batch_to_genus.py \
  --sample-root /scratch/zt1/project/jmeng-prj/user/nhu/202610/realmixes \
  --run-root /scratch/zt1/project/jmeng-prj/user/nhu/mrknowitall_test_20261013/realmix_runs \
  --repo /home/nhu/mrknowitall_test_20261013/MrKnowItAll \
  --sprout /home/nhu/mrknowitall_test_20261013/SPrOUT \
  --backbone /home/nhu/mrknowitall_test_20261013/backbone_corrected \
  --cache /scratch/zt1/project/jmeng-prj/user/nhu/mrknowitall_test_20261013/kew_cache \
  --max-parallel-samples 3 --exclude EMseq001
```

Each sample gets `order`, `family`, and `genus` results plus `RESULT.txt`, `RESULT.md`,
and `RESULT.pea.json.gz`. The summary
separately labels the leading genus by per-locus evidence and by SPrOUT's aggregate score;
they can disagree and neither is a calibrated presence probability. If more than 80 genera
are eligible, the controller screens the Kew species tree in successive rounds. Each
candidate family is divided into a few broad branches; a small anchor panel tests those
branches, and only supported branches and close neighbors enter the next round. The
default limit is five rounds (`--max-branch-rounds`). This is a coarse-to-fine search,
not a partition of every genus into jobs. The old full-family fallback is disabled.
If 81–160 genera remain, small Slurm panels shortlist them for a final joint panel of
at most 80 genera. Above 160, the sample is explicitly reported as unresolved instead
of launching a huge alignment. Configure the limits with `--genus-shard-size`,
`--max-parallel-shards`, `--max-final-genera`, and `--max-postscreen-genera`.
Anchor screens appear under `branch/round_XX/`.
Use a new `--run-root` when upgrading a sample that already has a legacy
`branch/plan.json`; the controller refuses to mix old and new genus screens.
