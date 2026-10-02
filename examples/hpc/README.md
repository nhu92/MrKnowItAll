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
- `ttu_kew_family_refinement.sbatch`: build family references and repeat SPrOUT.
- `ttu_kew_genus_refinement.sbatch`: build genus references and repeat SPrOUT.
- `ttu_kew_species_refinement.sbatch`: optional species-level experiment.

The separate refinement scripts are useful for debugging. For a normal sample, use
`ttu_mix_to_genus.sbatch`.
