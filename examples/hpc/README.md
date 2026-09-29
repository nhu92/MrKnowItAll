# UMD Zaratan SPrOUT smoke test

The Slurm script runs the public SPrOUT 50-gene order-to-family example in a new directory.
It samples paired reads deterministically and never overwrites an existing pipeline output.

From an existing SPrOUT checkout containing the example FASTQs:

```bash
run="$PWD/testrun_20260929"
mkdir -p "$run/logs"
curl -fsSL \
  https://raw.githubusercontent.com/nhu92/MrKnowItAll/main/examples/hpc/zaratan_sprout_testrun.sbatch \
  -o "$run/sprout_testrun.sbatch"
cd "$run"
bash -n sprout_testrun.sbatch
sbatch sprout_testrun.sbatch
```

Defaults: Zaratan `standard`, one task, 16 CPU, 64 GiB, six hours, 4% of the source read pairs,
and four concurrent downstream workers. `02_exon_trees.py` gives each MAFFT process the same
thread count as the Python worker pool, so 4 × 4 prevents the 64 × 64 oversubscription implied
by the original example.

Optional submission-time overrides include `SPROUT_SAMPLE_FRACTION`, `SPROUT_SAMPLE_SEED`,
`SPROUT_SOURCE`, `SPROUT_CONDA_ENV`, and `SPROUT_TREE_WORKERS`.
