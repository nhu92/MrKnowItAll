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

## TTU Nocona: Kew family refinement

`ttu_kew_family_refinement.sbatch` reuses the completed smoke test's extracted exons and order
call, generates a count-balanced Kew family panel with MrKnowItAll, then reruns SPrOUT tree,
distance, and family prediction stages. It never modifies the original smoke-test directory.

The full curated backbone is intentionally not stored in Git. Before submission, install the
provided alignment ZIP in the TTU checkout so that `data/backbone` contains 353 FASTA files:

```bash
git clone https://github.com/nhu92/MrKnowItAll.git \
  /lustre/scratch/nhu/202606/SPrOUT/MrKnowItAll
cd /lustre/scratch/nhu/202606/SPrOUT/MrKnowItAll
conda activate sprout
pip install -e .
sprout-ref install-backbone \
  --archive /path/on/ttu/angiosperms_353_v2_interim_targetfile_gene_alignments.zip \
  --destination data/backbone --overwrite
```

Then submit from a fresh run directory:

```bash
run=/lustre/scratch/nhu/202606/SPrOUT/kew_family_refinement_20260929
mkdir -p "$run/logs"
cp /lustre/scratch/nhu/202606/SPrOUT/MrKnowItAll/examples/hpc/ttu_kew_family_refinement.sbatch "$run/"
cd "$run"
bash -n ttu_kew_family_refinement.sbatch
sbatch ttu_kew_family_refinement.sbatch
```
