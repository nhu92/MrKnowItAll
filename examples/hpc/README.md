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

The full curated backbone is intentionally not stored in Git. The uploaded archive is expected
under `/lustre/scratch/nhu/202609/871ref`. First clone/update the repository, then submit the
preparation job. It runs the unit/lint tests on TTU, installs all 353 alignments, applies the
versioned taxonomy corrections and mislabel blacklist, and audits the correction hit counts:

```bash
git clone https://github.com/nhu92/MrKnowItAll.git \
  /lustre/scratch/nhu/202609/871ref/MrKnowItAll
run=/lustre/scratch/nhu/202609/871ref/prepare_backbone_run
mkdir -p "$run/logs"
cp /lustre/scratch/nhu/202609/871ref/MrKnowItAll/examples/hpc/ttu_prepare_backbone.sbatch "$run/"
cd "$run"
sbatch ttu_prepare_backbone.sbatch
```

After that job succeeds, submit family refinement from a fresh run directory:

```bash
run=/lustre/scratch/nhu/202606/SPrOUT/kew_family_refinement_20260929
mkdir -p "$run/logs"
cp /lustre/scratch/nhu/202609/871ref/MrKnowItAll/examples/hpc/ttu_kew_family_refinement.sbatch "$run/"
cd "$run"
bash -n ttu_kew_family_refinement.sbatch
sbatch ttu_kew_family_refinement.sbatch
```

The Kew reference builder runs the 353 independent locus-QC and MAFFT jobs with one process per
allocated CPU (`--threads 32` by default on this script). Each MAFFT process uses one thread, so
the job stays within its Slurm allocation. Override with `SPROUT_REF_BUILD_WORKERS`; progress is
printed after every ten completed loci.

`SPROUT_PROJECT` must equal the project token embedded in the reused query exon headers. The
default is `sprout_test`, matching `zaratan_sprout_testrun.sbatch`; the refinement script verifies
that token occurs in the resulting tree leaves before distance aggregation so a mismatch cannot
silently produce an all-zero ranking.

After family prediction, `ttu_kew_genus_refinement.sbatch` reads the selected families and builds
a balanced Kew genus panel. It requests 32 CPUs and 64 GiB, uses all 32 CPUs for independent
reference loci, and keeps the downstream SPrOUT tree pool at 4 × 4 threads to avoid nested
oversubscription. The defaults expect the project-token repair results produced during the first
TTU test; set `SPROUT_FAMILY_CANDIDATES` when using another family candidate file.

`ttu_kew_species_refinement.sbatch` performs the final genus-to-species step. Species are grouped
by the full `Genus_species` binomial (never by epithet alone), one QC-passing Kew recovery is
selected per species, and SPrOUT emits both the complete species ranking and the positive-z-score
candidate list.
