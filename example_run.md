# Example run: 11 real mixtures on Zaratan (14 October 2026)

This is the actual full-reads, no-upstream-reuse order → family → genus run.
The controller runs on a Zaratan login node only to download Kew data and submit
Slurm stages. FASTQ assembly, alignments, trees, and scoring run on compute nodes.
An explicit `unresolved_overflow` result is possible; it is not a genus call.

## Inputs and fixed versions

- Paired FASTQ.GZ: `EMseq001` through `EMseq011` under
  `/scratch/zt1/project/jmeng-prj/user/nhu/202610/realmixes/`.
- [MrKnowItAll](https://github.com/nhu92/MrKnowItAll) at `ccc26b6` for the controller
  and [SPrOUT](https://github.com/nhu92/SPrOUT) at `94e0a5f` initially;
  SPrOUT `4f16d37` fixes empty tree files for the five retried samples.
- SPrOUT's checked-in `sample_data/50targetfiles.fasta`, `gene.list.txt`, and
  `order106_refs_50genes/` for the initial 50-gene order search.
- A locally installed, corrected 353-locus backbone at
  `/home/nhu/mrknowitall_test_20261013/backbone_corrected`. Its source ZIP SHA-256 is
  `42b0c210961c13376c8e01f8a2673148cd4f6a330e5b7b2deca6ae51d3819931`.
  The ZIP is available as a
  [GitHub Release asset](https://github.com/nhu92/MrKnowItAll/releases/tag/backbone-source-2026-10-14).
  `sprout-ref install-backbone` applies the versioned 57-rule taxonomy correction file.
- Kew Tree of Life release 4.0, cached under
  `/scratch/zt1/project/jmeng-prj/user/nhu/mrknowitall_hierarchy_20261013/kew_cache`.
  Kew recoveries are assembled Angiosperms353 FASTA, not raw reads.

To recreate the backbone on another system, use the same MrKnowItAll commit
and its bundled correction TSV:

```bash
curl -fL -O https://github.com/nhu92/MrKnowItAll/releases/download/backbone-source-2026-10-14/angiosperms_353_v2_interim_targetfile_gene_alignments.zip
echo '42b0c210961c13376c8e01f8a2673148cd4f6a330e5b7b2deca6ae51d3819931  angiosperms_353_v2_interim_targetfile_gene_alignments.zip' | sha256sum -c -
sprout-ref install-backbone \
  --archive angiosperms_353_v2_interim_targetfile_gene_alignments.zip \
  --destination backbone_corrected
sprout-ref inspect-backbone backbone_corrected
```

## Validate, then submit

```bash
base=/scratch/zt1/project/jmeng-prj/user/nhu
repo="$base/mrknowitall_hierarchy_20261013/repo_v2"
sprout=/home/nhu/mrknowitall_test_20261013/SPrOUT
backbone=/home/nhu/mrknowitall_test_20261013/backbone_corrected
cache="$base/mrknowitall_hierarchy_20261013/kew_cache"
sample_root="$base/202610/realmixes"
run="$base/mrknowitall_full_20261014"
python=/home/nhu/miniconda3/envs/sprout/bin/python

# Use the repository package, not an older editable installation in the environment.
PYTHONPATH="$repo/src" "$python" -m pytest -q "$repo/tests"
"$python" "$repo/examples/hpc/zaratan_batch_to_genus.py" \
  --sample-root "$sample_root" --run-root "$run/realmix_runs" \
  --repo "$repo" --sprout "$sprout" --backbone "$backbone" --cache "$cache" \
  --max-parallel-samples 11 --max-parallel-shards 2 --poll-seconds 30 \
  --validate-only

# Start one durable controller. At submission time the root directory did not exist.
test ! -e "$run" || { echo "Run directory already exists" >&2; exit 2; }
mkdir -p "$run"
cd "$run"
nohup "$python" "$repo/examples/hpc/zaratan_batch_to_genus.py" \
  --sample-root "$sample_root" --run-root "$run/realmix_runs" \
  --repo "$repo" --sprout "$sprout" --backbone "$backbone" --cache "$cache" \
  --max-parallel-samples 11 --max-parallel-shards 2 --poll-seconds 30 \
  > "$run/controller.log" 2>&1 < /dev/null &
echo $! > "$run/CONTROLLER_PID"
```

The controller discovered exactly 11 paired samples. It submitted these order
jobs in parallel, then submitted family, branch-screen, shard, and final-genus
jobs only when each sample's preceding result was available. This is *one
controller plus many Slurm jobs*, not one giant compute-node job per sample.
The complete, stage-by-stage record is
[real11_job_manifest.tsv](examples/hpc/real11_job_manifest.tsv): 80 Slurm submissions,
75 successful stages and five preserved failed first-order attempts. The
replacement order jobs and all later stages succeeded.

| Sample | Mix | Initial order job |
| --- | --- | ---: |
| EMseq001 | brassica | 23276232 |
| EMseq002 | prunus | 23276234 |
| EMseq003 | citrus | 23276233 |
| EMseq004 | vaccinium | 23276236 |
| EMseq005 | mix1 | 23276235 |
| EMseq006 | mix2 | 23276237 |
| EMseq007 | mix3 | 23276238 |
| EMseq008 | mix5 | 23276240 |
| EMseq009 | mix6 | 23276239 |
| EMseq010 | mix7 | 23276241 |
| EMseq011 | mix8 | 23276242 |

The submitted scripts use Slurm account `jmeng-prj-aac`, partition `standard`,
32 CPUs per stage, 8-hour stage time limit, and 102400 MB local temporary space.
The controller requests 96 GB for order and small shards, and 192 GB for family,
branch, and final-genus stages (overriding the refinement script's 96 GB default).

### Effective options

| Option | Value | Meaning |
| --- | ---: | --- |
| `--max-parallel-samples` | 11 | Independent sample pipelines active at once |
| `--max-parallel-shards` | 2 | Concurrent small genus panels *per sample* |
| `--poll-seconds` | 30 | Slurm status polling interval |
| `--branch-threshold` | 120 | Begin tree-guided screening above this genus count |
| `--max-branch-rounds` | 5 | Maximum successive coarse screens |
| `--max-final-genera` | 80 | Maximum genera in final joint panel |
| `--max-postscreen-genera` | 160 | Stop unresolved above this count; no full-panel fallback |
| `--genus-shard-size` | 40 | Upper limit per small genus panel |
| `--order-memory` / `--refine-memory` / `--shard-memory` | 96G / 192G / 96G | Slurm memory requests |
| `--conda-env` | sprout | Compute-node Conda environment |

These are defaults unless explicitly shown in the command. The order stage uses
50 genes, SPrOUT's 106-reference order panel, Z ≥ 0.9, at least two supporting
loci, and at most eight order candidates. Family/genus stages use Z ≥ 0,
at most 12/20 candidates, and Kew panels with four/two representative
recoveries per taxon. Branch anchors use one representative per genus. Kew panel
construction requires at least 50 recovered loci, uses eight download workers
for prefetching, and up to 32 build threads on each compute node; SPrOUT tree
work uses four workers by default.

Other controller switches: repeat `--sample NAME` to select samples, use
`--exclude NAME` to omit one, `--submit-orders-only` to queue order work before
starting the controller, `--reuse-upstream-from OLD_RUN_ROOT` to reuse completed
order/family stages, or `--validate-only` to check inputs without submitting.
None of these selection/reuse switches was used in the full run above.

## An export failure and its restart

The first order attempts for five samples reached the report exporter, but
empty `.tre` files made the tree-rich exporter raise an exception. This was
fixed in SPrOUT `4f16d37`; no failed output was treated as a successful stage.
The failed `order/` directories were preserved as `order_failed_JOBID/` under
their respective samples. The five samples were then restarted in the **same**
run root, with these additional Slurm order jobs:

| Sample | Preserved failed job | Retry job |
| --- | ---: | ---: |
| EMseq001 | 23276232 | 23280759 |
| EMseq005 | 23276235 | 23280758 |
| EMseq007 | 23276238 | 23280757 |
| EMseq009 | 23276239 | 23280756 |
| EMseq010 | 23276241 | 23280760 |

The recovery controller used the same paths and limits, plus
`--max-parallel-samples 5` and repeated `--sample` for only these five. It wrote
`recovery.log`. Do not copy a failed `JOB_ID` into a new run: the controller
will resume that recorded job rather than submit a replacement.

```bash
# After preserving each failed order/ directory as order_failed_JOBID/:
cd "$run"
nohup "$python" "$repo/examples/hpc/zaratan_batch_to_genus.py" \
  --sample-root "$sample_root" --run-root "$run/realmix_runs" \
  --repo "$repo" --sprout "$sprout" --backbone "$backbone" --cache "$cache" \
  --max-parallel-samples 5 --max-parallel-shards 2 --poll-seconds 30 \
  --sample EMseq001 --sample EMseq005 --sample EMseq007 \
  --sample EMseq009 --sample EMseq010 \
  > "$run/recovery.log" 2>&1 < /dev/null &
echo $! > "$run/RECOVERY_PID"
```

## Inspect results and all stage jobs

```bash
squeue -u nhu
tail -n 50 "$run/controller.log"
tail -n 50 "$run/recovery.log"
find "$run/realmix_runs" -name JOB_ID -print | sort | while read -r file; do
  printf '%s\t%s\n' "${file#"$run/realmix_runs"/}" "$(cat "$file")"
done

# Each sample has RESULT.txt, RESULT.md, and RESULT.pea.json.gz if its
# order/family/branch/report stages completed. An unresolved genus is labelled.
for sample_dir in "$run"/realmix_runs/EMseq*; do
  printf '%s: ' "${sample_dir##*/}"
  grep -E '^(genus_status|genus_top_by_evidence|remaining_genera)=' \
    "$sample_dir/RESULT.txt" || true
done
```

`RESULT.md` is the readable report; `RESULT.pea.json.gz` opens in
[SPrOUT-PEA](https://nhu92.github.io/SPrOUT-PEA/). Per-stage rankings,
candidate evidence, QC, and reference FASTA remain under each sample directory.
The report contains pruned tree views, not FASTQ or full alignments.

### Final outputs in this run

All 11 samples produced `RESULT.txt`, `RESULT.md`, and a gzip-valid
`RESULT.pea.json.gz`. A completed workflow does not imply a reliable genus call:

| Sample | Mix | Genus-level outcome |
| --- | --- | --- |
| EMseq001 | brassica | Unresolved: 177 genera remain |
| EMseq002 | prunus | Prunus (evidence and SPrOUT score agree) |
| EMseq003 | citrus | Citrus (both agree) |
| EMseq004 | vaccinium | Vaccinium by evidence; Cavendishia by SPrOUT score, so inspect the disagreement |
| EMseq005 | mix1 | Unresolved: 747 genera remain |
| EMseq006 | mix2 | Unresolved: 391 genera remain |
| EMseq007 | mix3 | Zizania (both agree); Ginkgo is outside this angiosperm reference system |
| EMseq008 | mix5 | Unresolved: 257 genera remain |
| EMseq009 | mix6 | Unresolved: 302 genera remain |
| EMseq010 | mix7 | Unresolved: 339 genera remain |
| EMseq011 | mix8 | Unresolved: 315 genera remain |

`unresolved_overflow` means the tree-guided screen could not get under the
160-genus post-screen limit. These are not negative calls for the expected
genera. The workflow deliberately does not align every remaining genus at once
or force a false winner.

## Can every non-read input be obtained from GitHub?

| Requirement | Availability |
| --- | --- |
| MrKnowItAll code, Slurm scripts, taxonomy-correction TSV | Yes: this repository |
| SPrOUT code, 50-gene target file, gene list, 106-order reference panel | Yes: [SPrOUT repository](https://github.com/nhu92/SPrOUT/tree/main/sample_data) |
| Conda/package instructions | Yes: `environment.yml`, `pyproject.toml`, and SPrOUT setup instructions; exact Conda solve is not frozen |
| Kew recovery FASTA, manifest, and species tree | Public at [Kew Tree of Life](https://treeoflife.kew.org/data-access); fetched automatically by the code, **not stored in GitHub** |
| Original 353-locus interim alignment ZIP needed for this exact corrected backbone | Yes: [GitHub Release asset](https://github.com/nhu92/MrKnowItAll/releases/tag/backbone-source-2026-10-14), verified by the SHA-256 above; the 179 MB installed backbone is regenerated from the ZIP and bundled corrections |

No non-read input must be obtained privately: source code and the alignment ZIP
are on GitHub, while the program downloads Kew release 4.0 data directly from
Kew. A GitHub clone alone, **without network access to Kew**, is still not an
offline reproduction bundle. The installed backbone and Kew cache are derived
data, not tracked Git files. Do not substitute another backbone or Kew release
without recording it; reference choices may differ. For Kew data, follow its
[citation and release-number guidance](https://treeoflife.kew.org/cite-us).
