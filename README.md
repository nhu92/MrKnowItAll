# MrKnowItAll

MrKnowItAll is a small bioinformatics extension for
[SPrOUT](https://github.com/nhu92/SPrOUT). It helps identify the plant genera in a mixed
paired-end sequencing sample.

It runs SPrOUT first, downloads suitable Angiosperms353 sequences from the
[Kew Tree of Life](https://treeoflife.kew.org/), builds smaller reference sets, and repeats the
search from order to family to genus.

## What it does

```text
paired FASTQ reads
    → SPrOUT exon assembly
    → order search
    → Kew family references
    → family search
    → Kew genus references
    → genus result
```

- Accepts paired `FASTQ` or `FASTQ.GZ` files.
- Uses public Kew Angiosperms353 assembled sequences.
- Checks new sequences against a curated 353-gene backbone.
- Keeps references balanced between taxonomic groups.
- Reports both strong and low-abundance genus candidates.
- Runs without manual sequence editing.

MrKnowItAll does not assemble Kew raw reads. Kew already provides assembled Angiosperms353
recovery FASTA files. The program checks, filters, aligns, and formats those sequences for
SPrOUT. Only the user's sample reads are assembled with HybPiper.

## Install

```bash
conda env create -f environment.yml
conda activate sprout-refbuilder
pip install -e .
```

MAFFT is recommended. The mixed-read workflow also needs the SPrOUT environment with HybPiper,
fastp, BWA, trimAl, FastTree, and SeqKit.

## Prepare the backbone

The large alignment backbone is not stored in GitHub. Install it from the original ZIP:

```bash
sprout-ref install-backbone \
  --archive /path/to/angiosperms_353_gene_alignments.zip \
  --destination data/backbone
```

The installer checks all 353 alignments and applies the included taxonomy corrections.

## Run a mixed sample to genus

The easiest HPC entry point is:

```text
examples/hpc/ttu_mix_to_genus.sbatch
```

Example submission on TTU Nocona:

```bash
run=/path/to/new/run
mkdir -p "$run/logs"
cp examples/hpc/ttu_mix_to_genus.sbatch "$run/"
cd "$run"

sbatch --export=ALL,\
SPROUT_READ1=/path/to/sample_R1.fastq.gz,\
SPROUT_READ2=/path/to/sample_R2.fastq.gz,\
SPROUT_PROJECT=my_sample \
ttu_mix_to_genus.sbatch
```

Use `SPROUT_SAMPLE_FRACTION=0.04` for a quick test or `SPROUT_SAMPLE_FRACTION=1.0` for all reads.

## Main results

```text
FINAL_GENUS_CANDIDATES.txt   short genus list
FINAL_GENUS_RANKING.csv      complete SPrOUT genus ranking
FINAL_GENUS_EVIDENCE.csv     locus votes and scores
RUN_SUMMARY.txt              simple run summary
```

Intermediate order, family, reference, tree, and QC files remain in the run directory.

## Build a reference directly from Kew

```bash
sprout-ref sync
sprout-ref search "Abatia rugosa"

sprout-ref build \
  --taxon "Abatia rugosa" \
  --backbone data/backbone \
  --output runs/abatia-rugosa
```

The resulting `ref/` directory can be used with SPrOUT `02_exon_trees.py -r`.

## Web app

```bash
export SPROUT_REF_BACKBONE=/absolute/path/to/backbone
export SPROUT_REF_WORK_DIR=/absolute/path/to/runs
sprout-ref web --host 127.0.0.1 --port 8000
```

Open <http://127.0.0.1:8000>.

## Important limits

- This is a research workflow, not a final taxonomic authority.
- A species-level reference is created only when that species has real sequence data.
- Results depend on read depth, recovered genes, Kew coverage, and reference quality.
- Low-abundance mixture members may need more reads or more target genes.

See [How it works](docs/how-it-works.md) for a short explanation and
[HPC examples](examples/hpc/README.md) for the server scripts.

## Tests

```bash
pytest
ruff check src tests
```

## Citation and data

Code is released under the MIT License. Kew data, SPrOUT, Angiosperms353 targets, and the curated
backbone keep their original licenses and citation requirements.

- [Kew Tree of Life public releases](https://sftp.kew.org/pub/treeoflife/)
- [Baker et al. 2022](https://doi.org/10.1093/sysbio/syab035)
- [Johnson et al. 2019](https://doi.org/10.1093/sysbio/syy086)
- [SPrOUT](https://github.com/nhu92/SPrOUT)
