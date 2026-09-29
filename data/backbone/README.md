# Curated Angiosperms353 backbone

The supplied interim archive was audited as 353 aligned loci, 188,167 locus sequences,
821 specimen/source entries and 809 unique binomial names. The two CSV summaries and a
checksum manifest are versioned here; the approximately 189 MB of uncompressed FASTA is
kept outside ordinary Git history.

Install or reproduce the local backbone from the original ZIP:

```bash
sprout-ref install-backbone \
  --archive /path/to/angiosperms_353_v2_interim_targetfile_gene_alignments.zip \
  --destination data/backbone
```

Add `--overwrite` only when intentionally replacing an existing local installation.
The installer extracts only recognized FASTA/CSV files, flattens paths safely, checks
the 353-locus invariant and records the source SHA-256 in `backbone_install.json`. By default it
also applies the versioned taxonomy corrections bundled with MrKnowItAll, excludes specimens
flagged as phylogenetic mislabels, regenerates both summaries, and records correction hit counts.
Use `--no-taxonomy-corrections` only to reproduce the uncorrected archive exactly.

Validate the installed resource with:

```bash
sprout-ref inspect-backbone data/backbone
```

The original working description called this an “871 × 353” collection; the counts
above are the reproducible counts from the supplied archive and are used by this repo.
