# Curated 871-species backbone

Put the SPrOUT curated alignments here (one aligned FASTA per Angiosperms353 gene,
named `<gene_id>.fasta`). Large biological data are intentionally not copied into the
source distribution because the full 871 × 353 resource is not present in the public
SPrOUT repository and its redistribution terms cannot be inferred.

Validate the installed resource with:

```bash
sprout-ref inspect-backbone data/backbone
```

The public SPrOUT repository currently exposes demonstration panels for 50 genes, with
106 order-level and 298 family-level references. Those are useful for a smoke test but
must not be described as the full 871-species backbone.

