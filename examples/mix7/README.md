# mix7 / EM-10 controlled hierarchy demo

Hidden composition supplied for validation:

| Species | Expected order | Expected family |
|---|---|---|
| *Allium sativum* | Asparagales | Asparagaceae |
| *Asparagus officinalis* | Asparagales | Asparagaceae |
| *Brassica oleracea* | Brassicales | Brassicaceae |
| *Artocarpus heterophyllus* | Rosales | Moraceae |

Run after installing the full alignment ZIP into `data/backbone`:

```bash
sprout-ref demo-mix7 \
  --backbone data/backbone \
  --summary data/backbone/species_summary.csv \
  --output runs/mix7-demo
```

Outputs:

- `01_order/ref`: two high-coverage, distinct-family representatives per available order.
- `02_family/ref`: two high-coverage, distinct-genus representatives per family inside the
  simulated order calls `Asparagales, Brassicales, Rosales`.
- `03_genus/ref`: up to two species per represented genus inside the simulated family calls
  `Asparagaceae, Brassicaceae, Moraceae`, followed by exact Kew recovery integration for the
  four validation species.

This is an **oracle-controlled reference-construction demo**, not a blind identification run:
the screenshot contains no FASTQ data, so the correct order and family calls are simulated.
Ground truth is recorded in `demo_report.json`; it is not used to choose the order/family
skeletons, but it is deliberately integrated at the genus validation stage to test exact Kew
fetching, automatic QC, alignment projection, and coverage of all four components.

Observed with Kew Release 4.0 and the supplied 353-locus backbone:

| Species | Recovery | Accepted | Rejected | Missing |
|---|---|---:|---:|---:|
| *Allium sativum* | oneKP GJPF | 281 | 30 | 42 |
| *Asparagus officinalis* | INSDC GCF_001876935.1 | 347 | 2 | 4 |
| *Brassica oleracea* | INSDC GCF_000695525.1 | 350 | 0 | 3 |
| *Artocarpus heterophyllus* | INSDC GCA_025403435.1 | 349 | 2 | 2 |

All three output panels contain 353 equal-width FASTA alignments. Rejected or missing target
loci do not invalidate a panel: the curated skeleton for that locus is retained, and the reason
is recorded rather than repaired manually.
