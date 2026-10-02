# mix7 reference demo

This small demo builds order, family, and genus reference panels for a known four-species
mixture:

- *Allium sativum*
- *Asparagus officinalis*
- *Brassica oleracea*
- *Artocarpus heterophyllus*

After installing the 353-gene backbone, run:

```bash
sprout-ref demo-mix7 \
  --backbone data/backbone \
  --summary data/backbone/species_summary.csv \
  --output runs/mix7-demo
```

The demo checks reference construction and Kew sequence retrieval. It is not a blind read
identification test because the expected taxa are supplied in advance. Use the HPC mixed-read
workflow for an unknown FASTQ sample.
