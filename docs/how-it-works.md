# How MrKnowItAll works

MrKnowItAll adds automatic reference building to SPrOUT.

## Two kinds of sequence data

The user's mixed sample starts as paired reads. SPrOUT uses fastp and HybPiper to recover target
genes from those reads. These recovered genes are the unknown query sequences.

Kew provides assembled Angiosperms353 recovery FASTA files. MrKnowItAll downloads these assembled
genes, not the original Kew reads. They become the known reference sequences.

## Reference building

For each Kew sample, the program:

1. Reads its order, family, genus, and species from Kew metadata and the Kew species tree.
2. Separates its sequences by Angiosperms353 gene ID.
3. Rejects sequences that are very short, ambiguous, unusual, or inconsistent with the curated
   backbone.
4. Keeps the best sequence when a gene has more than one candidate.
5. Adds accepted sequences to the backbone alignment with MAFFT.
6. Writes one aligned FASTA file per gene for SPrOUT.

The curated backbone is used for checking and alignment. It is not copied into every final Kew
panel.

## Step-by-step identification

Searching every Kew plant at once would be slow and biased toward groups with more public data.
MrKnowItAll therefore uses three smaller searches:

```text
order → family → genus
```

Each family or genus receives the same requested number of reference samples. This makes the
SPrOUT scores easier to compare.

## Mixed-sample candidate selection

A strong component can make a weaker component's global z-score look small. MrKnowItAll combines:

- the main SPrOUT ranking;
- positive z-scores;
- the nearest reference for each query exon; and
- support from at least two separate locus trees.

The final evidence table shows why each genus was retained.

## Main output files

- `FINAL_GENUS_CANDIDATES.txt`: short result list.
- `FINAL_GENUS_RANKING.csv`: all genus scores.
- `FINAL_GENUS_EVIDENCE.csv`: locus votes, query votes, distances, and SPrOUT scores.
- `selection.tsv`: Kew samples used in a reference panel.
- `qc.json`: sequence checks and rejection reasons.
- `panel_report.json`: Kew release and build settings.

No rejected sequence is manually repaired. Missing genes remain missing.
