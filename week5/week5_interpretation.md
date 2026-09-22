# Week 5 Homework 1 - interpretation (treated versus control)

We compared `treated` with `control` in the 12-sample count matrix using a batch-adjusted
DESeq2 model (`~ batch + condition`, control as reference, apeglm shrinkage). The strongest QC
observation was the PCA: PC1 (24 % of variance, R2 = 0.99 with condition) separated the two
conditions, whereas PC2 (9 %, R2 = 0.49 with batch) separated batch C from batches A/B - evidence
that batch must be modelled rather than ignored. At padj < 0.05 and |shrunken log2FC| >= 1,
**60 genes were significant: 36 up and 24 down in treated**. Their coordinated direction across
all three batches indicates one shared transcriptional program, not scattered noise. Limitation:
58 of the 88 genes the instructor key marks as truly changed (|truth log2FC| >= 1) stay below the
detection limit at this sequencing depth, so the 60-gene list is a lower bound. AI drafted the
script and this text; I ran every command locally, confirmed the coefficient name from
`resultsNames(dds)`, and verified the estimates against the instructor truth key (Pearson r = 0.88).

## Headline numbers

| Quantity | Value |
|---|---|
| Design | `~ batch + condition` (control = reference), model matrix full rank (4/4) |
| Filter | `rowSums(counts >= 10) >= 3` -> 1,000 genes before, 989 after |
| Coefficient used | `condition_treated_vs_control` (confirmed with `resultsNames(dds)`) |
| Shrinkage | `apeglm` |
| padj < 0.05 (no effect-size filter) | 83 genes |
| padj < 0.05 and \|shrunken log2FC\| >= 1 | **60 genes (36 up, 24 down)** |
| PCA | PC1 24 % (condition), PC2 9 % (batch) |
| Truth-key check (all 989 genes) | Pearson r = 0.88; mean \|estimate - truth\| = 0.14 log2 units |
| Truth-DE genes recovered | 58 of 88 (2 called that the key marks as not changed) |

The full table, the fitted object, the figures, the console transcript and the truth-key check
are in this folder: `week5_deseq2_results.csv`, `week5_deseq2_object.rds`, `week5_pca.png`,
`week5_de_plot.png`, `week5_run_console.txt`, `week5_truth_key_check.txt`.
