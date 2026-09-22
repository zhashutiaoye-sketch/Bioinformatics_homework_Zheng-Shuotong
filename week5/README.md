# Week 5 - Transcriptomics: RNA-seq principles and differential expression (DESeq2)

Student: 郑烁曈 (SUAT24000155)

## Homework 1 - `homework1/`

From verified counts to an interpretable DESeq2 result.

| File | What it is |
|---|---|
| `week5_deseq2_analysis.R` | the analysis script (steps 1-12 of the assignment; started from the teacher's starter) |
| `week5_deseq2_results.csv` | full shrunken result table, all 989 filtered genes (padj, shrunken log2FC, direction) |
| `week5_pca.png` | PCA of the 12 libraries - PC1 24 % (condition), PC2 9 % (batch) |
| `week5_de_plot.png` | volcano plot with the \|log2FC\| = 1 and padj = 0.05 thresholds |
| `week5_interpretation.md` | 100-150 word interpretation + headline numbers |
| `week5_AI_verification_log.md` | AI use, checks performed, and the AI errors that had to be fixed |
| `week5_deseq2_object.rds` | fitted DESeq2 object |
| `session_info.txt` | `sessionInfo()` of the run |
| `week5_run_console.txt` | complete console transcript of the run |
| `week5_truth_key_check.txt` | supporting check against the instructor truth key (batch effect, shrinkage, accuracy) |

Headline result: `~ batch + condition`, control as reference, apeglm shrinkage ->
**60 genes significant (36 up / 24 down in treated)** at padj < 0.05 and |shrunken log2FC| >= 1;
estimates agree with the instructor truth key at Pearson r = 0.88.

The script reads the input CSVs from this folder or from the parent course folder.

## Homework 2 - `homework2/`

Complete RNA-seq analysis of `EasyMultiProfiler-Web/tests/RNAseq_output.csv` +
`RNAseq_mapping.csv` (24 mouse libraries, 24,394 genes, 6 groups x 4:
DMSO / T4400 / T3976 each +/- LIPUS) through the EasyMultiProfiler-Web app.

- `Week5_Homework_2_report.md` / `.pdf` - full report (11 pages, 7 figures)
- `figures/` - the app's own plots (`emp_*`) and an independent base-R re-analysis
  (`hw2_pca.png`, `hw2_panel_heatmap.png`, `hw2_libsize.png`)
- `tables/` - the five extra DESeq2 contrasts, sample QC, panel counts, PCA scores
- `emp_run_bundle/` - the app's one-click run bundle (plots, tables, log)
- `run_contrasts.py`, `verify_contrasts.py`, `hw2_stats.py`, `hw2_figures.R` - the scripts that were run

Headline results: DMSO vs T4400 = 238 genes (172 up / 66 down); DMSO vs T3976 = 1;
DMSO vs DMSO+LIPUS = 0; T4400 vs T4400+LIPUS = 0; T3976 vs T3976+LIPUS = 4. The T4400 arm shows
an innate-immune / acute-phase + ECM-remodelling program (`Saa3`, `Lcn2`, `Ccl3`, `Mmp13`, `Kng1`,
`Try5`, `Clec4e`, `Il1a`, `Cxcl10`, `Serpin*`), with one library per group carrying an amplified
version of the same program - the main QC caveat of this dataset.

## Teacher-provided inputs (unchanged, kept in this folder)

`Week5_Homework_1_Instructions.md`, `Week5_Homework_2_Instruction.md`,
`Week5_Homework_Verification_Checklist.md`, `Week5_Homework_Starter.R`,
`Week5_Homework_Count_Matrix.csv`, `Week5_Homework_Sample_Metadata.csv`,
`Week5_Homework_Gene_Annotation_Instructor_Key.csv`.

The lecture PDF is not redistributed here (copyright); it stays in the course folder.
