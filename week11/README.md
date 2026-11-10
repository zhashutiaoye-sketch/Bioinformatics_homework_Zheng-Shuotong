# Week 11 — Multi-Omics Paper Analysis

Assignment: choose **any three** of the four papers in `Week 11 Homework.docx` and answer all three
questions for each (analysis workflow / data interpretation or modelling / multi-omics integration),
citing figures, analytical methods and key results.

Papers answered here: **Paper 1, Paper 2, Paper 3** (Paper 4 left unattempted, as only three are required).

## Files

| File | What it is |
| --- | --- |
| `Week 11 Homework.docx` | Teacher template, unmodified |
| `Week 11 Homework - Papers 1-3 answered (EN).docx` | **Submission**: answers inserted under each question of a copy of the template |
| `Week 11 Homework - Papers 1-3 answered (EN).pdf` | Rendered copy of the above (13 pages) for quick reading |
| `Paper 1.pdf` | Klingelhuber et al., *Nat Metab* 2024 — spatiotemporal proteomic map of human adipogenesis |
| `Paper 2.pdf` | Oh et al., *Cell Metab* 2020 — gut-microbiome signature predicts cirrhosis |
| `Paper 3.pdf` | Priya et al., *Nat Microbiol* 2022 — shared and disease-specific host gene–microbiome associations |
| `Paper 4.pdf` | α-ketoglutarate signalling to chromatin (not answered) |

## Structure of the submission

- One-line answer-overview table (question → conclusion) directly under the Instructions.
- Per question: bold `Answer — <label>` heading, a `Short answer:` sentence, then 5–7 bullets whose
  bold labels carry the citation (e.g. `Fig. 1c,d; Extended Data Fig. 4; Methods 'PCP analysis'`)
  followed by the concrete numbers (n, median CoV, AUC ± CI, FDR, λ, counts of proteins/species/genes).
- Template text is preserved byte-for-byte; only answer blocks, the selection line, the overview
  table and a note under Paper 4 are added.

## Papers, in one line each

1. **Paper 1** — 1-h DIA LC–MS/MS time course across four human adipocyte models (days 0–14, n = 3),
   Spectronaut directDIA at 1% q-value, median CoV 12.3%; PCA + cross-model temporal correlation
   define 1,979 conserved proteins clustered into early/intermediate/late modules; transcriptome
   integration plus PCP/SVM spatial proteomics nominate C19orf12, validated by co-IP (TOM40/22/5),
   knockdown and Seahorse.
2. **Paper 2** — MiSeq shotgun metagenomes (MetaPhlAn2/HUMAnN2) and UPLC–qTOF metabolomes
   (MS-DIAL/MS-FINDER) → 310 species, 7,355 pathways, 435 metabolites; random forest + 10-fold CV gives
   19 species and 17 metabolites at AUC 0.91, frozen models validate externally at AUC 0.86–0.95 and
   separate cirrhosis from fibrosis at 0.84–0.94 after adding AST.
3. **Paper 3** — 208 paired mucosal RNA-seq/16S samples (CRC, IBD, IBS) processed to VST genes and
   CLR multi-rank taxa; sparse CCA gives pathway-level components, gene-wise lasso with stability
   selection gives gene–taxon pairs; shared pathways pair with disease-specific microbes and shared
   microbes with disease-specific host genes (e.g. *Streptococcus*–RIPK4 in CRC vs DPEP2 in IBS).
