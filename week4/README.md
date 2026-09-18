# Week 4 — Genomics: Sequencing, Variant Interpretation, and AI-Assisted Analysis

**Submission (Zheng Shuotong).** Deliverable in two formats:

- `Homework for week 4 - Genomics and variant interpretation (EN).md` — main report (markdown).
- `Homework for week 4 - Genomics and variant interpretation (EN).pdf` — same report as PDF (17 pages, figures embedded).

## Folder contents

| Path | What it is |
|---|---|
| `Week4_Homework.md` | the assignment sheet as downloaded (unchanged) |
| `Genomics — Sequencing, Variant Interpretation, and AI-Assisted Analysis-Week 4.pdf` | lecture slides (unchanged) |
| `Week_4_AI_Driven_Genomics_Deep_Analysis_WetLab_AI_FINAL.pdf` | Week-4 reading material (unchanged) |
| `data/` | the Week-4 student package: `variants_q4.tsv`, `sample_manifest.csv`, `demo_fastq/` (FASTQ + `fastqc_snapshot.tsv`), `README_data.md` |
| `figures/` | `Q1_workflow.png`, `Q2_workflow.png`, `Q3_locus_chain.png`, `Q4_prioritization.png` |
| `analysis/` | `q4_filter_Zheng-Shuotong.R` (Q4 base-R filters), `q4_filter_output.txt` (real console output), `q4_shortlist.tsv`, `q4_filter_funnel.txt`, `make_figures.R`, `q2_fastq_qc_notes.md`, `q4_filter_starter.R` |

## How to reproduce

```bash
cd analysis
"D:/R/R-4.4.3/bin/Rscript.exe" q4_filter_Zheng-Shuotong.R   # Q4 funnel + shortlist
"D:/R/R-4.4.3/bin/Rscript.exe" make_figures.R               # the four figures
# FastQC metrics were recomputed from data/demo_fastq/*.fastq.gz with Python (see the Appendix table)
```

## What the report contains

- **Q1** — assay strategy (ATAC-seq + H3K27ac CUT&Tag first, then sequence/variant, methylation, 3D contact, MPRA, CRISPR perturbation), assay comparison table, workflow figure, ~200-word explanation, required ending.
- **Q2** — FASTQ→interpretable-results workflow (GRCh38, QC gate, trimming, alignment, duplicate marking, BQSR, calling, filtering, VEP, visualisation, reproducibility), five FastQC modules re-derived from the raw demo reads (adapter 15.0 %, one template ×18, GC shoulder 10 %, p05 quality ≈ 11 after cycle 50, N 0 %), AI-audit table, workflow figure, required ending.
- **Q3** — layer-by-layer observations vs interpretations vs missing evidence, alternative explanation, CRISPRi/deletion functional test with nascent-transcription readout, integrated locus figure, ~200-word interpretation, required ending.
- **Q4** — my filtering logic, base-R funnel (12 → 4 → 2 candidates; TP53 `splice_acceptor_variant` first, KRAS codon-12 second), four verification passes against live data (Ensembl GRCh38 vs GRCh37 overlaps, ClinVar ID resolution, gnomAD r4 frequencies, PubMed), false-lead critique table, ~200-word interpretation, prioritisation figure, required ending.
- **References** — 32 PubMed IDs, each resolved through NCBI E-utilities before citation; all database queries (Ensembl REST, gnomAD r4 GraphQL, ClinVar) are listed in the Appendix with the commands used.

Note: `data/` files are teaching synthetics (see `data/README_data.md`); nothing in the report is a clinical finding.
