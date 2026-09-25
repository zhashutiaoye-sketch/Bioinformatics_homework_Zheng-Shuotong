# Week 6 — Microbiome / 16S bioinformatics

Student: 郑烁曈 (SUAT24000155)

| Item | Where |
|---|---|
| Teacher's assignment sheet | `Homework-Week6.docx` (image-only slide: "Homework for week 5/6 — use the 16S data files in the EMP `tests` folder and EMP-web to run the analysis; submit the final result through the EMP-web system; then generate a scientific hypothesis and state which parameters support it") |
| Lecture deck | `Week6_Microbiome_Bioinformatics.pdf` |
| Reading material (rubric, week-6 concepts) | `Week6_Microbiome_Reading_Material_Revised_2026.pdf` |
| **Homework 1 — 16S analysis with EMP-web + hypothesis** | `homework1/` |

## homework1/ — 16S microbiome: full EMP-web analysis, verification and hypothesis

Main report: `homework1/Week6_16S_Homework_report.md` (same content as `.pdf`).

| File / folder | What it contains |
|---|---|
| `homework1/Week6_16S_Homework_report.md` / `.pdf` | full report: data & parameters, the complete EMP-web procedure, the defects found in the app/data, the independent R verification, the statistical reading, the 100–150-word interpretation, and the scientific hypothesis with its supporting parameters |
| `homework1/AI_verification_log.md` | AI use: what was verified, what was wrong and corrected, what was rejected, and how each number was double-checked |
| `homework1/figures/` | app plots (`alpha_*`, `scatter_*`, `heatmap_*`, `barplot_top20`, `structure_top10`, `sankey_*`, `network_*`, `volcano`) and the independent base-R figures (`r_*`) |
| `homework1/tables/` | `r_summary_table.csv`, `r_permanova.csv`, `r_dispersion_tests.csv`, `r_alpha_tests.csv`, `r_diff_taxa_wilcoxon.csv`, `r_diff_taxa_edgeR_sensitivity.csv`, `r_pcoa_coords.csv`, plus the app's own contrast tables in `tables/app_contrasts/` |
| `homework1/app_run_bundle/run1_raw_import_132samples/` | the app's one-click bundle for the file **as provided** (132 libraries) — evidence that the differential step silently produces no table |
| `homework1/app_run_bundle/run2_analysis_130samples/` | the app's one-click bundle for the 130-library QC subset (used in the report): plots (PNG+PDF), tables, `summary.txt`, `bundle.zip` |
| `homework1/app/api_responses/` + `app/api_call_log.md` | every API request/parameter/response of the app run, including both GitHub sync responses |
| `homework1/scripts/` | `hw6_16s_verify.R` (independent verification, 25-point rubric item), the API drivers (`stepA_official.py`, `stepB_deep.py`, `stepC_appplots.py`, `prep_inputs.py`), `compare_app_r.py` |
| `homework1/console/` | R console transcripts (`hw6_16s_verify_console.txt`, `r_rerun_stdout.txt`), `app_vs_r_agreement.txt`, `interpretation_wordcount.txt` |
| `homework1/data/` | the derived 130-sample table, the metadata, and the untouched 132-sample original |

Headline results: 470 level-7 taxa × 130 libraries (2 libraries excluded for missing metadata);
Bray–Curtis PERMANOVA `~ Group` R² = 0.047, p = 0.001 with homogeneous dispersion (p = 0.449);
treatment shifts the UC community (paired, subject-blocked PERMANOVA R² = 0.026, p = 0.007) but not the IBS
community (R² = 0.011, p = 0.132); no taxon reaches FDR < 0.05 with Wilcoxon (edgeR finds 3–5, i.e. the
result is method-dependent); baseline responders and non-responders are statistically indistinguishable
(p = 0.87 / 0.98). Hypothesis: the microbiome response to treatment is a community-level, UC-specific
within-patient shift, so a response read-out must be a change metric — not a baseline classification.

Submission: synced through the EMP-web system (`POST /api/github/sync`) into this repository under
`EMP2026/Week_06/microbiome_16s/weekly/runs/2026-09-25T04-30-07-829Z-coi9gg/` (commit `bbc4000b`,
33 files) and additionally under `EMP2026/Week_05/...` because the assignment sheet is headed "week 5/6"
(commit `927e8cab`).
