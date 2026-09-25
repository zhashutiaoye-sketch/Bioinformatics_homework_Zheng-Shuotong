# AI verification log — Week 6 homework (16S microbiome with EMP-web)

What this file documents (required by the Week 6 reading material: "Record at least one AI suggestion that
you verified, corrected, or rejected"): which AI suggestions were accepted, which were **wrong and had to be
corrected**, which were **rejected**, and how each resulting number was checked against a second source.

Roles: the AI agent (Hermes, model deepseek-v4-flash) drove the EasyMultiProfiler-Web HTTP API and wrote the
first draft of `scripts/hw6_16s_verify.R` and of this report. The student ran the app, re-ran the R script
from a clean session, and inspected the outputs. Environment: R 4.4.3, vegan 2.7-3, edgeR 4.4.2,
EasyMultiProfiler (Web) v9.0.4.

## 1. Suggestions that were verified and kept

| AI suggestion | Verification performed | Outcome |
|---|---|---|
| Report PERMANOVA together with a dispersion test, and use subject-blocked permutations for the paired before/after contrast | `vegan::adonis2(..., permutations = how(nperm=999, blocks=subject))` and `betadisper` + `permutest` re-run by the student; both diagnostics appear in `tables/r_permanova.csv` and `r_dispersion_tests.csv` | kept — the paired UC result (p = 0.007) is only credible because dispersion is homogeneous (p = 0.799) |
| Rarefy to the smallest library depth before alpha/beta diversity, but keep raw counts for the count model | library sizes recomputed (median 10,197; min 1,624; max 20,283, 12.5× range) → rarefaction to 1,624 reads; `edgeR` run separately on raw counts | kept |
| Exclude the two libraries without metadata instead of letting `NA` propagate | the failure is reproducible through the app API: `POST /api/analyze/differential` → `"Column Group has beed deteced missing value"`; after excluding the two libraries the same endpoint returns a full table | kept and documented in the report (§3.1) |
| Check the app's alpha indices against an independent implementation rather than trusting the app | Pearson r between app and R: 0.87 (Shannon), 0.90 (Chao1), 0.91 (observed) over the same 130 samples | kept — with the added observation that the *absolute* values differ because the app computes them after the genus collapse |

## 2. AI suggestions that were WRONG and were corrected after checking

| Wrong AI suggestion | How it failed | Correction |
|---|---|---|
| `adonis2(as.formula("~ Group"), data = md)` inside a helper function | vegan ≥ 2.7 evaluates the response **inside** the formula: first error `找不到对象'Group'` ("object 'Group' not found"), then, after passing the distance object as `dd`, `vegdist(as.matrix(lhs)): input data must be numeric` — the name `dd` had already been used for the read-depth data frame, so vegan tried to compute a distance matrix from a data frame | the distance object is now assigned to the distinctive name `.hw6_dmat` and the formula is built as `.hw6_dmat ~ Group` (`scripts/hw6_16s_verify.R`), verified by re-running the whole script to exit code 0 |
| `sprintf("%d", median(depth))` for the depth summary | R aborted with `'%d'的格式无效` ("invalid format '%d'") three times before the whole QC section printed | all depth/richness values are formatted with `%.0f`; the script now runs end-to-end (see `console/r_rerun_stdout.txt`) |
| "The genus collapse yields genus-level labels" | false for this table: the collapsed features came back as `__`, `s__`, `s__ramosum`, `s__muciniphila` — the deepest *assigned* rank is used, and one feature is an aggregate of all unassigned reads | feature names are re-derived from the level-7 strings in R (genus when annotated, otherwise the deepest annotated rank or the species epithet) — the app's own label string is never quoted alone in the report |
| "A taxon-level hypothesis can be supported by the differential taxa table" | both the app (all 10 contrasts) and the independent Wilcoxon + BH analysis return **0 taxa at FDR < 0.05** | the hypothesis was reformulated at the community (centroid) level; a taxon-specific claim was explicitly rejected |
| "The app's default scatter shows the Bray–Curtis PCoA" | `ordination="auto"` falls back to `scatter_auto_Group.png`, titled **"PCA (assay matrix, Euclidean)"** (PC1 42.5 %) — Euclidean distance on raw counts, the data state the course material warns against for compositional data | the report uses the EMP Bray–Curtis PCoA (`scatter_emp_Group.png`) and the R PCoA; the auto-PCA plot is kept only as an example of the wrong data state |

## 3. Suggestions that were rejected

| Rejected suggestion | Reason |
|---|---|
| Push the analysis to the app's *Week 5* slot as well and describe it as a week-5 result | the assignment sheet is headed "Homework for week **5/6**" but the 16S track's week 5 slot is titled "解读与假设"; the run was synced to the 16S `week_06` slot as the unambiguous target, and the choice is stated in the report |
| Read a single "Firmicutes/Bacteroidetes ratio" as a response biomarker | rejected on the course material's own grounds ("not a universal biomarker"); the descriptive phylum bar plot is presented as composition only, with no test attached |
| Install `pheatmap` and treat the interactive heatmap as the primary heatmap figure | `pheatmap` was installed so the interface works, but the app's interactive heatmap has **no taxon labels and no colour legend**, so the one-click bundle heatmap (`04_top40_heatmap.png`, z-scored, labelled) is used as the heatmap figure instead |
| Report `edgeR` hits as "the differential taxa" | kept only as a **sensitivity check**: Wilcoxon gives 0 taxa, edgeR gives 3–5 for the same contrasts; the disagreement is reported as a limitation (method-dependent result), not as a finding |
| Present the 132-library import as equivalent to the 130-library import | rejected: the raw import reproduces a silent failure inside the one-click pipeline; both bundles are shipped so the difference is auditable |

## 4. Independent checks of the numbers quoted in the report

* Every statistic in the report exists twice: app (`tables/app_contrasts/`) and independent R (`tables/`).
  Agreement table produced by `scripts/compare_app_r.py` (output `console/app_vs_r_agreement.txt`).
* The interpretation block is 100–150 words by requirement; counted programmatically
  (`console/interpretation_wordcount.txt`, 133 words).
* `sessionInfo()` of the verification run is the tail of `console/hw6_16s_verify_console.txt`
  (R 4.4.3, vegan 2.7-3, edgeR 4.4.2, ggplot2 4.0.3, ragg 1.5.2).
* The app's bundle was unzipped and its files listed before use; each quoted plot/table file is present in
  `app_run_bundle/` with the same name as in `summary.txt`.
* Claimed app defects were reproduced, not inferred: the missing `05_diff_taxa` file in run 1, the explicit
  `NA`-in-group error message, the `pheatmap` requirement message, and the EMP-ordination cache message are
  all stored verbatim under `app/api_responses/` and `app_run_bundle/run1_raw_import_132samples/summary.txt`.
