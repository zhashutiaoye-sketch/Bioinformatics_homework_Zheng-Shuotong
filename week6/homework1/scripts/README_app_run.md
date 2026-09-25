# Scripts — Week 6 homework (16S / EasyMultiProfiler-Web)

Everything here was actually executed; the outputs in `../figures/`, `../tables/`, `../console/` and
`../app_run_bundle/` come from these scripts.

## A. Driving the app (the assignment's "complete analysis procedure")

The app's HTTP API is the same code path its buttons use. Start it first (port 8000 is occupied on this
machine by another local program, hence 8010; inputs are staged in an ASCII directory because
`/api/import/path` cannot read non-ASCII paths from this host's R build):

```bash
cd C:/Users/郑烁曈/EasyMultiProfiler-Web
R_LIBS_USER="D:/R/R-userlib" API_HOST=127.0.0.1 API_PORT=8010 \
  EMP_ALLOWED_ROOTS="C:/emp_hw6" EMP_ENABLE_USER_R=false EMP_CORS_ORIGIN='*' NO_PROXY='*' \
  "D:/R/R-4.4.3/bin/Rscript.exe" --vanilla webapp/backend/run_api.R
```

Then, with `D:/python/python.exe` (Python 3.11, stdlib only — `emp_api.py` is a small urllib helper):

| Script | What it does | Output |
|---|---|---|
| `prep_inputs.py` | reads the original `16S_level-7.csv` + `16S_mapping.csv`, drops the 2 libraries that have no metadata → `16S_level-7_mapped130.csv` | `../data/`, `qc_sample_overlap.json` |
| `stepA_official.py` | `POST /api/session` → `/api/import/path` → workflow `profile`/`validate` → the **one-click** `/api/workflows/microbiome_16s/run_all` (group_var=Group, taxonomy_level=Genus, alpha=shannon, beta=bray, ordination=PCoA) → polls the job → downloads the bundle | `../app_run_bundle/run2_analysis_130samples/`, `app/api_responses/10–15*.json` |
| `stepB_deep.py <Genus\|Family>` | imports a second experiment, collapses the taxonomy (`keep_top_n=0`, `drop_unassigned=false`), runs alpha (all 7 indices), dimension (PCA/PCoA/NMDS) and 10 differential contrasts (Group and Group_sub, pairwise + multi-group) | `../tables/app_contrasts/*.csv`, `app/api_responses/diff_*` |
| `stepC_appplots.py`, `stepC2.py` | the visualisation endpoints: alpha (8 calls), scatter (`auto`/`emp`/`assay_pca`), heatmap top-40/top-30, barplot top-20, structure top-10, volcano, sankey (Phylum→Genus), Spearman network; base64 PNG responses are decoded to files | `../figures/*.png` |
| `compare_app_r.py` | joins the app's alpha/ordination tables with the independent R tables and reports Pearson r / Spearman rho | `../console/app_vs_r_agreement.txt` |

`submit_sync.py` (kept in the working directory, not shipped) performed the EMP-web submission:
`POST /api/github/sync` → responses stored in `../app/api_responses/20_github_sync.json` and
`21_github_sync_week05.json`.

## B. Independent verification (the rubric's "R script", 25 points)

```bash
cd D:/大三课件/生信/week6/homework1/scripts
R_LIBS_USER="D:/R/R-userlib" "D:/R/R-4.4.3/bin/Rscript.exe" --vanilla hw6_16s_verify.R
```

`hw6_16s_verify.R` runs from a clean session: identity/orientation checks, read-depth QC, prevalence filter,
rarefaction (seed 20260925), alpha indices, Bray–Curtis distance, PCoA, PERMANOVA (999 permutations,
subject-blocked for paired contrasts), betadisper/permutest, Wilcoxon + BH and an edgeR sensitivity check,
six figures, the summary table, and `sessionInfo()`. Inputs are resolved from the script folder, its parent,
`../data/`, or the staging directory (`resolve_input()`); the console transcript is
`../console/hw6_16s_verify_console.txt` (a second run from this folder is `../console/r_rerun_stdout.txt`).
