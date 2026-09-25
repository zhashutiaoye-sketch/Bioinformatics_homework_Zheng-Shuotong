# EMP-web API call log — Week 6 homework (16S)

Everything the assignment calls "the EMP-web procedure" was executed through the app's own HTTP API
(Plumber, `webapp/backend/plumber.R`), which is the code path the interface buttons call. Session
`jlRMRNfHbzrDVVBk2MC7gUBf`; experiments `microbiome16s` (one-click) and `genus_full` (deep run).
Raw responses: `api_responses/*.json`.

## Launch (loopback ⇒ `emp_auth_required()` is FALSE, principal "local")

```bash
cd C:/Users/郑烁曈/EasyMultiProfiler-Web
R_LIBS_USER="D:/R/R-userlib" API_HOST=127.0.0.1 API_PORT=8010 \
  EMP_ALLOWED_ROOTS="C:/emp_hw6" EMP_ENABLE_USER_R=false EMP_CORS_ORIGIN='*' NO_PROXY='*' \
  "D:/R/R-4.4.3/bin/Rscript.exe" --vanilla webapp/backend/run_api.R
# health: GET http://127.0.0.1:8010/api/health -> {"status":"ok","version":"9.0.4"}
```

Port 8000 (the app default) is already held by another local program on this machine, hence 8010.
Non-ASCII paths break `/api/import/path` on this host, so the inputs were staged in `C:/emp_hw6/inputs`
and registered through `EMP_ALLOWED_ROOTS`.

## 1. Import and workflow validation (UI: Data → Import, Workflow → Validate)

| # | Call | Parameters | Result (`api_responses/`) |
|---|---|---|---|
| 1 | `POST /api/session` | — | session id `jlRMRNfHbzrDVVBk2MC7gUBf` |
| 2 | `POST /api/import/path` | `data_path=C:/emp_hw6/inputs/16S_level-7_mapped130.csv`, `metadata_path=…/16S_mapping.csv`, `data_type=tax`, `assay_name=counts`, `start_level=Species`, `tax_sep=;` | 130 samples × 470 features, `sample_overlap.matched=130/130`, no warnings (`10_import.json`) |
| 3 | `POST /api/workflows/microbiome_16s/profile` | `tax_sep=;` | max taxonomy depth 7, 470 features on every rank (`11_profile.json`) |
| 4 | `POST /api/workflows/microbiome_16s/validate` | `tax_sep=;` | `has_features`, `has_samples`, `taxonomy_depth_ok` = true (`12_validate.json`) |

Run 1 (on the file *as provided*, 132 libraries) used exactly this sequence and is kept in
`../app_run_bundle/run1_raw_import_132samples/` for comparison.

## 2. One-click 16S pipeline (UI: Workflow → 16S → Run all)

`POST /api/workflows/microbiome_16s/run_all`
→ `{group_var: "Group", taxonomy_level: "Genus", alpha_index: "shannon", beta_method: "bray", ord_method: "PCoA"}`
→ background job, polled with `GET /api/jobs/<id>`, result with `GET /api/jobs/<id>/result`
(`13_runall_job.json`, `14_runall_result.json`):

| Step | Status | Time |
|---|---|---|
| 1 Taxonomy prep (Genus, top 40) | ok | 0.55 s |
| 2 Alpha diversity (shannon, simpson, invsimpson, chao1, ace, observed, pielou) | ok | 5.26 s |
| 3 Beta diversity (PCoA, PCA, NMDS) | ok | 4.43 s |
| 4 Heatmap (top 40) + barplot (top 15) | ok | 1.89 s |
| 5 Differential taxa (Wilcoxon, `group_var=Group`) | ok | 18.96 s |
| 6 Metadata snapshot | ok | 0.01 s |

Bundle: `GET /api/bundles/<session_id>` → `20260925-121232.zip` (3,718 KB),
download `GET /api/bundles/<session_id>/<name>` (see `15_bundles.json`); contents in
`../app_run_bundle/run2_analysis_130samples/`.

## 3. Deep run (experiment `genus_full`)

| Call | Parameters | Result |
|---|---|---|
| `POST /api/import/path` | same files, `experiment_name=genus_full` | `import_genus_full.json` |
| `POST /api/workflows/microbiome_16s/prepare/taxonomy` | `collapse_level=Genus`, `keep_top_n=0`, `drop_unassigned=false`, `min_total_abundance=0` | 470 → 137 features, snapshot `20260925-121349_m16s_taxonomy_Genus` (`collapse_genus_full.json`) |
| `POST /api/analyze/alpha` | `method=shannon` (returns every index) | 130 rows × 7 indices (`alpha_genus_full.json`) |
| `POST /api/analyze/dimension` | `method` = PCA / PCoA / NMDS | coordinates (`dim_genus_full_*.json`); NMDS returned 0 rows through the API |
| `POST /api/analyze/differential` × 10 | `method=wilcox.test`, `group_var`=Group or Group_sub, `ref_group`/`test_group`, `subset_two_groups` true (pairwise) or false (multi-group), `filter_low=true` | `diff_genus_full_*.json`; used `data` payload serialised to `../tables/app_contrasts/*.csv` |

Contrasts: IBS_before↔IBS_after, UC_before↔UC_after, IBS_before↔UC_before, IBS_after↔UC_after,
IBS_before_great↔IBS_before_poor, IBS_after_great↔IBS_after_poor, UC_before_great↔UC_before_poor,
UC_after_great↔UC_after_poor, plus the multi-group runs on `Group` and `Group_sub`.

## 4. Visualisation (UI: the plot buttons)

| Call | Parameters | Output file in `../figures/` |
|---|---|---|
| `POST /api/visualize/alpha` × 8 | `group=Group` (7 indices) and `group=Group_sub`, `metric=shannon` | `alpha_<index>.png`, `alpha_shannon_Group_sub.png` |
| `POST /api/visualize/scatter` | `group=Group`, `ordination=auto` | `scatter_auto_Group.png` (PCA of the raw assay — see report §3.4) |
| `POST /api/visualize/scatter` | `group=Group`, `ordination=emp` (after `POST /api/analyze/dimension`) | `scatter_emp_Group.png`, `scatter_emp_Group_sub.png` |
| `POST /api/visualize/barplot` | `group=Group`, `mode=top20`, `top_n=20` | `barplot_top20.png` |
| `POST /api/visualize/heatmap` | `group=Group`, `top_n=40` (needs `pheatmap`) | `heatmap_top40.png` |
| `POST /api/visualize/structure` | `group=Group`, `top_n=10` | `structure_top10.png` |
| `POST /api/visualize/volcano` | `fc_cutoff=1.0`, `p_cutoff=0.05` | `volcano.png` |
| `POST /api/workflows/microbiome_16s/visualize/sankey` | `from_level=Phylum`, `to_level=Genus`, `top_n=15` | `sankey_phylum_genus.png` |
| `POST /api/workflows/microbiome_16s/visualize/network` | `method=spearman`, `cutoff=0.6`, `top_n=40` | `network_spearman.png` |

Plots are returned as base64 PNG and were decoded to files by `../scripts/stepC_appplots.py`.

## 5. Submission (UI: Export → Sync)

`POST /api/github/sync` — `{track_id: "microbiome_16s", assignment_id: "week_06",
session_id: "jlRMRNfHbzrDVVBk2MC7gUBf", experiment: "genus_full", include_rds: false}`
→ writes `EMP2026/Week_06/microbiome_16s/weekly/runs/<run_id>/` into the student's bound repository and
appends to `.local_run/data/students/SUAT24000155/sync_log.jsonl`; response in
`api_responses/20_github_sync.json`.
