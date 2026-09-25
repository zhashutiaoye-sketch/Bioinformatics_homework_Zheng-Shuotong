import json, os, csv, sys
from emp_api import post, get, dump, wait_job, wait_job, OUT

SID = open(os.path.join(OUT, "session_main.txt")).read().strip()
BASE = {
    "data_path": "C:/emp_hw6/inputs/16S_level-7_mapped130.csv",
    "metadata_path": "C:/emp_hw6/inputs/16S_mapping.csv",
    "data_type": "tax",
    "assay_name": "counts",
    "start_level": "Species",
    "tax_sep": ";",
}
RES = os.path.join(OUT, "results")
os.makedirs(RES, exist_ok=True)

def to_csv(rows, path):
    if not rows:
        open(path, "w").write("")
        return
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow(r)

def diff(exp, group_var, ref, test, tag, subset=True):
    r = post("/api/analyze/differential", {
        "session_id": SID, "experiment": exp, "method": "wilcox.test",
        "group_var": group_var, "ref_group": ref, "test_group": test,
        "filter_low": True, "subset_two_groups": subset, "cores": "auto"})
    dump("diff_%s.json" % tag, r, "app")
    data = r.get("data")
    rows = json.loads(data) if isinstance(data, str) else (data or [])
    if rows:
        to_csv(rows, os.path.join(RES, "diff_%s.csv" % tag))
    sig = [x for x in rows if (x.get("fdr") or 1) < 0.05]
    print("  %-34s n=%3d  fdr<0.05: %2d  %s" % (
        tag, len(rows), len(sig),
        ", ".join(sorted({x["feature"] for x in sig})[:6])), flush=True)
    return rows

def collapse(exp, level, keep_top_n=0, drop_unassigned=False):
    r = post("/api/workflows/microbiome_16s/prepare/taxonomy", {
        "session_id": SID, "experiment": exp, "collapse_level": level,
        "keep_top_n": keep_top_n, "drop_unassigned": drop_unassigned,
        "min_total_abundance": 0, "tax_sep": ";"})
    dump("collapse_%s.json" % exp, r, "app")
    print("  collapse %s -> %s: %s features (from %s)" % (
        exp, level, r.get("n_features_after"), r.get("n_features_before")), flush=True)
    return r

CONTRASTS = [
    ("Group", "IBS_before", "IBS_after", "IBS_before_vs_IBS_after"),
    ("Group", "UC_before", "UC_after", "UC_before_vs_UC_after"),
    ("Group", "IBS_before", "UC_before", "IBS_before_vs_UC_before"),
    ("Group", "IBS_after", "UC_after", "IBS_after_vs_UC_after"),
    ("Group_sub", "IBS_before_great", "IBS_before_poor", "IBS_before_great_vs_poor"),
    ("Group_sub", "IBS_after_great", "IBS_after_poor", "IBS_after_great_vs_poor"),
    ("Group_sub", "UC_before_great", "UC_before_poor", "UC_before_great_vs_poor"),
    ("Group_sub", "UC_after_great", "UC_after_poor", "UC_after_great_vs_poor"),
]

if __name__ == "__main__":
    level = sys.argv[1] if len(sys.argv) > 1 else "Genus"
    exp = "%s_full" % level.lower()
    print("=== import experiment %s ===" % exp, flush=True)
    imp = post("/api/import/path", dict(BASE, experiment_name=exp, session_id=SID))
    dump("import_%s.json" % exp, imp, "app")
    print("  samples=%s features=%s overlap=%s" % (imp.get("samples"), imp.get("features"),
                                                   (imp.get("sample_overlap") or {}).get("matched")), flush=True)
    collapse(exp, level)
    prof = post("/api/workflows/microbiome_16s/profile", {"session_id": SID, "experiment": exp, "tax_sep": ";"})
    dump("profile_%s.json" % exp, prof, "app")
    print("  profile:", json.dumps(prof.get("profile", {}), ensure_ascii=False)[:300], flush=True)

    # alpha diversity (all indices in one call)
    a = post("/api/analyze/alpha", {"session_id": SID, "experiment": exp, "method": "shannon"})
    dump("alpha_%s.json" % exp, a, "app")
    arows = json.loads(a["data"]) if isinstance(a.get("data"), str) else []
    to_csv(arows, os.path.join(RES, "alpha_%s.csv" % exp))
    print("  alpha rows=%d columns=%s" % (len(arows), a.get("columns")), flush=True)

    # ordinations
    for m in ["PCA", "PCoA", "NMDS"]:
        r = post("/api/analyze/dimension", {"session_id": SID, "experiment": exp, "method": m})
        dump("dim_%s_%s.json" % (exp, m), r, "app")
        rows = json.loads(r["data"]) if isinstance(r.get("data"), str) else []
        to_csv(rows, os.path.join(RES, "dim_%s_%s.csv" % (exp, m)))
        print("  dimension %s rows=%d cols=%s" % (m, len(rows), r.get("columns")), flush=True)

    print("=== differential contrasts (%s) ===" % exp, flush=True)
    for gv, ref, test, tag in CONTRASTS:
        diff(exp, gv, ref, test, "%s_%s" % (exp, tag))
    # multi-group (all levels at once)
    diff(exp, "Group", None, None, "%s_Group_all" % exp, subset=False)
    diff(exp, "Group_sub", None, None, "%s_Group_sub_all" % exp, subset=False)
