import json, os, sys
from emp_api import post, get, dump, wait_job, OUT

EXP = "microbiome16s"
SIDFILE = os.path.join(OUT, "session_main.txt")

IMPORT = {
    "data_path": "C:/emp_hw6/inputs/16S_level-7_mapped130.csv",
    "metadata_path": "C:/emp_hw6/inputs/16S_mapping.csv",
    "experiment_name": EXP,
    "data_type": "tax",
    "assay_name": "counts",
    "start_level": "Species",
    "tax_sep": ";",
}

if __name__ == "__main__":
    # 1. session + import
    s = post("/api/session", {})
    sid = s.get("session_id") or s.get("session", {}).get("session_id")
    print("SESSION", sid, flush=True)
    open(SIDFILE, "w").write(sid)
    imp = post("/api/import/path", dict(IMPORT, session_id=sid))
    dump("10_import.json", imp, "app")
    print("IMPORT", {k: imp.get(k) for k in ("samples", "features", "sample_overlap", "warnings")}, flush=True)

    # 2. workflow profile + validate
    for path, name in [("/api/workflows/microbiome_16s/profile", "11_profile.json"),
                       ("/api/workflows/microbiome_16s/validate", "12_validate.json")]:
        r = post(path, {"session_id": sid, "experiment": EXP, "tax_sep": ";"})
        dump(name, r, "app")
        print(name, json.dumps(r, ensure_ascii=False)[:500], flush=True)

    # 3. one-click 16S pipeline (the app's official "run all")
    r, j = post("/api/workflows/microbiome_16s/run_all", {
        "session_id": sid, "experiment": EXP, "group_var": "Group",
        "taxonomy_level": "Genus", "alpha_index": "shannon",
        "beta_method": "bray", "ord_method": "PCoA"}), None
    j = wait_job(r["job_id"], "run_all")
    dump("13_runall_job.json", j, "app")
    print("run_all:", j.get("status"), flush=True)
    res = get("/api/jobs/%s/result" % r["job_id"])
    dump("14_runall_result.json", res, "app")
    print("bundle:", res.get("zip_path"), flush=True)
    bl = get("/api/bundles/%s" % sid)
    dump("15_bundles.json", bl, "app")
    print("bundles:", [b["name"] for b in bl.get("bundles", [])], flush=True)

    # download the zip
    if bl.get("bundles"):
        name = bl["bundles"][-1]["name"]
        url = "%s/api/bundles/%s/%s" % ("http://127.0.0.1:8010", sid, name)
        import urllib.request
        dst = os.path.join(OUT, "bundle_app_run_main.zip")
        urllib.request.urlretrieve(url, dst)
        print("downloaded", dst, os.path.getsize(dst), flush=True)
