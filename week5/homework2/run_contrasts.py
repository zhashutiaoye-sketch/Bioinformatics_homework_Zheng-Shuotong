"""Run the remaining DESeq2 contrasts through the EasyMultiProfiler-Web API.

For each contrast: POST /api/analyze/differential (DESeq2, two-group subset),
then GET /api/analyze/diff_raw/<sid>/<experiment> to capture the full table.
Results are written next to this script as contrast_<ref>_vs_<test>.csv
plus a machine-readable summary in contrast_summary.json.
"""
import csv
import json
import sys
import urllib.request

API = "http://127.0.0.1:8010"
SESSION = sys.argv[1] if len(sys.argv) > 1 else open("sid.txt").read().strip()
EXPERIMENT = "rnaseq_week5"

CONTRASTS = [
    ("DMSO", "T3976"),
    ("DMSO", "DMSO+LIPUS"),
    ("T4400", "T4400+LIPUS"),
    ("T3976", "T3976+LIPUS"),
    ("DMSO+LIPUS", "T4400+LIPUS"),
]


def post(path, payload, timeout=600):
    req = urllib.request.Request(
        API + path,
        data=json.dumps(payload).encode("utf-8"),
        headers={"Content-Type": "application/json"},
        method="POST",
    )
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return json.loads(r.read().decode("utf-8"))


def get(path, timeout=600):
    with urllib.request.urlopen(API + path, timeout=timeout) as r:
        return json.loads(r.read().decode("utf-8"))


def analyse(ref, test):
    res = post("/api/analyze/differential", {
        "session_id": SESSION,
        "experiment": EXPERIMENT,
        "method": "DESeq2",
        "group_var": "Group",
        "ref_group": ref,
        "test_group": test,
        "filter_low": True,
        "subset_two_groups": True,
    })
    raw = get(f"/api/analyze/diff_raw/{SESSION}/{EXPERIMENT}")
    rows = raw.get("data") or []
    if isinstance(rows, str):
        rows = json.loads(rows)
    name = f"contrast_{ref}_vs_{test}.csv"
    with open(name, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    tested = [r for r in rows if r.get("padj") is not None]
    fdr = [r for r in tested if r["padj"] < 0.05]
    deg = [r for r in fdr if abs(r["log2FoldChange"]) >= 1]
    up = [r for r in deg if r["log2FoldChange"] > 0]
    down = [r for r in deg if r["log2FoldChange"] < 0]
    top = sorted(deg, key=lambda r: -abs(r["log2FoldChange"]))[:10]
    summary = {
        "ref": ref,
        "test": test,
        "n_features": len(rows),
        "n_tested_padj": len(tested),
        "n_padj_lt_0.05": len(fdr),
        "n_deg_padj_and_lfc": len(deg),
        "n_up": len(up),
        "n_down": len(down),
        "max_abs_lfc": round(max((abs(r["log2FoldChange"]) for r in deg), default=0), 3),
        "top10": [{k: r[k] for k in ("feature", "log2FoldChange", "padj")} for r in top],
    }
    print(json.dumps({k: v for k, v in summary.items() if k != "top10"},
                     ensure_ascii=False))
    print("   top:", ", ".join(f"{r['feature']}({r['log2FoldChange']:+.1f})"
                                for r in top[:6]))
    return summary


if __name__ == "__main__":
    out = []
    for ref, test in CONTRASTS:
        print(f"[{ref} vs {test}]")
        out.append(analyse(ref, test))
    with open("contrast_summary.json", "w", encoding="utf-8") as fh:
        json.dump(out, fh, ensure_ascii=False, indent=2)
    print("written: contrast_summary.json")
