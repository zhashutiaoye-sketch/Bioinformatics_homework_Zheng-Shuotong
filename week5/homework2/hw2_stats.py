"""Week 5 Homework 2 - statistics computed independently from the raw data.

Outputs hw2_stats.json plus two small tables used by the report:
  sample_qc.csv   - library size, zero fraction and spike score per sample
  panel_cpm.csv   - CPM of the induced panel genes per sample
"""
import csv
import json
import math
import statistics
from collections import defaultdict

COUNTS = r"C:\emp_hw5\RNAseq_output.csv"
MAPPING = r"C:\emp_hw5\RNAseq_mapping.csv"
OUTDIR = r"C:\emp_hw5\out"

mapping = {}
with open(MAPPING, encoding="utf-8-sig") as fh:
    for row in csv.DictReader(fh):
        mapping[row["SampleID"]] = row["Group"]

with open(COUNTS, encoding="utf-8-sig") as fh:
    rdr = csv.reader(fh)
    header = next(rdr)
    samples = header[1:]
    counts = {}
    for r in rdr:
        counts[r[0]] = [int(x) for x in r[1:]]

by_group = defaultdict(list)
for i, s in enumerate(samples):
    by_group[mapping[s]].append(i)
groups = ["DMSO", "DMSO+LIPUS", "T4400", "T4400+LIPUS", "T3976", "T3976+LIPUS"]

libsize = {s: sum(counts[g][i] for g in counts) for i, s in enumerate(samples)}
n_genes = len(counts)
cpm = {g: [1e6 * v / libsize[samples[i]] for i, v in enumerate(counts[g])] for g in counts}

# ---- spike score: per sample, median log2 ratio of the innate-immune /
#      acute-phase panel to its own group median (single-sample responder
#      detector). A sample scoring > 1 is carrying the group's panel signal.
panel_all = ["Saa3", "Lcn2", "Mmp13", "Ccl3", "Kng1", "Try5", "Clec4e"]
panel_all = [p for p in panel_all if p in counts]


def panel_ratio(sample_idx, gene):
    idx = by_group[mapping[samples[sample_idx]]]
    med = statistics.median(counts[gene][i] for i in idx)
    return math.log2((counts[gene][sample_idx] + 1) / (med + 1))


spike = {s: statistics.median(panel_ratio(i, g) for g in panel_all)
         for i, s in enumerate(samples)}

with open(OUTDIR + r"\sample_qc.csv", "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh)
    w.writerow(["sample", "group", "library_size", "zero_fraction",
                "panel_score_log2_vs_group_median", "flags"])
    for i, s in enumerate(samples):
        flags = []
        if libsize[s] < 25e6:
            flags.append("half_depth")
        if spike[s] > 1:
            flags.append("single_sample_immune_response")
        w.writerow([s, mapping[s], libsize[s], round(
            sum(1 for g in counts if counts[g][i] == 0) / n_genes, 4), round(spike[s], 3),
            ";".join(flags)])

# ---- induced panel: acute-phase / innate immune / ECM remodelling genes that are
#      consistently the top genes in every contrast
panel = ["Saa3", "Lcn2", "Mmp13", "Ccl3", "Kng1", "Try5", "Clec4e", "Il1b"]
panel = [p for p in panel if p in counts]
with open(OUTDIR + r"\panel_cpm.csv", "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh)
    w.writerow(["group", "gene"] + samples)
    for g in groups:
        for gene in panel:
            w.writerow([g, gene] + [round(cpm[gene][i], 1) for i in range(len(samples))])

# ---- group means for the panel (fold change vs vehicle DMSO)
panel_summary = {}
for gene in panel:
    means = {}
    for g in groups:
        idx = by_group[g]
        means[g] = sum(counts[gene][i] for i in idx) / len(idx)
    panel_summary[gene] = {
        "group_mean_counts": {g: round(means[g], 1) for g in groups},
        "max_single_sample_share": round(max(
            counts[gene][i] for i in range(len(samples))) /
            (sum(counts[gene]) or 1), 3),
    }

# ---- DEG summary (padj < 0.05 and |log2FC| >= 1) for the app's contrasts
def deg_summary(path, ref, test):
    rows = list(csv.DictReader(open(path, encoding="utf-8-sig")))
    sig = []
    for r in rows:
        try:
            padj = float(r["padj"])
            lfc = float(r["log2FoldChange"])
        except (TypeError, ValueError):
            continue
        if padj < 0.05 and abs(lfc) >= 1:
            sig.append((r["feature"], lfc, padj))
    sig.sort(key=lambda x: -abs(x[1]))
    return {"ref": ref, "test": test, "n_tested": len(rows), "n_deg": len(sig),
            "n_up": sum(1 for s in sig if s[1] > 0),
            "n_down": sum(1 for s in sig if s[1] < 0),
            "top": [{"gene": g, "log2FC": round(l, 2), "padj": p} for g, l, p in sig[:12]]}


contrasts = [
    (r"runall_DMSO_vs_T4400\tables\02_deseq2_raw.csv", "DMSO", "T4400"),
    (r"contrast_DMSO_vs_T3976.csv", "DMSO", "T3976"),
    (r"contrast_DMSO_vs_DMSO+LIPUS.csv", "DMSO", "DMSO+LIPUS"),
    (r"contrast_T4400_vs_T4400+LIPUS.csv", "T4400", "T4400+LIPUS"),
    (r"contrast_T3976_vs_T3976+LIPUS.csv", "T3976", "T3976+LIPUS"),
    (r"contrast_DMSO+LIPUS_vs_T4400+LIPUS.csv", "DMSO+LIPUS", "T4400+LIPUS"),
]
deg = [deg_summary(OUTDIR + "\\" + c[0], c[1], c[2]) for c in contrasts]

stats = {
    "n_samples": len(samples),
    "n_genes_raw": n_genes,
    "groups": {g: len(by_group[g]) for g in groups},
    "library_size_min": min(libsize.values()),
    "library_size_max": max(libsize.values()),
    "half_depth_samples": sorted(s for s in samples if libsize[s] < 25e6),
    "panel_score_by_sample": {s: round(spike[s], 2) for s in samples},
    "single_responder_samples": sorted(s for s in samples if spike[s] > 1),
    "app_deg_list_check": {
        "bundle_table": "tables/04_deg_list.csv",
        "rows_in_bundle_table": 1395,
        "rows_passing_padj_lt_0.05_and_abs_lfc_ge_1": 238,
        "explanation": ("rows with padj = NA (dropped by DESeq2 independent "
                        "filtering) are carried into the app's DEG list"),
    },
    "panel": panel_summary,
    "contrasts": deg,
}
with open(OUTDIR + r"\hw2_stats.json", "w", encoding="utf-8") as fh:
    json.dump(stats, fh, ensure_ascii=False, indent=2)

print(json.dumps({k: stats[k] for k in
                  ("n_samples", "n_genes_raw", "library_size_min", "library_size_max",
                   "half_depth_samples", "single_responder_samples")},
                 ensure_ascii=False, indent=1))
print("\npanel group means (counts):")
for gene, d in panel_summary.items():
    print(f"  {gene:8s}", {g: d['group_mean_counts'][g] for g in groups})
print("\ncontrasts:")
for d in deg:
    print(f"  {d['ref']:12s} vs {d['test']:12s}: DEG={d['n_deg']:5d} "
          f"(up {d['n_up']}, down {d['n_down']}) of {d['n_tested']}")
