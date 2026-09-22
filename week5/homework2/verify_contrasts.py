"""Independent sanity check of the EasyMultiProfiler DESeq2 contrasts.

Loads the raw count matrix + mapping, computes simple pseudo-bulk log2 fold
changes for two contrasts, and correlates them with the log2FoldChange column
the web app returned, to detect any mis-assigned contrast (e.g. group levels
containing '+').
"""
import csv
import math
import statistics
from collections import defaultdict

COUNTS = r"C:\emp_hw5\RNAseq_output.csv"
MAPPING = r"C:\emp_hw5\RNAseq_mapping.csv"
OUT = r"C:\emp_hw5\out"

groups = {}
with open(MAPPING, encoding="utf-8-sig") as fh:
    for row in csv.DictReader(fh):
        groups[row["SampleID"]] = row["Group"]

with open(COUNTS, encoding="utf-8-sig") as fh:
    rdr = csv.reader(fh)
    header = next(rdr)
    samples = header[1:]
    counts = {r[0]: [int(x) for x in r[1:]] for r in rdr}

by_group = defaultdict(list)
for i, s in enumerate(samples):
    by_group[groups[s]].append(i)
print("samples per group:", {k: len(v) for k, v in sorted(by_group.items())})


def mean_counts(gene, group):
    idx = by_group[group]
    v = counts[gene]
    return sum(v[i] for i in idx) / len(idx) if idx else 0.0


def pseudo_lfc(gene, ref, test):
    a, b = mean_counts(gene, ref) + 1, mean_counts(gene, test) + 1
    return math.log2(b / a)


def load_app_table(path):
    with open(path, encoding="utf-8-sig") as fh:
        return {r["feature"]: r for r in csv.DictReader(fh)}


def report(ref, test, app_csv):
    app = load_app_table(app_csv)
    common = [g for g in counts if g in app]
    mine = [pseudo_lfc(g, ref, test) for g in common]
    theirs = [float(app[g]["log2FoldChange"] or 0) for g in common]
    # pearson
    n = len(mine)
    mx, my = statistics.mean(mine), statistics.mean(theirs)
    cov = sum((a - mx) * (b - my) for a, b in zip(mine, theirs))
    sx = math.sqrt(sum((a - mx) ** 2 for a in mine))
    sy = math.sqrt(sum((b - my) ** 2 for b in theirs))
    r = cov / (sx * sy) if sx and sy else float("nan")
    big = sum(1 for g in common if abs(pseudo_lfc(g, ref, test)) >= 1)
    big_app = sum(1 for g in common if abs(float(app[g]["log2FoldChange"] or 0)) >= 1)
    print(f"{ref} vs {test}: genes={n} pearson(pseudo-bulk LFC, app LFC)={r:.3f} "
          f"| |LFC|>=1 pseudo-bulk={big} app={big_app}")
    for g in sorted(common, key=lambda g: -abs(pseudo_lfc(g, ref, test)))[:5]:
        print(f"    {g:12s} pseudo={pseudo_lfc(g, ref, test):+.2f} "
              f"app={float(app[g]['log2FoldChange'] or 0):+.2f} "
              f"padj={app[g]['padj']}")
    return r


print("\n--- contrast assignment check ---")
report("DMSO", "T4400", OUT + r"\runall_DMSO_vs_T4400\tables\02_deseq2_raw.csv")
report("DMSO+LIPUS", "T4400+LIPUS", OUT + r"\contrast_DMSO+LIPUS_vs_T4400+LIPUS.csv")
report("DMSO", "DMSO+LIPUS", OUT + r"\contrast_DMSO_vs_DMSO+LIPUS.csv")

# how different are the LIPUS groups by simple pseudo-bulk ratios?
for ref, test in [("DMSO", "DMSO+LIPUS"), ("T4400", "T4400+LIPUS"),
                  ("T3976", "T3976+LIPUS")]:
    d = [abs(pseudo_lfc(g, ref, test)) for g in counts]
    print(f"{ref} vs {test}: median |pseudo log2FC| = {statistics.median(d):.3f}; "
          f"genes with |LFC|>=1: {sum(1 for x in d if x >= 1)}")
