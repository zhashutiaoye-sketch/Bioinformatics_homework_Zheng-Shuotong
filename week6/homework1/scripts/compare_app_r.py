import csv, re, statistics
H = "D:/大三课件/生信/week6/homework1"

def rd(p):
    return list(csv.DictReader(open(p, encoding="utf-8-sig")))

def corr(a, b):
    n = len(a); ma = sum(a) / n; mb = sum(b) / n
    cov = sum((x - ma) * (y - mb) for x, y in zip(a, b))
    sa = sum((x - ma) ** 2 for x in a) ** .5
    sb = sum((y - mb) ** 2 for y in b) ** .5
    return cov / (sa * sb)

def rank(v):
    s = sorted(range(len(v)), key=lambda i: v[i]); rk = [0] * len(v)
    for pos, i in enumerate(s): rk[i] = pos
    return rk

# --- alpha agreement: app vs R (rarefied) ---
app = {r["primary"]: r for r in rd(H + "/tables/app_contrasts/alpha_genus_full.csv")}
rb = {x["SampleID"]: x for x in rd(H + "/tables/r_alpha_indices.csv")}
common = sorted(set(app).intersection(rb))
print("samples joined:", len(common))
for acol, rcol, name in [("shannon", "shannon", "Shannon"), ("simpson", "simpson", "Simpson"),
                         ("chao1", "chao1", "Chao1"), ("observerd_index", "observed", "Observed")]:
    ss = [s for s in common if app[s].get(acol) not in (None, "", "NA")
          and rb[s].get(rcol) not in (None, "", "NA")]
    a = [float(app[s][acol]) for s in ss]; b = [float(rb[s][rcol]) for s in ss]
    print("  %-9s app vs R  Pearson r = %.3f  (n=%d)  app median=%.3f  R median=%.3f"
          % (name, corr(a, b), len(ss), statistics.median(a), statistics.median(b)))

# --- beta agreement ---
apb = {x["primary"]: x for x in rd(H + "/tables/app_contrasts/dim_genus_full_PCoA.csv")}
rbc = {x["SampleID"]: x for x in rd(H + "/tables/r_pcoa_coords.csv")}
cc = sorted(set(apb).intersection(rbc))
for ax_a, ax_r in [("PCoA1", "Axis1"), ("PCoA2", "Axis2")]:
    a = [float(apb[s][ax_a]) for s in cc]; b = [float(rbc[s][ax_r]) for s in cc]
    print("  %s vs %s  Pearson r = %.3f" % (ax_a, ax_r, corr(a, b)))

# --- differential p-value agreement (matched taxon labels) ---
def norm(x):
    x = x.strip().lower()
    for p in ("k__", "p__", "c__", "o__", "f__", "g__", "s__"):
        if x.startswith(p):
            x = x[3:]
    return re.sub(r"[^a-z0-9]", "", x)

pairs = [("diff_genus_full_IBS_before_great_vs_poor.csv", "IBS_before_great vs IBS_before_poor"),
         ("diff_genus_full_IBS_before_vs_IBS_after.csv", "IBS_before vs IBS_after"),
         ("diff_genus_full_UC_before_vs_UC_after.csv", "UC_before vs UC_after"),
         ("diff_genus_full_IBS_before_vs_UC_before.csv", "IBS_before vs UC_before")]
rw = rd(H + "/tables/r_diff_taxa_wilcoxon.csv")
for f, vs in pairs:
    ap = rd(H + "/tables/app_contrasts/" + f)
    rob = [x for x in rw if x["vs"] == vs]
    m = {norm(x["taxon"].split(".")[0]): float(x["p"]) for x in rob}
    joined = []
    for x in ap:
        k = norm(re.split(r"[;(]", x["feature"])[0])
        if k in m and m[k] > 0:
            joined.append((float(x["pvalue"]), m[k]))
    if len(joined) > 4:
        a = rank([j[0] for j in joined]); b = rank([j[1] for j in joined])
        print("  %-42s matched=%3d  Spearman rho(p)=%.3f  min p app=%.4f R=%.4f"
              % (vs, len(joined), corr(a, b), min(j[0] for j in joined), min(j[1] for j in joined)))
    else:
        print("  %-42s matched=%3d (too few)" % (vs, len(joined)))

# --- edgeR vs app differential overlap (FDR<0.05) ---
eg = rd(H + "/tables/r_diff_taxa_edgeR_sensitivity.csv")
for vs in ["IBS_before vs IBS_after", "UC_before vs UC_after", "IBS_before vs UC_before"]:
    sig = [x for x in eg if x["comparison"] == vs and float(x["FDR"]) < .05]
    print("  edgeR %-24s FDR<0.05 taxa: %s" % (vs, ", ".join(x["taxon"] for x in sig) or "none"))
