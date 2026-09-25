import json, os, base64
from emp_api import post, dump, save_png, OUT

SID = open(os.path.join(OUT, "session_main.txt")).read().strip()
EXP = "genus_full"
PLOTS = os.path.join(OUT, "app_plots")
os.makedirs(PLOTS, exist_ok=True)

def viz(name, path, body):
    r = post(path, dict(body, session_id=SID, experiment=EXP))
    dump("viz_%s.json" % name, {k: (v[:80] + "..." if isinstance(v, str) and len(v) > 80 else v)
                                for k, v in r.items()}, "app")
    b64 = r.get("plot") or r.get("image")
    if isinstance(b64, str) and len(b64) > 200:
        p = save_png(b64, os.path.join(PLOTS, name + ".png"))
        print("  %-34s OK  %8d bytes" % (name, os.path.getsize(p)), flush=True)
    else:
        print("  %-34s FAIL %s" % (name, json.dumps(r, ensure_ascii=False)[:200]), flush=True)

if __name__ == "__main__":
    print("== alpha diversity plots (one per index, by Group) ==", flush=True)
    for m in ["shannon", "simpson", "invsimpson", "chao1", "ace", "observed", "pielou"]:
        viz("alpha_%s" % m, "/api/visualize/alpha", {"group": "Group", "metric": m})
    viz("alpha_shannon_Group_sub", "/api/visualize/alpha", {"group": "Group_sub", "metric": "shannon"})

    print("== ordination ==", flush=True)
    for ordi in ["auto", "emp", "assay_pca"]:
        viz("scatter_%s_Group" % ordi, "/api/visualize/scatter", {"group": "Group", "ordination": ordi})
    viz("scatter_auto_Group_sub", "/api/visualize/scatter", {"group": "Group_sub", "ordination": "auto"})

    print("== composition / abundance ==", flush=True)
    viz("heatmap_top40", "/api/visualize/heatmap", {"group": "Group", "top_n": 40})
    viz("barplot_top20", "/api/visualize/barplot", {"group": "Group", "mode": "top20", "top_n": 20})
    viz("structure_top10", "/api/visualize/structure", {"group": "Group", "top_n": 10})
    viz("sankey_phylum_genus", "/api/workflows/microbiome_16s/visualize/sankey",
        {"from_level": "Phylum", "to_level": "Genus", "top_n": 15})
    viz("network_spearman", "/api/workflows/microbiome_16s/visualize/network",
        {"method": "spearman", "cutoff": 0.6, "top_n": 40})
    viz("volcano", "/api/visualize/volcano", {"fc_cutoff": 1.0, "p_cutoff": 0.05})
    print("done ->", PLOTS, flush=True)
