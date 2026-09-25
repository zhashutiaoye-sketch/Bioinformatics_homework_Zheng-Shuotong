import json, os
from emp_api import post, dump, save_png, OUT

SID = open(os.path.join(OUT, "session_main.txt")).read().strip()
EXP = "genus_full"
PLOTS = os.path.join(OUT, "app_plots")

def viz(name, path, body):
    r = post(path, dict(body, session_id=SID, experiment=EXP))
    b64 = r.get("plot") or r.get("image")
    if isinstance(b64, str) and len(b64) > 200:
        p = save_png(b64, os.path.join(PLOTS, name + ".png"))
        print("  %-30s OK %8d bytes" % (name, os.path.getsize(p)), flush=True)
        return True
    print("  %-30s FAIL %s" % (name, json.dumps(r, ensure_ascii=False)[:220]), flush=True)
    return False

# EMP ordination: run the dimension analysis in the same breath as the plot
d = post("/api/analyze/dimension", {"session_id": SID, "experiment": EXP, "method": "PCoA"})
print("dimension PCoA:", d.get("n_rows"), flush=True)
viz("scatter_emp_Group", "/api/visualize/scatter", {"group": "Group", "ordination": "emp"})
viz("scatter_emp_Group_sub", "/api/visualize/scatter", {"group": "Group_sub", "ordination": "emp"})
viz("heatmap_top40", "/api/visualize/heatmap", {"group": "Group", "top_n": 40})
viz("heatmap_top30_Group_sub", "/api/visualize/heatmap", {"group": "Group_sub", "top_n": 30})
