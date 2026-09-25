import csv, json, os
SRC = "C:/Users/郑烁曈/EasyMultiProfiler-Web/tests/16S_level-7.csv"
MAP = "C:/Users/郑烁曈/EasyMultiProfiler-Web/tests/16S_mapping.csv"
OUTD = "C:/emp_hw6/inputs"
OUT = "C:/emp_hw6/out"

rows = list(csv.reader(open(SRC, encoding="utf-8-sig")))
hdr = rows[0][1:]
mapped = [r["SampleID"] for r in csv.DictReader(open(MAP, encoding="utf-8-sig"))]
keep = [s for s in hdr if s in set(mapped)]
drop = [s for s in hdr if s not in set(mapped)]
print("assay samples:", len(hdr), "metadata samples:", len(mapped))
print("kept:", len(keep))
print("dropped (no metadata):", drop)
with open(os.path.join(OUTD, "16S_level-7_mapped130.csv"), "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(["SampleID"] + keep)
    for row in rows[1:]:
        d = dict(zip(hdr, row[1:]))
        w.writerow([row[0]] + [d[s] for s in keep])
print("written", os.path.join(OUTD, "16S_level-7_mapped130.csv"))
# also a metadata-only copy in ASCII dir
import shutil
shutil.copy(MAP, os.path.join(OUTD, "16S_mapping.csv"))
json.dump({"assay_samples": len(hdr), "metadata_samples": len(mapped), "kept": len(keep), "dropped": drop},
          open(os.path.join(OUT, "qc_sample_overlap.json"), "w"), indent=1)
