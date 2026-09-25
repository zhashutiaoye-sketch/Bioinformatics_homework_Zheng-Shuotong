import re, os
H = "D:/大三课件/生信/week6/homework1"
md = open(H + "/Week6_16S_Homework_report.md", encoding="utf-8").read()
m = re.search(r"### Interpretation \(100–150 words.*?\n\n(.*?)\n\n---", md, re.S)
block = m.group(1)
# strip markdown blockquote markers
text = " ".join(l.lstrip("> ").strip() for l in block.strip().split("\n"))
words = text.split()
out = ("Interpretation word count check\n"
       "===============================\n"
       "requirement: 100-150 words (Week 6 reading material, 'One-Hour Homework' rubric, 20 points)\n"
       "counted:     %d words\n"
       "status:      %s\n\n" % (len(words), "OK" if 100 <= len(words) <= 150 else "OUT OF RANGE"))
open(H + "/console/interpretation_wordcount.txt", "w", encoding="utf-8").write(out + text + "\n")
print(out)
print("FIRST:", " ".join(words[:12]))
print("LAST :", " ".join(words[-12:]))
