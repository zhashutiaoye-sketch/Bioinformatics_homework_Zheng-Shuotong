# =============================================================================
# Week 4 homework figures (Q1-Q4). Base R graphics; ASCII labels only.
# Run:  Rscript make_figures.R   (from the analysis/ folder)
# Output: ../figures/Q1_workflow.png, Q2_workflow.png, Q3_locus_chain.png,
#         Q4_prioritization.png
# =============================================================================
outdir <- file.path("..", "figures")
dir.create(outdir, showWarnings = FALSE)

newfig <- function(name, w, h, res = 130) png(file.path(outdir, name), width = w, height = h, res = res)

# --- helpers ----------------------------------------------------------------
# Auto-fitted labelled box: shrinks font until the longest line fits the width
# and the number of lines fits the height.
fbox <- function(x, y, w, h, lines, fill = "white", border = "black",
                 cex_max = 0.78, cex_min = 0.48, lwd = 1.4, font = 2) {
  lines <- unlist(strsplit(lines, "\n"))
  n <- length(lines)
  cex_w <- 1.048 * w / max(nchar(lines))
  cex_h <- if (n > 1) h * 0.85 / (3.3 * (n - 1)) else cex_max
  cex <- max(cex_min, min(cex_max, cex_w, cex_h))
  rect(x - w / 2, y - h / 2, x + w / 2, y + h / 2, col = fill, border = border, lwd = lwd)
  ys <- if (n == 1) y else seq(y + 1.65 * cex * (n - 1), y - 1.65 * cex * (n - 1), length.out = n)
  text(x, ys, lines, cex = cex, font = font)
  invisible(cex)
}
conn <- function(x1, y1, x2, y2, col = "grey25", lwd = 1.6, lty = 1, code = 2)
  arrows(x1, y1, x2, y2, col = col, lwd = lwd, lty = lty, length = 0.1, code = code)
cap <- function(x, y, txt, cex = 0.65, col = "grey25", font = 3) text(x, y, txt, cex = cex, col = col, font = font)
frame <- function(title, sub = NULL) {
  par(mar = c(0, 0, if (is.null(sub)) 2.2 else 3.4, 0))
  plot.new(); plot.window(xlim = c(0, 100), ylim = c(0, 100))
  title(main = title, cex.main = 0.98, font.main = 2)
  if (!is.null(sub)) mtext(sub, side = 3, line = 0.25, cex = 0.7, font = 3)
}

# =============================================================== FIGURE 1 (Q1)
newfig("Q1_workflow.png", 1900, 1300)
frame("Q1 - Assay strategy: from \"Gene X is upregulated\" to a causal regulatory mechanism")
fbox(50, 96, 74, 6,
     "BIOLOGICAL QUESTION: Gene X is significantly upregulated in disease vs control. Which mechanism explains it?\nregulatory variant | chromatin accessibility | TF binding or histone modification | DNA methylation | enhancer-promoter contact",
     fill = "#E8F1FB", cex_max = 0.72)
text(25, 88.5, "A. STATE of the locus (what is happening)", cex = 0.85, font = 4)
text(75, 88.5, "B. SEQUENCE and CAUSALITY (why it happens)", cex = 0.85, font = 4)

fbox(25, 82, 42, 8, "STEP 1  RNA-seq (already available)\nmeasures: transcript level, isoforms, allele-specific expression\ncannot prove: which regulatory element, or the direction of causality",
     fill = "#F2F7F2")
conn(25, 78, 25, 74.5)
fbox(25, 69.5, 42, 9, "STEP 2  ATAC-seq + H3K27ac CUT&Tag (or CUT&RUN)\nmeasures: nucleosome-free accessibility and active-enhancer state\nwith input / IgG controls in the matched cell type and condition",
     fill = "#F2F7F2")
conn(25, 65, 25, 61.5)
fbox(25, 56.5, 42, 9, "STEP 3  WGBS / EM-seq (base-resolution 5mC)\nmeasures: DNA methylation of the candidate element\ncannot prove: causality; bulk 5mC hides cell-to-cell heterogeneity",
     fill = "#F2F7F2")
conn(25, 52, 25, 48.5)
fbox(25, 43.5, 42, 9, "STEP 4  Hi-C / Micro-C (or promoter-capture Hi-C)\nmeasures: physical element-promoter contact in 3D space\ncannot prove: that the contact is functional (co-accessibility is not causation)",
     fill = "#F2F7F2")

fbox(75, 82, 42, 8, "STEP 5  WGS / WES (or targeted sequencing of the element)\nmeasures: sequence variants at the element and at Gene X\ncannot prove: which variant changes activity (association only)",
     fill = "#FDF3EC")
conn(75, 78, 75, 74.5)
fbox(75, 69.5, 42, 9, "STEP 6  MPRA / lentiMPRA in the disease-relevant cell type\nmeasures: allele-specific regulatory activity of the element\nlimits: episomal / lentiviral context, not the native locus",
     fill = "#FDF3EC")
conn(75, 65, 75, 61.5)
fbox(75, 56.5, 42, 9, "STEP 7  CRISPRi / CRISPRa / enhancer deletion + RNA-seq\nmeasures: necessity and sufficiency of the element for Gene X\nlimits: perturbation spread, neighbouring genes, off-target effects",
     fill = "#FDF3EC")
conn(75, 52, 75, 48.5)
fbox(75, 43.5, 42, 9, "STEP 8  Same perturbation with ATAC / H3K27ac / Hi-C readout\nmeasures: whether the mechanism runs through accessibility, chromatin\nstate, methylation or 3D contact (mechanism, not only correlation)",
     fill = "#FDF3EC")

conn(46.5, 82, 53.5, 82, col = "#2E6FA8", lty = 2, code = 1); conn(46.5, 82, 53.5, 82, col = "#2E6FA8", lty = 2)
cap(50, 84.4, "variants inside\npeaks", cex = 0.58, col = "#2E6FA8")
conn(53.5, 69.5, 46.5, 69.5, col = "#B5651D", lty = 2); cap(50, 72, "validate in the\nnative locus", cex = 0.58, col = "#B5651D")
conn(53.5, 56.5, 46.5, 56.5, col = "#B5651D", lty = 2); cap(50, 59, "read out the\nchromatin state", cex = 0.58, col = "#B5651D")

fbox(50, 32, 60, 8.5, "DECISION NODE after Step 2\nAccessible + H3K27ac at the element -> proceed to sequence and causality (Steps 5-8)\nNo mark, methylated element, no contact -> test promoter-proximal or post-transcriptional\nmechanisms instead (CAGE / RAMPAGE, RNA stability, reporter assays)",
     fill = "#FFF8DC", cex_max = 0.68)
conn(25, 39, 34, 36.3); conn(75, 39, 66, 36.3)

fbox(50, 20, 80, 8, "INTEGRATION: one model per layer, then ask which layer the perturbation moves\nAccessibility -> chromatin state -> methylation -> 3D contact -> expression -> perturbation (see Q3)",
     fill = "#EFE6F7")
conn(50, 27.7, 50, 24)

fbox(50, 8, 94, 6.5, "THE BIOLOGICAL QUESTION CHOOSES THE ASSAY because the hypothesis must name a mechanism,\nand each assay can only eliminate one class of explanation.",
     fill = "#E4F5E4", lwd = 2, cex_max = 0.8)
conn(50, 16, 50, 11.3)
dev.off()

# =============================================================== FIGURE 2 (Q2)
newfig("Q2_workflow.png", 2000, 1150)
frame("Q2 - Human paired-end sequencing: raw FASTQ to interpretable genomic results",
      "Reference: GRCh38 / hg38, toplevel assembly (analysis set); every downstream coordinate is reported in this build")
fbox(13, 88, 25, 7.5, "1. RAW DATA\nS01_CTRL_WGS_R1/R2.fastq.gz\nrecord md5 + read group", fill = "#E8F1FB")
fbox(41, 88, 27, 7.5, "2. QC PER SAMPLE\nFastQC -> MultiQC\nper-base Q, GC, adapter,\nduplication, overrepresented", fill = "#E8F1FB")
fbox(70, 88, 25, 7.5, "3. TRIM IF REQUIRED\nfastp / Trimmomatic\nadapter + quality trim,\nthen re-run FastQC", fill = "#E8F1FB")
fbox(90.5, 88, 17, 7.5, "4. QC GATE\nfail -> stop\npass -> continue", fill = "#FFF8DC")
conn(25.5, 88, 27.5, 88); conn(54.5, 88, 57.5, 88); conn(82.5, 88, 82, 88)

fbox(13, 72, 25, 7.5, "5. ALIGN TO GRCh38\nbwa-mem2 / BWA-MEM\nread group stored,\nMAPQ recorded", fill = "#F2F7F2")
fbox(41, 72, 27, 7.5, "6. PROCESS MAPPED READS\nsamtools sort/index;\nMarkDuplicates; GATK BQSR;\nflag/MAPQ filtering (samtools view)", fill = "#F2F7F2")
fbox(70, 72, 25, 7.5, "7. CALL VARIANTS\nGATK HaplotypeCaller\nGVCF + joint genotyping\nper-site DP audit", fill = "#F2F7F2")
fbox(90.5, 72, 17, 7.5, "8. FILTER CALLS\nhard filters or VQSR\nkeep PASS only", fill = "#F2F7F2")
conn(13, 84, 13, 76); conn(13, 72, 27.5, 72); conn(54.5, 72, 57.5, 72); conn(82.5, 72, 82, 72)

fbox(24, 54, 33, 9, "9. ANNOTATE\nEnsembl VEP (MANE Select, GRCh38):\nconsequence in Sequence Ontology terms, gene, rsID, population AF\nSpliceAI / CADD / REVEL as supporting evidence, not proof",
     fill = "#FDF3EC")
fbox(66, 54, 36, 9, "10. INTERPRET AND PRIORITISE\nvariant QC (DP, GQ, allele balance) -> ancestry-aware AF filter ->\nconsequence class -> ClinVar / ClinGen / PubMed -> curated shortlist\nwith an explicit evidence ladder (see Q4)",
     fill = "#FDF3EC")
conn(90.5, 68, 90.5, 60, code = 1); arrows(90.5, 60, 84.5, 55.5, length = 0.1, lwd = 1.6)
conn(40.5, 54, 48, 54)

fbox(22, 39, 36, 8.5, "11. VISUALISE\nIGV read-level inspection of every prioritised site,\ncoverage tracks, karyoploteR / ggplot2 summaries", fill = "#EFE6F7")
fbox(70, 39, 36, 8.5, "12. REPRODUCIBILITY RECORD\nassembly + analysis set, tool versions, command log,\nthresholds, sample sheet, raw -> processed manifest", fill = "#EFE6F7")
conn(24, 49.5, 24, 43.5); conn(66, 49.5, 66, 43.5)

fbox(50, 25, 92, 9, "ASSAY-SPECIFIC BRANCH: same FASTQ -> BAM spine, different downstream analysis\nWGS / WES: germline and somatic SNV + indel, CNV, SV, ancestry-aware AF filters      RNA-seq: STAR / salmon -> featureCounts -> DESeq2 or edgeR\nATAC-seq: MACS3 BAMPE peaks -> consensus peaks -> DESeq2 -> TF motif enrichment      ChIP-seq / CUT&Tag: peak calling with input control -> H3K27ac state tracks\nWGBS / EM-seq: Bismark or bwa-meth -> per-CpG methylation -> DSS / methylKit DMRs",
     fill = "#EAF3FA", cex_max = 0.66)
conn(22, 34.5, 22, 29.5); conn(70, 34.5, 70, 29.5)

fbox(50, 12, 90, 7.5, "INTERPRETATION: technical quality -> statistical confidence -> biological meaning -> clinical actionability\nA PASS variant in a well-covered region is still a hypothesis until phenotype, segregation and functional data agree",
     fill = "#E4F5E4")
conn(50, 20.5, 50, 16)

fbox(50, 3.5, 90, 5, "THE ANALYST, NOT THE AI, IS RESPONSIBLE for the reference build, the thresholds, the quality gates and the final interpretation",
     fill = "#FFF0F0", lwd = 2, cex_max = 0.78)
conn(50, 8, 50, 6.1)
dev.off()

# =============================================================== FIGURE 3 (Q3)
newfig("Q3_locus_chain.png", 1900, 1250)
frame("Q3 - Integrated locus view: candidate element upstream of Gene Y",
      "Chain of evidence: accessibility -> chromatin state -> methylation -> 3D contact -> expression -> perturbation")

XA <- 15; XB <- 96   # track band start / end
x_el <- 35; x_pr <- 70

# coordinate axis
lines(c(XA, XB), c(93.5, 93.5), lwd = 1.8)
for (x in seq(XA, XB, by = 9)) lines(c(x, x), c(92.6, 93.9), lwd = 1)
text(c(XA, (XA + XB) / 2, XB), 91.3, c("locus start", "GRCh38 coordinates, one build for all tracks", "locus end"), cex = 0.62)

# vertical guides for element and promoter
segments(x_el, 41, x_el, 92, lty = 3, col = "darkgreen", lwd = 1.2)
segments(x_pr, 41, x_pr, 92, lty = 3, col = "darkorange", lwd = 1.2)

band <- function(ylo, yhi, name) {
  rect(XA, ylo, XB, yhi, col = "#F4F8FC", border = "grey80")
  text(XA - 1.2, (ylo + yhi) / 2, name, cex = 0.66, font = 2, adj = 1)
}
xs <- seq(XA, XB, length.out = 400)
band(80, 90, "ATAC-seq\nread pileup")
sig <- 8.2 * exp(-((xs - x_el) / 2.3)^2) + 5.2 * exp(-((xs - x_pr) / 1.8)^2) + 0.4
polygon(c(XA, xs, XB), c(80, 80 + pmin(sig, 9.4), 80), col = "#4C7FB8", border = NA)
text(XB - 1.5, 86.5, "accessible at the element (and at the promoter)", cex = 0.62, adj = 1, col = "#1F4E79")

band(66, 76, "H3K27ac\nCUT&Tag")
sig2 <- 8.0 * exp(-((xs - x_el) / 2.9)^2) + 7.6 * exp(-((xs - x_pr) / 2.1)^2) + 0.4
polygon(c(XA, xs, XB), c(66, 66 + pmin(sig2, 9.4), 66), col = "#D98C5F", border = NA)
text(XB - 1.5, 72.5, "active-enhancer mark at the element, active promoter at Gene Y", cex = 0.62, adj = 1, col = "#8A4A1E")

band(52, 62, "DNA methylation\n5mC fraction")
sig3 <- 9.0 - 7.0 * exp(-((xs - x_el) / 2.4)^2) - 4.0 * exp(-((xs - x_pr) / 1.6)^2)
polygon(c(XA, xs, XB), c(52.6, 52.6 + pmin(pmax(sig3, 0.5), 9.0), 52.6), col = "#C3B0D8", border = NA)
text(XB - 1.5, 58.5, "hypomethylated CpGs at the element, methylated elsewhere", cex = 0.62, adj = 1, col = "#5B4370")

band(38, 48, "Hi-C / Micro-C\ncontact arcs")
t <- seq(0, pi, length.out = 200)
lines(x_el + (x_pr - x_el) * (1 - cos(t)) / 2, 38.8 + 8.0 * sin(t), col = "#2E6FA8", lwd = 2.4)
text(x_el + 3, 40.6, "anchored contact: element <-> Gene Y promoter", cex = 0.62, font = 2, col = "#2E6FA8", adj = 0)

band(24, 34, "RNA-seq\nGene Y")
rect(x_pr - 1.6, 26.5, x_pr + 1.6, 31.5, col = "#FFE0A3", border = "darkorange", lwd = 1.3)
rect(x_pr + 6, 27.4, x_pr + 20, 30.6, col = "#C7D6EE", border = "navy", lwd = 1.3)
segments(x_pr + 6, 27.4, x_pr + 10, 30.6, col = "navy", lwd = 1.2); segments(x_pr + 12, 27.4, x_pr + 16, 30.6, col = "navy", lwd = 1.2)
arrows(x_pr + 22, 29, x_pr + 16, 29, length = 0.1, lwd = 1.5, col = "navy")
text(x_pr + 13, 25.4, "Gene Y (exons; arrow = direction of transcription)", cex = 0.6, col = "navy")
arrows(x_pr + 13, 32.6, x_pr + 13, 33.8, length = 0.09, lwd = 1.9, col = "#2F5F3F")
text(XA, 21.5, "RNA-seq (same cells, matched condition): Gene Y is higher in disease than control, and allele-level expression can show which allele is active",
     cex = 0.62, adj = 0, col = "#2F5F3F")

# element / promoter labels drawn last so the signal curves cannot hide them
rect(x_el - 9, 90.6, x_el + 9, 92.3, col = "white", border = "darkgreen", lwd = 1.2)
text(x_el, 91.45, "candidate element", cex = 0.62, font = 2, col = "darkgreen")
rect(x_pr - 8, 90.6, x_pr + 8, 92.3, col = "white", border = "darkorange", lwd = 1.2)
text(x_pr, 91.45, "Gene Y promoter", cex = 0.62, font = 2, col = "darkorange")

# evidence chain
cols <- c("#C9E4C5", "#D98C5F", "#C3B0D8", "#9BC2E6", "#A8D5A8", "#E8A0A0")
labs <- c("Accessibility\nATAC peak", "Chromatin state\nH3K27ac", "Methylation\nlow 5mC",
          "3D contact\nHi-C loop", "Expression\nRNA-seq up", "Perturbation\nCRISPRi / deletion")
for (i in 1:6) {
  x <- 4 + (i - 1) * 15.6
  fbox(x + 7, 13, 14.4, 9, labs[i], fill = cols[i], cex_max = 0.66)
  if (i < 6) arrows(x + 14.4, 13, x + 16.8, 13, length = 0.09, lwd = 1.7)
}
cap(50, 6, "Correlation chain: every layer is consistent with an enhancer of Gene Y, but none of them alone is causal - only the perturbation box tests necessity and sufficiency.", cex = 0.68)
cap(50, 3, "Missing evidence to be explicit about: cell-type-matched data, allele-level (AS-ATAC / ASE) information, cell-to-cell heterogeneity, and the neighbouring genes inside the same contact domain.", cex = 0.62)
dev.off()

# =============================================================== FIGURE 4 (Q4)
newfig("Q4_prioritization.png", 1900, 1300)
frame("Q4 - Variant prioritization: plan-first filtering with a data audit (12 synthetic variants -> 1-2 candidates)")

# left: step-wise funnel (vertical)
fy <- c(95, 90.5, 86, 81.5, 77, 72.5, 68)
fn <- c("00  all records: 12", "01  FILTER == PASS: 10", "02  DP >= 20: 9", "03  GQ >= 30: 9",
        "04  AF <= 0.001: 6", "05  gene symbol present: 5", "06  impact consequence: 4")
for (i in 1:7) {
  fbox(16, fy[i], 28, 3.9, fn[i], fill = ifelse(i == 7, "#C9E4C5", "#EAF3FA"), cex_max = 0.7)
  if (i < 7) arrows(16, fy[i] - 2.1, 16, fy[i + 1] + 2.1, length = 0.08, lwd = 1.6)
}
cap(16, 66.2, "per-step output:", cex = 0.6)
cap(16, 64.0, "analysis/q4_filter_output.txt", cex = 0.6)

# right: audit panel
fbox(65, 81.5, 62, 27, "AUDIT BEFORE ACCEPTING ANY RANKING - this is where AI output must be checked against data
1. Build consistency: Ensembl overlap GRCh38 vs GRCh37 for every POS -> the table mixes builds (TP53, BRCA2, HLA-A, ATM fit GRCh38; CFTR, LDLR, KRAS and chr4 fit GRCh37; chrX:153870000 sits inside MECP2 in neither build)
2. ClinVar identifiers: VCV000012345 does not return the TP53 record (ClinVar variation 12345 = TNFRSF1A c.295T>A) -> the IDs are placeholders and every ID must be resolved before use
3. Population frequency: gnomAD r4 reports rs6025 (F5 Leiden, GRCh38 1:169549811 C>T) at AF 0.0173 (genomes) / 0.0219 (exomes), not the 0.42 written in the teaching table
4. Annotation depth: the CONSEQUENCE column comes from a single annotator; a splice call needs SpliceAI-style support plus RT-PCR evidence, not the label alone",
     fill = "#FFF8DC", cex_max = 0.65)
conn(30, 81.5, 34, 81.5, lty = 2, code = 1); conn(30, 81.5, 34, 81.5, lty = 2)
cap(24, 84.6, "audit each", cex = 0.58)
cap(24, 82.6, "surviving set", cex = 0.58)

fbox(27, 50, 48, 16, "TOP CANDIDATE  chr17:7673803 G>A  TP53
splice_acceptor_variant | PASS | DP 80 | GQ 99 | AF 1e-05
Gene-level constraint: TP53 is loss-of-function intolerant (gnomAD v2/v4) and is annotated Pathogenic in this table
Verified: the GRCh38 coordinate lies inside TP53 (chr17:7,661,779 - 7,687,546); gnomAD r4 lists 17-7673803-G-A (genomes AF 0, exomes 4e-06)
Residual risk: single caller, no orthogonal confirmation, no segregation data, germline vs tumour status unknown",
     fill = "#E4F5E4", cex_max = 0.62)
fbox(77, 50, 42, 16, "SECOND CANDIDATE  chr12:25398284 C>A  KRAS
missense_variant (codon-12 hotspot) | PASS | DP 58 | GQ 91 | AF 1.5e-04
The table label is conflicting; resolved ClinVar records for the real rs121913530 alleles are germline Pathogenic / Likely pathogenic (c.34G>T p.Gly12Cys)
Build check: the table coordinate is GRCh37; the GRCh38 equivalent is chr12:25,245,351
Residual risk: codon-12 hotspot variants are usually somatic and treatment-relevant, so the tissue context changes the interpretation",
     fill = "#FDF3EC", cex_max = 0.62)
conn(16, 66.2, 16, 58); conn(16, 58, 27, 58, code = 1, lty = 2); arrows(34, 66.2, 45, 58.4, length = 0.1, lwd = 1.6)

fbox(50, 30, 94, 12, "FALSE-LEAD REVIEW - ask what would make the top variant wrong, then decide which concerns change the decision
Matters scientifically: (a) no orthogonal confirmation of the TP53 splice call -> inspect reads in IGV and re-call with a second caller; (b) mixed genome builds -> re-lift every coordinate before comparing tracks; (c) germline vs tumour context -> changes which ACMG/AMP rules apply and how the result is reported
Does not matter here: the synthetic AF values are teaching placeholders (the filtering logic does not change); the intergenic and intronic rows were excluded by design, so their removal is not a finding",
     fill = "#FFF0F0", cex_max = 0.65)
conn(50, 41.8, 50, 36.2)

fbox(50, 19, 94, 8.5, "KNOWN EVIDENCE -> COMPUTATIONAL INFERENCE -> SCIENTIFIC HYPOTHESIS -> REQUIRED EXPERIMENT
KNOWN: a rare PASS splice-acceptor candidate in a loss-of-function-intolerant tumour suppressor.  INFERRED: it may disrupt the canonical acceptor site and reduce functional p53.
HYPOTHESIS: the allele alters TP53 splicing in the disease tissue.  REQUIRED: orthogonal confirmation, RT-PCR or minigene splicing assay, and tissue-context or segregation data.",
     fill = "#E8F1FB", cex_max = 0.65)
conn(50, 24.4, 50, 23.4)

fbox(50, 7, 96, 7, "Variant chr17:7673803 G>A (TP53 splice acceptor) may influence TP53 transcript integrity and p53-dependent phenotypes by affecting canonical splice-acceptor use;\nthis can be tested by orthogonal variant confirmation plus an RT-PCR or minigene splicing assay in the relevant tissue.",
     fill = "#E4F5E4", lwd = 2, cex_max = 0.76)
conn(50, 14.6, 50, 10.6)
dev.off()

cat("figures written to", normalizePath(outdir), "\n")
print(list.files(outdir))
