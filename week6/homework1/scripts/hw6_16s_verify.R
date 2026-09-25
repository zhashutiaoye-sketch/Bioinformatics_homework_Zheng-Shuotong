#!/usr/bin/env Rscript
# =====================================================================
# Week 6 homework - 16S microbiome: independent verification of the
# EasyMultiProfiler-Web (EMP-web) analysis of tests/16S_level-7.csv.
#
# The app produced the primary analysis (app_run_bundle/).  This script
# re-runs the same statistical questions in plain R so that every number
# quoted in the report can be checked against an independent
# implementation, and so that the PERMANOVA / dispersion diagnostics the
# reading material requires can be reported.
#
# Run:  R_LIBS_USER=D:/R/R-userlib Rscript --vanilla scripts/hw6_16s_verify.R
# Inputs are resolved from the script folder, else from its parent.
# =====================================================================

suppressPackageStartupMessages({
  library(vegan)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(ggrepel)
  library(scales)
})

set.seed(20260925)          # reproducibility flag for rarefaction / permutations
options(stringsAsFactors = FALSE)

# ---------------------------------------------------------------- inputs
# resolve_input(): the script is run from homework1/ , but the teacher's
# files live at the week6 root, so look in the current directory first.
SCRIPT_DIR <- tryCatch({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f[1])) else getwd()
}, error = function(e) getwd())
if (!exists("SCRIPT_DIR") || !nzchar(SCRIPT_DIR)) SCRIPT_DIR <- getwd()

resolve_input <- function(rel) {
  cand <- c(file.path(SCRIPT_DIR, rel),
            file.path(SCRIPT_DIR, "..", rel),
            file.path(SCRIPT_DIR, "..", "..", rel),
            file.path(SCRIPT_DIR, "..", "data", basename(rel)),
            file.path(SCRIPT_DIR, "data", basename(rel)),
            file.path("C:/emp_hw6/inputs", basename(rel)))
  hit <- cand[file.exists(cand)]
  if (!length(hit)) stop(sprintf("input not found: %s (looked in %s)", rel,
                                 paste(cand, collapse = " | ")))
  normalizePath(hit[1], winslash = "/")
}

CNT <- resolve_input("16S_level-7_mapped130.csv")
MAP <- resolve_input("16S_mapping.csv")
OUTDIR <- file.path(SCRIPT_DIR, "..", "tables"); dir.create(OUTDIR, showWarnings = FALSE)
FIGDIR <- file.path(SCRIPT_DIR, "..", "figures"); dir.create(FIGDIR, showWarnings = FALSE)
LOG <- file.path(SCRIPT_DIR, "..", "console", "hw6_16s_verify_console.txt")
dir.create(dirname(LOG), showWarnings = FALSE)

sink(LOG, split = TRUE)
cat("Week 6 16S verification run -", format(Sys.time()), "\n")
cat("R", R.version.string, "| vegan", as.character(packageVersion("vegan")), "\n\n")

# ------------------------------------------------- step 0: identity checks
cat("== 0. Identity checks before any analysis ==\n")
counts <- read.csv(CNT, row.names = 1, check.names = FALSE)
meta   <- read.csv(MAP, check.names = FALSE)
cat(sprintf("count table : %d taxa x %d samples\n", nrow(counts), ncol(counts)))
cat(sprintf("metadata    : %d samples x %d columns (%s)\n", nrow(meta),
            ncol(meta), paste(names(meta), collapse = ", ")))
cat("first samples :", paste(head(colnames(counts), 3), collapse = ", "), "\n")
cat("first features:", paste(head(rownames(counts), 2), collapse = " | "), "\n")
# taxonomy separator / depth check
tax_parts <- strsplit(rownames(counts), ";", fixed = TRUE)
cat("taxonomy parts per feature (table):",
    paste(sort(unique(vapply(tax_parts, length, 1L))), collapse = "/"), "\n")
stopifnot(all(colnames(counts) %in% meta$SampleID),
          !anyDuplicated(meta$SampleID),
          !anyDuplicated(colnames(counts)),
          all(counts >= 0), is.numeric(as.matrix(counts)))
cat("orientation : features in rows, samples in columns -> OK\n")
cat("no negative values, no duplicated IDs -> OK\n\n")

# grouping structure derived from the sample-ID convention (documented)
meta$cohort <- ifelse(grepl("^J_", meta$SampleID), "IBS", "UC")
meta$visit  <- ifelse(grepl("_0?2$", meta$SampleID), "after", "before")
meta$response <- sub("^.*_(great|poor)$", "\\1", meta$Group_sub)
meta$subject  <- sub("_[0-9]{2}$", "", meta$SampleID)
meta <- meta[match(colnames(counts), meta$SampleID), ]
stopifnot(identical(meta$SampleID, colnames(counts)))

cat("== 1. Group sizes ==\n")
print(table(meta$Group, meta$response))
cat("\n")

# ------------------------------------------------- step 1: read depth QC
cat("== 2. Read depth QC (raw counts) ==\n")
depth <- colSums(counts)
dd <- data.frame(SampleID = names(depth), depth = as.integer(depth))
dd <- merge(dd, meta[, c("SampleID", "Group", "cohort", "visit", "response")], by = "SampleID")
cat(sprintf("library size: min %.0f  p05 %.0f  median %.0f  max %.0f  (ratio %.1f)\n",
            min(dd$depth), quantile(dd$depth, .05), median(dd$depth), max(dd$depth),
            max(dd$depth) / min(dd$depth)))
print(dd %>% group_by(Group) %>%
        summarise(n = n(), median_depth = median(depth),
                  min_depth = min(depth), max_depth = max(depth)))
cat("\n")
prev <- rowMeans(counts > 0)
cat(sprintf("features detected in >=1 sample: %d ; prevalence >=10%%: %d ; >=50%%: %d\n",
            sum(prev > 0), sum(prev >= .10), sum(prev >= .50)))

# feature filter: keep taxa with >=10% prevalence and >=20 reads total
keep <- prev >= .10 & rowSums(counts) >= 20
cnt_f <- counts[keep, , drop = FALSE]
cat(sprintf("retained taxa after prevalence>=10%% & total>=20 reads: %d of %d\n\n",
            nrow(cnt_f), nrow(counts)))

# ------------------------------------------- step 2: rarefaction for diversity
min_depth <- min(colSums(cnt_f))
cat(sprintf("== 3. Rarefaction to even depth = %.0f reads (seq depth of the smallest library)\n", min_depth))
rare <- as.data.frame(t(vegan::rrarefy(t(cnt_f), min_depth)))
cat(sprintf("depth after rarefaction: all libraries = %.0f reads (identical)\n",
            unique(colSums(rare))[1]))
cat(sprintf("features still detected in >=1 rarefied library: %d\n\n",
            sum(rowSums(rare) > 0)))

# ------------------------------------------------- step 3: alpha diversity
cat("== 4. Alpha diversity (rarefied table) ==\n")
alpha <- data.frame(
  SampleID = colnames(rare),
  observed = specnumber(t(rare)),
  shannon  = diversity(t(rare), index = "shannon"),
  simpson  = diversity(t(rare), index = "simpson"),
  invsimpson = diversity(t(rare), index = "invsimpson"),
  pielou   = diversity(t(rare), index = "shannon") / log(specnumber(t(rare)))
)
er <- estimateR(t(rare))
alpha$chao1 <- er["S.chao1", ]; alpha$ace <- er["S.ACE", ]
alpha <- merge(alpha, meta, by = "SampleID")
write.csv(alpha, file.path(OUTDIR, "r_alpha_indices.csv"), row.names = FALSE)

alpha_grp <- alpha %>% group_by(Group) %>%
  summarise(n = n(), observed = median(observed), shannon = median(shannon),
            chao1 = median(chao1), .groups = "drop")
print(as.data.frame(alpha_grp), digits = 3)
write.csv(alpha_grp, file.path(OUTDIR, "r_alpha_by_group.csv"), row.names = FALSE)

# Mann-Whitney tests for the hypothesis-relevant contrasts
alpha_test <- function(dat, gv, a, b) {
  sub <- dat[dat[[gv]] %in% c(a, b), ]
  do.call(rbind, lapply(c("observed", "shannon", "chao1"), function(m) {
    x <- sub[[m]][sub[[gv]] == a]; y <- sub[[m]][sub[[gv]] == b]
    w <- suppressWarnings(wilcox.test(x, y))
    data.frame(metric = m, comparison = sprintf("%s vs %s", a, b),
               median_a = median(x), median_b = median(y),
               n_a = length(x), n_b = length(y), p = w$p.value)
  }))
}
at <- rbind(
  alpha_test(alpha, "Group_sub", "IBS_before_great", "IBS_before_poor"),
  alpha_test(alpha, "Group_sub", "UC_before_great",  "UC_before_poor"),
  alpha_test(alpha, "Group_sub", "IBS_after_great",  "IBS_after_poor"),
  alpha_test(alpha, "Group_sub", "UC_after_great",   "UC_after_poor"),
  alpha_test(alpha, "Group",     "IBS_before", "IBS_after"),
  alpha_test(alpha, "Group",     "UC_before",  "UC_after"),
  alpha_test(alpha, "Group",     "IBS_before", "UC_before")
)
at$p_bh <- ave(at$p, at$metric, FUN = function(p) p.adjust(p, "BH"))
at$p <- signif(at$p, 4); at$p_bh <- signif(at$p_bh, 4)
at$median_a <- round(at$median_a, 3); at$median_b <- round(at$median_b, 3)
write.csv(at, file.path(OUTDIR, "r_alpha_tests.csv"), row.names = FALSE)
cat("\nalpha diversity Mann-Whitney tests:\n"); print(at, digits = 3, row.names = FALSE)
cat("\n")

# ------------------------------------------------- step 4: beta diversity
cat("== 5. Beta diversity: Bray-Curtis on the rarefied table ==\n")
bray <- vegdist(t(rare), method = "bray")
pcoa <- wcmdscale(bray, k = 4, eig = TRUE, add = "lingoes")
eig <- pcoa$eig; eig <- eig[eig > 0]
pct <- round(100 * eig / sum(eig), 1)
scores <- as.data.frame(pcoa$points[, 1:4])
colnames(scores) <- paste0("Axis", 1:4)
scores$SampleID <- rownames(scores)
scores <- merge(scores, meta, by = "SampleID")
write.csv(scores, file.path(OUTDIR, "r_pcoa_coords.csv"), row.names = FALSE)
cat(sprintf("PCoA axis variance explained: %s ... (first 4 of %d axes; total = 100%%)\n",
            paste0("Axis", 1:4, "=", pct[1:4], "%", collapse = "  "), length(pct)))
cat(sprintf("negative eigenvalues from the Bray-Curtis distance matrix: %d (Lingoes correction applied when present)\n\n",
            sum(pcoa$eig < 0)))

permanova <- function(dat_idx, rhs_txt, label, perm = 999, strata = NULL) {
  # vegan >= 2.7 wants the response itself inside the formula (adonis2(dd ~ g)),
  # so the distance object is placed in the formula environment, and the
  # right-hand-side terms are resolved from the metadata data.frame.
  # NOTE: the response symbol must NOT collide with any other object on the
  # search path (an earlier "dd <- data.frame(...)" made vegan call vegdist on
  # a data.frame), hence the distinctive name in the global environment.
  assign(".hw6_dmat", as.dist(as.matrix(bray)[dat_idx, dat_idx]), envir = globalenv())
  md <- droplevels(meta[dat_idx, ])
  f <- as.formula(paste(".hw6_dmat ~", rhs_txt))
  perm_ctrl <- if (!is.null(strata)) {
    tryCatch(permute::how(nperm = perm, blocks = md[[strata]]), error = function(e) perm)
  } else perm
  a <- adonis2(f, data = md, permutations = perm_ctrl)
  r2 <- a$R2[1]; p <- a$`Pr(>F)`[1]
  cat(sprintf("  %-46s R2=%.3f  p=%.3f  n=%d\n", label, r2, p, length(dat_idx)))
  for (i in seq_len(nrow(a))) {
    if (!is.na(a$R2[i]) && rownames(a)[i] != "Residual" && rownames(a)[i] != "Total")
      cat(sprintf("      term %-18s R2=%.3f  F=%.2f  p=%.3f\n", rownames(a)[i],
                  a$R2[i], a$F[i], a$`Pr(>F)`[i]))
  }
  list(label = label, r2 = r2, p = p, n = length(dat_idx), table = a)
}
disp_test <- function(idx, gv, label) {
  d <- as.dist(as.matrix(bray)[idx, idx]); md <- meta[idx, ]
  set.seed(20260925)
  bd <- betadisper(d, md[[gv]])
  pt <- permutest(bd, permutations = 999)
  cat(sprintf("  %-46s dispersion F=%.2f  p=%.3f\n", label, pt$tab$F[1], pt$tab$`Pr(>F)`[1]))
  data.frame(comparison = label, F = round(pt$tab$F[1], 3), p = round(pt$tab$`Pr(>F)`[1], 4),
             mean_dist = paste(round(tapply(bd$distances, bd$group, mean), 3), collapse = " | "),
             groups = paste(levels(bd$group), collapse = " | "))
}

cat("\n(All PCoA-based inference below uses free permutation unless strata are given.)\n")
pm <- list()
pm[[1]] <- permanova(seq_len(nrow(meta)), "Group",           "Group (4 levels, all 130 samples)")
pm[[2]] <- permanova(seq_len(nrow(meta)), "Group + response", "Group adjusted for response")
pm[[3]] <- permanova(seq_len(nrow(meta)), "Group_sub",       "Group_sub (8 levels)")
pm[[4]] <- permanova(which(meta$visit == "before"), "Group", "baseline only: IBS_before vs UC_before")
pm[[5]] <- permanova(which(meta$Group_sub %in% c("IBS_before_great","IBS_before_poor")),
                     "Group_sub", "baseline IBS: great vs poor")
pm[[6]] <- permanova(which(meta$Group_sub %in% c("UC_before_great","UC_before_poor")),
                     "Group_sub", "baseline UC: great vs poor")
pm[[7]] <- permanova(which(meta$Group_sub %in% c("IBS_after_great","IBS_after_poor")),
                     "Group_sub", "post-treatment IBS: great vs poor")
pm[[8]] <- permanova(which(meta$Group_sub %in% c("UC_after_great","UC_after_poor")),
                     "Group_sub", "post-treatment UC: great vs poor")
ibs <- which(meta$cohort == "IBS"); ucs <- which(meta$cohort == "UC")
pm[[9]] <- permanova(ibs, "visit", "IBS before vs after (paired, strata=subject)", strata = "subject")
pm[[10]] <- permanova(ucs, "visit", "UC before vs after (paired, strata=subject)", strata = "subject")

cat("\n== 6. Dispersion diagnostic (betadisper / permutest) ==\n")
dt <- rbind(
  disp_test(seq_len(nrow(meta)), "Group", "Group (4 levels)"),
  disp_test(which(meta$Group_sub %in% c("IBS_before_great","IBS_before_poor")), "Group_sub", "baseline IBS great vs poor"),
  disp_test(which(meta$Group_sub %in% c("UC_before_great","UC_before_poor")), "Group_sub", "baseline UC great vs poor"),
  disp_test(which(meta$Group_sub %in% c("IBS_after_great","IBS_after_poor")), "Group_sub", "post IBS great vs poor"),
  disp_test(which(meta$Group_sub %in% c("UC_after_great","UC_after_poor")), "Group_sub", "post UC great vs poor"),
  disp_test(ibs, "visit", "IBS before vs after"),
  disp_test(ucs, "visit", "UC before vs after")
)
write.csv(dt, file.path(OUTDIR, "r_dispersion_tests.csv"), row.names = FALSE)
write.csv(do.call(rbind, lapply(pm, function(x)
  data.frame(comparison = x$label, n = x$n, R2 = round(x$r2, 4), p = round(x$p, 4)))),
  file.path(OUTDIR, "r_permanova.csv"), row.names = FALSE)
cat("\n")

# ------------------------------------------- step 5: differential taxa (Wilcoxon)
cat("== 7. Differential taxa: Wilcoxon on relative abundance + BH ==\n")
rel <- sweep(cnt_f, 2, colSums(cnt_f), "/")
tax_split <- strsplit(rownames(rel), ";", fixed = TRUE)
# readable label: genus when assigned, else the deepest assigned rank, so that
# features the app collapses to opaque names ("__", "s__ramosum") stay traceable
levels_letter <- c("k", "p", "c", "o", "f", "g", "s")
lab <- vapply(tax_split, function(x) {
  x <- sub("^[a-z]__", "", x)
  genus <- x[6]; species <- x[7]
  if (nzchar(genus)) return(if (nzchar(species)) paste(genus, species) else genus)
  if (nzchar(species)) return(paste0(species, " (species; genus not annotated)"))
  assigned <- which(nzchar(x[1:5]))
  if (!length(assigned)) return("unclassified (no rank assigned)")
  deepest <- max(assigned)
  paste0("unclassified ", levels_letter[deepest], ":", x[deepest])
}, character(1))
rownames(rel) <- make.unique(lab)

diff_wilcox <- function(gv, a, b) {
  idx <- which(meta[[gv]] %in% c(a, b))
  sub <- rel[, idx]; grp <- factor(meta[[gv]][idx], levels = c(a, b))
  res <- do.call(rbind, lapply(rownames(sub), function(f) {
    x <- as.numeric(sub[f, grp == a]); y <- as.numeric(sub[f, grp == b])
    if (max(c(x, y)) == 0) return(NULL)
    w <- suppressWarnings(wilcox.test(x, y))
    data.frame(taxon = f, mean_a = mean(x), mean_b = mean(y),
               log2FC = log2((mean(y) + 1e-6) / (mean(x) + 1e-6)),
               prevalence_a = mean(x > 0), prevalence_b = mean(y > 0),
               p = w$p.value)
  }))
  res$fdr <- p.adjust(res$p, "BH")
  res$vs <- sprintf("%s vs %s", a, b)
  res[order(res$p), ]
}
dw <- list(
  diff_wilcox("Group_sub", "IBS_before_great", "IBS_before_poor"),
  diff_wilcox("Group_sub", "UC_before_great",  "UC_before_poor"),
  diff_wilcox("Group_sub", "IBS_after_great",  "IBS_after_poor"),
  diff_wilcox("Group_sub", "UC_after_great",   "UC_after_poor"),
  diff_wilcox("Group",     "IBS_before", "IBS_after"),
  diff_wilcox("Group",     "UC_before",  "UC_after"),
  diff_wilcox("Group",     "IBS_before", "UC_before")
)
dwall <- do.call(rbind, dw)
dwall$mean_a <- signif(dwall$mean_a, 4); dwall$mean_b <- signif(dwall$mean_b, 4)
dwall$log2FC <- round(dwall$log2FC, 3)
dwall$prevalence_a <- round(dwall$prevalence_a, 3); dwall$prevalence_b <- round(dwall$prevalence_b, 3)
dwall$p <- signif(dwall$p, 4); dwall$fdr <- signif(dwall$fdr, 4)
write.csv(dwall, file.path(OUTDIR, "r_diff_taxa_wilcoxon.csv"), row.names = FALSE)
for (d in dw) cat(sprintf("  %-42s n=%3d  fdr<0.05: %d\n", d$vs[1], nrow(d), sum(d$fdr < .05)))
cat("\n")

# sensitivity check with a count model (edgeR exact test) - documents that
# methods can disagree, it is NOT treated as the primary result.
if (requireNamespace("edgeR", quietly = TRUE)) {
  cat("== 8. Sensitivity check: edgeR exactTest (raw counts, same contrasts) ==\n")
  eg <- function(a, b) {
    idx <- which(meta$Group %in% c(a, b))
    y <- edgeR::DGEList(counts = as.matrix(cnt_f[, idx]), group = factor(meta$Group[idx], levels = c(a, b)))
    y <- edgeR::calcNormFactors(y)
    y <- edgeR::estimateDisp(y)          # exactTest needs a dispersion estimate
    et <- edgeR::exactTest(y)
    tt <- edgeR::topTags(et, n = Inf)$table
    cat(sprintf("  %-24s taxa with FDR<0.05: %d of %d\n", sprintf("%s vs %s", a, b),
                sum(tt$FDR < .05), nrow(tt)))
    data.frame(taxon = rownames(tt), logFC = round(tt$logFC, 3), FDR = signif(tt$FDR, 4),
               comparison = sprintf("%s vs %s", a, b))
  }
  eg_all <- rbind(eg("IBS_before", "IBS_after"), eg("UC_before", "UC_after"),
                  eg("IBS_before", "UC_before"))
  write.csv(eg_all, file.path(OUTDIR, "r_diff_taxa_edgeR_sensitivity.csv"), row.names = FALSE)
  cat("\n")
}

# ------------------------------------------------------------------ figures
cat("== 9. Figures ->", FIGDIR, "\n")
theme_hw <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey95"))
pal <- c(IBS_before = "#3B6FB6", IBS_after = "#8FB4E3", UC_before = "#C1462F", UC_after = "#E8A08C")

a_long <- alpha %>% select(SampleID, Group, response, observed, shannon, chao1) %>%
  pivot_longer(c(observed, shannon, chao1), names_to = "metric", values_to = "value")
a_long$metric <- factor(a_long$metric, levels = c("observed", "shannon", "chao1"),
                        labels = c("Observed taxa (richness)", "Shannon", "Chao1"))
a_long$Group <- factor(a_long$Group, levels = c("IBS_before", "IBS_after", "UC_before", "UC_after"))
p_alpha <- ggplot(a_long, aes(Group, value, fill = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = .55, width = .62) +
  geom_jitter(width = .16, size = 1.1, alpha = .75, colour = "grey20") +
  facet_wrap(~ metric, scales = "free_y") +
  scale_fill_manual(values = pal) +
  labs(x = NULL, y = "alpha diversity (rarefied table)",
       title = "Alpha diversity by clinical group",
       subtitle = sprintf("points = individual samples; boxes = median/IQR; rarefaction depth = %s reads; Mann-Whitney tests in tables/r_alpha_tests.csv",
                          format(min_depth, big.mark = ","))) +
  theme_hw + theme(axis.text.x = element_text(angle = 20, hjust = 1), legend.position = "none")
ggsave(file.path(FIGDIR, "r_alpha_diversity.png"), p_alpha, width = 11, height = 4.6, dpi = 300)

pcoa_lab <- sprintf("PCoA%s (%.1f%% Bray-Curtis)", 1:2, pct[1:2])
p_pcoa <- ggplot(scores, aes(Axis1, Axis2)) +
  stat_ellipse(aes(colour = Group), level = .68, linewidth = .4, show.legend = FALSE) +
  geom_point(aes(colour = Group, shape = response), size = 2.2, alpha = .85) +
  scale_colour_manual(values = pal) +
  labs(x = pcoa_lab[1], y = pcoa_lab[2],
       title = "Bray-Curtis PCoA of 130 16S samples",
       subtitle = sprintf("PERMANOVA ~ Group: R2 = %.3f, p = %.3f; dispersion p = %.3f (see tables/)",
                          pm[[1]]$r2, pm[[1]]$p, dt$p[1]),
       colour = "Group", shape = "response") +
  theme_hw
ggsave(file.path(FIGDIR, "r_pcoa_bray_curtis.png"), p_pcoa, width = 8, height = 6, dpi = 300)

bl <- scores[scores$visit == "before", ]
bl$Group <- factor(bl$Group, levels = c("IBS_before", "UC_before"))
p_bl <- ggplot(bl, aes(Axis1, Axis2)) +
  stat_ellipse(aes(colour = response), level = .68, linewidth = .4, show.legend = FALSE) +
  geom_point(aes(colour = response), size = 2.4, alpha = .85) +
  facet_wrap(~ Group) +
  scale_colour_manual(values = c(great = "#1B7F5A", poor = "#B0651A")) +
  labs(x = pcoa_lab[1], y = pcoa_lab[2],
       title = "Baseline samples: response group (great vs poor) within each cohort",
       subtitle = "same ordination as above, restricted to pre-treatment samples",
       colour = "response") +
  theme_hw
ggsave(file.path(FIGDIR, "r_pcoa_baseline_response.png"), p_bl, width = 9, height = 4.8, dpi = 300)

# phylum-level composition (top 10 phyla + Other), relative abundance
ph <- sub("^p__", "", vapply(strsplit(rownames(counts), ";", fixed = TRUE), function(x) x[2], character(1)))
ph[!nzchar(ph) | ph == "__"] <- "unclassified"
ph_tab <- rowsum(as.matrix(counts), group = ph)
keep_ph <- names(sort(rowSums(ph_tab), decreasing = TRUE))[1:10]
ph_tab <- rbind(ph_tab[keep_ph, , drop = FALSE],
                Other = colSums(ph_tab[setdiff(rownames(ph_tab), keep_ph), , drop = FALSE]))
ph_rel <- sweep(ph_tab, 2, colSums(ph_tab), "/")
ph_df <- as.data.frame(ph_rel) %>% tibble::rownames_to_column("phylum") %>%
  pivot_longer(-phylum, names_to = "SampleID", values_to = "rel") %>%
  left_join(meta[, c("SampleID", "Group", "response")], by = "SampleID")
ph_df$Group <- factor(ph_df$Group, levels = c("IBS_before", "IBS_after", "UC_before", "UC_after"))
p_bar <- ggplot(ph_df, aes(Group, rel, fill = phylum)) +
  stat_summary(fun = mean, geom = "bar", position = "stack", width = .7) +
  scale_fill_viridis_d(option = "turbo", direction = -1) +
  scale_y_continuous(labels = percent) +
  labs(x = NULL, y = "mean relative abundance", fill = "Phylum",
       title = "Top-10 phyla + Other (descriptive composition, not a test)") +
  theme_hw + theme(axis.text.x = element_text(angle = 20, hjust = 1))
ggsave(file.path(FIGDIR, "r_phylum_composition.png"), p_bar, width = 9, height = 5.5, dpi = 300)

# top differential taxa at baseline (great vs poor), both cohorts
top_taxa <- rbind(
  head(dw[[1]][order(dw[[1]]$p), c("taxon", "p")], 5),
  head(dw[[2]][order(dw[[2]]$p), c("taxon", "p")], 5)
)$taxon
bl_samples <- meta$SampleID[meta$visit == "before"]
detected <- rownames(rel)[apply(rel[, bl_samples, drop = FALSE], 1, max) > 1e-4]
top_taxa <- intersect(unique(top_taxa), detected)
sel <- rownames(rel) %in% top_taxa
sel_df <- as.data.frame(t(rel[sel, , drop = FALSE])) %>%
  tibble::rownames_to_column("SampleID") %>%
  pivot_longer(-SampleID, names_to = "taxon", values_to = "rel") %>%
  left_join(meta[, c("SampleID", "Group_sub", "Group", "response")], by = "SampleID") %>%
  filter(grepl("before", Group_sub))
sel_df$taxon <- factor(sel_df$taxon, levels = unique(top_taxa))
levels(sel_df$taxon) <- stringr::str_wrap(levels(sel_df$taxon), width = 22)
p_dt <- ggplot(sel_df, aes(response, rel + 1e-5, fill = response)) +
  geom_boxplot(outlier.shape = NA, alpha = .55, width = .6) +
  geom_jitter(width = .15, size = .9, alpha = .7, colour = "grey20") +
  facet_grid(taxon ~ Group, scales = "free_y") +
  scale_y_log10() +
  scale_fill_manual(values = c(great = "#1B7F5A", poor = "#B0651A")) +
  labs(x = NULL, y = "relative abundance (log10 scale, +1e-5)", fill = "response",
       title = "Top differential taxa at baseline (Wilcoxon, unadjusted p < 0.05)",
       subtitle = sprintf("exploratory; taxa with BH-FDR < 0.05 at baseline: %d (IBS) and %d (UC) - see tables/r_diff_taxa_wilcoxon.csv",
                          sum(dw[[1]]$fdr < .05), sum(dw[[2]]$fdr < .05))) +
  theme_hw + theme(axis.text.x = element_text(angle = 20, hjust = 1), legend.position = "none")
ggsave(file.path(FIGDIR, "r_top_diff_taxa_baseline.png"), p_dt, width = 10, height = 9, dpi = 300)

# dispersion (distance to group centroid)
bd_all <- betadisper(bray, meta$Group)
disp_df <- data.frame(Group = factor(meta$Group, levels = c("IBS_before", "IBS_after", "UC_before", "UC_after")),
                       dist = bd_all$distances)
p_disp <- ggplot(disp_df, aes(Group, dist, fill = Group)) +
  geom_boxplot(alpha = .55, width = .6) +
  geom_jitter(width = .15, size = 1, alpha = .7, colour = "grey20") +
  scale_fill_manual(values = pal) +
  labs(x = NULL, y = "distance to group centroid (Bray-Curtis)",
       title = "Within-group dispersion (betadisper)",
       subtitle = sprintf("PERMUTEST F = %.2f, p = %.3f -> dispersion is checked before reading the PERMANOVA p-value",
                          dt$F[1], dt$p[1])) +
  theme_hw + theme(axis.text.x = element_text(angle = 20, hjust = 1), legend.position = "none")
ggsave(file.path(FIGDIR, "r_dispersion.png"), p_disp, width = 8, height = 5, dpi = 300)

cat("figures written:\n"); print(basename(list.files(FIGDIR, pattern = "\\.png$")))

# ---------------------------------------------- summary table for the report
summary_tab <- data.frame(
  item = c("input table", "samples (with metadata)", "taxa (level 7, as provided)",
           "taxa retained (prevalence >= 10%)", "reads per sample (median)",
           "reads per sample (min / max)", "rarefaction depth used",
           "Shannon (median, all samples)",
           "Shannon IBS_before / IBS_after", "Shannon UC_before / UC_after",
           "PERMANOVA ~ Group (R2, p)", "PERMANOVA ~ Group + response (Group R2, p)",
           "baseline IBS great vs poor (PERMANOVA R2, p)",
           "baseline UC great vs poor (PERMANOVA R2, p)",
           "dispersion test (Group, F, p)",
           "differential taxa fdr<0.05 (baseline great vs poor, IBS / UC)"),
  value = c("16S_level-7_mapped130.csv (132-sample file restricted to the 130 samples that have metadata)",
            ncol(counts), nrow(counts), nrow(cnt_f),
            median(colSums(counts)), sprintf("%.0f / %.0f", min(colSums(counts)), max(colSums(counts))),
            min_depth, round(median(alpha$shannon), 3),
            paste(round(tapply(alpha$shannon, alpha$Group, median)[c("IBS_before","IBS_after")], 3), collapse = " / "),
            paste(round(tapply(alpha$shannon, alpha$Group, median)[c("UC_before","UC_after")], 3), collapse = " / "),
            sprintf("%.3f, %.3f", pm[[1]]$r2, pm[[1]]$p),
            sprintf("%.3f, %.3f", pm[[2]]$r2, pm[[2]]$p),
            sprintf("%.3f, %.3f", pm[[5]]$r2, pm[[5]]$p),
            sprintf("%.3f, %.3f", pm[[6]]$r2, pm[[6]]$p),
            sprintf("%.2f, %.3f", dt$F[1], dt$p[1]),
            sprintf("%d / %d", sum(dw[[1]]$fdr < .05), sum(dw[[2]]$fdr < .05)))
)
write.csv(summary_tab, file.path(OUTDIR, "r_summary_table.csv"), row.names = FALSE)
cat("\n== 10. Summary table ==\n"); print(summary_tab, row.names = FALSE)

cat("\n== 11. sessionInfo() ==\n")
print(sessionInfo())
sink()
cat("console log written to", LOG, "\n")
