# ================================================================
# Week 5 Homework 1 - Bulk RNA-seq differential expression with DESeq2
# Student : 郑烁曈  SUAT24000155
# Course  : Bioinformatics: From Multi-Omics Data to Discovery
# Inputs  : Week5_Homework_Count_Matrix.csv (1,000 genes x 12 samples)
#           Week5_Homework_Sample_Metadata.csv (condition + batch)
# Outputs : week5_deseq2_results.csv, week5_pca.png, week5_de_plot.png,
#           week5_deseq2_object.rds, session_info.txt,
#           week5_truth_key_check.txt (supporting verification),
#           week5_run_console.txt (console transcript, written by the shell)
#
# This script follows the required workflow 1-12 of the assignment sheet.
# The starter file Week5_Homework_Starter.R was used as the scaffold; every
# "TODO" it left is answered in place and marked with TODO-ANSWER.
# ================================================================

# ---------------------------------------------------------------
# 0. Working directory + library path (robust when run via Rscript)
# ---------------------------------------------------------------
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = FALSE))
} else if (!is.null(sys.frames()[[1]]$ofile)) {
  dirname(normalizePath(sys.frames()[[1]]$ofile, winslash = "/", mustWork = FALSE))
} else {
  getwd()
}
setwd(script_dir)
cat("Working directory:", getwd(), "\n")

# Packages live in a user library on this machine (D:/R/R-userlib).
if (!requireNamespace("DESeq2", quietly = TRUE)) {
  .libPaths(c("D:/R/R-userlib", .libPaths()))
}
suppressPackageStartupMessages({
  library(DESeq2)   # testing + shrinkage
  library(apeglm)   # required shrinkage estimator
  library(ggplot2)  # figures
  library(ggrepel)  # non-overlapping point labels
  library(dplyr)    # table handling
  library(tibble)   # rownames_to_column()
  library(Matrix)   # rankMatrix() for the full-rank check
})
cat("DESeq2:", as.character(packageVersion("DESeq2")),
    "| apeglm:", as.character(packageVersion("apeglm")),
    "| R:", R.version.string, "\n\n")

count_file    <- "Week5_Homework_Count_Matrix.csv"
metadata_file <- "Week5_Homework_Sample_Metadata.csv"

# ---------------------------------------------------------------
# 1. Import (counts stay as raw integers - no TPM/CPM/z-score)
# ---------------------------------------------------------------
counts <- read.csv(count_file, row.names = 1, check.names = FALSE)
coldata <- read.csv(metadata_file, row.names = 1, check.names = FALSE)

cat("Count matrix :", nrow(counts), "genes x", ncol(counts), "samples\n")
cat("Metadata     :", nrow(coldata), "samples x", ncol(coldata), "variables\n\n")

# ---------------------------------------------------------------
# 2. Mandatory validation (stopifnot: fail loudly, never silently)
# ---------------------------------------------------------------
stopifnot(ncol(counts) == nrow(coldata))                            # same sample count
stopifnot(identical(colnames(counts), rownames(coldata)))           # same IDs, same order
stopifnot(all(counts >= 0))                                          # non-negative
stopifnot(all(as.matrix(counts) == round(as.matrix(counts))))        # integers

# Extra checks demanded by the verification checklist
stopifnot(!any(duplicated(colnames(counts))))                        # no duplicated sample IDs
stopifnot(!any(duplicated(rownames(counts))))                        # no duplicated genes
stopifnot(all(is.finite(as.matrix(counts))))                         # no NA/Inf

coldata$condition <- relevel(factor(coldata$condition), ref = "control")
coldata$batch     <- factor(coldata$batch)
stopifnot(levels(coldata$condition)[1] == "control")                 # control is reference

cat("--- verification ---\n")
cat("samples identical to metadata row names, same order : TRUE\n")
cat("all counts non-negative integers                  : TRUE\n")
cat("duplicated sample IDs                             : ", sum(duplicated(colnames(counts))), "\n")
cat("condition levels (reference first)                : ", paste(levels(coldata$condition), collapse = ", "), "\n")
cat("batch levels                                      : ", paste(levels(coldata$batch), collapse = ", "), "\n")
cat("library sizes (M reads)                           : ",
    paste(round(colSums(counts) / 1e6, 3), collapse = ", "), "\n")
cat("design balance (batch x condition):\n")
print(table(coldata$batch, coldata$condition))
cat("\n")

# ---------------------------------------------------------------
# 3. Construct the DESeq2 object
# ---------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = coldata,
  design    = ~ batch + condition
)

# TODO-ANSWER (why batch belongs in the design):
# The three replicate batches A/B/C enter as a blocking factor. Each batch is a
# separate group of samples processed together, so it contributes its own shift
# in expression (library size + handling); the model estimates that shift as a
# nuisance coefficient and removes it, so the condition estimate is compared
# *within* batch instead of being confounded with it. The design is balanced
# (2 control + 2 treated per batch), therefore batch and condition are not
# collinear and the model matrix stays full rank. Section 12 quantifies what
# happens when batch is left out.
mm <- model.matrix(design(dds), colData(dds))
cat("model matrix columns :", paste(colnames(mm), collapse = ", "), "\n")
cat("model matrix rank    :", as.integer(rankMatrix(mm)[1]), "of", ncol(mm), "columns ->",
    if (rankMatrix(mm)[1] == ncol(mm)) "FULL RANK" else "RANK DEFICIENT", "\n\n")
stopifnot(as.integer(rankMatrix(mm)[1]) == ncol(mm))

# ---------------------------------------------------------------
# 4. Pre-filter: >= 10 counts in >= 3 samples
# ---------------------------------------------------------------
genes_before <- nrow(dds)
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]
cat("filter rule          : rowSums(counts >= 10) >= 3\n")
cat("genes before filtering:", genes_before, "\n")
cat("genes after filtering :", nrow(dds), "\n\n")

# ---------------------------------------------------------------
# 5. Fit the model and inspect the coefficient names
# ---------------------------------------------------------------
dds <- DESeq(dds)
coef_names <- resultsNames(dds)
cat("resultsNames(dds):\n"); print(coef_names)

# TODO-ANSWER (confirm the coefficient name instead of assuming it):
# resultsNames() is printed above. With condition releveled to "control" the
# treated-versus-control coefficient is literally "condition_treated_vs_control"
# (the other two entries are Intercept and the batch_B / batch_C block effects).
target_coef <- "condition_treated_vs_control"
if (!target_coef %in% coef_names) {
  stop("Expected coefficient was not found. Inspect resultsNames(dds) and update target_coef.")
}
cat("\nconfirmed coefficient:", target_coef,
    "->", if (grepl("_vs_", target_coef)) "treated versus control (positive = up in treated)" else "?",
    "\n\n")

# ---------------------------------------------------------------
# 6. Extract the contrast and shrink with apeglm
# ---------------------------------------------------------------
res <- results(dds, contrast = c("condition", "treated", "control"), alpha = 0.05)
res_shrunk <- lfcShrink(dds, coef = target_coef, type = "apeglm", res = res)

res_df <- as.data.frame(res_shrunk) |>
  rownames_to_column("gene_id") |>
  mutate(
    significant = !is.na(padj) & padj < 0.05 & abs(log2FoldChange) >= 1,
    direction = case_when(
      significant & log2FoldChange > 0 ~ "Up in treated",
      significant & log2FoldChange < 0 ~ "Down in treated",
      TRUE ~ "Not significant"
    )
  ) |>
  arrange(padj)

# ---------------------------------------------------------------
# 7. Export the complete shrunken result table (ALL genes kept)
# ---------------------------------------------------------------
write.csv(res_df, "week5_deseq2_results.csv", row.names = FALSE)

n_sig     <- sum(res_df$significant)
n_up      <- sum(res_df$direction == "Up in treated")
n_down    <- sum(res_df$direction == "Down in treated")
n_padj    <- sum(!is.na(res_df$padj) & res_df$padj < 0.05)
cat("thresholds      : padj < 0.05  AND  |shrunken log2FC| >= 1\n")
cat("genes tested (non-NA padj) :", sum(!is.na(res_df$padj)), "of", nrow(res_df), "\n")
cat("padj < 0.05 only           :", n_padj, "\n")
cat("significant (padj + effect):", n_sig, " (up in treated:", n_up, ", down in treated:", n_down, ")\n")
cat("complete table written to week5_deseq2_results.csv (", nrow(res_df), " rows )\n\n")

# ---------------------------------------------------------------
# 8. PCA on VST-transformed counts
# ---------------------------------------------------------------
# nsub defaults to 1000, which is more than the 989 genes that survive the
# filter; using all of them for the dispersion-mean fit is the recommended
# fallback (equivalent to varianceStabilizingTransformation) and keeps VST.
vsd <- vst(dds, blind = FALSE, nsub = nrow(dds))
pca_df <- plotPCA(vsd, intgroup = c("condition", "batch"), returnData = TRUE)
percent_var <- round(100 * attr(pca_df, "percentVar"))

# how much of each PC is explained by condition vs by batch (linear model R^2)
r2 <- function(y, g) summary(lm(y ~ g))$r.squared
cat("variance explained by the design (R^2 of lm on each PC):\n")
# plotPCA(returnData = TRUE) returns the first two principal components only.
for (pc in c("PC1", "PC2")) {
  cat(sprintf("  %s: condition %.2f | batch %.2f\n", pc,
              r2(pca_df[[pc]], pca_df$condition), r2(pca_df[[pc]], pca_df$batch)))
}
cat("\n")

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = condition, shape = batch, label = name)) +
  geom_point(size = 4) +
  geom_text_repel(size = 3, max.overlaps = Inf, seed = 1) +
  labs(
    title = "Week 5 RNA-seq PCA (VST, blind = FALSE)",
    subtitle = "colour = condition, shape = batch",
    x = paste0("PC1: ", percent_var[1], "% variance"),
    y = paste0("PC2: ", percent_var[2], "% variance")
  ) +
  theme_bw(base_size = 12)

ggsave("week5_pca.png", p_pca, width = 7, height = 5, dpi = 300)

# ---------------------------------------------------------------
# 9. Volcano plot (axes + thresholds labelled)
# ---------------------------------------------------------------
plot_df <- res_df |>
  mutate(neg_log10_padj = -log10(pmax(padj, 1e-300)))

p_volcano <- ggplot(plot_df, aes(x = log2FoldChange, y = neg_log10_padj, color = direction)) +
  geom_point(alpha = 0.7, size = 1.8) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  scale_color_manual(values = c("Up in treated" = "#C0392B",
                                "Down in treated" = "#2F6DB3",
                                "Not significant" = "grey70")) +
  labs(
    title = "Differential expression: treated versus control",
    subtitle = paste0("dashed lines: |log2FC| = 1 and padj = 0.05 ; ",
                      n_sig, " significant genes (", n_up, " up / ", n_down, " down)"),
    x = "Shrunken log2 fold change (apeglm)",
    y = "-log10 adjusted p value",
    color = NULL
  ) +
  theme_bw(base_size = 12)

ggsave("week5_de_plot.png", p_volcano, width = 7, height = 5, dpi = 300)

# ---------------------------------------------------------------
# 10. Save the fitted object and sessionInfo()
# ---------------------------------------------------------------
saveRDS(dds, "week5_deseq2_object.rds")
capture.output(sessionInfo(), file = "session_info.txt")
cat("saved: week5_deseq2_object.rds, session_info.txt\n\n")

# ---------------------------------------------------------------
# 11. Supporting diagnostics (not part of the required outputs)
#     (a) how many genes are called when batch is dropped
#     (b) agreement with the instructor truth log2 fold changes
# ---------------------------------------------------------------
no_batch <- DESeqDataSetFromMatrix(countData = counts, colData = coldata,
                                   design = ~ condition)
no_batch <- no_batch[rowSums(counts(no_batch) >= 10) >= 3, ]
no_batch <- DESeq(no_batch)
nb_shrunk <- lfcShrink(no_batch, coef = "condition_treated_vs_control", type = "apeglm")
nb_sig <- sum(!is.na(nb_shrunk$padj) & nb_shrunk$padj < 0.05 & abs(nb_shrunk$log2FoldChange) >= 1)
nb_sig_padj <- sum(!is.na(nb_shrunk$padj) & nb_shrunk$padj < 0.05)

unshrunk_sig <- sum(!is.na(res$padj) & res$padj < 0.05 & abs(res$log2FoldChange) >= 1)

key <- read.csv("Week5_Homework_Gene_Annotation_Instructor_Key.csv",
                stringsAsFactors = FALSE)
key$truth <- as.numeric(key$truth_log2FC_for_instructor)
cmp <- merge(res_df, key[, c("gene_id", "truth")], by = "gene_id")
truth_de <- abs(cmp$truth) >= 1
called   <- cmp$significant
agree_pearson  <- cor(cmp$log2FoldChange, cmp$truth, method = "pearson")
agree_spearman <- cor(cmp$log2FoldChange, cmp$truth, method = "spearman")
mae <- mean(abs(cmp$log2FoldChange - cmp$truth))

lines <- c(
  "Week 5 supporting verification - Week5_Homework_Gene_Annotation_Instructor_Key.csv",
  paste0("generated: ", format(Sys.time())),
  "",
  "A. Effect of including batch in the design (~ batch + condition vs ~ condition)",
  paste0("   significant genes (~ batch + condition): ", n_sig,
         "   [padj<0.05 and |shrunken log2FC|>=1]"),
  paste0("   significant genes (~ condition only)    : ", nb_sig),
  paste0("   padj<0.05 only      (~ batch + condition): ", n_padj),
  paste0("   padj<0.05 only      (~ condition only)  : ", nb_sig_padj),
  "",
  "B. Shrunken vs unshrunken effect-size threshold",
  paste0("   significant with apeglm-shrunken log2FC : ", n_sig),
  paste0("   significant with unshrunken log2FC      : ", unshrunk_sig),
  "",
  paste0("C. Agreement with the instructor truth log2 fold changes (all ",
         nrow(cmp), " genes)"),
  paste0("   Pearson  r = ", round(agree_pearson, 4)),
  paste0("   Spearman r = ", round(agree_spearman, 4)),
  paste0("   mean |estimate - truth| = ", round(mae, 4), " log2 units"),
  paste0("   truth-defined DE genes (|truth log2FC| >= 1): ", sum(truth_de)),
  paste0("   of those recovered by the pipeline          : ", sum(truth_de & called)),
  paste0("   called significant but not truth-DE         : ", sum(!truth_de & called)),
  "",
  "Note: the key is an instructor file shipped with the weekly folder; it is used",
  "here only as an independent check on the estimates, never as model input."
)
writeLines(lines, "week5_truth_key_check.txt")
cat(paste(lines, collapse = "\n"), "\n\n")

# ---------------------------------------------------------------
# 12. Interpretation input - numbers quoted in week5_interpretation.md
# ---------------------------------------------------------------
cat("=== numbers for the 100-150 word interpretation ===\n")
cat("comparison        : condition treated vs control (batch-adjusted)\n")
cat("PC1 / PC2 variance:", percent_var[1], "% /", percent_var[2], "%\n")
cat("significant genes :", n_sig, "(", n_up, "up,", n_down, "down )\n")
top <- head(res_df[res_df$significant, ], 5)
cat("top 5 significant genes:\n"); print(top[, c("gene_id", "log2FoldChange", "padj")])
cat("\nDONE\n")
