## Week 5 Homework 2 figures - base R (+ DESeq2 for the VST transform)
## hw2_pca.png           VST + PCA of all 24 libraries
## hw2_panel_heatmap.png A: induced panel ordered by group, B: 30 most variable genes
## hw2_libsize.png       library size per sample
## Inputs: C:/emp_hw5/RNAseq_output.csv, C:/emp_hw5/RNAseq_mapping.csv

.libPaths(c("D:/R/R-userlib", .libPaths()))
suppressPackageStartupMessages(library(DESeq2))

counts <- as.matrix(read.csv("C:/emp_hw5/RNAseq_output.csv", row.names = 1, check.names = FALSE))
meta <- read.csv("C:/emp_hw5/RNAseq_mapping.csv", stringsAsFactors = FALSE)
stopifnot(identical(colnames(counts), meta$SampleID))
levels_g <- c("DMSO", "DMSO+LIPUS", "T4400", "T4400+LIPUS", "T3976", "T3976+LIPUS")
group <- factor(meta$Group, levels = levels_g)
libsize <- colSums(counts)
cols <- c("DMSO" = "#C0392B", "DMSO+LIPUS" = "#E67E22", "T4400" = "#1F6FB2",
          "T4400+LIPUS" = "#5FA8DC", "T3976" = "#1E8449", "T3976+LIPUS" = "#58B368")

## ---- VST (blind) so that the PCA is not driven by the mean-variance trend
dds0 <- DESeqDataSetFromMatrix(countData = counts,
                               colData = data.frame(row.names = colnames(counts)),
                               design = ~ 1)
vsd <- assay(vst(dds0, blind = TRUE))
cpm <- t(t(counts) / libsize * 1e6)
logcpm <- log2(cpm + 1)

## ---- 1. PCA on the VST matrix
v <- apply(vsd, 1, var)
sel <- order(v, decreasing = TRUE)[1:2000]
pca <- prcomp(t(vsd[sel, ]), center = TRUE)
pv <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

png("hw2_pca.png", width = 2300, height = 1650, res = 200)
par(mar = c(4.5, 4.5, 3.5, 11))
plot(pca$x[, 1], pca$x[, 2], col = cols[group], pch = 19, cex = 1.7,
     xlab = paste0("PC1 (", pv[1], "% variance)"),
     ylab = paste0("PC2 (", pv[2], "% variance)"),
     main = "PCA of 24 RNA-seq libraries (DESeq2 VST, 2000 most variable genes)")
abline(h = 0, v = 0, col = "grey90")
text(pca$x[, 1], pca$x[, 2], labels = rownames(pca$x), pos = 4, cex = 0.72, col = "grey20")
legend("right", inset = c(-0.40, 0), legend = levels_g, pch = 19, col = cols,
       pt.cex = 1.5, bty = "n", title = "Group", cex = 0.9)
mtext(paste("PC1 is driven by the amplitude of the per-group innate-immune response",
            "(right-hand samples) - see heatmap A"), side = 1, line = 2.8, cex = 0.68)
dev.off()
write.csv(round(pca$x[, 1:4], 4), "hw2_pca_scores.csv")

## ---- 2. two-block heatmap -------------------------------------------------
panel <- c("Saa3", "Lcn2", "Mmp13", "Ccl3", "Kng1", "Try5", "Clec4e")
panel <- panel[panel %in% rownames(logcpm)]
means <- sapply(levels_g, function(g) rowMeans(logcpm[panel, group == g, drop = FALSE]))
score <- rowSums(t(logcpm[panel, ]) - t(means[, as.character(group)]))
ord_group <- order(as.numeric(group), score)
zrows <- function(genes) {
  z <- t(scale(t(logcpm[genes, , drop = FALSE])))
  z[z > 2] <- 2; z[z < -2] <- -2
  z
}
zp <- zrows(panel)[, ord_group]
ord_var <- setdiff(rownames(logcpm)[order(v, decreasing = TRUE)][1:30], panel)
mv <- logcpm[ord_var, ]
hc <- hclust(as.dist(1 - cor(mv, method = "spearman")), method = "average")
zv <- zrows(ord_var)[, hc$order]

bwr <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100)
draw_strip <- function(gvec, main, cex.main = 1.05, legend_here = FALSE) {
  n <- length(gvec)
  par(mar = c(0.6, 11, 3.2, 6))
  plot(NA, xlim = c(0.5, n + 0.5), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
  rect(seq_len(n) - 0.5, 0, seq_len(n) + 0.5, 1,
       col = cols[as.character(gvec)], border = NA)
  box(col = "grey70")
  title(main = main, cex.main = cex.main)
  if (legend_here) {
    legend(x = n + 0.8, y = 0.5, legend = levels_g, fill = cols, border = NA,
           bty = "n", cex = 0.78, xpd = NA, title = "Group")
  }
}
draw_heat <- function(zmat, col_labels, row_labels, cex.col = 0.62, cex.row = 0.7,
                      note = NULL) {
  nr <- nrow(zmat); nc <- ncol(zmat)
  par(mar = c(7.5, 11, 0.6, 6))
  ## this R build's image() wants nrow(z) = length(x) and ncol(z) = length(y)
  image(seq_len(nc), seq_len(nr), t(zmat[nr:1, , drop = FALSE]), col = bwr,
        axes = FALSE, xlab = "", ylab = "")
  box(col = "grey70")
  axis(1, at = seq_len(nc), labels = col_labels, las = 2, cex.axis = cex.col, tick = FALSE)
  axis(2, at = seq_len(nr), labels = rev(row_labels), las = 1, cex.axis = cex.row, tick = FALSE)
  if (!is.null(note)) mtext(note, side = 1, line = 6.4, cex = 0.62)
}

png("hw2_panel_heatmap.png", width = 2900, height = 2200, res = 190)
layout(matrix(1:4, nrow = 2, byrow = FALSE), heights = c(1.4, 7))
draw_strip(group[ord_group],
           "A. Induced innate-immune / acute-phase panel - samples ordered by group",
           legend_here = TRUE)
draw_heat(zp[, ord_group], colnames(zp)[ord_group], rownames(zp), cex.row = 0.82,
          note = "row z-score of log2 CPM (clipped at +/-2)")
draw_strip(group[hc$order],
           "B. 30 most variable genes - samples clustered by Spearman correlation")
draw_heat(zv[, hc$order],
          sprintf("%s (%.0fM)", colnames(zv)[hc$order], libsize[hc$order] / 1e6),
          rownames(zv), cex.col = 0.58, cex.row = 0.62,
          note = "library size in brackets; M = million reads")
dev.off()

## ---- 3. library size per sample
png("hw2_libsize.png", width = 2200, height = 1100, res = 190)
par(mar = c(8, 4.5, 3, 1))
bp <- barplot(libsize / 1e6, col = cols[group], border = NA, las = 2, cex.names = 0.7,
              ylab = "library size (million reads)", ylim = c(0, 48),
              main = "Library size per sample (dashed line = 25 M)")
abline(h = 25, col = "red", lty = 2)
legend("topright", legend = levels_g, fill = cols, border = NA, bty = "n", cex = 0.8)
dev.off()

cat("figures written\n")
cat("PCA variance PC1/PC2:", pv[1], "% /", pv[2], "%\n")
cat("panel order:", paste(colnames(zp), collapse = ", "), "\n")
cat("cluster order:", paste(colnames(zv), collapse = ", "), "\n")
