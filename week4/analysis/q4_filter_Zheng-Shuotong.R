# =============================================================================
# Week 4 - Q4 variant prioritization filter
# Zheng Shuotong (student). Base R only (no CRAN/Bioconductor dependency) so the
# script runs on any R >= 4.x installation, including the R 4.4.3 used here.
#
# Input : ../data/variants_q4.tsv  (synthetic teaching table, comment line '#')
# Output: q4_shortlist.tsv, q4_filter_funnel.txt, console audit log
#
# Every threshold below is a decision I made and justify in the report:
#   FILTER == PASS        -> only callers' PASS sites (drop LowQual / FAIL)
#   DP >= 20              -> depth floor: below this, a heterozygous site can be
#                            lost by chance; 20x gives ~>= 10 reads per allele
#   GQ >= 30              -> genotype posterior >= 0.999 for the called genotype
#   AF <= 0.001           -> rare in population (ACMG PM2-style evidence),
#                            using a deliberately conservative ceiling
#   gene symbol present   -> cannot interpret without a gene context
#   impact consequence    -> LoF and splice classes first, missense second tier
# Note: these are reviewer-defensible defaults, not universal rules; a clinical
# pipeline would replace the AF ceiling with ancestry-matched gnomAD filters.
# =============================================================================

data_path <- file.path("..", "data", "variants_q4.tsv")

raw <- read.delim(data_path, comment.char = "#", stringsAsFactors = FALSE,
                  check.names = FALSE, colClasses = "character")
num <- c("POS", "DP", "GQ", "AF")
for (nm in num) raw[[nm]] <- suppressWarnings(as.numeric(raw[[nm]]))

cat("Read", nrow(raw), "variant records from", data_path, "\n")
cat("Columns:", paste(names(raw), collapse = ", "), "\n\n")

# ---- My thresholds ----------------------------------------------------------
min_dp    <- 20
min_gq    <- 30
max_af    <- 0.001
pass_only <- TRUE

lof_csq      <- c("stop_gained", "frameshift_variant",
                  "splice_acceptor_variant", "splice_donor_variant",
                  "start_lost", "stop_lost")
missense_csq <- c("missense_variant")
impact_csq   <- c(lof_csq, missense_csq)

# ---- Stepwise funnel --------------------------------------------------------
funnel <- data.frame(step = character(), kept = integer(), stringsAsFactors = FALSE)
add <- function(label, dt) funnel[nrow(funnel) + 1, ] <<- list(label, nrow(dt))

v <- raw
add("00_all_records", v)

if (pass_only) { v <- v[v$FILTER == "PASS", ];  add("01_FILTER_PASS", v) }
v <- v[v$DP >= min_dp, ];    add(sprintf("02_DP>=%d", min_dp), v)
v <- v[v$GQ >= min_gq, ];    add(sprintf("03_GQ>=%d", min_gq), v)
v <- v[v$AF <= max_af, ];    add(sprintf("04_AF<=%g", max_af), v)
v <- v[v$GENE != "." & !is.na(v$GENE), ]; add("05_gene_symbol_present", v)
v <- v[v$CONSEQUENCE %in% impact_csq, ];  add("06_impact_consequence", v)

write.table(funnel, "q4_filter_funnel.txt", sep = "\t", row.names = FALSE, quote = FALSE)
cat("=== filter funnel ===\n"); print(funnel, row.names = FALSE)

# ---- Tiering of survivors --------------------------------------------------
tier <- function(csq) ifelse(csq %in% lof_csq, "tier1_putative_LoF_splice", "tier2_missense")
clin_rank <- c(Pathogenic = 1, Likely_pathogenic = 2,
               Conflicting_interpretations_of_pathogenicity = 3,
               Uncertain_significance = 4, Likely_benign = 5, Benign = 6)

v$TIER      <- tier(v$CONSEQUENCE)
v$CLIN_RANK <- ifelse(v$CLINVAR_SIG %in% names(clin_rank), clin_rank[v$CLINVAR_SIG], 9)
v <- v[order(v$TIER, v$CLIN_RANK, v$AF), ]

shortlist <- v[, c("CHROM", "POS", "REF", "ALT", "GENE", "CONSEQUENCE", "DP",
                   "GQ", "AF", "CLINVAR_SIG", "TIER", "NOTE")]
write.table(shortlist, "q4_shortlist.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n=== prioritised shortlist (tier -> ClinVar rank -> AF) ===\n")
print(shortlist, row.names = FALSE)

cat("\n=== excluded sites and why ===\n")
exc <- raw[!(paste(raw$CHROM, raw$POS) %in% paste(v$CHROM, v$POS)), ]
exc$reason <- with(exc, ifelse(FILTER != "PASS", paste0("FILTER=", FILTER),
                        ifelse(DP < min_dp, paste0("DP=", DP, " < ", min_dp),
                        ifelse(GQ < min_gq, paste0("GQ=", GQ, " < ", min_gq),
                        ifelse(AF > max_af, paste0("AF=", AF, " > ", max_af),
                        ifelse(GENE == ".", "no gene symbol",
                               paste0("consequence=", CONSEQUENCE)))))))
print(exc[, c("CHROM", "POS", "GENE", "CONSEQUENCE", "reason")], row.names = FALSE)

cat("\nNOTE: coordinate/build consistency and ClinVar-ID validity are NOT\n")
cat("checked by this script; they were audited separately against Ensembl\n")
cat("(GRCh38 vs GRCh37) and ClinVar before the ranking was accepted.\n")
