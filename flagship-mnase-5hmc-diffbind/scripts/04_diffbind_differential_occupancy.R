#!/usr/bin/env Rscript
#
# 04_diffbind_differential_occupancy.R — classify 5hmC peaks by which MNase fraction
# (mono/submono/sonicated) they're enriched in, using DiffBind.
#
# RECONSTRUCTED SCRIPT: the original DiffBind driver for this analysis was never saved to disk
# (only a 504MB .RData workspace and the downstream output tree survive). This script reproduces
# that output structure — one bed file per pairwise fraction comparison, per fold-change
# threshold, split into posFold / negFold / NS (non-significant) — from the surviving sample
# sheet schema and the standard DiffBind workflow. Treat FDR cutoff (0.05) and normalization
# method (DESeq2, library-size) as reasonable DiffBind defaults, not a recovered exact match.
#
# INPUTS:
#   argv[1]   Path to a DiffBind-format sample sheet CSV (see
#             config/diffbind_sampleSheet_template.csv). Required columns: SampleID, Tissue,
#             Factor, Condition, Treatment, Replicate, bamReads, ControlID, bamControl, Peaks,
#             PeakCaller. Treatment must be one of mono/submono/sonicated.
#   argv[2]   Output directory (created if absent)
#
# OUTPUTS (under argv[2]):
#   <comp>_full_report.tsv                          Full DiffBind report (all consensus peaks, unfiltered)
#   fold_<threshold>/posFold_<comp>_F<t>_results.bed  Sites enriched in fraction 2 (Fold >= log2(threshold)), FDR < 0.05
#   fold_<threshold>/negFold_<comp>_F<t>_results.bed  Sites enriched in fraction 1 (Fold <= -log2(threshold)), FDR < 0.05
#   fold_<threshold>/NS_<comp>_F<t>_results.bed       Non-significant sites (FDR >= 0.05), same for both thresholds
#
# Three pairwise comparisons are attempted: comp1 = mono vs submono, comp2 = mono vs sonicated,
# comp3 = submono vs sonicated. Any comparison referencing a fraction absent from the sample
# sheet (e.g. mono, when only the public submono/sonicated demo data is supplied) is skipped
# with a message rather than failing.
#
# REQUIRES: R packages DiffBind, GenomicRanges

suppressMessages(library(DiffBind))
suppressMessages(library(GenomicRanges))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript 04_diffbind_differential_occupancy.R <sampleSheet.csv> <outdir>")
}
sample_sheet_path <- args[1]
outdir <- args[2]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

fold_thresholds <- c(1.5, 2)
fdr_cutoff <- 0.05

contrast_defs <- list(
  comp1 = c("mono", "submono"),
  comp2 = c("mono", "sonicated"),
  comp3 = c("submono", "sonicated")
)

samples <- read.csv(sample_sheet_path, stringsAsFactors = FALSE)
dbObj <- dba(sampleSheet = samples)
dbObj <- dba.count(dbObj, bParallel = FALSE)
dbObj <- dba.normalize(dbObj)

present_fractions <- unique(dbObj$samples$Treatment)
added_contrasts <- character(0)

for (comp_name in names(contrast_defs)) {
  frac1 <- contrast_defs[[comp_name]][1]
  frac2 <- contrast_defs[[comp_name]][2]

  if (!(frac1 %in% present_fractions) || !(frac2 %in% present_fractions)) {
    message(sprintf("[%s] skipping %s vs %s: fraction not present in sample sheet",
                     comp_name, frac1, frac2))
    next
  }

  dbObj <- dba.contrast(dbObj,
                         group1 = dbObj$masks[[frac1]], name1 = frac1,
                         group2 = dbObj$masks[[frac2]], name2 = frac2)
  added_contrasts <- c(added_contrasts, comp_name)
}

if (length(added_contrasts) == 0) {
  stop("No contrasts could be built — sample sheet doesn't contain at least two matching fractions.")
}

dbObj <- dba.analyze(dbObj, method = DBA_DESEQ2)

write_sorted_bed <- function(df, path) {
  if (nrow(df) == 0) {
    file.create(path)
    return(invisible())
  }
  df <- df[order(df$seqnames, df$start), ]
  write.table(df[, c("seqnames", "start", "end", "Fold", "FDR")],
              file = path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

for (i in seq_along(added_contrasts)) {
  comp_name <- added_contrasts[i]

  full <- dba.report(dbObj, contrast = i, th = 1)
  full_df <- as.data.frame(full)

  write.table(full_df, file = file.path(outdir, sprintf("%s_full_report.tsv", comp_name)),
              sep = "\t", quote = FALSE, row.names = FALSE)

  sig <- full_df[!is.na(full_df$FDR) & full_df$FDR < fdr_cutoff, ]
  ns  <- full_df[is.na(full_df$FDR) | full_df$FDR >= fdr_cutoff, ]

  for (ft in fold_thresholds) {
    fold_dir <- file.path(outdir, sprintf("fold_%s", ft))
    dir.create(fold_dir, recursive = TRUE, showWarnings = FALSE)
    tag <- sprintf("F%s", sub("\\.0$", "", as.character(ft)))
    log2_thresh <- log2(ft)

    pos <- sig[sig$Fold >=  log2_thresh, ]
    neg <- sig[sig$Fold <= -log2_thresh, ]

    write_sorted_bed(pos, file.path(fold_dir, sprintf("posFold_%s_%s_results.bed", comp_name, tag)))
    write_sorted_bed(neg, file.path(fold_dir, sprintf("negFold_%s_%s_results.bed", comp_name, tag)))
    write_sorted_bed(ns,  file.path(fold_dir, sprintf("NS_%s_%s_results.bed", comp_name, tag)))

    message(sprintf("[%s] fold %s: %d posFold, %d negFold, %d NS",
                     comp_name, ft, nrow(pos), nrow(neg), nrow(ns)))
  }
}
