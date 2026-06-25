#!/usr/bin/env Rscript
# inspect_imwong_position.R — drill-down at the 18S alignment positions where
# the Imwong et al. published forward primers (PlasmoM_N1F, PlasmoM_N2F)
# bind. Settles whether the published forwards fail our entropy/coverage
# cutoff by a hair or by a mile.
#
# Provenance: Phase 5 Step 1 audit (HANDOFF.md). Chapter 6 of the v1 book
# located the published forwards at 18S aln 1683-1702 (N1F) and 1718-1740
# (N2F); collectively the forward-primer region is 1683-1740.
#
# Outputs a per-position TSV with: pos, coverage, entropy, snp_count, plus a
# yes/no flag for whether each of a set of (cov_thr, entropy_thr) pairs
# would have admitted that position into a candidate window.
#
# Usage:
#   inspect_imwong_position.R \
#     --alignment outputs/alignment/ma_18S.target.fasta \
#     --perpos    outputs/entropy/18S.per_position.tsv \
#     --windows   outputs/entropy/18S.windows.tsv \
#     --primers   "PlasmoM_N1F=ATGGCCGTTTTTAGTTCGTG,PlasmoM_N2F=GTTAATTCCGATAACGAACGAGA" \
#     --thresholds "0.90:0.20,0.85:0.25,0.85:0.30,0.80:0.30,0.80:0.40,0.70:0.50" \
#     --out       outputs/audit/imwong_forward_position.tsv

suppressPackageStartupMessages({
  library(Biostrings)
  library(tidyverse)
  library(data.table)
})

parse_args <- function(known) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) %% 2 != 0) stop("expected --flag value pairs")
  out <- known
  for (i in seq(1, length(args), by = 2)) {
    k <- sub("^--", "", args[i])
    if (!k %in% names(known)) stop(sprintf("unknown flag: %s", args[i]))
    out[[k]] <- args[i + 1]
  }
  out
}

opt <- parse_args(list(
  alignment    = NA_character_,
  perpos       = NA_character_,
  windows      = NA_character_,
  primers      = NA_character_,
  thresholds   = NA_character_,
  out          = NA_character_
))
for (k in names(opt)) if (is.na(opt[[k]])) stop(sprintf("--%s required", k))

# ---- re-derive primer alignment positions ----
aln    <- readDNAStringSet(opt$alignment)
unaln  <- DNAStringSet(gsub("-", "", as.character(aln)))
names(unaln) <- names(aln)

# parse "Name1=SEQ1,Name2=SEQ2"
prim_pairs <- strsplit(opt$primers, ",")[[1]]
primers <- tibble(
  name = sub("=.*$", "", prim_pairs),
  seq  = sub("^.*=", "", prim_pairs)
)

map_one <- function(name, seq) {
  p <- DNAString(seq)
  rows <- list()
  for (j in seq_along(unaln)) {
    fwd <- matchPattern(p, unaln[[j]], max.mismatch = 1)
    rev <- matchPattern(reverseComplement(p), unaln[[j]], max.mismatch = 1)
    strand <- if (length(fwd) > 0) "+" else if (length(rev) > 0) "-" else NA_character_
    if (is.na(strand)) next
    m <- if (strand == "+") fwd[1] else rev[1]
    letters <- strsplit(as.character(aln[[j]]), "")[[1]]
    aln_pos <- which(letters != "-")
    rows[[length(rows) + 1]] <- tibble(
      primer = name, ref = names(aln)[j], strand = strand,
      aln_start = aln_pos[start(m)], aln_end = aln_pos[end(m)]
    )
  }
  bind_rows(rows)
}

primer_positions <- pmap_dfr(primers, function(name, seq) map_one(name, seq))

# Region of interest = union of primer spans (use min start, max end).
roi_start <- min(primer_positions$aln_start, na.rm = TRUE)
roi_end   <- max(primer_positions$aln_end,   na.rm = TRUE)
cat(sprintf("[imwong] ROI: 18S aln %d..%d (%d bp)\n",
            roi_start, roi_end, roi_end - roi_start + 1))

# ---- pull per-position metrics for the ROI ----
pp <- fread(opt$perpos, sep = "\t")
roi <- pp[pos >= roi_start & pos <= roi_end]

# Pipeline windows that cover any pos in ROI (regardless of whether they're
# candidate-accepted).
win <- fread(opt$windows, sep = "\t")
overlapping_windows <- win[start <= roi_end & end >= roi_start]

# Pipeline candidate windows that cover any pos in ROI
pc_path <- file.path(dirname(opt$perpos),
                     sub("per_position", "primer_candidates", basename(opt$perpos)))
if (file.exists(pc_path)) {
  pc <- fread(pc_path, sep = "\t")
  pc_in_roi <- pc[start <= roi_end & end >= roi_start]
} else {
  pc_in_roi <- data.table()
}

# ---- per-threshold admission test ----
# Each (cov_thr, ent_thr) pair: would this position have been admitted to a
# candidate window if used as the threshold? Position-level proxy: position
# coverage ≥ cov_thr AND position entropy ≤ ent_thr. The pipeline's actual
# decision is window-level (mean over a 25-bp window), but the per-position
# flag here is what matters for whether *the position* survives.
thr_pairs <- strsplit(opt$thresholds, ",")[[1]]
thr_df <- tibble(
  pair  = thr_pairs,
  cov_thr = as.numeric(sub("^([0-9.]+):.*$", "\\1", thr_pairs)),
  ent_thr = as.numeric(sub("^.*:([0-9.]+)$", "\\1", thr_pairs))
)

admit_cols <- map(seq_len(nrow(thr_df)), function(i) {
  col_name <- sprintf("admit_%s", gsub("[.:]", "_", thr_df$pair[i]))
  v <- roi$coverage >= thr_df$cov_thr[i] & roi$entropy <= thr_df$ent_thr[i]
  tibble(!!col_name := v)
})
admit_tbl <- bind_cols(admit_cols)

# Distance from each position to the nearest window edge (window space)
# Useful for the "by a hair" verdict.
nearest_window_dist <- function(p) {
  if (nrow(overlapping_windows) == 0) return(c(NA, NA))
  # window starts before this position
  before <- overlapping_windows[end < p]
  after  <- overlapping_windows[start > p]
  d_before <- if (nrow(before) > 0) p - max(before$end)   else NA_integer_
  d_after  <- if (nrow(after)  > 0) min(after$start) - p  else NA_integer_
  c(d_before, d_after)
}

dist_v <- t(sapply(roi$pos, nearest_window_dist))
roi$dist_to_window_left  <- dist_v[, 1]
roi$dist_to_window_right <- dist_v[, 2]

# In-candidate-window flag
in_candidate <- map_lgl(roi$pos, function(p) {
  any(pc_in_roi$start <= p & pc_in_roi$end >= p)
})
roi$in_candidate_window <- in_candidate

# Bind admission columns
roi_out <- bind_cols(roi, admit_tbl)
fwrite(roi_out, opt$out, sep = "\t")

# Console summary
cat("\n[imwong] Threshold-admission summary (positions in ROI):\n")
n_pos <- nrow(roi)
for (col in names(admit_tbl)) {
  admitted <- sum(admit_tbl[[col]])
  cat(sprintf("  %s : %d / %d positions admitted (%.1f%%)\n",
              col, admitted, n_pos, 100 * admitted / n_pos))
}
cat(sprintf("\n[imwong] Per-position metrics in ROI:\n"))
cat(sprintf("  mean entropy  : %.4f\n", mean(roi$entropy)))
cat(sprintf("  median entropy: %.4f\n", median(roi$entropy)))
cat(sprintf("  max entropy   : %.4f\n", max(roi$entropy)))
cat(sprintf("  mean coverage : %.4f\n", mean(roi$coverage)))
cat(sprintf("  min coverage  : %.4f\n", min(roi$coverage)))
cat(sprintf("  positions in any candidate window: %d / %d\n",
            sum(roi$in_candidate_window), n_pos))
