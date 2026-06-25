#!/usr/bin/env Rscript
# length_classification.R — flag partial sequences relative to the longest
# sequence in a per-marker length table.
#
# Provenance: mit_similarity.Rmd lines 429-460. "is_partial" is true when
# `len < 0.90 * max(len)` for the marker.
#
# Usage:
#   length_classification.R --in outputs/qc/<marker>_lengths.tsv \
#                           --out outputs/qc/<marker>_length_classification.tsv \
#                           [--partial-frac 0.9]

suppressPackageStartupMessages({
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
  `in`            = NA_character_,
  out             = NA_character_,
  `partial-frac`  = "0.9"
))
stopifnot(!is.na(opt$`in`), !is.na(opt$out))
PARTIAL_FRAC <- as.numeric(opt[["partial-frac"]])

df <- fread(opt$`in`, sep = "\t", header = FALSE, col.names = c("id", "len"))

max_len <- max(df$len)

df <- df |> mutate(
  pct_of_max    = len / max_len,
  missing_bases = max_len - len,
  is_partial    = pct_of_max < PARTIAL_FRAC
)

dir.create(dirname(opt$out), showWarnings = FALSE, recursive = TRUE)
fwrite(df, opt$out, sep = "\t")

cat(sprintf("[length_classification] %d/%d partial (cutoff < %.2f * max=%d) -> %s\n",
            sum(df$is_partial), nrow(df), PARTIAL_FRAC, max_len, opt$out))
