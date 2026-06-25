#!/usr/bin/env Rscript
# cross_marker_resolution.R — build the v1 deliverable: pair × {MIT, 18S,
# 18S-loose, combined} resolution table. Rows are every unordered pair of
# the canonical target species.
#
# Provenance: NEW (Phase 4 acceptance criterion 5 in HANDOFF.md).
#
# Resolution semantics for one marker on one pair (A, B):
#   - NA    if either species has zero references for the marker.
#   - FALSE if the pair appears in the marker's strict suspicious table.
#   - TRUE  otherwise.
#
# `combined_resolved`:
#   TRUE if mit_resolved is TRUE OR 18S_resolved is TRUE.
#   FALSE if both are FALSE.
#   NA if both are NA.
#   (Mixed TRUE/NA collapses to TRUE; mixed FALSE/NA collapses to FALSE.)
#
# Usage:
#   cross_marker_resolution.R \
#     --species-targets falciparum,vivax,...,simium \
#     --mit-suspicious   outputs/cross_species/mit_suspicious.tsv \
#     --18s-suspicious   outputs/cross_species/18S_suspicious.tsv \
#     --18s-suspicious-loose outputs/cross_species/18S_suspicious_loose.tsv \
#     --mit-lengths      outputs/qc/mit_lengths.tsv \
#     --18s-coverage     outputs/cross_species/18S_species_coverage.tsv \
#     --out              outputs/cross_species/resolution_table.tsv

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
  `species-targets`         = NA_character_,
  `mit-suspicious`          = NA_character_,
  `18s-suspicious`          = NA_character_,
  `18s-suspicious-loose`    = NA_character_,
  `mit-lengths`             = NA_character_,
  `18s-coverage`            = NA_character_,
  out                       = NA_character_
))
for (k in names(opt)) if (is.na(opt[[k]])) stop(sprintf("--%s required", k))

# --- canonicalise species names ---
canon <- function(x) {
  x |> str_to_lower() |> str_replace("^p\\.?", "") |>
       str_replace("^cynomologi$", "cynomolgi")
}
unord_pair <- function(a, b) {
  paste(sort(c(canon(a), canon(b))), collapse = " vs ")
}

targets <- strsplit(opt[["species-targets"]], ",")[[1]] |> str_trim() |> canon()

# --- MIT references per species (from mit_lengths.tsv ID prefix) ---
mit_lens <- fread(opt[["mit-lengths"]], sep = "\t", header = FALSE,
                  col.names = c("id", "len"))
mit_counts <- mit_lens |>
  mutate(species = canon(str_extract(id, "^P[A-Za-z]+"))) |>
  count(species, name = "n_sequences")
mit_refs <- setNames(rep(0L, length(targets)), targets)
mit_refs[mit_counts$species] <- mit_counts$n_sequences

# --- 18S references per species (from species_coverage TSV) ---
s18_cov <- fread(opt[["18s-coverage"]], sep = "\t")
s18_refs <- setNames(rep(0L, length(targets)), targets)
s18_refs[canon(s18_cov$species)] <- s18_cov$n_sequences

# --- pair sets from suspicious tables (canonicalised) ---
load_pairs <- function(path) {
  d <- tryCatch(fread(path, sep = "\t"), error = function(e) data.table())
  if (nrow(d) == 0) return(character(0))
  if (!all(c("query_species", "subject_species") %in% names(d))) return(character(0))
  unique(map2_chr(canon(d$query_species), canon(d$subject_species), unord_pair))
}
mit_suspicious_pairs       <- load_pairs(opt[["mit-suspicious"]])
s18_suspicious_pairs       <- load_pairs(opt[["18s-suspicious"]])
s18_suspicious_loose_pairs <- load_pairs(opt[["18s-suspicious-loose"]])

# --- build the unordered pair grid ---
pairs <- expand.grid(a = targets, b = targets, stringsAsFactors = FALSE) |>
  filter(a < b) |>
  mutate(pair = paste(a, "vs", b))

resolve <- function(pair, sus_pairs, refs_a, refs_b) {
  if (refs_a == 0 || refs_b == 0) return(NA)
  !(pair %in% sus_pairs)
}
combine <- function(m, s) {
  if (is.na(m) && is.na(s)) return(NA)
  isTRUE(m) || isTRUE(s)
}

# Borderline notes — known facts that downstream readers want at-a-glance.
note_for <- function(a, b, mit_resolved, s18_refs_a, s18_refs_b) {
  notes <- character(0)
  if (identical(sort(c(a, b)), c("simium", "vivax")) && isFALSE(mit_resolved)) {
    notes <- c(notes,
      "MIT vivax/simium ambiguity (P. simium recently zoonotic from P. vivax)")
  }
  if (s18_refs_a <= 1 || s18_refs_b <= 1) {
    notes <- c(notes,
      sprintf("18S coverage low: %s=%d, %s=%d", a, s18_refs_a, b, s18_refs_b))
  }
  paste(notes, collapse = "; ")
}

res <- pairs |>
  rowwise() |>
  mutate(
    mit_resolved        = resolve(pair, mit_suspicious_pairs,
                                  mit_refs[[a]], mit_refs[[b]]),
    `18S_resolved`      = resolve(pair, s18_suspicious_pairs,
                                  s18_refs[[a]], s18_refs[[b]]),
    `18S_resolved_loose`= resolve(pair, s18_suspicious_loose_pairs,
                                  s18_refs[[a]], s18_refs[[b]]),
    combined_resolved   = combine(mit_resolved, `18S_resolved`),
    notes               = note_for(a, b, mit_resolved,
                                   s18_refs[[a]], s18_refs[[b]])
  ) |>
  ungroup() |>
  select(pair, mit_resolved, `18S_resolved`, `18S_resolved_loose`,
         combined_resolved, notes)

dir.create(dirname(opt$out), showWarnings = FALSE, recursive = TRUE)
fwrite(res, opt$out, sep = "\t")

# Summary that the report and the human running this can sanity-check.
n           <- nrow(res)
n_mit_TRUE  <- sum(res$mit_resolved        %in% TRUE)
n_mit_FALSE <- sum(res$mit_resolved        %in% FALSE)
n_mit_NA    <- sum(is.na(res$mit_resolved))
n_18_TRUE   <- sum(res$`18S_resolved`      %in% TRUE)
n_18_FALSE  <- sum(res$`18S_resolved`      %in% FALSE)
n_18_NA     <- sum(is.na(res$`18S_resolved`))
n_comb_TRUE <- sum(res$combined_resolved   %in% TRUE)
n_comb_FALSE<- sum(res$combined_resolved   %in% FALSE)
n_comb_NA   <- sum(is.na(res$combined_resolved))

cat(sprintf("[cross_marker_resolution] %d total pairs\n", n))
cat(sprintf("  mit_resolved:       T=%d  F=%d  NA=%d\n", n_mit_TRUE, n_mit_FALSE, n_mit_NA))
cat(sprintf("  18S_resolved:       T=%d  F=%d  NA=%d\n", n_18_TRUE,  n_18_FALSE,  n_18_NA))
cat(sprintf("  combined_resolved:  T=%d  F=%d  NA=%d\n", n_comb_TRUE,n_comb_FALSE,n_comb_NA))
cat(sprintf("-> %s\n", opt$out))
