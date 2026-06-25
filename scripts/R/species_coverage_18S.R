#!/usr/bin/env Rscript
# species_coverage_18S.R — for each canonical target species, count how many
# 18S references in the panel map to it. Surfaces species the assay's 18S
# leg cannot resolve at all (n_sequences = 0).
#
# Provenance: NEW (Phase 4 acceptance criterion 4 in HANDOFF.md).
#
# Usage:
#   species_coverage_18S.R --acc-map outputs/cross_species/18S_acc_to_header.tsv \
#     --species-targets falciparum,vivax,knowlesi,malariae,ovale,coatneyi,inui,fieldi,cynomolgi,simiovale,simium \
#     --out outputs/cross_species/18S_species_coverage.tsv

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
  `acc-map`         = NA_character_,
  `species-targets` = NA_character_,
  out               = NA_character_
))
stopifnot(!is.na(opt[["acc-map"]]),
          !is.na(opt[["species-targets"]]),
          !is.na(opt$out))

targets <- strsplit(opt[["species-targets"]], ",")[[1]] |> str_trim() |> str_to_lower()

m <- fread(opt[["acc-map"]], sep = "\t", header = TRUE,
           select = c("acc", "species"))
# species column is "p.<name>"; strip the prefix to canonicalise.
m <- m |> mutate(canon = str_to_lower(str_replace(species, "^p\\.?", "")))

cov <- tibble(species = targets) |>
  rowwise() |>
  mutate(n_sequences = sum(m$canon == species, na.rm = TRUE)) |>
  ungroup() |>
  mutate(is_covered = n_sequences >= 1)

dir.create(dirname(opt$out), showWarnings = FALSE, recursive = TRUE)
fwrite(cov, opt$out, sep = "\t")
cat(sprintf("[species_coverage_18S] %d/%d targets covered (≥1 ref) -> %s\n",
            sum(cov$is_covered), nrow(cov), opt$out))
