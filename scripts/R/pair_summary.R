#!/usr/bin/env Rscript
# pair_summary.R — collapse a marker's cross-species BLAST hits to one row
# per unordered species pair, with strict-ambiguity flag and ladder class.
#
# Provenance: NEW (Phase 4 extension; not in the legacy Rmd). Mirrors the
# brief's "richer per-pair summary" intent. Ladder definitions come from
# CONTEXT-analysis.md.
#
# Usage:
#   pair_summary.R --marker {mit|18S} \
#     --blast outputs/blast/self_blast_{marker}.tsv \
#     [--acc-map outputs/cross_species/18S_acc_to_header.tsv] \
#     --pident-min 99 --cov-shorter-min 0.99 --shorter-len-min 10 \
#     --top-n-hsps 5 \
#     --out outputs/cross_species/{marker}_pair_summary.tsv

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
  marker             = NA_character_,
  blast              = NA_character_,
  `acc-map`          = NA_character_,
  `pident-min`       = "99",
  `cov-shorter-min`  = "0.99",
  `shorter-len-min`  = "10",
  `top-n-hsps`       = "5",
  out                = NA_character_
))
stopifnot(opt$marker %in% c("mit", "18S"),
          !is.na(opt$blast), !is.na(opt$out))
if (opt$marker == "18S" && is.na(opt[["acc-map"]])) {
  stop("--acc-map is required when --marker 18S")
}

PIDENT_STRICT  <- as.numeric(opt[["pident-min"]])
COVS_STRICT    <- as.numeric(opt[["cov-shorter-min"]])
SHORTER_MIN    <- as.integer(opt[["shorter-len-min"]])
TOP_N          <- as.integer(opt[["top-n-hsps"]])

SPECIES_PAT <- "(vivax|falciparum|knowlesi|malariae|ovale|coat|inui|fieldi|cyno|simiovale|simium)"
blast_cols  <- c("query","subject","pident","aln_len","mismatch","gapopen",
                 "qlen","slen","evalue","bitscore")

# Canonicalise species names: lowercase, strip leading "p." or "p" (so that
# "psimium", "p.simium", and "simium" all collapse to "simium"). Then merge
# legacy typos.
canon_species <- function(x) {
  x |>
    str_to_lower() |>
    str_replace("^p\\.?", "") |>
    str_replace("^cynomologi$", "cynomolgi")
}

# Same scaffold as cross_species_filter.R (read BLAST, optionally join 18S
# acc map, filter to panel species, compute coverage metrics). We keep this
# logic local rather than calling the other script so pair_summary can be
# run end-to-end from a single BLAST tsv.
df <- fread(opt$blast, sep = "\t", header = FALSE, col.names = blast_cols)

if (opt$marker == "18S") {
  m <- fread(opt[["acc-map"]], sep = "\t", header = TRUE,
             select = c("acc", "species"))
  df <- df |>
    left_join(m |> rename(query   = acc, query_species   = species), by = "query") |>
    left_join(m |> rename(subject = acc, subject_species = species), by = "subject") |>
    mutate(
      query   = paste0(query,   "_", coalesce(query_species,   "unknown")),
      subject = paste0(subject, "_", coalesce(subject_species, "unknown"))
    ) |>
    select(-query_species, -subject_species)
}

df <- df |>
  filter(str_detect(query,   regex(SPECIES_PAT, ignore_case = TRUE)),
         str_detect(subject, regex(SPECIES_PAT, ignore_case = TRUE)))

if (opt$marker == "mit") {
  df <- df |> mutate(
    query_species   = str_extract(query,   "^[^_]+"),
    subject_species = str_extract(subject, "^[^_]+")
  )
} else {
  df <- df |> mutate(
    query_species   = sub("^[^_]+_", "", query),
    subject_species = sub("^[^_]+_", "", subject)
  )
}
df <- df |> mutate(
  query_species   = canon_species(query_species),
  subject_species = canon_species(subject_species)
)

df <- df |> filter(query != subject, query_species != subject_species)

df <- df |> mutate(
  qcov        = aln_len / qlen,
  scov        = aln_len / slen,
  shorter_len = pmin(qlen, slen),
  cov_shorter = aln_len / shorter_len,
  mincov      = pmin(qcov, scov)
)

# Ladder class per HSP. Priority ladder (highest priority first):
#   high_confidence_possible_ambiguity  (mincov ≥ 0.85 & pident ≥ 98)
#   moderate_possible_ambiguity         (mincov ≥ 0.70 & pident ≥ 98)
#   subject_partial                     (qcov ≥ 0.85 & scov < 0.50)
#   query_partial                       (scov ≥ 0.85 & qcov < 0.50)
#   short_high_identity_local_match     (mincov < 0.50 & pident ≥ 98)
#   lower_concern                       (otherwise)
ladder_one <- function(pident, qcov, scov, mincov) {
  if (mincov >= 0.85 & pident >= 98) return("high_confidence_possible_ambiguity")
  if (mincov >= 0.70 & pident >= 98) return("moderate_possible_ambiguity")
  if (qcov   >= 0.85 & scov  <  0.50) return("subject_partial")
  if (scov   >= 0.85 & qcov  <  0.50) return("query_partial")
  if (mincov <  0.50 & pident >= 98) return("short_high_identity_local_match")
  "lower_concern"
}
LADDER_PRIORITY <- c(
  high_confidence_possible_ambiguity = 1L,
  moderate_possible_ambiguity        = 2L,
  subject_partial                    = 3L,
  query_partial                      = 4L,
  short_high_identity_local_match    = 5L,
  lower_concern                      = 6L
)

df <- df |> mutate(
  ladder_class = pmap_chr(list(pident, qcov, scov, mincov), ladder_one),
  ladder_rank  = LADDER_PRIORITY[ladder_class]
)

# Top-N HSPs per query–subject pair (matches the strict pipeline so the
# strict_ambiguity column is consistent with `<marker>_suspicious.tsv`).
df_topN <- df |>
  arrange(query, subject, desc(bitscore), desc(aln_len)) |>
  group_by(query, subject) |>
  slice_head(n = TOP_N) |>
  ungroup()

df_topN <- df_topN |> mutate(
  pair = pmap_chr(list(query_species, subject_species),
                  \(a, b) paste(sort(c(a, b)), collapse = " vs "))
)

pair_summary <- df_topN |>
  group_by(pair) |>
  summarise(
    n_hits          = n(),
    best_pident     = max(pident),
    max_aln_len     = max(aln_len),
    max_cov_shorter = max(cov_shorter),
    max_mincov      = max(mincov),
    min_mismatch    = min(mismatch),
    min_gapopen     = min(gapopen),
    strict_ambiguity = any(pident      >= PIDENT_STRICT &
                           cov_shorter >= COVS_STRICT  &
                           shorter_len >= SHORTER_MIN),
    # Pair-level ladder = highest-priority class achieved by any HSP in the pair.
    ladder_class    = ladder_class[which.min(ladder_rank)],
    .groups = "drop"
  ) |>
  arrange(desc(strict_ambiguity), desc(best_pident))

dir.create(dirname(opt$out), showWarnings = FALSE, recursive = TRUE)
fwrite(pair_summary, opt$out, sep = "\t")
cat(sprintf("[pair_summary:%s] %d cross-species pairs (strict_ambiguity=TRUE: %d) -> %s\n",
            opt$marker, nrow(pair_summary), sum(pair_summary$strict_ambiguity),
            opt$out))
