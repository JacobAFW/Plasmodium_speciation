#!/usr/bin/env Rscript
# cross_species_filter.R — apply the canonical "suspicious cross-species
# near-identity" filter and the per-sequence "closest non-self match"
# summary to a self-BLAST tabular output.
#
# Provenance: mit_similarity.Rmd lines 25-139 (MIT) and 242-391 (18S).
# The MIT path extracts species from the self-describing query/subject IDs
# (regex `^([^_]+)`); the 18S path joins on a parsed acc->species mapping
# (`18S_acc_to_header.tsv`) before computing species.
#
# Usage:
#   cross_species_filter.R \
#     --marker {mit|18S} \
#     --blast outputs/blast/self_blast_<marker>.tsv \
#     [--acc-map outputs/cross_species/18S_acc_to_header.tsv]   # required for 18S
#     --pident-min 99 --cov-shorter-min 0.99 --shorter-len-min 10 --top-n-hsps 5 \
#     --suspicious-out outputs/cross_species/<marker>_suspicious.tsv \
#     --closest-out    outputs/cross_species/<marker>_closest_per_seq.tsv

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

# Minimal --flag value CLI parser. Unknown flags error out so typos don't
# silently get ignored.
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
  `cov-min`          = "0.99",
  `cov-metric`       = "cov_shorter",        # cov_shorter (default) or mincov
  `shorter-len-min`  = "10",
  `top-n-hsps`       = "5",
  `suspicious-out`   = NA_character_,
  `closest-out`      = NA_character_
))

stopifnot(opt$marker %in% c("mit", "18S"),
          opt[["cov-metric"]] %in% c("cov_shorter", "mincov"),
          !is.na(opt$blast),
          !is.na(opt[["suspicious-out"]]),
          !is.na(opt[["closest-out"]]))
if (opt$marker == "18S" && is.na(opt[["acc-map"]])) {
  stop("--acc-map is required when --marker 18S")
}

PIDENT_MIN      <- as.numeric(opt[["pident-min"]])
COV_MIN         <- as.numeric(opt[["cov-min"]])
COV_METRIC      <- opt[["cov-metric"]]
SHORTER_LEN_MIN <- as.integer(opt[["shorter-len-min"]])
TOP_N_HSPS      <- as.integer(opt[["top-n-hsps"]])

# Same regex the legacy chunks use; identifies the panel-relevant species.
SPECIES_PAT <- "(vivax|falciparum|knowlesi|malariae|ovale|coat|inui|fieldi|cyno|simiovale|simium)"

blast_cols <- c("query","subject","pident","aln_len","mismatch","gapopen",
                "qlen","slen","evalue","bitscore")

df <- fread(opt$blast, sep = "\t", header = FALSE, col.names = blast_cols)

# --- 18S path: append _<species> to query/subject IDs via the acc map ---
# The legacy logic appends species to the IDs themselves (rather than carrying
# a separate column) so the rest of the pipeline can stay marker-agnostic.
if (opt$marker == "18S") {
  m <- fread(opt[["acc-map"]], sep = "\t", header = TRUE,
             select = c("acc", "species"))
  df <- df |>
    left_join(m |> rename(query   = acc, query_species   = species),   by = "query") |>
    left_join(m |> rename(subject = acc, subject_species = species),   by = "subject") |>
    mutate(
      query   = paste0(query,   "_", coalesce(query_species,   "unknown")),
      subject = paste0(subject, "_", coalesce(subject_species, "unknown"))
    ) |>
    select(-query_species, -subject_species)
}

# Species filter applies to query AND subject (both must be panel-relevant).
df <- df |>
  filter(str_detect(query,   regex(SPECIES_PAT, ignore_case = TRUE)),
         str_detect(subject, regex(SPECIES_PAT, ignore_case = TRUE)))

# Species extraction: MIT pulls from `^([^_]+)`; 18S splits on the FIRST `_`
# and takes the suffix (which is the appended species we built above).
if (opt$marker == "mit") {
  df <- df |> mutate(
    query_species   = str_to_lower(str_extract(query,   "^[^_]+")),
    subject_species = str_to_lower(str_extract(subject, "^[^_]+"))
  )
} else {
  # Legacy uses Python's str.split("_", n=1) which splits ONCE on the first
  # underscore and keeps the rest. R's str_split_i(., "_", 2) splits on ALL
  # underscores and takes the 2nd piece, which loses everything after the
  # second `_`. For NCBI predicted-RNA accessions like XR_002198256.1, the
  # ID after species-append is `XR_002198256.1_p.coatneyi` and the legacy
  # species column ends up as `002198256.1_p.coatneyi` (the legacy ladder
  # has been documented; we reproduce it for byte-equivalence).
  df <- df |> mutate(
    query_species   = str_to_lower(sub("^[^_]+_", "", query)),
    subject_species = str_to_lower(sub("^[^_]+_", "", subject))
  )
}

# Legacy: typo "pcynomologi" got normalised to "pcynomolgi" inline. Keep.
df <- df |> mutate(
  query_species   = str_replace(query_species,   "pcynomologi", "pcynomolgi"),
  subject_species = str_replace(subject_species, "pcynomologi", "pcynomolgi")
)

# Drop self-hits and within-species hits.
df <- df |> filter(query != subject, query_species != subject_species)

# Coverage metrics. mincov is the more conservative metric (always
# ≤ cov_shorter); used by the looser ladder.
df <- df |> mutate(
  qcov         = aln_len / qlen,
  scov         = aln_len / slen,
  shorter_len  = pmin(qlen, slen),
  cov_shorter  = aln_len / shorter_len,
  mincov       = pmin(qcov, scov)
)

# Top-N HSPs per query–subject pair (legacy: 5).
df_topN <- df |>
  arrange(query, subject, desc(bitscore), desc(aln_len)) |>
  group_by(query, subject) |>
  slice_head(n = TOP_N_HSPS) |>
  ungroup()

# --- Suspicious cross-species near-identity table ---
# species_pair is appended to match the legacy MIT output (the Rmd as
# checked-in omits it; replication parity documented in MEMORY.md).
suspicious <- df_topN |>
  filter(pident                >= PIDENT_MIN,
         .data[[COV_METRIC]]   >= COV_MIN,
         shorter_len           >= SHORTER_LEN_MIN) |>
  mutate(species_pair = pmap_chr(
    list(query_species, subject_species),
    \(q, s) paste(sort(c(q, s)), collapse = " vs ")
  ))

dir.create(dirname(opt[["suspicious-out"]]), showWarnings = FALSE, recursive = TRUE)
fwrite(suspicious, opt[["suspicious-out"]], sep = "\t")
cat(sprintf("[cross_species_filter:%s] %d suspicious cross-species hits -> %s\n",
            opt$marker, nrow(suspicious), opt[["suspicious-out"]]))

# --- Closest non-self match per sequence ---
best <- df |>
  arrange(query, desc(bitscore), desc(aln_len)) |>
  group_by(query) |>
  slice_head(n = 1) |>
  ungroup() |>
  mutate(pct_covered    = cov_shorter * 100,
         divergence_pct = 100 - pident) |>
  select(query, query_species, subject, subject_species, pident, divergence_pct,
         aln_len, pct_covered, mismatch, gapopen, bitscore) |>
  arrange(desc(pident))

dir.create(dirname(opt[["closest-out"]]), showWarnings = FALSE, recursive = TRUE)
fwrite(best, opt[["closest-out"]], sep = "\t")
cat(sprintf("[cross_species_filter:%s] %d closest-per-sequence rows -> %s\n",
            opt$marker, nrow(best), opt[["closest-out"]]))
