#!/usr/bin/env Rscript
# msp1_amplicon_specificity.R — in-silico PCR of Swangsri et al. 2025
# species-specific msp1 primers against MSP1 references for all 11 target
# species in the panel. Answers: do the published primers extend cleanly to
# other species, or do they collide / fail / mis-prime?
#
# This is the *amplicon-level* question. The full-gene MSP1 alignment is a
# separate script (msp1_full_gene_alignment.R) that asks whether MSP1 should
# join MIT and 18S as a third marker in the long-read panel.
#
# Inputs (required):
#   --refdir   data/reference/msp1/        # one *_msp1.fasta per species,
#                                          # output of fetch_msp1_refs.py
#   --out      outputs/msp1/amplicon_specificity.tsv
#
# Optional:
#   --max-mismatch  3  (per primer; default 3 — generous, surfaces near-misses)
#
# Output TSV columns:
#   primer_pair        — Pk / Pcy / Pin (Swangsri species-specific set)
#   target_species     — species the references in this row belong to
#   ref_id             — FASTA header accession
#   ref_len_bp         — reference sequence length
#   fwd_strand         — + / - / NA (NA = no match within max-mismatch)
#   fwd_start          — 1-based start of forward primer match on the ref
#   fwd_mismatches     — # mismatches in forward primer (NA if no match)
#   rev_strand         — same for reverse primer (reverse-comp tested)
#   rev_start          — start of reverse primer match
#   rev_mismatches     — mismatches in reverse primer
#   amplicon_len       — predicted amplicon length if both primers bind in
#                        the right orientation; NA otherwise
#   amplicon_gc        — predicted %GC of the amplicon (NA if no amplicon)
#   predicted_tm_proxy — rough Tm proxy = 4*(G+C) + 2*(A+T) over the amplicon
#                        (Wallace rule). Not a direct competitor to SYBR-melt
#                        Tm but flags whether the amplicon would even *yield*
#                        a melt peak in the Swangsri-relevant range
#   on_target          — TRUE if target_species == primer's nominal target
#                        (Pk/knowlesi, Pcy/cynomolgi, Pin/inui)
#
# Provenance: companion to chapter 7 of the speciation panel report.
# Reproduces what colleagues validating the Swangsri 2025 assay would do by
# hand in BioEdit, but at panel scale.

suppressPackageStartupMessages({
  library(Biostrings)
  library(tidyverse)
})

# ---- args ----------------------------------------------------------------
parse_args <- function() {
  a <- commandArgs(trailingOnly = TRUE)
  defaults <- list(
    refdir       = "data/reference/msp1",
    out          = "outputs/msp1/amplicon_specificity.tsv",
    max_mismatch = "3"
  )
  if (length(a) %% 2 != 0) stop("expected --flag value pairs")
  out <- defaults
  for (i in seq(1, length(a), by = 2)) {
    k <- gsub("-", "_", sub("^--", "", a[i]))
    if (!k %in% names(defaults)) stop(sprintf("unknown flag: %s", a[i]))
    out[[k]] <- a[i + 1]
  }
  out$max_mismatch <- as.integer(out$max_mismatch)
  out
}
opt <- parse_args()

# ---- Swangsri 2025 primer table -----------------------------------------
# Sequences taken verbatim from Swangsri et al. 2025 (Sci Rep 15:22872),
# Methods §"PCR and melting curve". Each primer pair is the species-specific
# pair the paper designed against P. knowlesi / P. cynomolgi / P. inui msp1.
primers <- tribble(
  ~primer_pair, ~nominal_target, ~fwd_name,     ~fwd_seq,                    ~rev_name,     ~rev_seq,
  "Pk",         "knowlesi",      "Pk_MSP1_F626","TCAACGGGGTTAATGTCACCG",     "Pk_MSP1_R816","TGTAGAAGATGCTGCAGGGG",
  "Pcy",        "cynomolgi",     "Pcy_MSP1_F600","ACTACGGAGAATGGTAAAAGGAA",  "Pcy_MSP1_R699","AGCTTCCGTACTGCCTATCG",
  "Pin",        "inui",          "Pin_MSP1_F702","GACTCCTACTGTTTCGGGTG",     "Pin_MSP1_R820","CCTTCTCGTAACTTCCATCTTC"
)

# ---- helpers -------------------------------------------------------------
species_from_fname <- function(fname) {
  sub("_msp1\\.fasta$", "", basename(fname))
}

# Find best match for `primer` in `ref` allowing up to max_mm mismatches.
# Tests forward then reverse-complement of primer. Returns 1-row tibble:
# strand, start, end, mm (or all NA if no match within threshold).
best_match <- function(primer, ref, max_mm) {
  p     <- DNAString(primer)
  prc   <- reverseComplement(p)

  fwd <- vmatchPattern(p,   DNAStringSet(ref), max.mismatch = max_mm)[[1]]
  rev <- vmatchPattern(prc, DNAStringSet(ref), max.mismatch = max_mm)[[1]]

  # mismatch count for a given match: aligned subseq vs primer
  mm_count <- function(strand, s, e) {
    sub <- subseq(ref, s, e)
    if (strand == "+") {
      neditAt(sub, p, fixed = FALSE)
    } else {
      neditAt(sub, prc, fixed = FALSE)
    }
  }

  candidates <- bind_rows(
    if (length(fwd) > 0) tibble(strand = "+",
                                start = start(fwd), end = end(fwd))
    else tibble(),
    if (length(rev) > 0) tibble(strand = "-",
                                start = start(rev), end = end(rev))
    else tibble()
  )
  if (nrow(candidates) == 0) {
    return(tibble(strand = NA_character_, start = NA_integer_,
                  end = NA_integer_, mm = NA_integer_))
  }
  candidates <- candidates |>
    rowwise() |>
    mutate(mm = mm_count(strand, start, end)) |>
    ungroup() |>
    arrange(mm, start) |>
    slice(1)
  candidates
}

gc_pct <- function(seq) {
  letters_v <- strsplit(as.character(seq), "")[[1]]
  100 * sum(letters_v %in% c("G", "C")) / length(letters_v)
}

wallace_tm <- function(seq) {
  letters_v <- strsplit(toupper(as.character(seq)), "")[[1]]
  4 * sum(letters_v %in% c("G", "C")) +
  2 * sum(letters_v %in% c("A", "T"))
}

# ---- gather refs ---------------------------------------------------------
refdir <- normalizePath(opt$refdir, mustWork = TRUE)
ref_files <- list.files(refdir, pattern = "_msp1\\.fasta$", full.names = TRUE)
if (length(ref_files) == 0) {
  stop(sprintf("no *_msp1.fasta files in %s — has fetch_msp1_refs.py been run?",
               refdir))
}

# ---- main loop -----------------------------------------------------------
rows <- list()
for (f in ref_files) {
  target_species <- species_from_fname(f)
  refs <- readDNAStringSet(f)
  if (length(refs) == 0) next
  for (i in seq_along(refs)) {
    ref     <- refs[[i]]
    ref_id  <- names(refs)[i]
    ref_len <- length(ref)
    for (pp in seq_len(nrow(primers))) {
      pr  <- primers[pp, ]
      fwd <- best_match(pr$fwd_seq, ref, opt$max_mismatch)
      rev <- best_match(pr$rev_seq, ref, opt$max_mismatch)

      # Productive amplicon: fwd on + strand 5' of rev on - strand,
      # within a sensible distance (< 2000 bp for a "PCR-able" amplicon).
      amp_len <- NA_integer_
      amp_gc  <- NA_real_
      amp_tm  <- NA_real_
      if (!is.na(fwd$strand) && !is.na(rev$strand) &&
          fwd$strand == "+"  && rev$strand == "-" &&
          fwd$start  <  rev$end &&
          (rev$end - fwd$start + 1) < 2000) {
        amp_seq <- subseq(ref, fwd$start, rev$end)
        amp_len <- length(amp_seq)
        amp_gc  <- gc_pct(amp_seq)
        amp_tm  <- wallace_tm(amp_seq)
      }

      rows[[length(rows) + 1]] <- tibble(
        primer_pair        = pr$primer_pair,
        target_species     = target_species,
        ref_id             = ref_id,
        ref_len_bp         = ref_len,
        fwd_strand         = fwd$strand,
        fwd_start          = fwd$start,
        fwd_mismatches     = fwd$mm,
        rev_strand         = rev$strand,
        rev_start          = rev$start,
        rev_mismatches     = rev$mm,
        amplicon_len       = amp_len,
        amplicon_gc        = amp_gc,
        predicted_tm_proxy = amp_tm,
        on_target          = pr$nominal_target == target_species
      )
    }
  }
}

result <- bind_rows(rows)
out_path <- normalizePath(opt$out, mustWork = FALSE)
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
write_tsv(result, out_path)
cat(sprintf("[msp1] wrote %d rows to %s\n", nrow(result), out_path))

# ---- console summary -----------------------------------------------------
cat("\n[msp1] cross-species amplification summary (any ref with productive amplicon):\n")
result |>
  group_by(primer_pair, target_species) |>
  summarise(
    refs_total       = n(),
    refs_amplifying  = sum(!is.na(amplicon_len)),
    min_fwd_mm       = suppressWarnings(min(fwd_mismatches, na.rm = TRUE)),
    min_rev_mm       = suppressWarnings(min(rev_mismatches, na.rm = TRUE)),
    .groups          = "drop"
  ) |>
  mutate(across(starts_with("min_"), ~ ifelse(is.infinite(.x), NA_integer_, .x))) |>
  arrange(primer_pair, target_species) |>
  print(n = Inf)
