#!/usr/bin/env Rscript
# pcy_binding_site_check.R — locate the Swangsri Pcy primer binding windows
# on the 13 full-CDS P. cynomolgi MSP1 references, extract the per-position
# mismatch pattern, and report position-from-3' (the axis that matters for
# PCR specificity).
#
# Background: the Phase A in-silico amplicon check (outputs/msp1/amplicon_
# specificity_strict.tsv) showed only 1/25 cyno refs amplify under Swangsri
# strict tolerance — the lone amplifier is AB266195.1 (Pt2/Berok, Swangsri's
# design ref). The tree at outputs/msp1/tree/msp1_tree.nwk shows that all 13
# full-CDS cyno records cluster monophyletically but split into a Pt2/Berok
# sub-clade (Pt2/Berok + Gombak) and a larger ATCC-30146 sub-clade (11 refs).
# This script answers: where in the primer do the failing-clade mismatches
# land — at the 3' end (real failure) or 5' end (wet-lab-tolerable)?
#
# Primer sequences are kept in sync with
# scripts/R/msp1_amplicon_specificity.R (the single source of truth for
# Swangsri 2025 primer sequences in this repo). If you edit the primers
# there, edit them here too.
#
# Usage:
#   pcy_binding_site_check.R \
#     [--refs        data/reference/msp1/cynomolgi_msp1.fasta] \
#     [--min-len     3000] \
#     [--max-mm      5] \
#     [--out-tsv     outputs/msp1/binding_sites/pcy_binding_windows.tsv] \
#     [--out-fwd-aln outputs/msp1/binding_sites/pcy_fwd_alignment.txt] \
#     [--out-rev-aln outputs/msp1/binding_sites/pcy_rev_alignment.txt]

suppressPackageStartupMessages({
  library(Biostrings)
  library(tidyverse)
})

# ---- args ----------------------------------------------------------------
parse_args <- function(known) {
  a <- commandArgs(trailingOnly = TRUE)
  if (length(a) %% 2 != 0) stop("expected --flag value pairs")
  out <- known
  if (length(a) >= 2) {
    for (i in seq(1, length(a), by = 2)) {
      k <- gsub("-", "_", sub("^--", "", a[i]))
      if (!k %in% names(known)) stop(sprintf("unknown flag: %s", a[i]))
      out[[k]] <- a[i + 1]
    }
  }
  out
}
opt <- parse_args(list(
  refs        = "data/reference/msp1/cynomolgi_msp1.fasta",
  min_len     = "3000",
  max_mm      = "5",
  out_tsv     = "outputs/msp1/binding_sites/pcy_binding_windows.tsv",
  out_fwd_aln = "outputs/msp1/binding_sites/pcy_fwd_alignment.txt",
  out_rev_aln = "outputs/msp1/binding_sites/pcy_rev_alignment.txt"
))
min_len <- as.integer(opt$min_len)
max_mm  <- as.integer(opt$max_mm)

# ---- primers (in sync with scripts/R/msp1_amplicon_specificity.R) -------
primers <- tribble(
  ~primer,  ~name,            ~seq,
  "fwd",    "Pcy_MSP1_F600",  "ACTACGGAGAATGGTAAAAGGAA",
  "rev",    "Pcy_MSP1_R699",  "AGCTTCCGTACTGCCTATCG"
)

# ---- sub-clade assignment (from drafts/msp1_tree_summary.md) -------------
# pt2_berok  = Pt2/Berok (AB266195.1) + Gombak (AB266191.1), the 2-tip sister.
# atcc_30146 = the other 11 full-CDS cyno refs (everything else).
SUB_PT2 <- c("AB266195.1", "AB266191.1")
sub_clade <- function(acc) ifelse(acc %in% SUB_PT2, "pt2_berok", "atcc_30146")

parse_acc <- function(headers) sub("\\s.*$", "", headers)

# ---- locate primer binding window ---------------------------------------
# Searches both strands of `ref` with up to `max_mm` mismatches. Returns the
# best (lowest-mismatch, 5'-most tie-break) hit as a 1-row tibble
# (strand, start, end, mm). NULL if no hit within tolerance.
locate_primer <- function(primer, ref, max_mm) {
  p    <- DNAString(primer)
  prc  <- reverseComplement(p)
  hf   <- matchPattern(p,   ref, max.mismatch = max_mm, fixed = FALSE)
  hr   <- matchPattern(prc, ref, max.mismatch = max_mm, fixed = FALSE)

  prim_chars <- strsplit(primer, "")[[1]]

  cands <- list()
  if (length(hf) > 0) {
    for (j in seq_along(hf)) {
      s <- start(hf)[j]; e <- end(hf)[j]
      obs_chars <- strsplit(as.character(subseq(ref, s, e)), "")[[1]]
      mm <- sum(prim_chars != obs_chars)
      cands[[length(cands) + 1]] <- tibble(strand = "+", start = s, end = e, mm = mm)
    }
  }
  if (length(hr) > 0) {
    for (j in seq_along(hr)) {
      s <- start(hr)[j]; e <- end(hr)[j]
      obs_5p <- as.character(reverseComplement(subseq(ref, s, e)))
      obs_chars <- strsplit(obs_5p, "")[[1]]
      mm <- sum(prim_chars != obs_chars)
      cands[[length(cands) + 1]] <- tibble(strand = "-", start = s, end = e, mm = mm)
    }
  }
  if (length(cands) == 0) return(NULL)
  bind_rows(cands) |> arrange(mm, start) |> slice(1)
}

# Positions (1-indexed from primer 5' end) where primer and observed window
# differ. If strand == "-", the observed window is reverse-complemented first
# so positions are in the primer's 5'->3' frame.
mismatch_positions <- function(primer_seq, observed_window, strand) {
  obs <- if (strand == "+") observed_window
         else as.character(reverseComplement(DNAString(observed_window)))
  p_chars <- strsplit(primer_seq, "")[[1]]
  o_chars <- strsplit(obs, "")[[1]]
  which(p_chars != o_chars)
}

# ---- load refs ----------------------------------------------------------
all_refs <- readDNAStringSet(opt$refs)
keep <- width(all_refs) >= min_len
refs <- all_refs[keep]
cat(sprintf("[pcy] %d refs total, %d at length >= %d bp\n",
            length(all_refs), length(refs), min_len))

# ---- iterate ------------------------------------------------------------
rows <- list()
for (i in seq_along(refs)) {
  ref   <- refs[[i]]
  hdr   <- names(refs)[i]
  acc   <- parse_acc(hdr)
  clade <- sub_clade(acc)
  for (pp in seq_len(nrow(primers))) {
    pr  <- primers[pp, ]
    hit <- locate_primer(pr$seq, ref, max_mm)
    if (is.null(hit)) {
      rows[[length(rows) + 1]] <- tibble(
        ref_id          = acc,
        sub_clade       = clade,
        primer          = pr$primer,
        primer_name     = pr$name,
        primer_seq      = pr$seq,
        primer_len      = nchar(pr$seq),
        hit_found       = FALSE,
        strand          = NA_character_,
        start           = NA_integer_,
        end             = NA_integer_,
        n_mismatches    = NA_integer_,
        mm_pos_from_5   = NA_character_,
        mm_pos_from_3   = NA_character_,
        mm_in_3pr_3     = NA_integer_,
        mm_in_3pr_5     = NA_integer_,
        observed_window = NA_character_
      )
      next
    }
    obs <- as.character(subseq(ref, hit$start, hit$end))
    mm5 <- mismatch_positions(pr$seq, obs, hit$strand)
    mm3 <- nchar(pr$seq) - mm5 + 1
    rows[[length(rows) + 1]] <- tibble(
      ref_id          = acc,
      sub_clade       = clade,
      primer          = pr$primer,
      primer_name     = pr$name,
      primer_seq      = pr$seq,
      primer_len      = nchar(pr$seq),
      hit_found       = TRUE,
      strand          = hit$strand,
      start           = hit$start,
      end             = hit$end,
      n_mismatches    = length(mm5),
      mm_pos_from_5   = paste(mm5, collapse = ","),
      mm_pos_from_3   = paste(mm3, collapse = ","),
      mm_in_3pr_3     = sum(mm3 <= 3L),
      mm_in_3pr_5     = sum(mm3 <= 5L),
      observed_window = obs
    )
  }
}
result <- bind_rows(rows)
dir.create(dirname(opt$out_tsv), recursive = TRUE, showWarnings = FALSE)
write_tsv(result, opt$out_tsv)
cat(sprintf("[pcy] wrote %d rows to %s\n", nrow(result), opt$out_tsv))

# ---- per-primer alignment files -----------------------------------------
write_alignment <- function(prim_filter, out_path) {
  rows_p <- result |> filter(primer == prim_filter, hit_found)
  if (nrow(rows_p) == 0) { writeLines("(no hits)", out_path); return(invisible(NULL)) }
  primer_seq  <- rows_p$primer_seq[1]
  primer_name <- rows_p$primer_name[1]
  L <- nchar(primer_seq)

  hdr_from5 <- paste0(formatC(seq_len(L) %% 10, width = 1), collapse = "")
  from3 <- rep(" ", L)
  for (k in 1:5) from3[L - k + 1] <- as.character(k)
  hdr_from3 <- paste(from3, collapse = "")

  display <- function(primer_seq, obs, strand) {
    obs_5p <- if (strand == "+") obs
              else as.character(reverseComplement(DNAString(obs)))
    p_chars <- strsplit(primer_seq, "")[[1]]
    o_chars <- strsplit(obs_5p, "")[[1]]
    paste(ifelse(p_chars == o_chars, ".", o_chars), collapse = "")
  }

  pad <- 32
  lines <- c(
    sprintf("# %s (primer length %d bp). 5' end on the left, 3' end on the right.",
            primer_name, L),
    "# Top header: digit mod 10 (5' counting). Bottom header: position from 3' end (1-5).",
    "# Matches = '.'; mismatches shown as the variant base. Strand flagged per row.",
    "# Sub-clade labels from drafts/msp1_tree_summary.md.",
    "",
    paste0(strrep(" ", pad), hdr_from5, "   (pos from 5', mod 10)"),
    paste0(strrep(" ", pad), hdr_from3, "   (pos from 3', 1-5)"),
    paste0(sprintf("%-30s  ", "PRIMER"), primer_seq, "   <- reference (5'->3')"),
    ""
  )
  rows_p <- rows_p |> arrange(sub_clade, ref_id)
  for (r in seq_len(nrow(rows_p))) {
    tag  <- sprintf("[%s] %s", rows_p$sub_clade[r], rows_p$ref_id[r])
    disp <- display(rows_p$primer_seq[r], rows_p$observed_window[r],
                    rows_p$strand[r])
    info <- sprintf("   strand=%s  mm=%d  3'≤3:%d  3'≤5:%d  posFrom3=[%s]",
                    rows_p$strand[r], rows_p$n_mismatches[r],
                    rows_p$mm_in_3pr_3[r], rows_p$mm_in_3pr_5[r],
                    rows_p$mm_pos_from_3[r])
    lines <- c(lines, sprintf("%-30s  %s%s", tag, disp, info))
  }
  writeLines(lines, out_path)
}
write_alignment("fwd", opt$out_fwd_aln)
write_alignment("rev", opt$out_rev_aln)
cat(sprintf("[pcy] wrote %s\n[pcy] wrote %s\n", opt$out_fwd_aln, opt$out_rev_aln))

# ---- console summary ----------------------------------------------------
cat("\n[pcy] productive-hit summary:\n")
print(result |> count(primer, hit_found))
cat("\n[pcy] per-sub-clade × primer mismatch totals:\n")
print(result |>
        filter(hit_found) |>
        group_by(sub_clade, primer) |>
        summarise(
          n           = n(),
          total_mm    = sum(n_mismatches),
          mm_at_3pr_3 = sum(mm_in_3pr_3),
          mm_at_3pr_5 = sum(mm_in_3pr_5),
          .groups     = "drop"
        ) |>
        arrange(primer, sub_clade))

cat("\n[pcy] per-ref breakdown:\n")
print(result |>
        filter(hit_found) |>
        select(sub_clade, ref_id, primer, strand, n_mismatches,
               mm_in_3pr_3, mm_in_3pr_5, mm_pos_from_3) |>
        arrange(primer, sub_clade, ref_id), n = Inf)
