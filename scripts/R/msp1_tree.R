#!/usr/bin/env Rscript
# msp1_tree.R — first-look NJ tree from MSP1 reference panel.
#
# Two stages, selected by --stage:
#   --stage profile   reads an aligned FASTA, reports per-species and per-record
#                     completeness statistics, recommends a length cutoff. No
#                     tree built. STOP-checkpoint 1 in the prompt.
#   --stage tree      applies the cutoff, optionally re-aligns, builds an NJ tree
#                     (TN93 distance, midpoint root), writes Newick + two figures
#                     (full + cynomolgi subtree, PDF + SVG). STOP-checkpoint 2.
#
# This is a first-look tree, not a publishable phylogeny: NJ, no bootstrap,
# no alignment trimming. Companion to chapter 8 (06d_msp1_swangsri_validation.qmd).
#
# Usage:
#   msp1_tree.R --stage profile \
#               --aln    outputs/msp1/tree/msp1_aln_all.fasta \
#               --out    outputs/msp1/tree/profile.tsv \
#               --report outputs/msp1/tree/profile_report.md
#
#   msp1_tree.R --stage tree \
#               --aln          outputs/msp1/tree/msp1_aln_all.fasta \
#               --min-len      3000 \
#               --filtered-aln outputs/msp1/tree/msp1_aln_filtered.fasta \
#               --realign      true \
#               --nwk          outputs/msp1/tree/msp1_tree.nwk \
#               --fig-full     reports/figures/msp1_tree_full \
#               --fig-cyno     reports/figures/msp1_tree_cyno
#
# Highlight set (bold + marker on tree tips):
#   AB266195.1 — Pcynomolgi Pt2 (Berok)  — Swangsri Pcy primer design reference
#   AB444063.1 — Pcynomolgi ATCC 30146   — lab strain used for plasmid clone
#
# Plot palette: Okabe-Ito (categorical, colour-blind safe).

suppressPackageStartupMessages({
  library(Biostrings)
  library(tidyverse)
})

# ---- args ----------------------------------------------------------------
parse_args <- function(known) {
  a <- commandArgs(trailingOnly = TRUE)
  if (length(a) %% 2 != 0) stop("expected --flag value pairs")
  out <- known
  for (i in seq(1, length(a), by = 2)) {
    k <- sub("^--", "", a[i])
    if (!k %in% names(known)) stop(sprintf("unknown flag: %s", a[i]))
    out[[k]] <- a[i + 1]
  }
  out
}

opt <- parse_args(list(
  stage           = NA_character_,
  aln             = NA_character_,
  out             = NA_character_,
  report          = NA_character_,
  `min-len`       = "3000",
  `filtered-aln`  = NA_character_,
  realign         = "true",
  nwk             = NA_character_,
  `fig-full`      = NA_character_,
  `fig-cyno`      = NA_character_,
  mafft           = "vvg-box/opt/umamba/envs/vvg-box/bin/mafft"
))
if (is.na(opt$stage) || !opt$stage %in% c("profile", "tree"))
  stop("--stage must be 'profile' or 'tree'")
if (is.na(opt$aln)) stop("--aln required")

# ---- shared helpers ------------------------------------------------------
parse_species <- function(headers) {
  # Headers are normalised by fetch_msp1_refs.py to:
  #   "ACC.V Pspecies msp1 | <descr>"
  m <- regmatches(headers, regexec("^\\S+\\s+P([A-Za-z]+)\\s+msp1", headers))
  out <- vapply(m, function(x) if (length(x) == 2) x[2] else NA_character_, "")
  out
}
parse_acc <- function(headers) sub("\\s.*$", "", headers)

# Inferno colormap from viridis. Sequential by design, but evenly-spaced
# samples skipping the near-black and pure-yellow endpoints give a usable
# categorical palette for ~5–8 species. The cyno subtree palette is recomputed
# from just the species present so the two colours used (cynomolgi + sister
# context, typically fieldi) sit far apart in the colormap rather than
# adjacent in a 7-species spread.
palette_for <- function(n) {
  if (n == 1) return("#57106E")  # inferno(0.25)
  # Skip the near-black low end (begin >= 0.25 starts at clean purple,
  # not dark-grey/black).
  viridisLite::inferno(n, begin = 0.25, end = 0.95)
}

# Highlights — accessions singled out in the chapter framing.
HIGHLIGHT <- c("AB266195.1", "AB444063.1")
HIGHLIGHT_LBL <- c(
  "AB266195.1" = "Pt2/Berok (design ref)",
  "AB444063.1" = "ATCC 30146 (lab strain)"
)

ungap_widths <- function(dss) {
  # Per-record non-gap count.
  m <- as.matrix(dss)
  rowSums(m != "-")
}

# ============================================================
# stage = profile
# ============================================================
if (opt$stage == "profile") {
  if (is.na(opt$out))    stop("--out required for stage=profile")
  if (is.na(opt$report)) stop("--report required for stage=profile")

  aln <- readDNAStringSet(opt$aln)
  n_seq <- length(aln)
  aln_w <- unique(width(aln))
  if (length(aln_w) != 1)
    stop(sprintf("alignment widths not uniform (%s)", paste(aln_w, collapse=",")))

  hdr <- names(aln)
  sp  <- parse_species(hdr)
  acc <- parse_acc(hdr)
  ungapped <- ungap_widths(aln)

  per_rec <- tibble(
    accession   = acc,
    species     = sp,
    aln_width   = aln_w,
    ungapped_bp = ungapped,
    coverage    = ungapped / aln_w
  )
  write_tsv(per_rec, opt$out)

  # Per-species summary
  per_sp <- per_rec |>
    group_by(species) |>
    summarise(
      n            = n(),
      median_bp    = as.integer(median(ungapped_bp)),
      max_bp       = max(ungapped_bp),
      n_ge_3kb     = sum(ungapped_bp >= 3000L),
      n_ge_5kb     = sum(ungapped_bp >= 5000L),
      median_cov   = round(median(coverage), 3),
      .groups      = "drop"
    ) |>
    arrange(desc(n))

  # Coverage breakpoints
  brk <- c(0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9, 1.0001)
  cov_hist <- per_rec |>
    mutate(bin = cut(coverage, breaks = brk, right = FALSE,
                     include.lowest = TRUE)) |>
    count(bin)

  # Highlight presence
  hl <- per_rec |>
    filter(accession %in% HIGHLIGHT) |>
    arrange(match(accession, HIGHLIGHT)) |>
    select(accession, species, ungapped_bp, coverage)

  # Cynomolgi-specific drilldown
  cyn <- per_rec |>
    filter(species == "cynomolgi") |>
    arrange(desc(ungapped_bp))

  # Markdown report
  lines <- c(
    "# MSP1 alignment profile (Step 1 — STOP-CHECKPOINT 1)",
    "",
    sprintf("- Input alignment: `%s`", opt$aln),
    sprintf("- Sequences: **%d**", n_seq),
    sprintf("- Alignment width: **%d** columns", aln_w),
    sprintf("- Overall mean per-record coverage: **%.3f** (gap fraction %.3f)",
            mean(per_rec$coverage), 1 - mean(per_rec$coverage)),
    "",
    "## Per-species record completeness",
    "",
    "| Species | n | median ungapped bp | max bp | n ≥ 3 kb | n ≥ 5 kb | median cov |",
    "|---|---:|---:|---:|---:|---:|---:|",
    paste(sprintf("| %s | %d | %d | %d | %d | %d | %.3f |",
                  per_sp$species, per_sp$n, per_sp$median_bp, per_sp$max_bp,
                  per_sp$n_ge_3kb, per_sp$n_ge_5kb, per_sp$median_cov),
          collapse = "\n"),
    "",
    "## Coverage distribution (per-record fraction non-gap)",
    "",
    "| Coverage bin | n records |",
    "|---|---:|",
    paste(sprintf("| %s | %d |", as.character(cov_hist$bin), cov_hist$n),
          collapse = "\n"),
    "",
    "## Cynomolgi drill-down (n=25 expected)",
    "",
    "| Accession | ungapped bp | coverage |",
    "|---|---:|---:|",
    paste(sprintf("| %s | %d | %.3f |",
                  cyn$accession, cyn$ungapped_bp, cyn$coverage),
          collapse = "\n"),
    "",
    "## Highlight set",
    "",
    "| Accession | Species | ungapped bp | coverage | Note |",
    "|---|---|---:|---:|---|",
    paste(sprintf("| %s | %s | %d | %.3f | %s |",
                  hl$accession, hl$species, hl$ungapped_bp, hl$coverage,
                  HIGHLIGHT_LBL[hl$accession]),
          collapse = "\n")
  )
  writeLines(lines, opt$report)

  # Console summary
  cat(sprintf("\n[profile] %d sequences, alignment width %d cols\n",
              n_seq, aln_w))
  cat("\n[profile] per-species summary:\n")
  print(per_sp, n = Inf)
  cat("\n[profile] coverage distribution:\n")
  print(cov_hist)
  cat("\n[profile] highlight accessions:\n")
  print(hl)
  cat("\n[profile] cynomolgi drilldown (top 10 by length):\n")
  print(head(cyn, 10))
  cat(sprintf("\n[profile] wrote per-record TSV: %s\n", opt$out))
  cat(sprintf("[profile] wrote report:        %s\n", opt$report))
  quit(status = 0)
}

# ============================================================
# stage = tree
# ============================================================
if (opt$stage == "tree") {
  if (is.na(opt$`filtered-aln`)) stop("--filtered-aln required for stage=tree")
  if (is.na(opt$nwk))            stop("--nwk required for stage=tree")
  if (is.na(opt$`fig-full`))     stop("--fig-full required for stage=tree")
  if (is.na(opt$`fig-cyno`))     stop("--fig-cyno required for stage=tree")

  suppressPackageStartupMessages({
    library(ape)
    library(phangorn)
    library(ggtree)
    library(treeio)
    library(tidytree)
    library(viridisLite)
  })

  min_len <- as.integer(opt$`min-len`)
  realign <- tolower(opt$realign) %in% c("true", "yes", "1")

  aln <- readDNAStringSet(opt$aln)
  ungapped <- ungap_widths(aln)
  keep <- ungapped >= min_len
  cat(sprintf("[tree] filter: ungapped >= %d bp -> %d / %d kept\n",
              min_len, sum(keep), length(aln)))
  filt_unaligned <- DNAStringSet(gsub("-", "", as.character(aln[keep])))
  names(filt_unaligned) <- names(aln)[keep]

  if (realign) {
    cat("[tree] re-aligning filtered subset with mafft --auto...\n")
    tmp_in  <- tempfile(fileext = ".fasta")
    writeXStringSet(filt_unaligned, tmp_in)
    cmd <- sprintf("%s --auto --thread 4 %s > %s 2> %s",
                   shQuote(opt$mafft), shQuote(tmp_in),
                   shQuote(opt$`filtered-aln`),
                   shQuote(paste0(opt$`filtered-aln`, ".log")))
    rc <- system(cmd)
    if (rc != 0) stop("mafft re-alignment failed")
    file.remove(tmp_in)
  } else {
    writeXStringSet(aln[keep], opt$`filtered-aln`)
  }

  filt_aln <- readDNAStringSet(opt$`filtered-aln`)
  cat(sprintf("[tree] filtered alignment: n=%d width=%d\n",
              length(filt_aln), unique(width(filt_aln))))

  # ---- distance + NJ tree ----
  dna <- as.DNAbin(filt_aln)
  cat("[tree] computing TN93 distance...\n")
  d <- dist.dna(dna, model = "TN93", pairwise.deletion = TRUE,
                as.matrix = FALSE)
  if (any(is.nan(d)) || any(is.na(d))) {
    n_bad <- sum(is.nan(d) | is.na(d))
    cat(sprintf("[tree] NOTE: %d NA/NaN distances under TN93; falling back to raw\n",
                n_bad))
    d <- dist.dna(dna, model = "raw", pairwise.deletion = TRUE,
                  as.matrix = FALSE)
  }
  cat("[tree] building NJ tree...\n")
  tr <- nj(d)
  tr <- midpoint(tr)
  write.tree(tr, opt$nwk)
  cat(sprintf("[tree] wrote Newick: %s\n", opt$nwk))

  # ---- annotate tips ----
  hdr <- tr$tip.label
  sp  <- parse_species(hdr)
  acc <- parse_acc(hdr)

  # Strain extraction from descriptive remainder for tip labels.
  strain_re <- ".*\\bstrain[: ]\\s*([^,|]+).*"
  rest <- sub("^\\S+\\s+P[A-Za-z]+\\s+msp1\\s*\\|\\s*", "", hdr)
  strain <- ifelse(grepl(strain_re, rest, ignore.case = TRUE),
                   sub(strain_re, "\\1", rest, ignore.case = TRUE),
                   NA_character_)
  isolate_re <- ".*\\bisolate\\s+([A-Za-z0-9._-]+).*"
  isolate <- ifelse(grepl(isolate_re, rest, ignore.case = TRUE),
                    sub(isolate_re, "\\1", rest, ignore.case = TRUE),
                    NA_character_)
  short_strain <- coalesce(strain, isolate)

  meta <- tibble(label = hdr, accession = acc, species = sp,
                 strain = short_strain,
                 highlight = acc %in% HIGHLIGHT,
                 highlight_note = HIGHLIGHT_LBL[acc])

  # Tip label: "ACC strain" (if strain present), else "ACC".
  meta <- meta |>
    mutate(tip_label = ifelse(is.na(strain) | strain == "",
                              accession,
                              paste0(accession, "  ", strain)))

  # ---- full tree figure ----
  # Alphabetical ordering -> deterministic species->colour mapping across runs.
  species_levels <- sort(unique(meta$species))
  pal <- setNames(palette_for(length(species_levels)), species_levels)

  cat(sprintf("[tree] drawing full tree (n=%d tips)...\n", length(tr$tip.label)))
  p_full <- ggtree(tr, layout = "rectangular") %<+% meta +
    geom_tippoint(aes(color = species), size = 1.0, na.rm = TRUE) +
    geom_tiplab(aes(label = ifelse(highlight, tip_label, ""),
                    fontface = "bold"),
                size = 2.4, color = "black", offset = 0.001) +
    geom_tippoint(aes(subset = highlight),
                  shape = 21, size = 2.6, fill = "yellow",
                  color = "black", stroke = 0.5) +
    scale_color_manual(values = pal, name = "species",
                       na.value = "grey60") +
    theme_tree2() +
    labs(caption = sprintf(
      "First-look NJ tree (TN93 distance, midpoint root). n=%d sequences, msp1 ungapped >= %d bp. No bootstrap.",
      length(tr$tip.label), min_len))

  ggsave(paste0(opt$`fig-full`, ".pdf"), p_full,
         width = 10, height = 14, device = cairo_pdf)
  ggsave(paste0(opt$`fig-full`, ".svg"), p_full,
         width = 10, height = 14)
  cat(sprintf("[tree] wrote %s.{pdf,svg}\n", opt$`fig-full`))

  # ---- cynomolgi subtree ----
  cyno_tips <- meta$label[meta$species == "cynomolgi"]
  if (length(cyno_tips) >= 2) {
    # Walk up to MRCA of cynomolgi tips, then take its descendants plus the
    # parent clade for context.
    cyno_mrca <- ape::getMRCA(tr, cyno_tips)
    parent <- tr$edge[tr$edge[, 2] == cyno_mrca, 1]
    if (length(parent) == 0) parent <- cyno_mrca  # already at root
    sub_tips <- tr$tip.label[phangorn::Descendants(tr, parent, "tips")[[1]]]
    sub_tr <- ape::keep.tip(tr, sub_tips)

    cat(sprintf("[tree] cyno subtree: %d cynomolgi tips + %d context tips\n",
                length(cyno_tips), length(sub_tips) - length(cyno_tips)))

    sub_meta <- meta |> filter(label %in% sub_tips)
    # Recompute the palette from just the species in the subtree so the
    # cynomolgi vs sister-context colours are far apart in inferno rather
    # than adjacent (which they would be on the 7-species spread).
    sub_levels <- sort(unique(sub_meta$species))
    sub_pal <- setNames(palette_for(length(sub_levels)), sub_levels)

    p_cyno <- ggtree(sub_tr, layout = "rectangular") %<+% sub_meta +
      geom_tippoint(aes(color = species), size = 1.4, na.rm = TRUE) +
      geom_tiplab(aes(label = tip_label,
                      fontface = ifelse(highlight, "bold", "plain")),
                  size = 2.6, offset = 0.0005, color = "black") +
      geom_tippoint(aes(subset = highlight),
                    shape = 21, size = 3.2, fill = "yellow",
                    color = "black", stroke = 0.6) +
      scale_color_manual(values = sub_pal, name = "species",
                         na.value = "grey60") +
      theme_tree2() +
      labs(caption = sprintf(
        "Cynomolgi clade + immediate context. NJ/TN93/midpoint, n=%d. Bold labels: highlight set.",
        length(sub_tips)))

    cyno_h <- max(6, 0.18 * length(sub_tips) + 2)
    ggsave(paste0(opt$`fig-cyno`, ".pdf"), p_cyno,
           width = 11, height = cyno_h, device = cairo_pdf)
    ggsave(paste0(opt$`fig-cyno`, ".svg"), p_cyno,
           width = 11, height = cyno_h)
    cat(sprintf("[tree] wrote %s.{pdf,svg}\n", opt$`fig-cyno`))
  } else {
    cat("[tree] insufficient cynomolgi tips for subtree (skipped)\n")
  }

  # ---- diagnostic summary on cynomolgi cluster integrity ----
  if (length(cyno_tips) >= 2) {
    cyno_idx <- which(tr$tip.label %in% cyno_tips)
    cyno_mrca <- ape::getMRCA(tr, cyno_tips)
    desc_of_mrca <- phangorn::Descendants(tr, cyno_mrca, "tips")[[1]]
    desc_labels <- tr$tip.label[desc_of_mrca]
    desc_meta <- meta |> filter(label %in% desc_labels)
    cyno_in_clade <- desc_meta |> filter(species == "cynomolgi") |> nrow()
    noncyno_in_clade <- desc_meta |> filter(species != "cynomolgi") |> nrow()
    cyno_off_clade <- length(cyno_tips) - cyno_in_clade
    cat(sprintf("\n[tree] cynomolgi clustering diagnostic:\n"))
    cat(sprintf("  cynomolgi tips total:                 %d\n", length(cyno_tips)))
    cat(sprintf("  cynomolgi tips inside MRCA clade:     %d\n", cyno_in_clade))
    cat(sprintf("  non-cynomolgi tips inside MRCA clade: %d\n", noncyno_in_clade))
    cat(sprintf("  cynomolgi tips outside MRCA clade:    %d\n", cyno_off_clade))

    # Pt2/Berok and ATCC 30146 placement
    for (a in HIGHLIGHT) {
      ix <- which(meta$accession == a)
      if (length(ix)) {
        in_mrca <- meta$label[ix] %in% desc_labels
        cat(sprintf("  %s (%s) -- inside cyno MRCA clade: %s\n",
                    a, HIGHLIGHT_LBL[a], in_mrca))
      } else {
        cat(sprintf("  %s -- NOT in filtered set\n", a))
      }
    }
  }

  quit(status = 0)
}
