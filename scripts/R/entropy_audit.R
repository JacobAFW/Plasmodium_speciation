#!/usr/bin/env Rscript
# entropy_audit.R — independent recomputation of per-position entropy +
# coverage, diff against the pipeline's output, coverage spot-check at five
# positions, top/bottom candidate-window eyeball, and (18S only) a bootstrap
# CI band over the small reference panel.
#
# Provenance: Phase 5 Step 1 audit (HANDOFF.md). Re-derives the same
# Shannon-entropy-over-A/C/G/T-frequencies metric the legacy seq_entropy.py
# emits, but via Biostrings::consensusMatrix and a hand-rolled R entropy
# helper. Diff tolerance: |diff| < 1e-9 (rounding only).
#
# Usage:
#   entropy_audit.R --marker {mit|18S} \
#     --alignment outputs/alignment/ma_<marker>.target.fasta \
#     --pipeline-perpos     outputs/entropy/<marker>.per_position.tsv \
#     --pipeline-windows    outputs/entropy/<marker>.windows.tsv \
#     --pipeline-candidates outputs/entropy/<marker>.primer_candidates.tsv \
#     --recompute-out  outputs/audit/<marker>_entropy_recompute.tsv \
#     --diff-out       outputs/audit/<marker>_entropy_diff.tsv \
#     --spotcheck-out  outputs/audit/<marker>_coverage_spotcheck.tsv \
#     --eyeball-tsv    outputs/audit/<marker>_window_eyeball.tsv \
#     --eyeball-png    reports/figures/<marker>_window_eyeball.png \
#     --eyeball-svg    reports/figures/<marker>_window_eyeball.svg \
#     [--bootstrap-tsv outputs/audit/<marker>_entropy_bootstrap.tsv \
#      --bootstrap-png reports/figures/<marker>_entropy_bootstrap.png \
#      --bootstrap-svg reports/figures/<marker>_entropy_bootstrap.svg \
#      --bootstrap-n 100 --seed 42]

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
  marker                  = NA_character_,
  alignment               = NA_character_,
  `pipeline-perpos`       = NA_character_,
  `pipeline-windows`      = NA_character_,
  `pipeline-candidates`   = NA_character_,
  `recompute-out`         = NA_character_,
  `diff-out`              = NA_character_,
  `spotcheck-out`         = NA_character_,
  `eyeball-tsv`           = NA_character_,
  `eyeball-png`           = NA_character_,
  `eyeball-svg`           = NA_character_,
  `bootstrap-tsv`         = NA_character_,
  `bootstrap-png`         = NA_character_,
  `bootstrap-svg`         = NA_character_,
  `bootstrap-n`           = "100",
  seed                    = "42"
))
stopifnot(opt$marker %in% c("mit", "18S"))
set.seed(as.integer(opt$seed))
for (out in c(opt[["recompute-out"]], opt[["diff-out"]], opt[["spotcheck-out"]],
              opt[["eyeball-tsv"]], opt[["eyeball-png"]], opt[["eyeball-svg"]],
              opt[["bootstrap-tsv"]], opt[["bootstrap-png"]], opt[["bootstrap-svg"]])) {
  if (!is.na(out)) dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
}

# ---------------------------------------------------------------------------
# 1. Independent per-position recomputation
# ---------------------------------------------------------------------------
aln    <- readDNAStringSet(opt$alignment)
n_seqs <- length(aln)
aln_w  <- unique(width(aln))
stopifnot(length(aln_w) == 1)

# consensusMatrix counts every IUPAC character row + "-"/"+"/"." rows.
# Use `as.prob=FALSE` to get raw counts. Rows by name: A, C, G, T, others, "-".
cm <- consensusMatrix(aln, as.prob = FALSE)

# Identify gap rows (legacy seq_entropy.py treats "-" and "." as gaps).
gap_rows <- intersect(c("-", "."), rownames(cm))
gap_counts <- if (length(gap_rows) == 0) {
  rep(0L, ncol(cm))
} else if (length(gap_rows) == 1) {
  cm[gap_rows, ]
} else {
  colSums(cm[gap_rows, , drop = FALSE])
}

acgt_rows <- intersect(c("A", "C", "G", "T"), rownames(cm))
acgt_mat  <- cm[acgt_rows, , drop = FALSE]
acgt_totals <- colSums(acgt_mat)

shannon_acgt <- function(col) {
  # col is a length-4 vector of A,C,G,T counts.
  total <- sum(col)
  if (total == 0) return(0)
  p <- col[col > 0] / total
  -sum(p * log2(p))
}

entropy_v <- apply(acgt_mat, 2, shannon_acgt)
snp_v     <- apply(acgt_mat, 2, function(col) max(0L, sum(col > 0) - 1L))
coverage_v <- (n_seqs - gap_counts) / n_seqs

recompute <- tibble(
  pos       = seq_len(aln_w),
  coverage  = as.numeric(coverage_v),
  snp_count = as.integer(snp_v),
  entropy   = as.numeric(entropy_v)
)
fwrite(recompute, opt[["recompute-out"]], sep = "\t")

# ---------------------------------------------------------------------------
# 2. Diff against the pipeline
# ---------------------------------------------------------------------------
pipe_pp <- fread(opt[["pipeline-perpos"]], sep = "\t")
stopifnot(nrow(pipe_pp) == nrow(recompute))

diff <- pipe_pp |>
  rename(pipeline_coverage = coverage,
         pipeline_entropy  = entropy,
         pipeline_snp      = snp_count) |>
  left_join(
    recompute |>
      rename(recompute_coverage = coverage,
             recompute_entropy  = entropy,
             recompute_snp      = snp_count),
    by = "pos"
  ) |>
  mutate(
    diff_coverage = recompute_coverage - pipeline_coverage,
    diff_entropy  = recompute_entropy  - pipeline_entropy,
    diff_snp      = recompute_snp      - pipeline_snp
  )
fwrite(diff, opt[["diff-out"]], sep = "\t")

max_cov_diff <- max(abs(diff$diff_coverage), na.rm = TRUE)
max_ent_diff <- max(abs(diff$diff_entropy),  na.rm = TRUE)
max_snp_diff <- max(abs(diff$diff_snp),      na.rm = TRUE)

cat(sprintf("[entropy_audit:%s] |diff_coverage|_max = %.3e\n", opt$marker, max_cov_diff))
cat(sprintf("[entropy_audit:%s] |diff_entropy|_max  = %.3e\n", opt$marker, max_ent_diff))
cat(sprintf("[entropy_audit:%s] |diff_snp|_max      = %d\n",   opt$marker, as.integer(max_snp_diff)))
if (max_cov_diff > 1e-9 || max_ent_diff > 1e-9 || max_snp_diff > 0L) {
  cat(sprintf("[entropy_audit:%s] WARN diff exceeds tolerance — investigate\n", opt$marker))
}

# ---------------------------------------------------------------------------
# 3. Coverage spot-check at 5 positions
#    Pick: lowest entropy, highest entropy, median entropy, lowest coverage,
#    highest coverage. Count non-gap bases manually from the alignment.
# ---------------------------------------------------------------------------
pick_pos <- function(metric, fn) {
  v <- recompute[[metric]]
  recompute$pos[which(fn(v))][1]
}
spotcheck_positions <- unique(c(
  pick_pos("entropy",  function(v) v == min(v[v > 0], na.rm = TRUE)), # min nonzero
  pick_pos("entropy",  function(v) v == max(v, na.rm = TRUE)),        # max
  pick_pos("entropy",  function(v) abs(v - median(v[v > 0], na.rm = TRUE)) == min(abs(v - median(v[v > 0], na.rm = TRUE)))),
  pick_pos("coverage", function(v) v == min(v, na.rm = TRUE)),
  pick_pos("coverage", function(v) v == max(v, na.rm = TRUE))
))
spotcheck_positions <- na.omit(spotcheck_positions)

manual_counts <- function(p) {
  col <- substring(as.character(aln), p, p)
  tab <- table(col)
  list(
    A = unname(tab["A"]),
    C = unname(tab["C"]),
    G = unname(tab["G"]),
    T = unname(tab["T"]),
    gap = unname(tab["-"]) %||% 0L,
    other = sum(!names(tab) %in% c("A","C","G","T","-"))
  )
}
`%||%` <- function(a, b) if (is.null(a) || is.na(a)) b else a

sc_rows <- lapply(spotcheck_positions, function(p) {
  mc <- manual_counts(p)
  pipe_row <- pipe_pp[pos == p]
  manual_nongap <- sum(unlist(mc[c("A","C","G","T")]), na.rm = TRUE) +
                   (mc$other %||% 0L)
  manual_coverage <- manual_nongap / n_seqs
  tibble(
    pos                   = p,
    manual_A              = mc$A %||% 0L,
    manual_C              = mc$C %||% 0L,
    manual_G              = mc$G %||% 0L,
    manual_T              = mc$T %||% 0L,
    manual_gap            = mc$gap %||% 0L,
    manual_other          = mc$other %||% 0L,
    manual_coverage       = manual_coverage,
    pipeline_coverage     = pipe_row$coverage,
    diff_coverage         = manual_coverage - pipe_row$coverage,
    pipeline_entropy      = pipe_row$entropy,
    pipeline_snp_count    = pipe_row$snp_count
  )
})
spotcheck <- bind_rows(sc_rows)
fwrite(spotcheck, opt[["spotcheck-out"]], sep = "\t")
cat(sprintf("[entropy_audit:%s] spot-check: max |diff_coverage| = %.3e\n",
            opt$marker, max(abs(spotcheck$diff_coverage), na.rm = TRUE)))

# ---------------------------------------------------------------------------
# 4. Top / bottom candidate-window eyeball
#    Top 3 = lowest mean_entropy among accepted candidates (primer_candidates).
#    Bottom 3 = highest mean_entropy among windows with mean_coverage ≥ 0.9
#               that didn't pass (in `windows.tsv` but not `primer_candidates`).
#    Render: tile heatmap (rows = sequences, cols = alignment positions, fill = base).
# ---------------------------------------------------------------------------
pipe_win <- fread(opt[["pipeline-windows"]], sep = "\t")
pipe_pc  <- fread(opt[["pipeline-candidates"]], sep = "\t")

# Accepted: present in primer_candidates. Use w=25 to keep the eyeball plot
# narrow; if no w=25 candidates exist (unlikely), fall back to all.
top_w <- pipe_pc[window == 25][order(mean_entropy)][1:3]
if (nrow(top_w) < 3) top_w <- pipe_pc[order(mean_entropy)][1:3]

# Rejected: in windows.tsv at w=25 with mean_coverage ≥ 0.9 but NOT in pc.
candidate_keys <- paste(pipe_pc$window, pipe_pc$start, pipe_pc$end, sep = "_")
rej <- pipe_win[window == 25 & mean_coverage >= 0.9]
rej$key <- paste(rej$window, rej$start, rej$end, sep = "_")
rej <- rej[!key %in% candidate_keys]
bot_w <- rej[order(-mean_entropy)][1:3]

eyeball <- list()
build_slice <- function(label, win_row) {
  s <- as.integer(win_row$start); e <- as.integer(win_row$end)
  slice_mat <- as.matrix(aln, use.names = TRUE)  # n_seqs x aln_len char matrix
  slice <- as.data.frame(slice_mat[, s:e, drop = FALSE])
  slice$ref_id <- names(aln)
  slice <- slice |>
    pivot_longer(-ref_id, names_to = "col_offset", values_to = "base") |>
    mutate(pos = s + as.integer(sub("V", "", col_offset)) - 1L,
           label = label,
           win_start = s, win_end = e,
           win_mean_entropy = as.numeric(win_row$mean_entropy)) |>
    select(label, win_start, win_end, win_mean_entropy, ref_id, pos, base)
  slice
}

for (i in seq_len(nrow(top_w))) {
  eyeball[[length(eyeball) + 1]] <- build_slice(
    sprintf("TOP_%d (mean_entropy=%.4f)", i, top_w$mean_entropy[i]),
    top_w[i])
}
for (i in seq_len(nrow(bot_w))) {
  eyeball[[length(eyeball) + 1]] <- build_slice(
    sprintf("BOT_%d (mean_entropy=%.4f)", i, bot_w$mean_entropy[i]),
    bot_w[i])
}
eyeball_df <- bind_rows(eyeball)
fwrite(eyeball_df, opt[["eyeball-tsv"]], sep = "\t")

# Tile heatmap
eyeball_df <- eyeball_df |>
  mutate(label_f = factor(label, levels = unique(label)),
         base_f  = factor(base, levels = c("A","C","G","T","-",
                                           setdiff(unique(base),
                                                   c("A","C","G","T","-")))))

eyeball_palette <- c("A" = "#1b9e77", "C" = "#377eb8", "G" = "#e7298a",
                     "T" = "#d95f02", "-" = "#888888")
p_eye <- ggplot(eyeball_df, aes(x = pos, y = ref_id, fill = base_f)) +
  geom_tile() +
  facet_wrap(~ label_f, scales = "free_x", ncol = 1) +
  scale_fill_manual(values = eyeball_palette, na.value = "white",
                    name = "base") +
  labs(x = "Alignment position", y = "Reference",
       title = sprintf("%s — top/bottom candidate-window eyeball", opt$marker)) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 6),
        strip.text  = element_text(size = 8, hjust = 0))

ggsave(opt[["eyeball-png"]], p_eye, width = 10, height = 12, dpi = 300)
ggsave(opt[["eyeball-svg"]], p_eye, width = 10, height = 12)
cat(sprintf("[entropy_audit:%s] eyeball figures written\n", opt$marker))

# ---------------------------------------------------------------------------
# 5. Bootstrap CI per position (18S only)
# ---------------------------------------------------------------------------
if (!is.na(opt[["bootstrap-tsv"]])) {
  B  <- as.integer(opt[["bootstrap-n"]])
  cat(sprintf("[entropy_audit:%s] bootstrap n=%d resamples, seed=%s\n",
              opt$marker, B, opt$seed))
  boot_mat <- matrix(NA_real_, nrow = aln_w, ncol = B)
  for (b in seq_len(B)) {
    idx <- sample.int(n_seqs, n_seqs, replace = TRUE)
    sub <- aln[idx]
    cm_b <- consensusMatrix(sub, as.prob = FALSE)
    acgt_b <- cm_b[intersect(c("A","C","G","T"), rownames(cm_b)), , drop = FALSE]
    boot_mat[, b] <- apply(acgt_b, 2, shannon_acgt)
  }
  bootstrap <- tibble(
    pos        = seq_len(aln_w),
    entropy_mean    = rowMeans(boot_mat, na.rm = TRUE),
    entropy_lower95 = apply(boot_mat, 1, quantile, probs = 0.025, na.rm = TRUE),
    entropy_upper95 = apply(boot_mat, 1, quantile, probs = 0.975, na.rm = TRUE),
    entropy_point   = recompute$entropy
  )
  fwrite(bootstrap, opt[["bootstrap-tsv"]], sep = "\t")

  p_boot <- ggplot(bootstrap, aes(x = pos)) +
    geom_ribbon(aes(ymin = entropy_lower95, ymax = entropy_upper95),
                alpha = 0.25, fill = "steelblue") +
    geom_line(aes(y = entropy_point), linewidth = 0.4) +
    labs(x = "Alignment position", y = "Shannon entropy (A/C/G/T)",
         title = sprintf("%s — bootstrap 95%% CI band (B = %d resamples)",
                         opt$marker, B)) +
    theme_bw()
  ggsave(opt[["bootstrap-png"]], p_boot, width = 10, height = 4, dpi = 300)
  ggsave(opt[["bootstrap-svg"]], p_boot, width = 10, height = 4)
  cat(sprintf("[entropy_audit:%s] bootstrap figures written\n", opt$marker))
}

cat(sprintf("[entropy_audit:%s] DONE\n", opt$marker))
