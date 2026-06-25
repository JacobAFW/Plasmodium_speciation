#!/usr/bin/env Rscript
# plot_entropy.R — render entropy-vs-position and coverage-vs-position panels
# for one marker, with primer-candidate windows and internal high-entropy
# blocks highlighted.
#
# Provenance: refactor of scripts/legacy/speciation_long/scripts/plot_entropy.R
# into a standalone, argv-driven script. The plotting helpers are unchanged
# in spirit; the caller now supplies all input/output paths and parameters.
#
# Usage:
#   plot_entropy.R --marker mit \
#                  --per-position outputs/entropy/mit.per_position.tsv \
#                  --windows      outputs/entropy/mit.windows.tsv \
#                  --primers      outputs/entropy/mit.primer_candidates.tsv \
#                  --entropy-png  reports/figures/mit_entropy_vs_position.png \
#                  --entropy-svg  reports/figures/mit_entropy_vs_position.svg \
#                  --coverage-png reports/figures/mit_coverage_vs_position.png \
#                  --coverage-svg reports/figures/mit_coverage_vs_position.svg \
#                  [--internal-window 50 --internal-thr 0.6 --internal-min-run 3]

suppressPackageStartupMessages({
  library(tidyverse)
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
  marker                = NA_character_,
  `per-position`        = NA_character_,
  windows               = NA_character_,
  primers               = NA_character_,
  `entropy-png`         = NA_character_,
  `entropy-svg`         = NA_character_,
  `coverage-png`        = NA_character_,
  `coverage-svg`        = NA_character_,
  `internal-window`     = "50",
  `internal-thr`        = "0.6",
  `internal-min-run`    = "3"
))
required <- c("marker","per-position","windows","primers",
              "entropy-png","entropy-svg","coverage-png","coverage-svg")
for (k in required) if (is.na(opt[[k]])) stop(sprintf("--%s required", k))

INT_WIN  <- as.integer(opt[["internal-window"]])
INT_THR  <- as.numeric(opt[["internal-thr"]])
INT_MINR <- as.integer(opt[["internal-min-run"]])

# ---- helpers (1:1 from legacy plot_entropy.R) ----
read_perpos <- function(path) {
  read_tsv(path, show_col_types = FALSE) |> mutate(pos = as.integer(pos))
}
read_windows <- function(path) {
  read_tsv(path, show_col_types = FALSE) |>
    mutate(start = as.integer(start),
           end   = as.integer(end),
           mid   = (start + end) / 2)
}
merge_blocks <- function(df, gap = 0L) {
  if (nrow(df) == 0) return(df)
  df <- df |> arrange(start, end)
  blocks <- list(); cur_s <- df$start[1]; cur_e <- df$end[1]
  for (i in 2:nrow(df)) {
    s <- df$start[i]; e <- df$end[i]
    if (s <= (cur_e + gap + 1)) {
      cur_e <- max(cur_e, e)
    } else {
      blocks[[length(blocks) + 1]] <- tibble(start = cur_s, end = cur_e)
      cur_s <- s; cur_e <- e
    }
  }
  blocks[[length(blocks) + 1]] <- tibble(start = cur_s, end = cur_e)
  bind_rows(blocks)
}
internal_entropy_blocks <- function(win, window_size, entropy_thr, min_consecutive) {
  w <- win |> filter(window == window_size, end_region == "none") |> arrange(start)
  if (nrow(w) == 0) return(tibble(start = integer(), end = integer()))
  w <- w |> mutate(is_high = mean_entropy >= entropy_thr)
  idx <- which(w$is_high)
  if (length(idx) == 0) return(tibble(start = integer(), end = integer()))
  grp <- cumsum(c(TRUE, diff(idx) != 1))
  runs <- tibble(i = idx, grp = grp) |>
    group_by(grp) |>
    summarise(n = n(),
              start = min(w$start[i]),
              end   = max(w$end[i]),
              .groups = "drop") |>
    filter(n >= min_consecutive) |>
    select(start, end)
  merge_blocks(runs, gap = 0L)
}
plot_metric <- function(perpos, primer_blocks, internal_blocks, ycol, title) {
  ggplot(perpos, aes(x = pos, y = .data[[ycol]])) +
    geom_rect(data = internal_blocks,
              aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
              inherit.aes = FALSE, alpha = 0.15) +
    geom_rect(data = primer_blocks,
              aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
              inherit.aes = FALSE, alpha = 0.25) +
    geom_line(linewidth = 0.4) +
    labs(x = "Alignment position", y = ycol, title = title) +
    theme_bw()
}

# ---- main ----
perpos <- read_perpos(opt[["per-position"]])
win    <- read_windows(opt$windows)
primer <- read_windows(opt$primers)

primer_blocks <- primer |> select(start, end) |> distinct() |> merge_blocks(gap = 0L)
internal_blocks <- internal_entropy_blocks(win,
                                           window_size     = INT_WIN,
                                           entropy_thr     = INT_THR,
                                           min_consecutive = INT_MINR)

p_entropy <- plot_metric(perpos, primer_blocks, internal_blocks, "entropy",
                         paste0(opt$marker, ": entropy vs position"))
p_cov     <- plot_metric(perpos, primer_blocks, internal_blocks, "coverage",
                         paste0(opt$marker, ": coverage vs position"))

# Save both PNG (for embedding) and SVG (for editable handoff). 300 dpi PNG.
for (out in c(opt[["entropy-png"]], opt[["entropy-svg"]],
              opt[["coverage-png"]], opt[["coverage-svg"]])) {
  dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
}
ggsave(opt[["entropy-png"]],  p_entropy, width = 10, height = 4, dpi = 300)
ggsave(opt[["entropy-svg"]],  p_entropy, width = 10, height = 4)
ggsave(opt[["coverage-png"]], p_cov,     width = 10, height = 4, dpi = 300)
ggsave(opt[["coverage-svg"]], p_cov,     width = 10, height = 4)

cat(sprintf("[plot_entropy:%s] wrote 4 figures to %s\n",
            opt$marker, dirname(opt[["entropy-png"]])))
