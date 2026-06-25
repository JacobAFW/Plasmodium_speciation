# 07_threshold_sensitivity — calibrate the interior-side candidate-window thresholds

One rule:

- `interior_threshold_sweep` — calls `scripts/python/interior_threshold_sweep.py` over a 4×4 `(cov_thr, min_run)` grid per marker (32 cells total). Each cell invokes `seq_entropy.py` to a temp directory, reads the resulting `interior_candidates.tsv`, and records `n_windows`, `n_blocks`, `median_entropy`, `median_coverage`, plus `imwong_recovered` and `imwong_n_windows` for 18S. Output: `outputs/sensitivity/interior_threshold_sweep.tsv`.

## Why this rule file exists

Phase 5 Step 1's audit established that the published Imwong forward primers (18S aln 1683–1740) are not rejected by entropy or coverage — they have full coverage and near-zero entropy in every reference. They are rejected by `seq_entropy.py`'s `end_region` filter, which only admits windows within 300 bp of the alignment core's edges. Step 2 added a parallel interior-side output (`interior_candidates.tsv`) that drops the `end_region` gate; this rule calibrates its thresholds.

## What the sweep found

| Finding | Detail |
|---|---|
| `min_run` is a no-op for the interior output | `min_run` only feeds `find_core_bounds`; the interior output does not depend on the core. At every (marker, cov_thr) row in the grid, varying `min_run ∈ {15, 18, 21, 25}` leaves `n_windows`, `n_blocks`, `median_entropy`, `median_coverage`, and `imwong_recovered` exactly equal to the value at `min_run = 15`. |
| Imwong recovery is unanimous on 18S | 16/16 18S cells show `imwong_recovered = TRUE` with 180 windows overlapping the published-forward footprint. The region clears entropy + coverage at every threshold in the grid. |
| `cov_thr` is the only meaningful axis | n_windows decreases monotonically with `cov_thr` (MIT: 11557 → 11200; 18S: 2438 → 1906). Median entropy stays flat (~0.05 for 18S, ~0.07 for MIT) — the marginal windows that drop out as `cov_thr` rises are low-coverage stragglers, not high-entropy noise. |

## Chosen defaults

`cov_thr = 0.90`, `entropy_thr = 0.20`, `min_run = 50`, `windows = 25,50` — **identical to the design side.** The constraint that hid the Imwong region was structural (end-anchoring), not threshold strictness. Once the structural gate is removed, v1's existing thresholds surface the region cleanly with no need to relax anything.

Headline counts at the chosen defaults:
- 18S: 1,906 windows / 30 contiguous blocks. Imwong region recovered.
- MIT: 11,200 windows / 30 contiguous blocks. (No published-primer target on MIT.)

## Outputs

```
outputs/sensitivity/
└── interior_threshold_sweep.tsv         # 32 rows (2 markers × 4 cov_thr × 4 min_run)
```

## How to run

```bash
source envs/activate.sh
snakemake --snakefile workflow/Snakefile --cores 1 \
  outputs/sensitivity/interior_threshold_sweep.tsv
```

The sweep is single-threaded by design — 32 cells × ~2 s/cell = ~70 s. Snakemake parallelism over a single rule's body wouldn't help.

## Methods

- Grid: `cov_thr ∈ {0.70, 0.80, 0.85, 0.90}` × `min_run ∈ {15, 18, 21, 25}` × marker ∈ {`mit`, `18S`}, with `entropy_thr = 0.20` and `windows = 25,50` held fixed.
- Each cell invokes `python scripts/python/seq_entropy.py` with the cell's parameters to a temp `--outprefix`, reads the resulting `interior_candidates.tsv` and `core_bounds.tsv`, and computes the summary metrics.
- "Imwong recovery" = any interior window overlaps the published forward footprint at 18S aln 1683–1740 (passed as `--imwong-start 1683 --imwong-end 1740`).
- "Blocks" = number of distinct contiguous runs of consecutive admitted starts, computed per window size. A run of 25 consecutive passing 25-bp windows is one block.

## Provenance

`HANDOFF.md` Phase 5 Step 2. Step 1 (`06_entropy_audit.README.md`) ruled out the entropy-method-bug hypothesis; this step ruled out the threshold-too-tight hypothesis. The remaining hypothesis — that the published forwards sit where the design-side analysis doesn't look — is confirmed by the recovery in the interior output.
