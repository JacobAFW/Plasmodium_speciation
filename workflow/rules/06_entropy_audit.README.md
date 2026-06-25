# 06_entropy_audit — audit the entropy pipeline before trusting downstream Phase 5 work

Three rules:

- `entropy_audit_mit` — independent recomputation of per-position Shannon entropy + coverage from the MAFFT alignment via `Biostrings::consensusMatrix`, diff against the pipeline's `mit.per_position.tsv`, coverage spot-check at five positions (lowest entropy, highest entropy, median entropy, lowest coverage, highest coverage), and a top-3 / bottom-3 candidate-window eyeball heatmap.
- `entropy_audit_18S` — identical to the MIT rule, plus a 100-bootstrap CI band over the 18S panel (the panel has only 17 target-species references after filtering; bootstrap surfaces how unstable entropy estimates are at that n).
- `inspect_imwong_position` — drill-down at the 18S alignment positions where the published Imwong PlasmoM forward primers bind. Outputs per-position entropy + coverage + threshold-admission flags for six (cov_thr, entropy_thr) pairs.

## Why this rule file exists

Chapter 6 of the v1 book surfaced that the two published forward primers (`PlasmoM_N1F`, `PlasmoM_N2F`) sit at 18S aln 1683–1740 with **zero** overlap with v1's entropy-derived candidate windows, even though the published reverse primer overlaps 70 of them. Before threshold sweeping in `07_threshold_sensitivity.smk` we audit whether v1's entropy method is correct at all (acceptance criteria 1–3), and we drill down on the specific mismatch (acceptance criteria 4–5).

## Acceptance criteria

| Check | Tolerance | Result on v1 |
|---|---|---|
| Independent recomputation matches pipeline (MIT) | `\|diff\| < 1e-9` | `\|diff_entropy\|_max = 4.44e-16` ✓ |
| Independent recomputation matches pipeline (18S) | `\|diff\| < 1e-9` | `\|diff_entropy\|_max = 6.66e-16` ✓ |
| Coverage spot-check (MIT, 5 positions) | exact to 4 dp | all positions exact ✓ |
| Coverage spot-check (18S, 5 positions) | exact to 4 dp | all positions exact (one at 1e-16 ε) ✓ |
| 18S bootstrap stability | CI width visualisable | median CI width 0.285 — substantial uncertainty at n=17 ✓ |
| Imwong position drill-down | "by a hair or by a mile" | **structural reject, not a threshold reject** (see below) |
| Eyeball top/bottom windows | visual sanity | figures produced; top windows are conserved blocks, bottom are gappy/variable regions ✓ |

## What the Imwong drill-down found (headline)

The 58-bp 18S region where `PlasmoM_N1F` and `PlasmoM_N2F` bind (aln 1683–1740) has:

- coverage = 1.0 across every position
- mean entropy = 0.011, median = 0, max = 0.32
- 56/58 positions admit at v1's default thresholds (cov ≥ 0.90, ent ≤ 0.20)
- **0/58 positions are in any candidate window**

The position is essentially optimal for primer binding by entropy + coverage. It is rejected because the pipeline's `seq_entropy.py` only admits windows whose `end_region ∈ {"left", "right"}` — i.e., within 300 bp of the core edges. The published forward region sits 1366 bp from the left core edge (317) and 478 bp from the right core edge (2218), so it lands in `end_region = "none"` and is unconditionally excluded regardless of threshold.

**This is a structural decision in the pipeline, not a bug, and not a threshold issue.** Phase 5 Step 2 (threshold sensitivity) won't surface this position because no `(cov_thr, min_run)` relaxation touches the `end_bp` constraint. The right Phase 5 action is to add an interior-windows path (lift the `end_region` restriction or widen `end_bp`) before primer3 selection — covered in Step 2's design.

## Outputs

```
outputs/audit/
├── mit_entropy_recompute.tsv         # independent per-position recomputation
├── mit_entropy_diff.tsv              # per-position diff vs pipeline
├── mit_coverage_spotcheck.tsv        # 5-position manual count vs pipeline
├── mit_window_eyeball.tsv            # alignment slice for top/bottom windows
├── 18S_entropy_recompute.tsv         # same for 18S
├── 18S_entropy_diff.tsv
├── 18S_coverage_spotcheck.tsv
├── 18S_window_eyeball.tsv
├── 18S_entropy_bootstrap.tsv         # mean + 95% CI per position
└── imwong_forward_position.tsv       # ROI drill-down with threshold flags

reports/figures/
├── mit_window_eyeball.{png,svg}      # top/bottom window heatmap
├── 18S_window_eyeball.{png,svg}
└── 18S_entropy_bootstrap.{png,svg}   # CI band overlay on point entropy
```

## How to run

```bash
source envs/activate.sh
snakemake --snakefile workflow/Snakefile --cores 4 \
  outputs/audit/mit_entropy_diff.tsv \
  outputs/audit/18S_entropy_diff.tsv \
  outputs/audit/imwong_forward_position.tsv
```

Or call any single rule by output filename.

## Methods

- Independent entropy: `Biostrings::consensusMatrix(aln, as.prob = FALSE)` → per-position A/C/G/T counts; gap-row counts excluded from coverage denominator; Shannon entropy computed over A/C/G/T frequencies (matches the legacy `seq_entropy.py` definition exactly).
- Bootstrap (18S): 100 resamples with replacement of the 17 reference sequences; entropy recomputed per resample; 2.5th/97.5th percentiles per position give the 95% CI band.
- Imwong drill-down: `Biostrings::matchPattern` with 1-mismatch tolerance against the unaligned references; ungapped → aligned coordinate translation via cumulative non-gap count.

## Provenance

`HANDOFF.md` Phase 5 Step 1. Audit acceptance criteria documented at the rule-level table above.
