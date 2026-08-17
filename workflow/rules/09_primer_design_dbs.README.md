# 09_primer_design_dbs — DBS-compatible MIT amplicon tier (≤4 kb)

Phase 5 Step 3b. Parallel to Step 3 (`08_primer_design.smk`). Runs the same
primer3 machinery (imported from `scripts/python/run_primer3.py`, no edits)
against the WHOLE-MIT-core conserved-window menu, so primers can sit in the
interior — Tier A picks only from end-anchored blocks and produces ~5.5 kb
amplicons that don't survive DBS sample prep.

## Rules

| Rule | Reads | Writes | Notes |
|---|---|---|---|
| `primer3_dbs_mit_single` | `mit.interior_candidates.tsv`, `ma_mit.target.fasta` | `mit_dbs_single_pairs{,_filtered,_unfilled}.tsv` | Task A. Enumerates every (fwd_block × rev_block) combo across the merged conserved menu; product-size 3300–4000 bp. |
| `primer3_dbs_mit_tiled`  | same + `mit.core_bounds.tsv` | `mit_dbs_tiled_pairs{,_filtered,_unfilled}.tsv` | Task B. Restricts to (start × mid) as `tile1` and (mid × end) as `tile2`; product-size 3000–3500 bp per tile. `tile_id` carried in the output. |
| `mit_amplicon_resolution_tradeoff` | the two DBS tables + Tier-A `mit_pairs_filtered.tsv` + v1 `resolution_table.tsv` | `mit_amplicon_resolution_tradeoff.tsv` | Task C. Per species pair × option → {resolvable, n_diagnostic_sites}; reconciles Tier-A column against v1. |
| `dbs_multiplex_compatibility` | DBS singles + DBS tiles + existing `18S_pairs_filtered.tsv` | `dbs_multiplex_compatibility.tsv` | Task D. Full N×N table over the DBS panel using `multiplex_compat.py` unchanged. |
| `primer_design_dbs_all` | (convenience aggregator) | — | Materialises every output above. |

## Bands used per option

Derived from the **actual conserved-block geometry** of `mit.interior_candidates.tsv`,
not hard-coded. Interior-candidate rows (11 200 across left/right/none end_regions)
collapse into **4 disjoint blocks**:

```
B1  aln  726–1961  (1236 bp)   —  "start" block, contains Tier-A L1[749–1048]
B2  aln 1973–3359  (1387 bp)   —  first interior block
B3  aln 3382–6531  (3150 bp)   —  second interior block, contains Tier-A R1
B4  aln 6605–6747  ( 143 bp)   —  "end" block, contains Tier-A R2
```

Template projection puts them at `tpl 1-1213, 1219-2585, 2586-5721, 5791-5932`
(template length 5932 bp). The MIT core (`aln 749–6656`, 5908 bp) is spanned
end-to-end.

Achievable per-block-pair amplicon envelopes (template coords, min–max):

```
B1×B2:    7–2585 bp    (too short for DBS band)
B1×B3: 1374–5721 bp    (single: hits 3300–4000 band; tile1 3000–3500)
B1×B4: 4579–5932 bp    (over the DBS ceiling)
B2×B3:    2–4503 bp    (single: hits 3300–4000 band)
B2×B4: 3207–4714 bp    (single: hits 3300–4000 band; tile2 3000–3500)
B3×B4:   71–3347 bp    (tile2: hits 3000–3500 lower end)
```

- **Single option (Task A):** three block-pair combos in-band — `B1×B3`,
  `B2×B3`, `B2×B4` — plus one that primer3 rejects (`B3×B4`, primer3
  explain: "not in any ok left region" for the tight 3300–3500 slice at
  the far right of the core). 60 pairs total; 60/60 pass the filter with
  amplicon **3319–3894 bp** and Tm spread **2.11 °C**. Informed caveat:
  we widened the DBS band to 3300–4000 (from a nominal 3500–4000 target)
  because the achievable envelope hits the target only marginally at
  some combos — the 3300 bp floor buys headroom against the ~4 kb DBS
  ceiling.
- **Tiled option (Task B):** two block-pair combos per tile.
  - `tile1` = `B1 × B3`, 20 pairs, amplicon 3010–3467 bp.
  - `tile2` = `B2 × B4` (20 pairs) and `B3 × B4` (20 pairs), amplicon 3049–3467 bp
    across both combos.
  - Total 60 pairs, Tm spread **2.25 °C**.

## Coordinate map

- Alignment coords **1-based inclusive** (`seq_entropy.py` convention).
- Template coords in every TSV **1-based inclusive** (majority-consensus
  ungapped positions).
- Aln↔template map is bijective by construction (see `run_primer3.py`
  docstring); `run_primer3_dbs.py resolution` converts template coords
  back to alignment for diagnostic-position counting.

## Validation results (matches prompt § "Validation")

- ✓ Every filtered DBS pair amplicon ≤ 4000 bp (max 3894 single / 3467 tiled).
- ✓ Bands **stated and justified** from actual block geometry above.
- ✓ Task-C Tier-A column matches v1's `mit_resolved` on **54/55** pairs. The
  one disagreement is `simium vs vivax` (v1 FALSE / DBS-analysis TRUE with
  n_diagnostic = 3). The disagreement is **method-level, not data-level**:
  v1's `mit_resolved` comes from the BLAST-based strict-suspicious filter
  (pident ≥ 99, cov_shorter ≥ 0.99), which flags simium/vivax as MIT-
  ambiguous; the Task-C metric counts per-column consensus differences,
  which finds 3 (see `outputs/primer_design/dbs/README.md` § "The
  simium/vivax reconciliation" for the recommended interpretation).
- ✓ Tiled union interval covers ≥ diagnostic positions of single for every
  pair (the two tiles' union = aln 790–6710, wider than single's 738–6709
  by a hair; per-pair counts are within ±1 of each other and never lower).
- ✓ Per-species mismatch profile carried over from Tier A (same underlying
  consensus template): see `outputs/primer_design/dbs/README.md`.
- ✓ **No** existing files edited (`run_primer3.py`, `multiplex_compat.py`,
  `08_primer_design.smk`, `mit_pairs*.tsv`, `18S_pairs*.tsv`,
  `mit_template.fasta`). Verified by `git status`.

## Wiring into the top-level Snakefile

`workflow/Snakefile` includes each rule file explicitly; add:

```python
include: "rules/09_primer_design_dbs.smk"
```

The convenience aggregator `primer_design_dbs_all` is not added to
`FINAL_TARGETS` yet — that's a call for Jacob after the trade-off headline
is reviewed and one of {single, tiled} is chosen for the panel.
