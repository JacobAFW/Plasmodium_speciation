# 08_primer_design — Phase 5 Step 3: primer3-driven primer-pair selection

Turns the end-anchored candidate windows in
`outputs/entropy/{mit,18S}.primer_candidates.tsv` into concrete primer pairs,
using primer3 constrained to those windows. This is the step that has been
open since May; it now proceeds
against the end-anchored candidates only (Imwong-forward chase is settled).

## What the rule does

Two `primer3_*` rules (one per marker) plus a `multiplex_compatibility` rule
and a convenience `primer_design_all` target.

Per marker:

1. **Template.** Build an **ungapped majority-base consensus** of
   `outputs/alignment/ma_{marker}.target.fasta`. Every alignment column
   that is *not* majority-gap becomes exactly one template position — the
   aln↔template map is bijective by construction, so a candidate window at
   aln `[a, b]` projects onto template positions without single-reference
   insertion drift.

2. **Blocks.** Load `primer_candidates.tsv`; collapse the heavily-overlapping
   sliding-window rows for each `end_region` into disjoint (aln_start,
   aln_end) blocks. `left` blocks become forward-primer zones; `right`
   blocks become reverse-primer zones. Project each block onto template
   coordinates.

3. **primer3.** For every (left_block × right_block) combination whose
   achievable amplicon span is inside the per-marker product-size range,
   invoke `primer3-py.bindings.design_primers(...)` with
   `SEQUENCE_PRIMER_PAIR_OK_REGION_LIST` set to
   `[[left_start, left_len, right_start, right_len]]` so primer3 is
   constrained to that combo's two zones. Ask for `PRIMER_NUM_RETURN=20`
   per combo.

4. **Thermodynamics.** For each returned pair, compute:
   - `f_hairpin_dG`, `r_hairpin_dG` via `primer3.calc_hairpin(seq).dg`
   - `pair_dimer_dG` via `primer3.calc_heterodimer(f_seq, r_seq).dg`
   All in kcal/mol (primer3-py returns cal/mol; we scale ÷1000).

5. **Emit** three per-marker TSVs:
   - `{marker}_pairs.tsv` — every returned pair.
   - `{marker}_pairs_filtered.tsv` — pairs whose hairpin ΔG ≥ −6 kcal/mol
     on both primers, heterodimer ΔG ≥ −8 kcal/mol on the pair, and Tm
     inside [58, 62] °C. Panel-rule-of-thumb cutoffs.
   - `{marker}_pairs_unfilled.tsv` — one row per (left_block, right_block)
     combo that yielded zero pairs, with primer3's own EXPLAIN strings for
     `PRIMER_LEFT_EXPLAIN`, `PRIMER_RIGHT_EXPLAIN`, `PRIMER_PAIR_EXPLAIN`.
     Satisfies HANDOFF Step 3's "if a window yields none, log the failing
     primer3 reason" acceptance check.

6. **Multiplex.** `multiplex_compatibility.tsv` — one row per
   (pair_a, pair_b) across BOTH markers' filtered pairs. Columns:
   `tm_spread` (max − min of the four primers' Tm), `cross_dimer_dG`
   (worst ΔG among the 4 cross-pair oligo pairings), `amplicon_len_diff`.

## Why python + primer3-py rather than an R + CLI wrapper

The prompt permits either. Python + `primer3-py` was picked because the
schema requires per-primer hairpin ΔG and per-pair heterodimer ΔG —
`primer3-py` exposes these as one-liner calls
(`calc_hairpin(seq).dg`, `calc_heterodimer(a,b).dg`) that share primer3's
thermodynamic core with the design call. Doing the same via the
primer3 CLI would need a separate `PRIMER_TASK=check_primers` pass per
primer (~4× the calls) and BOULDER-IO parsing for values that the Python
bindings return as attributes. The CLI is still installed via
`envs/install.sh` for Step 5's `mfeprimer` neighbourhood.

## Per-marker product-size range (informed caveats)

I derived the achievable product-size envelope from the actual left/right
candidate blocks *on the template* before picking a range, per the prompt's
"informed caveat" instruction. Numbers below are from the dry-run block
projection on the current `outputs/entropy/*.primer_candidates.tsv`.

- **MIT.** `[5000, 5900]` bp. The HANDOFF nominates 4000–5000 bp, but the
  actual end-anchored blocks — left aln 749–1048 (single block) and right
  aln 6357–6531 + 6605–6656 (two blocks) — give an achievable amplicon
  envelope of **5230–5824 bp**. Setting 4000–5000 would return zero pairs.
  Widened to 5000–5900 to cover the envelope. The "long-amplicon flanking
  primers spanning the MIT core" intent is preserved; the 4–5 kb figure
  from the HANDOFF is best read as an order-of-magnitude target, not a
  hard bound. **This is a HANDOFF-noted deviation and belongs in the
  the project session log.**

- **18S.** `[1200, 1700]` bp. Core is aln 317–2218 (`core_len = 1902`).
  Achievable envelope from the two left blocks (aln 317–368, 487–616) and
  two right blocks (aln 1988–2128, 2194–2218) is **1209–1682 bp**. Set to
  1200–1700 to cover it. Corresponds to "the full 18S gene flanked by
  its two conserved-end anchors".

Ranges are locked in `_PRODUCT_MIN` / `_PRODUCT_MAX` inside the `.smk`.
Overriding for either marker requires editing those two lines.

## Coordinate conventions

- The candidate TSV emits **1-based inclusive** alignment coords (matches
  `seq_entropy.py`).
- The template is **1-based inclusive**; `f_start`/`f_end`/`r_start`/`r_end`
  in the pairs TSVs are template coordinates.
- primer3-internal coords are 0-based; conversions happen inside
  `run_primer3.py` and are not surfaced.

## Outputs

Under `outputs/primer_design/`:

| File                                | Schema                                                   |
|-------------------------------------|----------------------------------------------------------|
| `{marker}_pairs.tsv`                | `pair_id, marker, window_id, pair_idx, f_seq, r_seq, f_start, f_end, r_start, r_end, amplicon_len, f_tm, r_tm, f_gc, r_gc, f_hairpin_dG, r_hairpin_dG, pair_dimer_dG, pair_penalty` |
| `{marker}_pairs_filtered.tsv`       | subset of the above                                     |
| `{marker}_pairs_unfilled.tsv`       | `window_id, reason, n_pairs, left_explain, right_explain, pair_explain` |
| `{marker}_template.fasta`           | the majority-consensus template used (single record)    |
| `multiplex_compatibility.tsv`       | `pair_a, pair_b, marker_a, marker_b, tm_spread, cross_dimer_dG, amplicon_len_diff` |

## Not done here

- **Step 4** (off-target genome manifest) and **Step 5** (MFEprimer specificity)
  live in later rule files. `mfeprimer` is installed in this pass but not
  invoked here.
- Down-selection to a *panel* (small set of chosen pairs) — the filtered TSV
  is the design candidate set. Panel selection is a separate step.
- Cross-checking that every `combined_resolved = TRUE` pair in
  `outputs/cross_species/resolution_table.tsv` remains resolvable by the
  chosen panel — spot-checked in the driver's stderr log; the systematic
  check waits on Step 5 in-silico PCR.

## How to run

```bash
source envs/activate.sh
snakemake --snakefile workflow/Snakefile --cores 4 primer_design_all
```

Or drive the script directly (per-marker):

```bash
python scripts/python/run_primer3.py \
    --marker mit \
    --aln outputs/alignment/ma_mit.target.fasta \
    --candidates outputs/entropy/mit.primer_candidates.tsv \
    --pairs-out outputs/primer_design/mit_pairs.tsv \
    --filtered-out outputs/primer_design/mit_pairs_filtered.tsv \
    --unfilled-out outputs/primer_design/mit_pairs_unfilled.tsv \
    --template-out outputs/primer_design/mit_template.fasta \
    --product-min 4000 --product-max 5000
```
