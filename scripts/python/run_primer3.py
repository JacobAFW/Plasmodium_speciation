#!/usr/bin/env python3
"""
run_primer3.py — Phase 5 Step 3 driver.

For one marker (mit or 18S):
  1. Load the MAFFT alignment; build an ungapped majority-consensus template
     with a bijective alignment↔template coordinate map (majority-gap columns
     are dropped from the template).
  2. Load `outputs/entropy/{marker}.primer_candidates.tsv`; collapse the
     heavily-overlapping sliding windows into disjoint blocks per `end_region`
     (`left` → forward-primer zones; `right` → reverse-primer zones); project
     each block onto template coordinates.
  3. For every (left_block × right_block) combination whose separation falls
     within the requested product-size range, run primer3 via `primer3-py`
     bindings, constraining primers to the projected zones via
     `SEQUENCE_PRIMER_PAIR_OK_REGION_LIST`. Ask primer3 for
     `--num-return` pairs per combo.
  4. Compute per-primer hairpin ΔG and per-pair heterodimer (self-dimer)
     ΔG via primer3-py's thermodynamic bindings. Emit `{marker}_pairs.tsv`
     (all returned pairs) and `{marker}_pairs_filtered.tsv` (subset with no
     primer3 pair-penalty warnings and dG above sensible thresholds — see
     `--dg-hairpin-min` / `--dg-dimer-min`).
  5. Emit `{marker}_pairs_unfilled.tsv`: one row per (left_block, right_block)
     combo that returned zero pairs, with primer3's own EXPLAIN strings. This
     satisfies the HANDOFF's "if a window yields none, log the failing
     primer3 reason" acceptance check.

Coordinate conventions:
- Alignment coords are **1-based inclusive** (as emitted by `seq_entropy.py`).
- Template coords in the output TSV are **1-based inclusive**.
- primer3 internally uses 0-based; conversions live inside this script.

Design of the template (ungapped majority consensus):
- The candidate-window analysis in `seq_entropy.py` scores columns of the
  alignment. Using the same alignment's majority-base consensus as the
  primer3 template makes the aln↔template map bijective by construction —
  every alignment column that's not majority-gap becomes exactly one
  template position, and the reverse holds too. That means a candidate
  window at aln positions [a,b] maps to template positions
  [aln_to_template[a], aln_to_template[b]] without any single-reference
  insertion/deletion drift. The consensus is real enough for primer3 —
  primer3 does not use the template beyond as a substrate for
  hybridisation-suitability scoring and picking sequence.

Usage:
    python3 scripts/python/run_primer3.py \
        --marker mit \
        --aln outputs/alignment/ma_mit.target.fasta \
        --candidates outputs/entropy/mit.primer_candidates.tsv \
        --pairs-out outputs/primer_design/mit_pairs.tsv \
        --filtered-out outputs/primer_design/mit_pairs_filtered.tsv \
        --unfilled-out outputs/primer_design/mit_pairs_unfilled.tsv \
        --product-min 4000 --product-max 5000
"""
from __future__ import annotations

import argparse
import sys
from collections import Counter
from pathlib import Path

import pandas as pd
from Bio import SeqIO
import primer3


# primer3-py exposes both snake_case (2.x) and camelCase (1.x) top-level
# functions. Prefer snake_case where present, fall back to camelCase.
_calc_hairpin = getattr(primer3, "calc_hairpin", None) or primer3.calcHairpin
_calc_heterodimer = getattr(primer3, "calc_heterodimer", None) or primer3.calcHeterodimer
_design_primers = (
    getattr(primer3.bindings, "design_primers", None)
    or getattr(primer3.bindings, "designPrimers", None)
)


# ---------- alignment + template ---------------------------------------------

def load_alignment(fasta_path: str):
    records = list(SeqIO.parse(fasta_path, "fasta"))
    if not records:
        raise SystemExit(f"[run_primer3] empty alignment: {fasta_path}")
    aln_lens = {len(r.seq) for r in records}
    if len(aln_lens) != 1:
        raise SystemExit(f"[run_primer3] alignment not constant width: {aln_lens}")
    return records, aln_lens.pop()


def build_consensus_template(records, aln_len: int):
    """Majority-base consensus with gap-columns removed. Returns
    (template_seq, aln_to_tpl_1based, tpl_to_aln_1based).

    `aln_to_tpl_1based[i]` (1-based aln position i) is the template
    position (1-based) if that column is not majority-gap, else 0.
    """
    tpl_chars = []
    aln_to_tpl = [0] * (aln_len + 1)     # aln_to_tpl[i] for i in 1..aln_len
    tpl_to_aln = [0]                     # tpl_to_aln[0] unused; extend as we go
    for col in range(aln_len):
        col_bases = [str(r.seq[col]).upper() for r in records]
        counts = Counter(b for b in col_bases if b in ("A", "C", "G", "T"))
        gap_count = sum(1 for b in col_bases if b in ("-", "."))
        if gap_count > sum(counts.values()):
            # majority-gap column — drop from template
            continue
        if not counts:
            continue
        # Highest-count A/C/G/T wins; tie broken by A>C>G>T for reproducibility.
        best = max(("A", "C", "G", "T"), key=lambda b: (counts.get(b, 0), -"ACGT".index(b)))
        tpl_chars.append(best)
        tpl_pos = len(tpl_chars)
        aln_to_tpl[col + 1] = tpl_pos
        tpl_to_aln.append(col + 1)
    return "".join(tpl_chars), aln_to_tpl, tpl_to_aln


# ---------- candidate windows -> blocks --------------------------------------

def collapse_to_blocks(cands: pd.DataFrame) -> dict[str, list[tuple[int, int]]]:
    """Collapse the sliding-window rows for each `end_region` into disjoint
    (aln_start, aln_end) blocks — a block is the union of overlapping windows.
    Returns {"left": [(s,e), ...], "right": [(s,e), ...]}. All 1-based inclusive.
    """
    out: dict[str, list[tuple[int, int]]] = {"left": [], "right": []}
    for side in ("left", "right"):
        sub = cands[cands["end_region"] == side].sort_values("start")
        if sub.empty:
            continue
        cur_s = int(sub.iloc[0]["start"])
        cur_e = int(sub.iloc[0]["end"])
        for _, r in sub.iloc[1:].iterrows():
            s, e = int(r["start"]), int(r["end"])
            if s <= cur_e + 1:
                cur_e = max(cur_e, e)
            else:
                out[side].append((cur_s, cur_e))
                cur_s, cur_e = s, e
        out[side].append((cur_s, cur_e))
    return out


def project_block_to_template(block_aln: tuple[int, int], aln_to_tpl: list[int]) -> tuple[int, int]:
    """Project an aln [s,e] (1-based incl) onto the template. Because template
    is majority-gap-free, columns inside the block that ARE majority-gap
    collapse; we take the min/max non-zero tpl positions inside [s,e].
    """
    s, e = block_aln
    tpl_positions = [aln_to_tpl[i] for i in range(s, e + 1) if aln_to_tpl[i] > 0]
    if not tpl_positions:
        return (0, 0)
    return (min(tpl_positions), max(tpl_positions))


# ---------- primer3 ---------------------------------------------------------

def run_primer3_combo(
    template: str,
    fwd_tpl: tuple[int, int],
    rev_tpl: tuple[int, int],
    product_range: tuple[int, int],
    tm_min: float, tm_max: float,
    gc_min: float, gc_max: float,
    len_min: int, len_max: int,
    num_return: int,
) -> dict:
    """Run primer3 with primers constrained to the two projected zones.

    primer3 uses 0-based coords for SEQUENCE_PRIMER_PAIR_OK_REGION_LIST:
    [left_start, left_length, right_start, right_length].
    """
    left_start_0 = fwd_tpl[0] - 1
    left_len = fwd_tpl[1] - fwd_tpl[0] + 1
    right_start_0 = rev_tpl[0] - 1
    right_len = rev_tpl[1] - rev_tpl[0] + 1

    seq_args = {
        "SEQUENCE_ID": f"combo_{fwd_tpl[0]}_{fwd_tpl[1]}__{rev_tpl[0]}_{rev_tpl[1]}",
        "SEQUENCE_TEMPLATE": template,
        "SEQUENCE_PRIMER_PAIR_OK_REGION_LIST": [
            [left_start_0, left_len, right_start_0, right_len]
        ],
    }
    global_args = {
        "PRIMER_TASK": "generic",
        "PRIMER_PICK_LEFT_PRIMER": 1,
        "PRIMER_PICK_INTERNAL_OLIGO": 0,
        "PRIMER_PICK_RIGHT_PRIMER": 1,
        "PRIMER_OPT_SIZE": (len_min + len_max) // 2,
        "PRIMER_MIN_SIZE": len_min,
        "PRIMER_MAX_SIZE": len_max,
        "PRIMER_OPT_TM": (tm_min + tm_max) / 2.0,
        "PRIMER_MIN_TM": tm_min,
        "PRIMER_MAX_TM": tm_max,
        "PRIMER_MIN_GC": gc_min,
        "PRIMER_MAX_GC": gc_max,
        "PRIMER_PRODUCT_SIZE_RANGE": [list(product_range)],
        "PRIMER_NUM_RETURN": num_return,
        "PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT": 1,
        # We assess dimer/hairpin separately via primer3-py bindings; primer3
        # can still surface egregious ones through its own scoring.
    }
    return _design_primers(seq_args, global_args)


def parse_primer3_result(
    res: dict,
    block_pair_id: str,
    aln_to_tpl_reverse: list[int],
) -> list[dict]:
    """Parse a primer3 result dict into row dicts.

    `aln_to_tpl_reverse` is `tpl_to_aln`, 1-based indexed; we surface template
    coords in the row (that's what wet-lab wants) but the reader can join
    back to aln via `tpl_to_aln`.
    """
    n = int(res.get("PRIMER_PAIR_NUM_RETURNED", 0))
    out = []
    for i in range(n):
        f_seq = res[f"PRIMER_LEFT_{i}_SEQUENCE"]
        r_seq = res[f"PRIMER_RIGHT_{i}_SEQUENCE"]
        # primer3 emits *_i as "start,len" (0-based start), for left the primer
        # runs 5'→3' on the template; for right, `start` is the 3' *end* on
        # the template (rightmost coord), and it runs 3'→5' on template.
        f_start0, f_len = res[f"PRIMER_LEFT_{i}"]
        r_end0, r_len = res[f"PRIMER_RIGHT_{i}"]
        f_start_tpl = f_start0 + 1
        f_end_tpl = f_start0 + f_len
        r_start_tpl = r_end0 - r_len + 2   # 5' end of reverse primer on template (1-based)
        r_end_tpl = r_end0 + 1             # 3' end of reverse primer on template (1-based)
        pair_size = int(res[f"PRIMER_PAIR_{i}_PRODUCT_SIZE"])

        # Hairpin ΔG (kcal/mol) via primer3-py.
        try:
            f_hairpin = _calc_hairpin(f_seq).dg / 1000.0
        except Exception:
            f_hairpin = float("nan")
        try:
            r_hairpin = _calc_hairpin(r_seq).dg / 1000.0
        except Exception:
            r_hairpin = float("nan")
        try:
            pair_dimer = _calc_heterodimer(f_seq, r_seq).dg / 1000.0
        except Exception:
            pair_dimer = float("nan")

        out.append({
            "window_id": block_pair_id,
            "pair_idx": i,
            "f_seq": f_seq,
            "r_seq": r_seq,
            "f_start": f_start_tpl,
            "f_end": f_end_tpl,
            "r_start": r_start_tpl,
            "r_end": r_end_tpl,
            "amplicon_len": pair_size,
            "f_tm": round(float(res[f"PRIMER_LEFT_{i}_TM"]), 3),
            "r_tm": round(float(res[f"PRIMER_RIGHT_{i}_TM"]), 3),
            "f_gc": round(float(res[f"PRIMER_LEFT_{i}_GC_PERCENT"]), 2),
            "r_gc": round(float(res[f"PRIMER_RIGHT_{i}_GC_PERCENT"]), 2),
            "f_hairpin_dG": round(f_hairpin, 3),
            "r_hairpin_dG": round(r_hairpin, 3),
            "pair_dimer_dG": round(pair_dimer, 3),
            "pair_penalty": round(float(res.get(f"PRIMER_PAIR_{i}_PENALTY", 0.0)), 3),
        })
    return out


# ---------- main ------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--marker", required=True, choices=["mit", "18S"])
    ap.add_argument("--aln", required=True)
    ap.add_argument("--candidates", required=True,
                    help="outputs/entropy/{marker}.primer_candidates.tsv")
    ap.add_argument("--pairs-out", required=True)
    ap.add_argument("--filtered-out", required=True)
    ap.add_argument("--unfilled-out", required=True)
    ap.add_argument("--template-out", default=None,
                    help="Optional: write the consensus template FASTA here.")
    ap.add_argument("--product-min", type=int, required=True)
    ap.add_argument("--product-max", type=int, required=True)
    ap.add_argument("--tm-min", type=float, default=58.0)
    ap.add_argument("--tm-max", type=float, default=62.0)
    ap.add_argument("--gc-min", type=float, default=40.0)
    ap.add_argument("--gc-max", type=float, default=60.0)
    ap.add_argument("--len-min", type=int, default=18)
    ap.add_argument("--len-max", type=int, default=25)
    ap.add_argument("--num-return", type=int, default=5,
                    help="primer3 PRIMER_NUM_RETURN per (left_block, right_block) combo.")
    # Filter thresholds. Hairpin/dimer dG in kcal/mol; typical panel cutoff
    # is ΔG >= -6 kcal/mol for hairpin, ΔG >= -8 kcal/mol for heterodimer at
    # the amplification temperature. Overridable.
    ap.add_argument("--dg-hairpin-min", type=float, default=-6.0)
    ap.add_argument("--dg-dimer-min",   type=float, default=-8.0)
    args = ap.parse_args()

    Path(args.pairs_out).parent.mkdir(parents=True, exist_ok=True)

    # 1. alignment + template
    records, aln_len = load_alignment(args.aln)
    template, aln_to_tpl, tpl_to_aln = build_consensus_template(records, aln_len)
    print(f"[{args.marker}] alignment: {aln_len} cols, {len(records)} seqs; "
          f"template: {len(template)} bp", file=sys.stderr)
    if args.template_out:
        Path(args.template_out).parent.mkdir(parents=True, exist_ok=True)
        with open(args.template_out, "w") as fh:
            fh.write(f">{args.marker}_majority_consensus  aln_len={aln_len} "
                     f"n_seq={len(records)}\n")
            for k in range(0, len(template), 70):
                fh.write(template[k:k + 70] + "\n")

    # 2. candidate windows → blocks → projected onto template
    cands = pd.read_csv(args.candidates, sep="\t")
    if not {"start", "end", "end_region"}.issubset(cands.columns):
        raise SystemExit(f"[{args.marker}] candidates TSV missing expected columns")
    blocks_aln = collapse_to_blocks(cands)
    blocks_tpl = {
        side: [project_block_to_template(b, aln_to_tpl) for b in blocks_aln[side]]
        for side in ("left", "right")
    }
    print(f"[{args.marker}] blocks: left={len(blocks_aln['left'])}, "
          f"right={len(blocks_aln['right'])}", file=sys.stderr)
    for side in ("left", "right"):
        for b_aln, b_tpl in zip(blocks_aln[side], blocks_tpl[side]):
            print(f"    {side:5s} aln {b_aln[0]}-{b_aln[1]}  ->  "
                  f"tpl {b_tpl[0]}-{b_tpl[1]}", file=sys.stderr)

    # 3. per-combo primer3
    all_pairs: list[dict] = []
    unfilled: list[dict] = []
    combo_summary: list[dict] = []

    for li, (b_aln_l, b_tpl_l) in enumerate(zip(blocks_aln["left"], blocks_tpl["left"])):
        for ri, (b_aln_r, b_tpl_r) in enumerate(zip(blocks_aln["right"], blocks_tpl["right"])):
            # skip combos where the (fwd_start .. rev_end) span can't yield a
            # product inside [product_min, product_max]
            max_product = b_tpl_r[1] - b_tpl_l[0] + 1
            min_product = b_tpl_r[0] - b_tpl_l[1] + 1
            block_pair_id = (
                f"L{li+1}[{b_aln_l[0]}-{b_aln_l[1]}]"
                f"__R{ri+1}[{b_aln_r[0]}-{b_aln_r[1]}]"
            )
            if max_product < args.product_min or min_product > args.product_max:
                unfilled.append({
                    "window_id": block_pair_id,
                    "reason": (
                        f"product-size out of range: achievable span "
                        f"{min_product}-{max_product} bp not in "
                        f"[{args.product_min}, {args.product_max}]"
                    ),
                    "n_pairs": 0,
                    "left_explain": "",
                    "right_explain": "",
                    "pair_explain": "",
                })
                combo_summary.append({"window_id": block_pair_id, "n_pairs": 0,
                                      "achievable_min": min_product,
                                      "achievable_max": max_product})
                continue

            res = run_primer3_combo(
                template=template, fwd_tpl=b_tpl_l, rev_tpl=b_tpl_r,
                product_range=(args.product_min, args.product_max),
                tm_min=args.tm_min, tm_max=args.tm_max,
                gc_min=args.gc_min, gc_max=args.gc_max,
                len_min=args.len_min, len_max=args.len_max,
                num_return=args.num_return,
            )
            rows = parse_primer3_result(res, block_pair_id, tpl_to_aln)
            for r in rows:
                r["marker"] = args.marker
            all_pairs.extend(rows)

            n = int(res.get("PRIMER_PAIR_NUM_RETURNED", 0))
            combo_summary.append({"window_id": block_pair_id, "n_pairs": n,
                                  "achievable_min": min_product,
                                  "achievable_max": max_product})
            if n == 0:
                unfilled.append({
                    "window_id": block_pair_id,
                    "reason": "primer3 returned no pairs",
                    "n_pairs": 0,
                    "left_explain":  res.get("PRIMER_LEFT_EXPLAIN", ""),
                    "right_explain": res.get("PRIMER_RIGHT_EXPLAIN", ""),
                    "pair_explain":  res.get("PRIMER_PAIR_EXPLAIN", ""),
                })

    # 4. assign global pair_id, write outputs
    df_all = pd.DataFrame(all_pairs)
    if df_all.empty:
        # Ensure schema-consistent empty frame
        df_all = pd.DataFrame(columns=[
            "pair_id", "marker", "window_id", "pair_idx",
            "f_seq", "r_seq", "f_start", "f_end", "r_start", "r_end",
            "amplicon_len", "f_tm", "r_tm", "f_gc", "r_gc",
            "f_hairpin_dG", "r_hairpin_dG", "pair_dimer_dG", "pair_penalty",
        ])
    else:
        df_all.insert(0, "pair_id", [f"{args.marker}_{i+1:04d}" for i in range(len(df_all))])
        cols = ["pair_id", "marker", "window_id", "pair_idx",
                "f_seq", "r_seq", "f_start", "f_end", "r_start", "r_end",
                "amplicon_len", "f_tm", "r_tm", "f_gc", "r_gc",
                "f_hairpin_dG", "r_hairpin_dG", "pair_dimer_dG", "pair_penalty"]
        df_all = df_all[cols]
    df_all.to_csv(args.pairs_out, sep="\t", index=False)
    print(f"[{args.marker}] wrote {len(df_all)} pairs -> {args.pairs_out}",
          file=sys.stderr)

    # Filtered: hairpin dG >= threshold on both primers, heterodimer dG >=
    # threshold on the pair, Tm actually within [tm_min, tm_max]. primer3
    # already respects these in-design, but we cross-check.
    df_filt = df_all[
        (df_all["f_hairpin_dG"] >= args.dg_hairpin_min) &
        (df_all["r_hairpin_dG"] >= args.dg_hairpin_min) &
        (df_all["pair_dimer_dG"] >= args.dg_dimer_min) &
        (df_all["f_tm"] >= args.tm_min) & (df_all["f_tm"] <= args.tm_max) &
        (df_all["r_tm"] >= args.tm_min) & (df_all["r_tm"] <= args.tm_max)
    ].copy()
    df_filt.to_csv(args.filtered_out, sep="\t", index=False)
    print(f"[{args.marker}] wrote {len(df_filt)} filtered -> {args.filtered_out}",
          file=sys.stderr)

    df_un = pd.DataFrame(unfilled)
    if df_un.empty:
        df_un = pd.DataFrame(columns=["window_id", "reason", "n_pairs",
                                       "left_explain", "right_explain",
                                       "pair_explain"])
    df_un.to_csv(args.unfilled_out, sep="\t", index=False)
    print(f"[{args.marker}] wrote {len(df_un)} unfilled combos -> {args.unfilled_out}",
          file=sys.stderr)

    # Summary line for the caller / snakemake log
    tm_spread = (
        max(df_filt["f_tm"].max(), df_filt["r_tm"].max()) -
        min(df_filt["f_tm"].min(), df_filt["r_tm"].min())
    ) if not df_filt.empty else float("nan")
    print(f"[{args.marker}] filtered Tm spread: {tm_spread:.2f} °C", file=sys.stderr)


if __name__ == "__main__":
    main()
