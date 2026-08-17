#!/usr/bin/env python3
"""
run_primer3_dbs.py — Phase 5 Step 3b: DBS-compatible MIT amplicon tier (≤4 kb).

Parallel driver to run_primer3.py, for the MIT marker only, targeting
amplicons that fit through the dried-blood-spot (DBS) sample prep, which
constrains ONT amplicons to <=4 kb.

Three subcommands (all read-only w.r.t. the Tier-A ~5.5 kb pipeline):

  single      Task A. Enumerate every (fwd_block × rev_block) combo across
              the WHOLE MIT core (interior + end conserved windows collapsed
              into disjoint blocks), keep those whose achievable amplicon
              overlaps a DBS-safe product-size band (default 3300-4000 bp,
              target 3500-4000 with ~200 bp headroom below), and run
              primer3-py per combo. Emits `_pairs.tsv`, `_filtered.tsv`,
              `_unfilled.tsv`.

  tiled       Task B. Same block menu as `single`, but restricted to two
              block-pair COMBOS that jointly tile the MIT core with a
              modest overlap:
                * tile1  = start-of-core block  × mid-of-core block
                * tile2  = mid-of-core block    × end-of-core block
              Each tile targets ~3000-3500 bp (widened to the achievable
              envelope where the block geometry demands it — see the
              per-tile note logged to stderr). Adds a `tile_id` column.

  resolution  Task C. For each species pair × option
              ({tier_a, single, tiled}), count fixed diagnostic positions
              — alignment columns where species A's well-covered consensus
              base ≠ species B's — falling inside the amplicon interval
              (fwd `start` → rev `end` in alignment coords; for `tiled`,
              the union of both tiles' intervals). Emit the trade-off
              table. Reconcile the `tier_a` column against
              `outputs/cross_species/resolution_table.tsv`.

Reuse: this script imports (does NOT edit) the block-collapse, template
build, block-projection, primer3-call, and result-parse machinery from
`run_primer3.py`. A separate verifier is reading that file concurrently,
so we deliberately keep it untouched.

Coordinate conventions (same as run_primer3.py):
- Alignment coords are 1-based inclusive.
- Template coords in the output TSVs are 1-based inclusive.
- primer3 uses 0-based internally; conversions live in run_primer3.py.
"""
from __future__ import annotations

import argparse
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd
from Bio import SeqIO

# Reuse — do NOT edit run_primer3.py.
from run_primer3 import (
    load_alignment,
    build_consensus_template,
    project_block_to_template,
    run_primer3_combo,
    parse_primer3_result,
)


# ---------- shared block-collapse (across ALL end_regions) --------------------

def collapse_all_to_blocks(cands: pd.DataFrame) -> list[tuple[int, int]]:
    """Collapse EVERY candidate-window row into disjoint (aln_start, aln_end)
    blocks — a block is the union of overlapping windows regardless of
    `end_region`. This is what makes the DBS design differ from Tier A:
    Tier A collapses per-side, this collapses the whole conserved menu so
    interior primers are eligible.
    """
    if cands.empty:
        return []
    sub = cands.sort_values("start").reset_index(drop=True)
    blocks = []
    cur_s = int(sub.iloc[0]["start"])
    cur_e = int(sub.iloc[0]["end"])
    for _, r in sub.iloc[1:].iterrows():
        s, e = int(r["start"]), int(r["end"])
        if s <= cur_e + 1:
            cur_e = max(cur_e, e)
        else:
            blocks.append((cur_s, cur_e))
            cur_s, cur_e = s, e
    blocks.append((cur_s, cur_e))
    return blocks


def load_candidates(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    if not {"start", "end", "end_region"}.issubset(df.columns):
        raise SystemExit(f"[dbs] candidates TSV missing expected columns: {path}")
    return df


def _print_blocks(tag: str, blocks_aln, blocks_tpl, file=sys.stderr):
    for i, (b_aln, b_tpl) in enumerate(zip(blocks_aln, blocks_tpl)):
        print(f"[{tag}] block B{i+1}: aln {b_aln[0]}-{b_aln[1]}  "
              f"({b_aln[1]-b_aln[0]+1} bp)  ->  tpl {b_tpl[0]}-{b_tpl[1]}",
              file=file)


# ---------- primer3 driver for a set of block-pair combos --------------------

def _combo_rows(
    template, aln_to_tpl, tpl_to_aln,
    combos,                     # list of (fwd_block_aln, rev_block_aln, tile_id_or_None)
    product_min, product_max,
    tm_min, tm_max, gc_min, gc_max, len_min, len_max, num_return,
    marker="mit",
):
    all_pairs, unfilled = [], []
    for fwd_aln, rev_aln, tile_id in combos:
        fwd_tpl = project_block_to_template(fwd_aln, aln_to_tpl)
        rev_tpl = project_block_to_template(rev_aln, aln_to_tpl)
        block_pair_id = (
            f"L[{fwd_aln[0]}-{fwd_aln[1]}]__R[{rev_aln[0]}-{rev_aln[1]}]"
        )
        max_product = rev_tpl[1] - fwd_tpl[0] + 1
        min_product = rev_tpl[0] - fwd_tpl[1] + 1
        if max_product < product_min or min_product > product_max:
            unfilled.append({
                "window_id": block_pair_id,
                "tile_id": tile_id or "",
                "reason": (
                    f"product-size out of range: achievable span "
                    f"{min_product}-{max_product} bp not in "
                    f"[{product_min}, {product_max}]"
                ),
                "n_pairs": 0,
                "left_explain": "", "right_explain": "", "pair_explain": "",
            })
            continue

        res = run_primer3_combo(
            template=template, fwd_tpl=fwd_tpl, rev_tpl=rev_tpl,
            product_range=(product_min, product_max),
            tm_min=tm_min, tm_max=tm_max, gc_min=gc_min, gc_max=gc_max,
            len_min=len_min, len_max=len_max, num_return=num_return,
        )
        rows = parse_primer3_result(res, block_pair_id, tpl_to_aln)
        for r in rows:
            r["marker"] = marker
            r["tile_id"] = tile_id or ""
        all_pairs.extend(rows)

        n = int(res.get("PRIMER_PAIR_NUM_RETURNED", 0))
        if n == 0:
            unfilled.append({
                "window_id": block_pair_id,
                "tile_id": tile_id or "",
                "reason": "primer3 returned no pairs",
                "n_pairs": 0,
                "left_explain":  res.get("PRIMER_LEFT_EXPLAIN", ""),
                "right_explain": res.get("PRIMER_RIGHT_EXPLAIN", ""),
                "pair_explain":  res.get("PRIMER_PAIR_EXPLAIN", ""),
            })
    return all_pairs, unfilled


_PAIR_COLS = [
    "pair_id", "marker", "tile_id", "window_id", "pair_idx",
    "f_seq", "r_seq", "f_start", "f_end", "r_start", "r_end",
    "amplicon_len", "f_tm", "r_tm", "f_gc", "r_gc",
    "f_hairpin_dG", "r_hairpin_dG", "pair_dimer_dG", "pair_penalty",
]


def _emit(df_all, unfilled, args, id_prefix):
    Path(args.pairs_out).parent.mkdir(parents=True, exist_ok=True)
    if df_all.empty:
        df_all = pd.DataFrame(columns=_PAIR_COLS)
    else:
        df_all.insert(0, "pair_id",
                      [f"{id_prefix}_{i+1:04d}" for i in range(len(df_all))])
        df_all = df_all[_PAIR_COLS]
    df_all.to_csv(args.pairs_out, sep="\t", index=False)
    print(f"[dbs] wrote {len(df_all)} pairs -> {args.pairs_out}", file=sys.stderr)

    # Filtered: hairpin/dimer/dG + Tm envelope + amplicon <= 4000 (DBS ceiling).
    if df_all.empty:
        df_filt = df_all.copy()
    else:
        df_filt = df_all[
            (df_all["f_hairpin_dG"] >= args.dg_hairpin_min) &
            (df_all["r_hairpin_dG"] >= args.dg_hairpin_min) &
            (df_all["pair_dimer_dG"] >= args.dg_dimer_min) &
            (df_all["f_tm"] >= args.tm_min) & (df_all["f_tm"] <= args.tm_max) &
            (df_all["r_tm"] >= args.tm_min) & (df_all["r_tm"] <= args.tm_max) &
            (df_all["amplicon_len"] <= 4000)
        ].copy()
    df_filt.to_csv(args.filtered_out, sep="\t", index=False)
    print(f"[dbs] wrote {len(df_filt)} filtered -> {args.filtered_out}",
          file=sys.stderr)

    un_cols = ["window_id", "tile_id", "reason", "n_pairs",
               "left_explain", "right_explain", "pair_explain"]
    df_un = pd.DataFrame(unfilled, columns=un_cols) if unfilled \
        else pd.DataFrame(columns=un_cols)
    df_un.to_csv(args.unfilled_out, sep="\t", index=False)
    print(f"[dbs] wrote {len(df_un)} unfilled combos -> {args.unfilled_out}",
          file=sys.stderr)

    if not df_filt.empty:
        tm_spread = (
            max(df_filt["f_tm"].max(), df_filt["r_tm"].max()) -
            min(df_filt["f_tm"].min(), df_filt["r_tm"].min())
        )
        amp_lo, amp_hi = int(df_filt["amplicon_len"].min()), int(df_filt["amplicon_len"].max())
        print(f"[dbs] filtered Tm spread: {tm_spread:.2f} °C  "
              f"amplicon range: {amp_lo}-{amp_hi} bp", file=sys.stderr)


# ---------- subcommand: single -----------------------------------------------

def cmd_single(args):
    records, aln_len = load_alignment(args.aln)
    template, aln_to_tpl, tpl_to_aln = build_consensus_template(records, aln_len)
    print(f"[dbs/single] alignment: {aln_len} cols, {len(records)} seqs; "
          f"template: {len(template)} bp", file=sys.stderr)

    cands = load_candidates(args.candidates)
    blocks_aln = collapse_all_to_blocks(cands)
    blocks_tpl = [project_block_to_template(b, aln_to_tpl) for b in blocks_aln]
    print(f"[dbs/single] merged conserved blocks: {len(blocks_aln)}",
          file=sys.stderr)
    _print_blocks("dbs/single", blocks_aln, blocks_tpl)

    # Log the achievable envelope per block pair before any primer3 call.
    print(f"[dbs/single] achievable amplicon per (fwd,rev) block pair "
          f"(tpl coords, min-max):", file=sys.stderr)
    for i in range(len(blocks_aln)):
        for j in range(i + 1, len(blocks_aln)):
            lo = blocks_tpl[j][0] - blocks_tpl[i][1] + 1
            hi = blocks_tpl[j][1] - blocks_tpl[i][0] + 1
            print(f"    B{i+1}xB{j+1}: {lo}-{hi} bp", file=sys.stderr)

    combos = [(blocks_aln[i], blocks_aln[j], None)
              for i in range(len(blocks_aln))
              for j in range(i + 1, len(blocks_aln))]
    all_pairs, unfilled = _combo_rows(
        template, aln_to_tpl, tpl_to_aln, combos,
        product_min=args.product_min, product_max=args.product_max,
        tm_min=args.tm_min, tm_max=args.tm_max,
        gc_min=args.gc_min, gc_max=args.gc_max,
        len_min=args.len_min, len_max=args.len_max,
        num_return=args.num_return,
        marker="mit",
    )
    _emit(pd.DataFrame(all_pairs), unfilled, args, id_prefix="mit_dbs")


# ---------- subcommand: tiled ------------------------------------------------

def _pick_tile_combos(blocks_aln, core_left, core_right):
    """Assign tile IDs from the actual block geometry.

    tile1 = start-block × any mid-block that yields ~3000-3500 amplicon.
    tile2 = any mid-block × end-block that yields ~3000-3500 amplicon.

    "Start" = the block containing (or nearest to) core_left; "end" = the
    block containing (or nearest to) core_right; "mid" = everything else.
    We emit BOTH tile1 and tile2 candidates and let primer3 filter by
    product-size range.
    """
    n = len(blocks_aln)
    # start block: closest to core_left
    start_idx = min(range(n), key=lambda i: abs(blocks_aln[i][0] - core_left))
    # end block: closest to core_right, must be > start_idx
    end_idx = max(range(n), key=lambda i: -abs(blocks_aln[i][1] - core_right))
    if end_idx <= start_idx:
        end_idx = n - 1
    mid_idxs = [i for i in range(n) if i not in (start_idx, end_idx)]

    combos = []
    # tile1: start × any block strictly downstream (mid or end)
    for j in range(start_idx + 1, end_idx + 1):
        combos.append((blocks_aln[start_idx], blocks_aln[j], "tile1"))
    # tile2: any mid block × end
    for i in mid_idxs:
        if i < end_idx:
            combos.append((blocks_aln[i], blocks_aln[end_idx], "tile2"))
    return combos, start_idx, mid_idxs, end_idx


def cmd_tiled(args):
    records, aln_len = load_alignment(args.aln)
    template, aln_to_tpl, tpl_to_aln = build_consensus_template(records, aln_len)
    print(f"[dbs/tiled] alignment: {aln_len} cols, {len(records)} seqs; "
          f"template: {len(template)} bp", file=sys.stderr)

    core = pd.read_csv(args.core_bounds, sep="\t")
    core_left = int(core.iloc[0]["core_left"])
    core_right = int(core.iloc[0]["core_right"])
    print(f"[dbs/tiled] MIT core: aln {core_left}-{core_right} "
          f"({core_right-core_left+1} bp)", file=sys.stderr)

    cands = load_candidates(args.candidates)
    blocks_aln = collapse_all_to_blocks(cands)
    blocks_tpl = [project_block_to_template(b, aln_to_tpl) for b in blocks_aln]
    print(f"[dbs/tiled] merged conserved blocks: {len(blocks_aln)}",
          file=sys.stderr)
    _print_blocks("dbs/tiled", blocks_aln, blocks_tpl)

    combos, s_i, m_is, e_i = _pick_tile_combos(blocks_aln, core_left, core_right)
    print(f"[dbs/tiled] start block = B{s_i+1}, end block = B{e_i+1}, "
          f"mid blocks = {['B'+str(i+1) for i in m_is]}", file=sys.stderr)
    for fwd_aln, rev_aln, tile_id in combos:
        print(f"    {tile_id}: fwd aln {fwd_aln} × rev aln {rev_aln}",
              file=sys.stderr)

    all_pairs, unfilled = _combo_rows(
        template, aln_to_tpl, tpl_to_aln, combos,
        product_min=args.product_min, product_max=args.product_max,
        tm_min=args.tm_min, tm_max=args.tm_max,
        gc_min=args.gc_min, gc_max=args.gc_max,
        len_min=args.len_min, len_max=args.len_max,
        num_return=args.num_return,
        marker="mit",
    )
    _emit(pd.DataFrame(all_pairs), unfilled, args, id_prefix="mit_dbstile")


# ---------- subcommand: resolution -------------------------------------------

_SPECIES_ALIASES = {
    # The MIT alignment carries both spellings; treat as one species.
    "cynomologi": "cynomolgi",
}

def _species_from_id(rec_id: str) -> str | None:
    m = re.match(r"^P([a-zA-Z]+)_", rec_id)
    if not m:
        return None
    sp = m.group(1).lower()
    return _SPECIES_ALIASES.get(sp, sp)


def _per_species_consensus(records, aln_len: int, min_cov: float = 0.5):
    """Build a per-species per-column consensus base.

    Returns {species: [base_or_None] * (aln_len+1)} using 1-based indexing.
    A column is only assigned a base for species S if:
      - fraction of S's sequences with an ACGT at that column ≥ min_cov, AND
      - one base is strictly the plurality among ACGT.
    Otherwise the column is None (untypable → not diagnostic on that side).
    """
    by_sp: dict[str, list] = defaultdict(list)
    for r in records:
        sp = _species_from_id(r.id)
        if sp is None:
            continue
        by_sp[sp].append(str(r.seq).upper())

    consensus: dict[str, list] = {}
    for sp, seqs in by_sp.items():
        n = len(seqs)
        col_bases: list = [None] * (aln_len + 1)  # 1-based
        for col in range(aln_len):
            bases = [s[col] for s in seqs]
            acgt = [b for b in bases if b in ("A", "C", "G", "T")]
            if len(acgt) < min_cov * n:
                continue
            cnt = Counter(acgt)
            top2 = cnt.most_common(2)
            if len(top2) == 1 or top2[0][1] > top2[1][1]:
                col_bases[col + 1] = top2[0][0]
        consensus[sp] = col_bases
    return consensus


def _diag_count(cons_a, cons_b, aln_lo: int, aln_hi: int) -> int:
    """Count columns in [aln_lo, aln_hi] where both consensuses are ACGT
    and differ."""
    n = 0
    for i in range(aln_lo, aln_hi + 1):
        a = cons_a[i]; b = cons_b[i]
        if a is not None and b is not None and a != b:
            n += 1
    return n


def _amplicon_interval_aln(pairs_tsv: str, tpl_to_aln: list[int]) -> tuple[int, int] | None:
    """Outer envelope of all filtered pairs in an option's TSV, converted
    from template coords to alignment coords (1-based inclusive).
    """
    df = pd.read_csv(pairs_tsv, sep="\t")
    if df.empty:
        return None
    tpl_lo = int(df["f_start"].min())
    tpl_hi = int(df["r_end"].max())
    return (tpl_to_aln[tpl_lo], tpl_to_aln[tpl_hi])


def _tiled_intervals_aln(pairs_tsv: str, tpl_to_aln: list[int]) -> list[tuple[int, int]]:
    """One (aln_lo, aln_hi) per tile_id. Union used for diagnostic count."""
    df = pd.read_csv(pairs_tsv, sep="\t")
    if df.empty or "tile_id" not in df.columns:
        return []
    out = []
    for tile, sub in df.groupby("tile_id"):
        tpl_lo = int(sub["f_start"].min())
        tpl_hi = int(sub["r_end"].max())
        out.append((tpl_to_aln[tpl_lo], tpl_to_aln[tpl_hi]))
    return out


def _diag_count_union(cons_a, cons_b, intervals: list[tuple[int, int]]) -> int:
    """Diagnostic-position count over the union of possibly-overlapping intervals."""
    if not intervals:
        return 0
    # merge intervals
    ivs = sorted(intervals)
    merged = [list(ivs[0])]
    for lo, hi in ivs[1:]:
        if lo <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], hi)
        else:
            merged.append([lo, hi])
    return sum(_diag_count(cons_a, cons_b, lo, hi) for lo, hi in merged)


def cmd_resolution(args):
    records, aln_len = load_alignment(args.aln)
    template, aln_to_tpl, tpl_to_aln = build_consensus_template(records, aln_len)
    print(f"[dbs/resolution] alignment: {aln_len} cols, {len(records)} seqs",
          file=sys.stderr)

    # per-species per-column consensus
    cons = _per_species_consensus(records, aln_len)
    sp_list = sorted(cons.keys())
    print(f"[dbs/resolution] species in MIT alignment: {sp_list}",
          file=sys.stderr)

    # option → alignment interval(s)
    tier_a_iv = _amplicon_interval_aln(args.tier_a_pairs, tpl_to_aln)
    single_iv = _amplicon_interval_aln(args.single_pairs,  tpl_to_aln)
    tiled_ivs = _tiled_intervals_aln( args.tiled_pairs,   tpl_to_aln)
    print(f"[dbs/resolution] tier_a interval (aln): {tier_a_iv}", file=sys.stderr)
    print(f"[dbs/resolution] single interval (aln): {single_iv}", file=sys.stderr)
    print(f"[dbs/resolution] tiled intervals (aln): {tiled_ivs}", file=sys.stderr)

    # unordered species pairs
    rows = []
    for i, a in enumerate(sp_list):
        for b in sp_list[i + 1:]:
            n_tierA  = _diag_count(cons[a], cons[b], *tier_a_iv) if tier_a_iv else 0
            n_single = _diag_count(cons[a], cons[b], *single_iv) if single_iv else 0
            n_tiled  = _diag_count_union(cons[a], cons[b], tiled_ivs)
            rows.append({
                "pair":                     f"{a} vs {b}",
                "tier_a_n_diagnostic":      n_tierA,
                "tier_a_resolvable":        bool(n_tierA >= 1),
                "single_n_diagnostic":      n_single,
                "single_resolvable":        bool(n_single >= 1),
                "tiled_n_diagnostic":       n_tiled,
                "tiled_resolvable":         bool(n_tiled >= 1),
            })
    df = pd.DataFrame(rows)

    # Reconciliation with the v1 resolution table.
    disagreements = []
    if args.resolution_table:
        rt = pd.read_csv(args.resolution_table, sep="\t")
        # rt uses "a vs b" text; our pairs are alphabetic (a < b). Canonicalise both.
        def canon(pair: str):
            parts = [p.strip() for p in pair.split("vs")]
            return " vs ".join(sorted(parts))
        rt["_canon"] = rt["pair"].map(canon)
        df["_canon"] = df["pair"].map(canon)
        merged = df.merge(
            rt[["_canon", "mit_resolved"]], on="_canon", how="left",
        )
        merged["tier_a_reconciled_with_v1"] = merged.apply(
            lambda r: (
                "match" if pd.notna(r["mit_resolved"]) and
                    (bool(r["tier_a_resolvable"]) == bool(r["mit_resolved"]))
                else ("v1_missing_row" if pd.isna(r["mit_resolved"]) else "DISAGREE")
            ),
            axis=1,
        )
        for _, r in merged.iterrows():
            if r["tier_a_reconciled_with_v1"] == "DISAGREE":
                disagreements.append(
                    (r["pair"], bool(r["tier_a_resolvable"]),
                     bool(r["mit_resolved"]),
                     int(r["tier_a_n_diagnostic"]))
                )
        df = merged.drop(columns=["_canon"])
    else:
        df["tier_a_reconciled_with_v1"] = "n/a"

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, sep="\t", index=False)
    print(f"[dbs/resolution] wrote {len(df)} pair rows -> {args.out}",
          file=sys.stderr)

    # summary lines
    def _resolvable(col):
        return int(df[col].astype(bool).sum())
    print(f"[dbs/resolution] resolvable counts (of {len(df)}): "
          f"tier_a={_resolvable('tier_a_resolvable')}, "
          f"single={_resolvable('single_resolvable')}, "
          f"tiled={_resolvable('tiled_resolvable')}", file=sys.stderr)

    lost_single = df[df["tier_a_resolvable"] & ~df["single_resolvable"]]
    lost_tiled  = df[df["tier_a_resolvable"] & ~df["tiled_resolvable"]]
    recovered   = df[df["tier_a_resolvable"] & ~df["single_resolvable"]
                     & df["tiled_resolvable"]]
    print(f"[dbs/resolution] pairs lost by SINGLE vs tier_a: {len(lost_single)}: "
          f"{list(lost_single['pair'])}", file=sys.stderr)
    print(f"[dbs/resolution] pairs lost by TILED  vs tier_a: {len(lost_tiled)}: "
          f"{list(lost_tiled['pair'])}",  file=sys.stderr)
    print(f"[dbs/resolution] pairs single loses but tiled recovers: {len(recovered)}: "
          f"{list(recovered['pair'])}", file=sys.stderr)

    if disagreements:
        print("[dbs/resolution] ⚠ Tier-A reconciliation DISAGREEMENTS "
              "(pair, computed_resolvable, v1_mit_resolved, n_diagnostic):",
              file=sys.stderr)
        for row in disagreements:
            print(f"    {row}", file=sys.stderr)


# ---------- CLI --------------------------------------------------------------

def _add_primer3_args(p):
    p.add_argument("--aln", required=True)
    p.add_argument("--candidates", required=True,
                   help="outputs/entropy/mit.interior_candidates.tsv "
                        "(has no end_region gate — covers whole MIT core).")
    p.add_argument("--pairs-out", required=True)
    p.add_argument("--filtered-out", required=True)
    p.add_argument("--unfilled-out", required=True)
    p.add_argument("--product-min", type=int, required=True)
    p.add_argument("--product-max", type=int, required=True,
                   help="DBS ceiling is 4000 bp; keep this ≤ 4000.")
    p.add_argument("--tm-min", type=float, default=58.0)
    p.add_argument("--tm-max", type=float, default=62.0)
    p.add_argument("--gc-min", type=float, default=40.0)
    p.add_argument("--gc-max", type=float, default=60.0)
    p.add_argument("--len-min", type=int, default=18)
    p.add_argument("--len-max", type=int, default=25)
    p.add_argument("--num-return", type=int, default=10)
    p.add_argument("--dg-hairpin-min", type=float, default=-6.0)
    p.add_argument("--dg-dimer-min",   type=float, default=-8.0)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p_single = sub.add_parser("single", help="Task A: one ≤4 kb MIT amplicon")
    _add_primer3_args(p_single)
    p_single.set_defaults(func=cmd_single)

    p_tiled = sub.add_parser("tiled", help="Task B: two tiled ~3-3.5 kb amplicons")
    _add_primer3_args(p_tiled)
    p_tiled.add_argument("--core-bounds", required=True,
                         help="outputs/entropy/mit.core_bounds.tsv")
    p_tiled.set_defaults(func=cmd_tiled)

    p_res = sub.add_parser("resolution", help="Task C: species-pair resolution trade-off")
    p_res.add_argument("--aln", required=True)
    p_res.add_argument("--tier-a-pairs", required=True)
    p_res.add_argument("--single-pairs", required=True)
    p_res.add_argument("--tiled-pairs",  required=True)
    p_res.add_argument("--resolution-table", default=None,
                       help="outputs/cross_species/resolution_table.tsv "
                            "(for Tier-A reconciliation).")
    p_res.add_argument("--out", required=True)
    p_res.set_defaults(func=cmd_resolution)

    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
