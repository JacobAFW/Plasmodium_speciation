#!/usr/bin/env python3
"""
multiplex_compat.py — Phase 5 Step 3 multiplex-compatibility table.

Takes every filtered primer pair across markers and emits a pairwise
compatibility table over all (pair_A, pair_B) combinations. Metrics:

- `tm_spread`     — |min(Tm)_A − min(Tm)_B| ; max over the four primers of
                    max(Tm) − min(Tm). Multiplex rule of thumb: ≤ 4 °C.
- `cross_dimer_dG` — worst (most negative) heterodimer ΔG among the 4
                    cross-pair oligo pairings (A_f × B_f, A_f × B_r,
                    A_r × B_f, A_r × B_r). Kcal/mol.
- `amplicon_len_diff` — |amplicon_len_A − amplicon_len_B|. For length-based
                    gel discrimination in a multiplex read-out.

Diagonal (pair_A == pair_B) is emitted for completeness with self-metrics
so `every pair × pair cell` in the HANDOFF acceptance check is present.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd
import primer3


_calc_heterodimer = getattr(primer3, "calc_heterodimer", None) or primer3.calcHeterodimer


def cross_dimer_dg(f1, r1, f2, r2) -> float:
    """Worst (most negative) heterodimer ΔG (kcal/mol) among the 4 cross
    pairings between two primer pairs."""
    dgs = []
    for a, b in ((f1, f2), (f1, r2), (r1, f2), (r1, r2)):
        try:
            dgs.append(_calc_heterodimer(a, b).dg / 1000.0)
        except Exception:
            dgs.append(float("nan"))
    return float(min(dgs))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs", nargs="+", required=True,
                    help="One or more marker _pairs_filtered.tsv files.")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    frames = []
    for p in args.pairs:
        df = pd.read_csv(p, sep="\t")
        if df.empty:
            print(f"[multiplex] warning: {p} is empty", file=sys.stderr)
        frames.append(df)
    all_pairs = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()

    if all_pairs.empty:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(columns=["pair_a", "pair_b", "marker_a", "marker_b",
                              "tm_spread", "cross_dimer_dG",
                              "amplicon_len_diff"]).to_csv(args.out, sep="\t", index=False)
        print(f"[multiplex] wrote empty table -> {args.out}", file=sys.stderr)
        return

    rows = []
    for _, a in all_pairs.iterrows():
        for _, b in all_pairs.iterrows():
            tms = [a["f_tm"], a["r_tm"], b["f_tm"], b["r_tm"]]
            spread = float(max(tms) - min(tms))
            xd = cross_dimer_dg(a["f_seq"], a["r_seq"], b["f_seq"], b["r_seq"])
            rows.append({
                "pair_a": a["pair_id"], "pair_b": b["pair_id"],
                "marker_a": a["marker"], "marker_b": b["marker"],
                "tm_spread": round(spread, 3),
                "cross_dimer_dG": round(xd, 3),
                "amplicon_len_diff": int(abs(a["amplicon_len"] - b["amplicon_len"])),
            })
    out = pd.DataFrame(rows)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, sep="\t", index=False)
    print(f"[multiplex] wrote {len(out)} rows -> {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
