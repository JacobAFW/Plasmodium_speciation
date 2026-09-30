#!/usr/bin/env python3
"""
offtarget_summary.py — roll the per-genome in-silico PCR evidence up to
HANDOFF Step 5's off-target grain: one row per
(pair_id, off_target_species, subject).

WHY THIS IS NOT `insilico_pcr.py --summary-out`
-----------------------------------------------
That summary is a full (pair x subject SEQUENCE) grid with zero rows kept,
which is exactly right for the target half ("did every target amplify?") and
exactly wrong here: a staged genome is thousands of contigs, so the grid would
be 180 pairs x 8090 A. gambiae contigs = 1.5M rows of "no, this primer pair
does not amplify contig 7,412 either". The off-target question is asked per
GENOME, and zero is the expected answer, so the zero-fill happens once per
(pair, genome) instead.

The off-target species label also comes from `MANIFEST.tsv`, not from the
panel-species ladder in `insilico_pcr.py` — the ladder exists to map a
reference header onto one of the 11 PANEL species and by construction returns
None for a host, a vector or a non-panel Plasmodium. Any subject sequence that
DOES hit the ladder is reported as a warning, because it would mean an
off-target genome carrying a header that reads like a panel species.

Usage:
  offtarget_summary.py --amplicons outputs/specificity/offtarget/*.amplicons.tsv \\
      --pairs outputs/primer_design/mit_pairs_filtered.tsv [...] \\
      --manifest data/staged/MANIFEST.tsv \\
      --efficiency-threshold 0.5 \\
      --out outputs/specificity/offtarget_amplification.tsv
"""
import argparse
import sys
from pathlib import Path

import pandas as pd


def slug(filename: str) -> str:
    """`plasmodium_yoelii.fna` -> `plasmodium_yoelii` (the subject key)."""
    return Path(filename).name.rsplit(".", 1)[0]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--amplicons", nargs="+", required=True,
                    help="per-genome amplicon TSVs from insilico_pcr.py; the "
                         "file stem is the subject key")
    ap.add_argument("--pairs", nargs="+", required=True,
                    help="the same primer-pair TSVs the screen was run with")
    ap.add_argument("--manifest", required=True, help="data/staged/MANIFEST.tsv")
    ap.add_argument("--efficiency-threshold", type=float, default=0.5)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    pairs = pd.concat(
        [pd.read_csv(p, sep="\t", dtype={"pair_id": str}).assign(
            pair_source=Path(p).name) for p in args.pairs],
        ignore_index=True)[["pair_id", "pair_source"]]

    man = pd.read_csv(args.manifest, sep="\t")
    man["subject"] = man["expected_filename"].map(slug)
    meta = man.set_index("subject")[["genus_species", "purpose"]]

    subjects = [slug(a).replace(".amplicons", "") for a in args.amplicons]
    unknown = sorted(set(subjects) - set(meta.index))
    if unknown:
        sys.exit(f"[offtarget_summary] ERROR: no MANIFEST.tsv row for {unknown}")

    frames = []
    for path, subj in zip(args.amplicons, subjects):
        df = pd.read_csv(path, sep="\t", dtype={"pair_id": str})
        df["subject"] = subj
        frames.append(df)
    amps = pd.concat(frames, ignore_index=True)

    # The panel ladder should never fire on an off-target genome. If it does,
    # say so loudly rather than silently carrying a bogus panel-species label.
    stray = amps.loc[amps["species"].notna(), ["subject", "reference", "species"]]
    if len(stray):
        print("[offtarget_summary] WARNING: panel-species ladder matched "
              f"{len(stray)} off-target amplicon rows (header looks like a panel "
              f"species):\n{stray.drop_duplicates().to_string(index=False)}",
              flush=True)

    grid = pd.MultiIndex.from_product(
        [pairs["pair_id"].tolist(), subjects],
        names=["pair_id", "subject"]).to_frame(index=False)

    if len(amps):
        agg = (amps.groupby(["pair_id", "subject"], as_index=False)
               .agg(amplicon_count=("amplicon_len", "size"),
                    amplicon_lengths=("amplicon_len",
                                      lambda s: ",".join(map(str, sorted(s)))),
                    max_efficiency=("efficiency", "max"),
                    n_contigs_hit=("reference", "nunique")))
    else:
        agg = pd.DataFrame(columns=["pair_id", "subject", "amplicon_count",
                                    "amplicon_lengths", "max_efficiency",
                                    "n_contigs_hit"])

    out = grid.merge(agg, on=["pair_id", "subject"], how="left")
    for col, fill in (("amplicon_count", 0), ("n_contigs_hit", 0)):
        out[col] = out[col].fillna(fill).astype(int)
    out["amplicon_lengths"] = out["amplicon_lengths"].fillna("")
    out["max_efficiency"] = out["max_efficiency"].fillna(0.0).round(4)

    out = out.merge(pairs, on="pair_id", how="left")
    out["off_target_species"] = out["subject"].map(meta["genus_species"])
    out["purpose"] = out["subject"].map(meta["purpose"])
    # HANDOFF Step 5's flag, verbatim: a product exists AND it would actually
    # prime at the annealing temperature.
    out["flag"] = ((out["amplicon_count"] > 0) &
                   (out["max_efficiency"] > args.efficiency_threshold))

    out = out[["pair_id", "pair_source", "off_target_species", "purpose",
               "subject", "amplicon_count", "n_contigs_hit",
               "amplicon_lengths", "max_efficiency", "flag"]]
    out = out.sort_values(["pair_id", "subject"])
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, sep="\t", index=False)

    print(f"[offtarget_summary] {len(pairs)} pairs x {len(subjects)} genomes "
          f"= {len(out)} rows -> {args.out}", flush=True)
    print(f"[offtarget_summary] {int((out.amplicon_count > 0).sum())} cells with "
          f"any product; {int(out.flag.sum())} flagged "
          f"(> {args.efficiency_threshold}); "
          f"{out.loc[out.flag, 'pair_id'].nunique()} distinct pairs flagged",
          flush=True)


if __name__ == "__main__":
    main()
