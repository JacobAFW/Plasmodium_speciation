#!/usr/bin/env python3
"""
panel_verdict.py — combine Step 5's target and off-target halves into one
per-pair verdict: `outputs/specificity/panel_pass_fail.tsv`.

OFF-TARGET MEANS NON-PLASMODIUM (corrected 2026-08-28)
-----------------------------------------------------
This panel is **genus-conserved by design**: the primers are meant to amplify
any *Plasmodium*, and speciation happens downstream by reading the amplicon,
not by which primer fired. So a pair that amplifies *P. berghei* or
*P. gallinaceum* is doing exactly what it was designed to do. Scoring that as a
specificity failure — which this script did until 2026-08-28 — inverts the
design intent and rejected 34 pairs that are in fact clean.

The failure that matters is priming on something that is **not Plasmodium** and
**is** in the tube: host (human, macaque), vector (Anopheles) or contaminant
(E. coli, M. tuberculosis). That is the pass/fail axis now. Congener
amplification is still measured and reported in full — `n_congener_flagged`,
`congener_flagged_species` — it just annotates rather than fails.

Which manifest `purpose` counts as a congener is `config.specificity.
congener_purpose`, not a literal, so a manifest that grows new categories does
not silently change a verdict.

VERDICTS
--------
  pass            amplifies every target species AND no NON-Plasmodium
                  off-target above the efficiency threshold
  fail_no_target  misses >= 1 target species
  fail_offtarget  >= 1 NON-Plasmodium genome (host / vector / contaminant) with
                  amplicon_count > 0 AND max_efficiency > threshold
  borderline      a NON-Plasmodium product exists but stays below the
                  threshold — review, do not auto-reject

Precedence is fail_no_target > fail_offtarget > borderline > pass, and it is
only a display order: `target_ok`, `n_offtarget_flagged`,
`n_offtarget_flagged_nonplasmodium`, `n_congener_flagged` and `n_offtarget_any`
are all kept as their own columns, so nothing a verdict outranks is hidden.

WHAT COUNTS AS "AMPLIFIES A TARGET"
-----------------------------------
Two readings exist and BOTH are reported, because they disagree by a lot:

  strict  (`target_ok`, drives the verdict) — >= 1 on-size amplicon AND
          max_efficiency > threshold, for every target species. This is the
          definition Session 14's target-half table was built on (Tier A
          26/40, DBS single 8/60, DBS tiled 0/60, 18S 0/20); using anything
          else here would silently restate that result.
  lenient (`target_ok_amplicon_only`) — HANDOFF's literal wording, >= 1
          on-size amplicon per target species, efficiency ignored.

The gap between them IS the finding, not an accounting detail: a pair that
puts an on-size product on every species but does it with a duplex that melts
below the annealing temperature has not amplified anything in a real tube.
`11_specificity.README.md` calls this the interesting failure mode, so the
species involved are carried as `n_species_low_eff` / `low_eff_species`.

The expected species set is per MARKER, taken from what that pair was actually
tested against, not from `config["species_targets"]`. The 18S reference DB has
no simiovale, so scoring 18S pairs out of 11 would fail every one of them for a
gap in the reference set rather than a defect in the primer.

Specificity is not resolution: a `pass` here says the pair amplifies the right
things, never that it can tell them apart. Read against
`outputs/cross_species/resolution_table.tsv`.

Usage:
  panel_verdict.py --target-coverage outputs/specificity/target_species_coverage.tsv \\
      --offtarget outputs/specificity/offtarget_amplification.tsv \\
      --manifest data/staged/MANIFEST.tsv \\
      --efficiency-threshold 0.5 \\
      --out outputs/specificity/panel_pass_fail.tsv
"""
import argparse
from pathlib import Path

import pandas as pd


def tier_of(pair_source: str) -> str:
    """`mit_dbs_single_pairs_filtered.tsv` -> `mit_dbs_single`."""
    return Path(pair_source).name.split(".")[0].replace("_pairs_filtered", "")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--target-coverage", required=True)
    ap.add_argument("--offtarget", required=True)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--efficiency-threshold", type=float, default=0.5)
    ap.add_argument("--congener-purpose", default="non_target_plasmodium",
                    help="manifest `purpose` that marks a non-panel Plasmodium; "
                         "these are annotated, never failed")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cov = pd.read_csv(args.target_coverage, sep="\t", dtype={"pair_id": str})
    off = pd.read_csv(args.offtarget, sep="\t", dtype={"pair_id": str})
    man = pd.read_csv(args.manifest, sep="\t")

    hosts = set(man.loc[man["purpose"] == "host", "genus_species"])

    cov["tier"] = cov["pair_source"].map(tier_of)
    cov["amplified"] = cov["n_amplified"] > 0
    cov["low_eff"] = cov["max_efficiency"] <= args.efficiency_threshold
    cov["amplified_strict"] = cov["amplified"] & ~cov["low_eff"]

    tgt = (cov.groupby(["pair_id", "pair_source", "tier"], as_index=False)
           .agg(n_target_species=("species", "nunique"),
                n_species_amplified=("amplified", "sum"),
                n_species_amplified_strict=("amplified_strict", "sum"),
                n_species_low_eff=("low_eff", "sum")))
    miss = (cov[~cov["amplified_strict"]].groupby("pair_id")["species"]
            .apply(lambda s: ",".join(sorted(s))).rename("missing_species"))
    lowe = (cov[cov["low_eff"]].groupby("pair_id")["species"]
            .apply(lambda s: ",".join(sorted(s))).rename("low_eff_species"))
    tgt = tgt.merge(miss, on="pair_id", how="left").merge(lowe, on="pair_id",
                                                          how="left")
    tgt["target_ok"] = (tgt["n_species_amplified_strict"]
                        == tgt["n_target_species"])
    tgt["target_ok_amplicon_only"] = (tgt["n_species_amplified"]
                                      == tgt["n_target_species"])

    hit = off[off["amplicon_count"] > 0]
    flagged = off[off["flag"]]
    ot = (off.groupby("pair_id", as_index=False)
          .agg(n_offtarget_tested=("subject", "nunique")))
    ot["n_offtarget_any"] = ot["pair_id"].map(
        hit.groupby("pair_id")["subject"].nunique()).fillna(0).astype(int)
    ot["n_offtarget_flagged"] = ot["pair_id"].map(
        flagged.groupby("pair_id")["subject"].nunique()).fillna(0).astype(int)
    ot["max_offtarget_efficiency"] = ot["pair_id"].map(
        off.groupby("pair_id")["max_efficiency"].max()).fillna(0.0).round(4)
    for name, src in (("offtarget_flagged_species", flagged),
                      ("offtarget_any_species", hit)):
        ot[name] = ot["pair_id"].map(
            src.groupby("pair_id")["off_target_species"]
            .apply(lambda s: ",".join(sorted(set(s))))).fillna("")

    # Split the flagged off-targets by manifest `purpose`. MIT and 18S are
    # genus-conserved BY DESIGN, so cross-amplifying a non-panel *Plasmodium*
    # is the intended behaviour, not a defect — the amplicon is what speciates,
    # not the primer. Only the NON-Plasmodium half is a specificity failure.
    ot["offtarget_flagged_purposes"] = ot["pair_id"].map(
        flagged.groupby("pair_id")["purpose"]
        .apply(lambda s: ",".join(sorted(set(s))))).fillna("")

    cong = flagged[flagged["purpose"] == args.congener_purpose]
    nonpl = flagged[flagged["purpose"] != args.congener_purpose]
    nonpl_any = hit[hit["purpose"] != args.congener_purpose]

    ot["n_offtarget_flagged_nonplasmodium"] = ot["pair_id"].map(
        nonpl.groupby("pair_id")["subject"].nunique()).fillna(0).astype(int)
    ot["n_offtarget_any_nonplasmodium"] = ot["pair_id"].map(
        nonpl_any.groupby("pair_id")["subject"].nunique()).fillna(0).astype(int)
    ot["offtarget_flagged_nonplasmodium_species"] = ot["pair_id"].map(
        nonpl.groupby("pair_id")["off_target_species"]
        .apply(lambda s: ",".join(sorted(set(s))))).fillna("")

    # Congeners: measured in full, reported in full, never a fail.
    ot["n_congener_flagged"] = ot["pair_id"].map(
        cong.groupby("pair_id")["subject"].nunique()).fillna(0).astype(int)
    ot["congener_flagged_species"] = ot["pair_id"].map(
        cong.groupby("pair_id")["off_target_species"]
        .apply(lambda s: ",".join(sorted(set(s))))).fillna("")
    ot["congener_amplification"] = ot["n_congener_flagged"] > 0

    # The headline specificity failure: host DNA is what a real sample is made
    # of, so a host product outranks every other off-target as a problem.
    host_hit = hit[hit["off_target_species"].isin(hosts)]
    host_flag = flagged[flagged["off_target_species"].isin(hosts)]
    ot["host_any_species"] = ot["pair_id"].map(
        host_hit.groupby("pair_id")["off_target_species"]
        .apply(lambda s: ",".join(sorted(set(s))))).fillna("")
    ot["host_flagged_species"] = ot["pair_id"].map(
        host_flag.groupby("pair_id")["off_target_species"]
        .apply(lambda s: ",".join(sorted(set(s))))).fillna("")
    ot["host_cross_amplification"] = ot["host_flagged_species"] != ""

    v = tgt.merge(ot, on="pair_id", how="left")
    missing_off = v["n_offtarget_tested"].isna().sum()
    if missing_off:
        raise SystemExit(f"[panel_verdict] ERROR: {missing_off} pairs have no "
                         "off-target rows — the screen did not cover every pair")

    # Only NON-Plasmodium products move the verdict. Congener amplification is
    # the design intent and is carried in its own columns.
    v["verdict"] = "pass"
    v.loc[v["n_offtarget_any_nonplasmodium"] > 0, "verdict"] = "borderline"
    v.loc[v["n_offtarget_flagged_nonplasmodium"] > 0, "verdict"] = "fail_offtarget"
    v.loc[~v["target_ok"], "verdict"] = "fail_no_target"

    v["missing_species"] = v["missing_species"].fillna("")
    v["low_eff_species"] = v["low_eff_species"].fillna("")
    v = v[["pair_id", "pair_source", "tier", "verdict", "target_ok",
           "target_ok_amplicon_only", "n_target_species",
           "n_species_amplified_strict", "n_species_amplified",
           "missing_species", "n_species_low_eff", "low_eff_species",
           "n_offtarget_tested", "n_offtarget_any", "n_offtarget_flagged",
           "n_offtarget_any_nonplasmodium",
           "n_offtarget_flagged_nonplasmodium",
           "offtarget_flagged_nonplasmodium_species",
           "congener_amplification", "n_congener_flagged",
           "congener_flagged_species", "max_offtarget_efficiency",
           "offtarget_any_species", "offtarget_flagged_species",
           "offtarget_flagged_purposes", "host_cross_amplification",
           "host_any_species", "host_flagged_species"]].sort_values("pair_id")

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    v.to_csv(args.out, sep="\t", index=False)

    print(f"[panel_verdict] {len(v)} pairs -> {args.out}", flush=True)
    print(v.groupby(["tier", "verdict"]).size().to_string(), flush=True)
    print("[panel_verdict] target_ok strict vs amplicon-only per tier:", flush=True)
    print(v.groupby("tier")[["target_ok", "target_ok_amplicon_only"]].sum()
          .to_string(), flush=True)
    print(f"[panel_verdict] host cross-amplification: "
          f"{int(v['host_cross_amplification'].sum())} pairs", flush=True)
    print(f"[panel_verdict] congener (expected, not a fail) amplification: "
          f"{int(v['congener_amplification'].sum())} pairs; "
          f"non-Plasmodium flagged: "
          f"{int((v['n_offtarget_flagged_nonplasmodium'] > 0).sum())} pairs",
          flush=True)


if __name__ == "__main__":
    main()
