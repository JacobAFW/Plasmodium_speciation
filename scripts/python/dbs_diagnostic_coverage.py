#!/usr/bin/env python3
"""
dbs_diagnostic_coverage.py — for every candidate primer pair, does the amplicon
still carry the fixed sites that separate *P. fieldi* from *P. simiovale*, and
does it still separate them in practice?

THE CONSTRAINT THIS TESTS, AND WHY IT CHANGED
---------------------------------------------
The original simiovale/fieldi investigation concluded that "all discriminating
power for this pair is in `AB434920` 4780-5928", i.e. in **cytb**, and that
became the constraint on the DBS <=4 kb tier: a short amplicon had to reach cytb
or the pair was lost.

That premise was an artefact. The cox1 locus set had been built by pooling
cytb-only fragments with cox1-only deposits, leaving 8 mutually covered columns
out of 6676 — a scan that could only ever return zero sites. Rebuilt properly,
cox1 carries **7** fixed sites (3378-4394) on *better* sampling than cytb's
**5** (4881-5835). See `simiovale_fieldi_report_corrected.md` sections 3b/6.

So the real constraint is **cox1 OR cytb**, which is materially looser: a short
amplicon trimmed from the 3' end that loses cytb but retains cox1 still
separates the pair. This script measures that, per pair.

TWO READ-OUTS, DELIBERATELY
---------------------------
  site coverage  how many of the 7 cox1 / 5 cytb fixed sites fall inside the
                 amplicon. Cheap, and it is the quantity the constraint is
                 phrased in.
  pair margin    the actual best-hit separation of fieldi from simiovale when
                 only this pair's amplicon region is visible: the smallest
                 own-species-minus-other-species %id margin over every
                 reference of both species. This is the decision-relevant
                 number, because a handful of fixed sites scattered in ~3.4 kb
                 does not automatically outvote thousands of shared positions.

They can disagree, and where they do the margin is what a real read-out does.

Coordinates come from the anchor reference named in
`config["diagnostic_windows"]["anchor"]`, projected anchor -> alignment column
-> consensus template via the same map the primer design used.

Usage:
  dbs_diagnostic_coverage.py --config workflow/config.yaml \
      --alignment outputs/alignment/ma_mit.target.fasta \
      --pairs MIT_DBS_single=outputs/primer_design/dbs/mit_dbs_single_pairs_filtered.tsv \
      --out outputs/primer_design/dbs/dbs_diagnostic_coverage.tsv
"""
import argparse
import sys
import tempfile
from pathlib import Path

import pandas as pd
import yaml
from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_primer3 import build_consensus_template            # noqa: E402
from resolution_margin import (blast_all_vs_all, best_hits,  # noqa: E402
                               species_from_id)


def anchor_to_aln(records, anchor_id):
    """1-based anchor position -> 1-based alignment column."""
    rec = next((r for r in records if r.id == anchor_id), None)
    if rec is None:
        raise SystemExit(f"[dbs_diagnostic] anchor {anchor_id!r} not in alignment")
    out, n = {}, 0
    for col, ch in enumerate(str(rec.seq), start=1):
        if ch not in "-.":
            n += 1
            out[n] = col
    return out


def fieldi_simiovale_margin(records, labels, interval, a, b):
    """Smallest own-minus-other %id margin across both species' references."""
    seqs = {}
    for rec in records:
        s = str(rec.seq)[interval[0] - 1:interval[1]]
        seqs[rec.id] = (s.replace("-", "").replace(".", "").upper(), True)
    with tempfile.TemporaryDirectory() as td:
        hits = best_hits(blast_all_vs_all(seqs, Path(td)))
    if not len(hits):
        return None, 0, 0
    hits = hits.assign(sp=hits["sseqid"].map(labels))
    margins, n_a, n_b = [], 0, 0
    for focal, other in ((a, b), (b, a)):
        for q in [k for k, v in labels.items() if v == focal]:
            h = hits[hits["qseqid"] == q]
            own, oth = h[h["sp"] == focal]["pident"], h[h["sp"] == other]["pident"]
            if own.empty or oth.empty:
                continue
            margins.append(float(own.max()) - float(oth.max()))
            n_a, n_b = (n_a + 1, n_b) if focal == a else (n_a, n_b + 1)
    return (min(margins) if margins else None), n_a, n_b


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True)
    ap.add_argument("--alignment", required=True)
    ap.add_argument("--pairs", action="append", required=True, metavar="TIER=TSV")
    ap.add_argument("--out", required=True)
    ap.add_argument("--no-margin", action="store_true",
                    help="site coverage only; skip the per-pair BLAST margin")
    args = ap.parse_args()

    cfg = yaml.safe_load(Path(args.config).read_text())["diagnostic_windows"]
    sp_a, sp_b = cfg["pair"].split("_vs_")

    records = list(SeqIO.parse(args.alignment, "fasta"))
    aln_len = len(records[0].seq)
    _, _, tpl_to_aln = build_consensus_template(records, aln_len)
    a2a = anchor_to_aln(records, cfg["anchor"])
    labels = {r.id: species_from_id(r.id) for r in records}

    # The unverified ATCC 30157 record is held out: it carries the simiovale
    # allele at 12/12 sites, so leaving it labelled as anything would corrupt
    # exactly the margin this script measures.
    keep = [r for r in records if labels.get(r.id) in (sp_a, sp_b)]
    print(f"[dbs_diagnostic] {sp_a} n={sum(1 for r in keep if labels[r.id]==sp_a)}, "
          f"{sp_b} n={sum(1 for r in keep if labels[r.id]==sp_b)}", flush=True)

    win = cfg["windows"]
    rows = []
    for spec in args.pairs:
        tier, tsv = spec.split("=", 1)
        df = pd.read_csv(tsv, sep="\t")
        for _, r in df.iterrows():
            lo, hi = tpl_to_aln[int(r["f_start"])], tpl_to_aln[int(r["r_end"])]
            row = {"tier": tier, "pair_id": r["pair_id"],
                   "amplicon_len_bp": int(r["amplicon_len"]),
                   "aln_start": lo, "aln_end": hi,
                   "anchor_start": min((k for k, v in a2a.items() if v >= lo),
                                       default=None),
                   "anchor_end": max((k for k, v in a2a.items() if v <= hi),
                                     default=None)}
            got = {}
            for name, w in win.items():
                cols = [a2a[s] for s in w["sites"] if s in a2a]
                got[name] = sum(1 for c in cols if lo <= c <= hi)
                row[f"n_{name}_sites"] = got[name]
                row[f"n_{name}_sites_total"] = len(w["sites"])
                ws, we = a2a.get(w["start"]), a2a.get(w["end"])
                row[f"{name}_window_fully_covered"] = bool(
                    ws and we and lo <= ws and hi >= we)
            full = [n for n, w in win.items() if got[n] == len(w["sites"])]
            any_ = [n for n in win if got[n] > 0]
            row["n_diagnostic_sites"] = sum(got.values())
            row["windows_complete"] = "+".join(sorted(full)) or "none"
            row["windows_partial"] = "+".join(sorted(set(any_) - set(full))) or "none"
            row["site_verdict"] = ("complete_window" if full else
                                   "partial_only" if any_ else "no_sites")
            rows.append(row)

    if not args.no_margin:
        cache = {}
        for row in rows:
            key = (row["aln_start"], row["aln_end"])
            if key not in cache:
                cache[key] = fieldi_simiovale_margin(keep, labels, key, sp_a, sp_b)
            m, n_a, n_b = cache[key]
            row["pair_margin_pct"] = None if m is None else round(m, 4)
            row["n_margins_a"], row["n_margins_b"] = n_a, n_b
            row["separates_pair"] = (
                "untestable_sole_reference" if m is None else
                "no" if m <= 0 else
                "yes" if (n_a and n_b) else "yes_one_side_only")
        print(f"[dbs_diagnostic] {len(cache)} distinct amplicon intervals BLASTed",
              flush=True)

    out = pd.DataFrame(rows)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, sep="\t", index=False)
    print(f"[dbs_diagnostic] {len(out)} pairs -> {args.out}\n", flush=True)

    cols = ["site_verdict", "windows_complete"] + \
           ([] if args.no_margin else ["separates_pair"])
    print(out.groupby(["tier"] + cols).size().to_string(), flush=True)


if __name__ == "__main__":
    main()
