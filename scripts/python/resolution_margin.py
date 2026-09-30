#!/usr/bin/env python3
"""
resolution_margin.py — per-pair species resolution measured as an ASSIGNMENT
MARGIN, not as a diagnostic-site threshold.

WHY THIS EXISTS
---------------
The project's older resolution criterion is "does this species pair carry >= 1
fixed diagnostic site inside the amplicon". On this panel that criterion returns
55/55 under every amplicon option, so it cannot discriminate between options and
cannot fail a pair — it is not decision-relevant. It is also optimistic in the
one place it matters: `fieldi vs simiovale` has 5-7 diagnostic sites and was
scored resolvable, while a real best-hit read-out called both *P. fieldi*
references *P. simiovale*.

This script measures what the assay actually does: extract each reference's
amplicon region, BLAST it against every other reference, and take the species of
the best hit. Two numbers come out, and neither is a threshold:

  loo_margin   (per reference) best own-species %id  -  best other-species %id,
               over the whole panel. Negative => this isolate would be called as
               some other species.
  pair_margin  (per species pair) the SMALLEST own-minus-other margin over every
               reference of both species, considering only those two species.
               Negative or zero => at least one reference is pulled across that
               boundary, so no best-hit rule can separate the pair.

A species with one reference has no conspecific to match, so its margin is
undefined rather than zero. That is reported as `untestable_sole_reference`, NOT
as a failure — the distinction is load-bearing after the 2026-08-28 reference
curation, which left *P. fieldi* with a single MIT reference.

REFERENCE CURATION THIS ASSUMES
-------------------------------
Run against the corrected reference set (see the reference-curation note in
`workflow/rules/03_alignment.README.md`):

  * `Pfieldi_AB444133.1_strain_ATCC_30157` relabelled `Punverified_ATCC-30157_AB444133`.
    It carries the *P. simiovale* allele at 12/12 fixed mitochondrial diagnostic
    sites, so it is not scored as any panel species. By default it is also held
    OUT of the BLAST database (`--unverified-in-db exclude`), matching the
    "cleaned DB" of the corrected simiovale/fieldi analysis: a record of
    unresolved identity cannot be a legitimate species call. `--unverified-in-db
    include` reproduces the as-shipped behaviour for comparison.
  * `Psimiovale_AY800109_1` dropped from the panel set as an exact duplicate of
    `Psimiovale_AB434920_1` (0 mismatches / 5809 mutually covered columns), so
    one simiovale haplotype is not double-weighted in a best-hit margin.

COORDINATES
-----------
Amplicon intervals come from the filtered primer tables, whose coordinates are
in the ungapped majority-consensus template. The template and its bijective
alignment<->template map are rebuilt with `run_primer3.build_consensus_template`,
the same construction the primer design used, so the projection is drift-free.
Dropping the duplicate reference does not move a single column: `tpl_to_aln` is
identical with and without it (verified 2026-08-28), which is why the existing
primer coordinates remain valid against the curated alignment.

Usage:
  resolution_margin.py \
      --alignment outputs/alignment/ma_mit.target.fasta \
      --marker mit \
      --pairs MIT_TierA=outputs/primer_design/mit_pairs_filtered.tsv \
      --pairs MIT_DBS_single=outputs/primer_design/dbs/mit_dbs_single_pairs_filtered.tsv \
      --tiled MIT_DBS_tiled=outputs/primer_design/dbs/mit_dbs_tiled_pairs_filtered.tsv \
      --species-targets falciparum,vivax,... \
      --loo-out outputs/cross_species/loo_assignment_mit.tsv \
      --pair-out outputs/cross_species/resolution_table_margin_mit.tsv
"""
import argparse
import itertools
import re
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd
from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_primer3 import build_consensus_template  # noqa: E402

# Species label from a self-describing MIT id (`Psimiovale_AB434920_1`). This is
# the prefix-strip route the R scripts use, not the `case_when` ladder — the two
# agree on every id in the panel, and the prefix is unambiguous by construction.
_ALIASES = {"cynomologi": "cynomolgi"}
UNVERIFIED = "unverified"


def species_from_id(rec_id: str):
    m = re.match(r"^P([a-zA-Z]+)_", rec_id)
    if not m:
        return None
    sp = m.group(1).lower()
    return _ALIASES.get(sp, sp)


# ---------------------------------------------------------------------------
# Amplicon intervals
# ---------------------------------------------------------------------------

def representative_interval(pairs_tsv: Path, tpl_to_aln, tiled: bool):
    """Alignment interval(s) of the pair you would actually order.

    Primary = lowest `pair_penalty`. For a tiled option that is the best pair
    per tile, and the amplicon region is the union of the tiles.
    """
    df = pd.read_csv(pairs_tsv, sep="\t")
    groups = ([g for _, g in df.groupby("tile_id")] if tiled and df["tile_id"].notna().any()
              else [df])
    ivs, chosen = [], []
    for g in groups:
        r = g.loc[g["pair_penalty"].idxmin()]
        ivs.append((tpl_to_aln[int(r["f_start"])], tpl_to_aln[int(r["r_end"])]))
        chosen.append(str(r["pair_id"]))
    # Tiles overlap by design. Merge before extraction, or the overlap is
    # concatenated twice and that shared stretch gets double the weight in
    # every identity the BLAST comparison computes. The UNMERGED endpoints are
    # still what the spanning test probes: a tiled assay needs every tile's
    # own primer footprints present, not just the outer two.
    merged = []
    for a, b in sorted(ivs):
        if merged and a <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(merged[-1][1], b))
        else:
            merged.append((a, b))
    return merged, sorted(ivs), chosen


def extract(records, intervals, probe_intervals, flank=30):
    """Degapped amplicon-region sequence per record, plus a spans flag.

    `intervals` are merged (what gets extracted); `probe_intervals` are the
    per-tile originals whose endpoints the spanning test probes.

    `spans` mirrors the speciation-power audit: the record must carry real
    sequence at EVERY primer footprint (+/- `flank` alignment columns), i.e. a
    real PCR product would exist. Non-spanning records are kept in the database
    (a real read can still hit them) but their own call is reported as such.
    """
    out = {}
    for rec in records:
        s = str(rec.seq)
        seq = "".join(s[a - 1:b] for a, b in intervals)
        seq = seq.replace("-", "").replace(".", "").upper()
        ok = True
        for a, b in probe_intervals:
            for pos in (a, b):
                lo, hi = max(0, pos - 1 - flank), min(len(s), pos + flank)
                if not any(c.upper() in "ACGT" for c in s[lo:hi]):
                    ok = False
        out[rec.id] = (seq, ok)
    return out


# ---------------------------------------------------------------------------
# Leave-one-out BLAST
# ---------------------------------------------------------------------------

def blast_all_vs_all(seqs, tmp: Path):
    fa = tmp / "amplicons.fasta"  # HARDCODE(scratch file inside a TemporaryDirectory, not a reference input; 2026-08-28)
    fa.write_text("".join(f">{k}\n{v}\n" for k, (v, _) in seqs.items() if v))
    subprocess.run(["makeblastdb", "-in", str(fa), "-dbtype", "nucl",
                    "-out", str(tmp / "db")], check=True, capture_output=True)
    cols = "qseqid sseqid pident length qstart qend qlen bitscore"
    res = subprocess.run(
        ["blastn", "-task", "blastn", "-dust", "no", "-evalue", "1e-5",
         "-max_target_seqs", "500", "-query", str(fa), "-db", str(tmp / "db"),
         "-outfmt", f"6 {cols}"],
        check=True, capture_output=True, text=True)
    if not res.stdout.strip():
        return pd.DataFrame(columns=cols.split())
    df = pd.DataFrame([l.split("\t") for l in res.stdout.strip().split("\n")],
                      columns=cols.split())
    for c in ("pident", "bitscore"):
        df[c] = df[c].astype(float)
    for c in ("length", "qstart", "qend", "qlen"):
        df[c] = df[c].astype(int)
    return df[df["qseqid"] != df["sseqid"]]           # self-hits never count


def best_hits(df, min_qcov=0.80):
    """Collapse HSPs to one row per (query, subject).

    Identity is taken from the top-bitscore HSP; coverage is the union of that
    pair's HSP query intervals, so a stack of short spurious local matches
    cannot manufacture a high-coverage hit.
    """
    rows = []
    for (q, s), g in df.groupby(["qseqid", "sseqid"], sort=False):
        iv, cov = sorted(zip(g["qstart"], g["qend"])), 0
        cur_s, cur_e = iv[0]
        for a, b in iv[1:]:
            if a <= cur_e + 1:
                cur_e = max(cur_e, b)
            else:
                cov += cur_e - cur_s + 1
                cur_s, cur_e = a, b
        cov += cur_e - cur_s + 1
        top = g.loc[g["bitscore"].idxmax()]
        rows.append({"qseqid": q, "sseqid": s, "pident": float(top["pident"]),
                     "bitscore": float(top["bitscore"]),
                     "qcov": cov / int(top["qlen"])})
    out = pd.DataFrame(rows)
    return out[out["qcov"] >= min_qcov] if len(out) else out


# ---------------------------------------------------------------------------
# Scoring
# ---------------------------------------------------------------------------

def score(hits, labels, targets, option, spans):
    """Per-reference leave-one-out call + margin, over the whole panel."""
    rows = []
    for q, sp in labels.items():
        if sp not in targets:
            continue
        h = hits[hits["qseqid"] == q].copy()
        h["sp"] = h["sseqid"].map(labels)
        h = h[h["sp"].notna()]
        if h.empty:
            rows.append({"option": option, "reference": q, "true_species": sp,
                         "called_species": None, "best_pident": None,
                         "runnerup_species": None, "runnerup_pident": None,
                         "loo_margin": None, "n_conspecific": 0,
                         "spans_amplicon": spans.get(q, False),
                         "basis": "no_hit"})
            continue
        top = h.loc[h["bitscore"].idxmax()]
        own = h[h["sp"] == sp]
        oth = h[h["sp"] != sp]
        n_con = len(own)
        best_own = own["pident"].max() if n_con else None
        best_oth = oth["pident"].max() if len(oth) else None
        ru = (oth.loc[oth["pident"].idxmax()] if len(oth) else None)
        rows.append({
            "option": option, "reference": q, "true_species": sp,
            "called_species": top["sp"], "best_pident": round(float(top["pident"]), 4),
            "runnerup_species": (None if ru is None else ru["sp"]),
            "runnerup_pident": (None if ru is None else round(float(ru["pident"]), 4)),
            "loo_margin": (None if best_own is None or best_oth is None
                           else round(float(best_own) - float(best_oth), 4)),
            "n_conspecific": n_con,
            "spans_amplicon": spans.get(q, False),
            "basis": ("untestable_sole_reference" if n_con == 0
                      else "correct" if top["sp"] == sp else "misassigned"),
        })
    return rows


def pair_margins(hits, labels, targets, option):
    """Per species pair: the smallest own-minus-other margin on either side."""
    by_sp = {t: [k for k, v in labels.items() if v == t] for t in targets}
    rows = []
    for a, b in itertools.combinations(sorted(targets), 2):
        margins, n_a, n_b = [], 0, 0
        for focal, other in ((a, b), (b, a)):
            for q in by_sp[focal]:
                h = hits[hits["qseqid"] == q].copy()
                h["sp"] = h["sseqid"].map(labels)
                own = h[h["sp"] == focal]["pident"]
                oth = h[h["sp"] == other]["pident"]
                if own.empty or oth.empty:
                    continue
                margins.append(float(own.max()) - float(oth.max()))
                if focal == a:
                    n_a += 1
                else:
                    n_b += 1
        # Sidedness is not cosmetic. After the 2026-08-28 curation *P. fieldi*
        # has ONE MIT reference, so every fieldi pair is measured from the
        # other species' references only: "no simiovale reference is pulled to
        # fieldi" is a real result, but it is silent on whether the single
        # fieldi reference is pulled to simiovale. Calling that plain
        # `resolved` would repeat exactly the optimism this metric replaces.
        if not margins:
            verdict, m = "untestable_sole_reference", None
        else:
            m = min(margins)
            if m <= 0:
                verdict = "CONFUSABLE"
            elif n_a and n_b:
                verdict = "resolved"
            else:
                verdict = "resolved_one_side_only"
        rows.append({"option": option, "pair": f"{a} vs {b}",
                     "n_refs_a": len(by_sp[a]), "n_refs_b": len(by_sp[b]),
                     "n_margins_a": n_a, "n_margins_b": n_b,
                     "pair_margin_pct": (None if m is None else round(m, 4)),
                     "verdict": verdict})
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--alignment", required=True)
    ap.add_argument("--pairs", action="append", default=[],
                    metavar="OPTION=TSV", help="untiled option")
    ap.add_argument("--tiled", action="append", default=[],
                    metavar="OPTION=TSV", help="tiled option (union of tiles)")
    ap.add_argument("--species-targets", required=True)
    ap.add_argument("--unverified-in-db", choices=["exclude", "include"],
                    default="exclude")
    ap.add_argument("--min-qcov", type=float, default=0.80)
    ap.add_argument("--loo-out", required=True)
    ap.add_argument("--pair-out", required=True)
    args = ap.parse_args()

    targets = [s.strip() for s in args.species_targets.split(",") if s.strip()]
    records = list(SeqIO.parse(args.alignment, "fasta"))
    aln_len = len(records[0].seq)
    _, _, tpl_to_aln = build_consensus_template(records, aln_len)

    labels = {r.id: species_from_id(r.id) for r in records}
    unver = [k for k, v in labels.items() if v == UNVERIFIED]
    if args.unverified_in_db == "exclude" and unver:
        records = [r for r in records if r.id not in unver]
        for k in unver:
            labels.pop(k)
    print(f"[resolution_margin] {len(records)} references in the database "
          f"({args.unverified_in_db} unverified: {unver or 'none'})", flush=True)

    opts = ([(s, False) for s in args.pairs] + [(s, True) for s in args.tiled])
    loo_rows, pair_rows = [], []
    for spec, tiled in opts:
        option, tsv = spec.split("=", 1)
        ivs, probes, chosen = representative_interval(Path(tsv), tpl_to_aln, tiled)
        bp = sum(b - a + 1 for a, b in ivs)
        print(f"[resolution_margin] {option}: {chosen} tiles {probes} "
              f"-> extracted {ivs} = {bp} cols", flush=True)
        seqs = extract(records, ivs, probes)
        spans = {k: v[1] for k, v in seqs.items()}
        with tempfile.TemporaryDirectory() as td:
            hits = best_hits(blast_all_vs_all(seqs, Path(td)), args.min_qcov)
        loo_rows += score(hits, labels, targets, option, spans)
        pair_rows += pair_margins(hits, labels, targets, option)

    for rows, out in ((loo_rows, args.loo_out), (pair_rows, args.pair_out)):
        Path(out).parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(rows).to_csv(out, sep="\t", index=False)
        print(f"[resolution_margin] {len(rows)} rows -> {out}", flush=True)

    p = pd.DataFrame(pair_rows)
    print("\n[resolution_margin] pair verdicts by option:", flush=True)
    print(p.groupby(["option", "verdict"]).size().to_string(), flush=True)
    bad = p[p["verdict"] != "resolved"]
    if len(bad):
        print("\n[resolution_margin] pairs not cleanly resolved:", flush=True)
        print(bad.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
