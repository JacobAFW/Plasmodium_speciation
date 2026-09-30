#!/usr/bin/env python3
"""
insilico_pcr.py — in-silico PCR of primer pairs against a subject FASTA,
using blastn to locate candidate binding sites and primer3 thermodynamics
to score them.

WHY blastn AND NOT MFEprimer
----------------------------
`CONTEXT-software.md` and `HANDOFF.md` (Step 5) both specify
`bioconda::mfeprimer`. That package does not exist — a repoquery against
bioconda for osx-arm64, osx-64, noarch and linux-64 returns no entries on
any platform. Jacob chose (2026-08-17) to build Step 5 on the already-
installed blastn rather than fetch an MFEprimer release binary from GitHub,
so this script is the substitute engine. See
`workflow/rules/11_specificity.README.md` for the efficiency model and how
it maps onto MFEprimer's `max_efficiency` threshold.

METHOD
------
1. Every primer of every pair is written as a query and blastn'd
   (`-task blastn-short`) against the subject to LOCATE candidate binding
   sites. blast is used only as a fast, sensitive locator — none of its
   scores enter the verdict.
2. Each located site is re-derived exactly: the hit's subject coordinates
   are extended to the primer's full footprint, that window is pulled from
   the subject, reverse-complemented if the primer binds the minus strand,
   and compared base-by-base with the primer. This gives exact mismatch
   counts (including at the 3' terminus) and the exact duplex to score.
   Gapped binding is not modelled — an indel inside a ~20-mer duplex does
   not prime in practice.
3. Sites are paired into amplicons: a plus-strand binder upstream of a
   minus-strand binder, product length = (minus-strand site 5'-most subject
   coord) - (plus-strand site 5'-most subject coord) + 1. Both primer
   orientations are tried, so a reference stored in either orientation is
   found.
4. Each amplicon is scored (see `efficiency()`).

Usage:
  insilico_pcr.py --pairs outputs/primer_design/mit_pairs_filtered.tsv \\
      --subject data/reference/mit_all.target.fasta \\
      --subject-kind reference \\
      --amplicons-out outputs/specificity/mit_target_amplicons.tsv \\
      --summary-out   outputs/specificity/mit_target_amplification.tsv
"""
import argparse
import math
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd
import primer3
from Bio import SeqIO
from Bio.Seq import Seq

# --------------------------------------------------------------------------
# Species normalisation — the canonical ladder from CONTEXT-data.md, applied
# to the FASTA header (not to a pre-parsed token), so that headers like
# "Plasmodium cf. malariae type 2" resolve to `malariae` rather than to the
# `p.cf` that a naive genus-word regex yields.
# --------------------------------------------------------------------------
# ORDER IS LOAD-BEARING: "simiovale" contains the substring "ovale" and the
# first pattern that matches wins, so simiovale (and simium, for symmetry) sit
# ABOVE ovale. This now matches CONTEXT-data.md, which was itself corrected on
# 2026-08-28 — until then this file carried a local prepend to work around the
# documented order being wrong.
_LADDER = [
    ("falciparum", "falciparum"),
    ("vivax",      "vivax"),
    ("knowlesi",   "knowlesi"),
    ("malariae",   "malariae"),
    ("simiovale",  "simiovale"),
    ("simium",     "simium"),
    ("ovale",      "ovale"),
    ("coat",       "coatneyi"),
    ("inui",       "inui"),
    ("fieldi",     "fieldi"),
    ("cyno|cynomolgi|pcynomolgi", "cynomolgi"),
]


def normalise_species(header: str):
    """Map a FASTA header to a panel species, or None if it is off-panel.

    Mirrors `normalise_species()` in CONTEXT-data.md, including the
    `pcynomologi -> pcynomolgi` repair applied before the ladder. Order
    matters: `simiovale` must be tested before `ovale` would otherwise
    claim it, which is why the ladder is a list and not a dict.

    Returns None for any off-panel header. That includes
    `Punverified_ATCC-30157_AB444133` (see the reference-curation note in
    workflow/rules/03_alignment.README.md), which is deliberate: a record of
    unresolved identity must not be scored as one of the 11 panel species.
    """
    x = header.lower().replace("pcynomologi", "pcynomolgi")
    for pattern, species in _LADDER:
        if re.search(pattern, x):
            return species
    return None


def revcomp(s: str) -> str:
    return str(Seq(s).reverse_complement())


# --------------------------------------------------------------------------
# Binding-site model
# --------------------------------------------------------------------------
class Site:
    """One exact primer footprint on the subject.

    `start`/`end` are 1-based inclusive subject coordinates of the footprint,
    always start <= end. `strand` is '+' when the primer sequence itself
    reads along the subject plus strand (so the primer's 3' end is at `end`),
    '-' when the primer binds the minus strand (3' end at `start`).
    """

    __slots__ = ("role", "start", "end", "strand", "n_mm", "n_mm3",
                 "tm", "dg", "eff")

    def __init__(self, role, start, end, strand, n_mm, n_mm3, tm, dg, eff):
        self.role, self.start, self.end, self.strand = role, start, end, strand
        self.n_mm, self.n_mm3 = n_mm, n_mm3
        self.tm, self.dg, self.eff = tm, dg, eff

    @property
    def five_prime(self) -> int:
        """Subject coordinate of the primer's 5' end — the amplicon boundary."""
        return self.start if self.strand == "+" else self.end


def count_mismatches(primer: str, site_seq: str, three_prime_window: int):
    """Exact (total, 3'-terminal-window) mismatch counts for one footprint.

    Split out of `efficiency()` so callers can reject a site on the cheap
    string comparison before paying for primer3's thermodynamic alignment.
    blast is a deliberately noisy locator (`-word_size 7`, `-evalue 1000`), so
    the overwhelming majority of located sites are junk that `--max-mismatch`
    discards; scoring them first is what made a genome-scale subject
    intractable. No verdict changes — the sites this lets us skip were being
    thrown away anyway.
    """
    n = min(len(primer), len(site_seq))
    comp = revcomp(site_seq)
    n_mm = sum(1 for a, b in zip(primer[:n], comp[:n]) if a != b)
    tail = min(three_prime_window, n)
    n_mm3 = sum(1 for a, b in zip(primer[n - tail:n], comp[n - tail:n]) if a != b)
    return n_mm, n_mm3


def efficiency(primer: str, site_seq: str, anneal_temp: float,
               three_prime_window: int, k: float):
    """Score one primer:template duplex in [0, 1].

    Substitute for MFEprimer's `max_efficiency`. Two terms, both standard
    PCR priming rules:

      * A mismatch anywhere in the 3'-terminal `three_prime_window` bases
        returns 0 — Taq cannot extend from an unpaired 3' end, so such a
        site does not prime regardless of how well the rest anneals.
      * Otherwise the score is a logistic in the duplex melting temperature
        relative to the annealing temperature:
            eff = 1 / (1 + exp(-(Tm_duplex - anneal_temp) / k))
        so eff = 0.5 exactly when the duplex melts at the annealing
        temperature. That is what makes HANDOFF's `max_efficiency > 0.5`
        threshold carry over unchanged and still mean something physical.

    `site_seq` is the subject strand the primer anneals to, 5'->3', already
    oriented so that it is the reverse complement of a perfectly-matching
    primer. Tm and dG come from primer3's thermodynamic alignment of the
    actual (possibly mismatched) duplex, so mismatches depress Tm rather
    than being counted heuristically.
    """
    n_mm, n_mm3 = count_mismatches(primer, site_seq, three_prime_window)

    res = primer3.calc_heterodimer(primer, site_seq)
    tm = res.tm if res.structure_found else float("-inf")
    dg = res.dg if res.structure_found else float("nan")

    if n_mm3 > 0 or not res.structure_found:
        eff = 0.0
    else:
        eff = 1.0 / (1.0 + math.exp(-(tm - anneal_temp) / k))
    return n_mm, n_mm3, tm, dg, eff


# --------------------------------------------------------------------------
# blastn locator
# --------------------------------------------------------------------------
_BLAST_FIELDS = ["qseqid", "sseqid", "qstart", "qend", "sstart", "send", "qlen", "slen"]


def locate_sites(primers: dict, subject: str, workdir: Path, threads: int,
                 word_size: int, db_dir: Path = None) -> pd.DataFrame:
    """blastn every primer against the subject; return raw hit coordinates."""
    # HARDCODE(scratch filename in a per-call TemporaryDirectory, not a project path; 2026-08-17): nothing outside this function reads it.
    qfa = workdir / "primers.fasta"
    with open(qfa, "w") as fh:
        for name, seq in primers.items():
            fh.write(f">{name}\n{seq}\n")

    # `--db-dir` keeps the blast database between runs. Building one over a
    # 3 Gbp host genome is ~20 min of the wall clock; a re-run of the off-target
    # half should not pay it twice. Default (None) keeps the old throwaway.
    if db_dir is None:
        dbdir = workdir / "db"
        dbdir.mkdir(exist_ok=True)
        db = dbdir / "subject"
    else:
        db_dir.mkdir(parents=True, exist_ok=True)
        db = db_dir / Path(subject).name
    if not Path(str(db) + ".nin").exists() and not Path(str(db) + ".nal").exists():
        subprocess.run(
            ["makeblastdb", "-in", subject, "-dbtype", "nucl", "-out", str(db)],
            check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE,
        )

    out = workdir / "hits.tsv"
    subprocess.run(
        ["blastn", "-task", "blastn-short",
         "-query", str(qfa), "-db", str(db),
         "-word_size", str(word_size),
         "-evalue", "1000",
         "-dust", "no", "-soft_masking", "false",
         "-max_target_seqs", "1000000",
         "-num_threads", str(threads),
         "-outfmt", "6 " + " ".join(_BLAST_FIELDS),
         "-out", str(out)],
        check=True, stderr=subprocess.PIPE,
    )

    if out.stat().st_size == 0:
        return pd.DataFrame(columns=_BLAST_FIELDS)
    # Narrow dtypes: a Gbp-scale subject yields tens of millions of hit rows,
    # and the default int64 + object columns are what would run the machine out
    # of memory long before the science gets interesting.
    return pd.read_csv(out, sep="\t", names=_BLAST_FIELDS,
                       dtype={"qseqid": "category", "sseqid": "category",
                              "qstart": "int32", "qend": "int32",
                              "sstart": "int32", "send": "int32",
                              "qlen": "int32", "slen": "int32"})


def exact_sites(hits: pd.DataFrame, primers: dict, subjects: dict,
                anneal_temp: float, three_prime_window: int, k: float,
                max_mismatch: int) -> dict:
    """Turn blast hits into exact, deduplicated, scored footprints.

    Keyed by (pair_id, subject_id) -> list[Site].
    """
    by_key: dict = {}
    seen: set = set()

    for h in hits.itertuples(index=False):
        primer = primers[h.qseqid]
        pair_id, role = h.qseqid.rsplit("|", 1)
        subj = subjects.get(h.sseqid)
        if subj is None:
            continue

        strand = "+" if h.sstart <= h.send else "-"
        # Extend the hit to the primer's full footprint. blast reports the
        # aligned segment only; the primer's unaligned flanks still occupy
        # subject bases and must be compared.
        if strand == "+":
            start = h.sstart - (h.qstart - 1)
            end = h.send + (h.qlen - h.qend)
        else:
            start = h.send - (h.qlen - h.qend)
            end = h.sstart + (h.qstart - 1)
        if start < 1 or end > len(subj):
            continue  # footprint runs off the end of the reference

        key = (pair_id, h.sseqid, role, start, end, strand)
        if key in seen:
            continue
        seen.add(key)

        window = subj[start - 1:end]
        # The strand the primer anneals to, 5'->3'.
        site_seq = revcomp(window) if strand == "+" else window
        # Cheap reject BEFORE primer3 — same discard set as scoring first.
        if count_mismatches(primer, site_seq, three_prime_window)[0] > max_mismatch:
            continue
        n_mm, n_mm3, tm, dg, eff = efficiency(
            primer, site_seq, anneal_temp, three_prime_window, k)

        by_key.setdefault((pair_id, h.sseqid), []).append(
            Site(role, start, end, strand, n_mm, n_mm3, tm, dg, eff))
    return by_key


def build_amplicons(sites: list, max_product: int):
    """Pair a plus-strand binder with a downstream minus-strand binder.

    Both primer-role assignments are tried (F upstream / R downstream, and
    R upstream / F downstream) so a reference stored in either orientation
    yields the amplicon.
    """
    plus = [s for s in sites if s.strand == "+"]
    minus = [s for s in sites if s.strand == "-"]
    out = []
    for p in plus:
        for m in minus:
            if m.role == p.role:
                continue  # both ends must be different primers of the pair
            length = m.five_prime - p.five_prime + 1
            if length <= 0 or length > max_product:
                continue
            out.append({
                "amplicon_len":  length,
                "amplicon_start": p.five_prime,
                "amplicon_end":   m.five_prime,
                "fwd_role":      p.role,
                "fwd_mm":        p.n_mm, "fwd_mm3": p.n_mm3,
                "fwd_tm":        round(p.tm, 2), "fwd_dg": round(p.dg, 1),
                "fwd_eff":       p.eff,
                "rev_role":      m.role,
                "rev_mm":        m.n_mm, "rev_mm3": m.n_mm3,
                "rev_tm":        round(m.tm, 2), "rev_dg": round(m.dg, 1),
                "rev_eff":       m.eff,
                "efficiency":    p.eff * m.eff,
            })
    return out


# --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pairs", nargs="+", required=True,
                    help="primer-pair TSV(s); need pair_id, f_seq, r_seq, amplicon_len")
    ap.add_argument("--subject", required=True, help="subject multi-FASTA")
    ap.add_argument("--subject-kind", default="reference",
                    help="free-text label recorded in the output (e.g. reference, offtarget)")
    ap.add_argument("--amplicons-out", required=True,
                    help="per-amplicon TSV (the evidence)")
    ap.add_argument("--summary-out",
                    help="optional per (pair_id, reference) TSV — the HANDOFF "
                         "Step 5 target-half shape. It is a FULL grid, so it is "
                         "worth writing for a reference set of tens of sequences "
                         "('did every target amplify?') and not for a genome of "
                         "thousands of contigs, where the off-target rollup is "
                         "the right grain and zero rows are the expectation")
    ap.add_argument("--coverage-out",
                    help="optional per (pair_id, species) rollup")
    ap.add_argument("--max-product", type=int, default=20000,
                    help="discard amplicons longer than this (bp)")
    ap.add_argument("--size-tolerance", type=float, default=0.25,
                    help="fraction of the designed amplicon_len an observed "
                         "amplicon may deviate by and still count as on-size")
    ap.add_argument("--three-prime-window", type=int, default=5,
                    help="a mismatch within this many 3'-terminal bases zeroes the primer")
    ap.add_argument("--max-mismatch", type=int, default=5,
                    help="discard binding sites with more total mismatches than this")
    ap.add_argument("--anneal-temp", type=float, default=55.0,
                    help="annealing temperature (C) the efficiency logistic is centred on")
    ap.add_argument("--logistic-k", type=float, default=3.0,
                    help="steepness (C) of the efficiency logistic")
    ap.add_argument("--word-size", type=int, default=7, help="blastn word size")
    ap.add_argument("--db-dir",
                    help="reuse/persist the blastn database here instead of a "
                         "throwaway temp dir (worth it for Gbp-scale subjects)")
    ap.add_argument("--threads", type=int, default=4)
    args = ap.parse_args()

    for tool in ("makeblastdb", "blastn"):
        if shutil.which(tool) is None:
            sys.exit(f"[insilico_pcr] ERROR: {tool} not on PATH — activate the env first")

    # ---- primer pairs -----------------------------------------------------
    frames = []
    for p in args.pairs:
        df = pd.read_csv(p, sep="\t", dtype={"pair_id": str})
        missing = {"pair_id", "f_seq", "r_seq", "amplicon_len"} - set(df.columns)
        if missing:
            sys.exit(f"[insilico_pcr] ERROR: {p} lacks columns {sorted(missing)}")
        df["pair_source"] = Path(p).name
        frames.append(df)
    pairs = pd.concat(frames, ignore_index=True)
    dupes = pairs["pair_id"][pairs["pair_id"].duplicated()].unique()
    if len(dupes):
        sys.exit(f"[insilico_pcr] ERROR: duplicate pair_id across inputs: {list(dupes)}")

    primers = {}
    for r in pairs.itertuples(index=False):
        primers[f"{r.pair_id}|F"] = r.f_seq.upper()
        primers[f"{r.pair_id}|R"] = r.r_seq.upper()
    print(f"[insilico_pcr] {len(pairs)} pairs / {len(primers)} primers", flush=True)

    # ---- subject ----------------------------------------------------------
    subjects, headers = {}, {}
    for rec in SeqIO.parse(args.subject, "fasta"):
        subjects[rec.id] = str(rec.seq).upper()
        headers[rec.id] = rec.description
    print(f"[insilico_pcr] {len(subjects)} subject sequences from {args.subject}",
          flush=True)

    species = {sid: normalise_species(h) for sid, h in headers.items()}
    n_unmapped = sum(1 for v in species.values() if v is None)
    if n_unmapped:
        print(f"[insilico_pcr] {n_unmapped} subject sequences did not map to a "
              f"panel species (recorded as species=NA, not dropped)", flush=True)

    # ---- locate + score ---------------------------------------------------
    with tempfile.TemporaryDirectory() as td:
        hits = locate_sites(primers, args.subject, Path(td), args.threads,
                            args.word_size,
                            Path(args.db_dir) if args.db_dir else None)
        print(f"[insilico_pcr] {len(hits)} raw blastn hits", flush=True)
        sites = exact_sites(hits, primers, subjects, args.anneal_temp,
                            args.three_prime_window, args.logistic_k,
                            args.max_mismatch)

    expected = dict(zip(pairs["pair_id"], pairs["amplicon_len"]))
    source = dict(zip(pairs["pair_id"], pairs["pair_source"]))

    amp_rows = []
    for (pair_id, sid), sl in sites.items():
        for a in build_amplicons(sl, args.max_product):
            exp = expected[pair_id]
            a.update({
                "pair_id":      pair_id,
                "pair_source":  source[pair_id],
                "reference":    sid,
                "species":      species[sid],
                "subject_kind": args.subject_kind,
                "expected_len": exp,
                "onsize":       abs(a["amplicon_len"] - exp) <= args.size_tolerance * exp,
            })
            amp_rows.append(a)

    amp_cols = ["pair_id", "pair_source", "reference", "species", "subject_kind",
                "amplicon_len", "expected_len", "onsize", "efficiency",
                "amplicon_start", "amplicon_end",
                "fwd_role", "fwd_mm", "fwd_mm3", "fwd_tm", "fwd_dg", "fwd_eff",
                "rev_role", "rev_mm", "rev_mm3", "rev_tm", "rev_dg", "rev_eff"]
    amps = (pd.DataFrame(amp_rows, columns=amp_cols)
            .sort_values(["pair_id", "reference", "amplicon_len"]))
    Path(args.amplicons_out).parent.mkdir(parents=True, exist_ok=True)
    amps.to_csv(args.amplicons_out, sep="\t", index=False)
    print(f"[insilico_pcr] wrote {len(amps)} amplicons to {args.amplicons_out}",
          flush=True)

    if not (args.summary_out or args.coverage_out):
        return

    # ---- summary: full (pair x reference) grid, zeros included -------------
    grid = pd.MultiIndex.from_product(
        [pairs["pair_id"].tolist(), list(subjects.keys())],
        names=["pair_id", "reference"]).to_frame(index=False)

    if len(amps):
        agg = (amps.groupby(["pair_id", "reference"], as_index=False)
               .agg(amplicon_count=("amplicon_len", "size"),
                    amplicon_lengths=("amplicon_len",
                                      lambda s: ",".join(map(str, sorted(s)))),
                    max_efficiency=("efficiency", "max"),
                    onsize_count=("onsize", "sum")))
    else:
        agg = pd.DataFrame(columns=["pair_id", "reference", "amplicon_count",
                                    "amplicon_lengths", "max_efficiency",
                                    "onsize_count"])

    summary = grid.merge(agg, on=["pair_id", "reference"], how="left")
    summary["amplicon_count"] = summary["amplicon_count"].fillna(0).astype(int)
    summary["onsize_count"] = summary["onsize_count"].fillna(0).astype(int)
    summary["amplicon_lengths"] = summary["amplicon_lengths"].fillna("")
    summary["max_efficiency"] = summary["max_efficiency"].fillna(0.0).round(4)
    summary["species"] = summary["reference"].map(species)
    summary["pair_source"] = summary["pair_id"].map(source)
    summary["subject_kind"] = args.subject_kind
    summary = summary[["pair_id", "pair_source", "species", "reference",
                       "subject_kind", "amplicon_count", "onsize_count",
                       "amplicon_lengths", "max_efficiency"]]
    if args.summary_out:
        summary.to_csv(args.summary_out, sep="\t", index=False)
        print(f"[insilico_pcr] wrote {len(summary)} (pair x reference) rows to "
              f"{args.summary_out}", flush=True)

    # ---- optional per-species rollup --------------------------------------
    if args.coverage_out:
        cov = summary[summary["species"].notna()].copy()
        cov["amplified"] = cov["onsize_count"] > 0
        roll = (cov.groupby(["pair_id", "species"], as_index=False)
                .agg(n_refs=("reference", "size"),
                     n_amplified=("amplified", "sum"),
                     max_efficiency=("max_efficiency", "max")))
        roll["frac_refs_amplified"] = (roll["n_amplified"] / roll["n_refs"]).round(3)
        roll["pair_source"] = roll["pair_id"].map(source)
        roll = roll[["pair_id", "pair_source", "species", "n_refs",
                     "n_amplified", "frac_refs_amplified", "max_efficiency"]]
        roll.to_csv(args.coverage_out, sep="\t", index=False)
        print(f"[insilico_pcr] wrote {len(roll)} (pair x species) rows to "
              f"{args.coverage_out}", flush=True)


if __name__ == "__main__":
    main()
