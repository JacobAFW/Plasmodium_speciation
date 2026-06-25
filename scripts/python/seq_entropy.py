#!/usr/bin/env python3
import argparse
import math
from collections import Counter
from pathlib import Path

import pandas as pd
from Bio import SeqIO


def shannon_entropy_basefreq(counts: Counter) -> float:
    """Shannon entropy over base frequencies (A,C,G,T only)."""
    total = sum(counts.values())
    if total == 0:
        return 0.0
    ent = 0.0
    for _, c in counts.items():
        p = c / total
        ent -= p * math.log2(p)
    return ent


def load_alignment(fasta_path: str):
    records = list(SeqIO.parse(fasta_path, "fasta"))
    if not records:
        raise ValueError(f"No sequences found in {fasta_path}")
    lengths = {len(r.seq) for r in records}
    if len(lengths) != 1:
        raise ValueError(
            f"Alignment is not a fixed length. Found lengths: {sorted(lengths)[:10]} ..."
        )
    aln_len = next(iter(lengths))
    return records, aln_len


def per_position_metrics(records, aln_len: int) -> pd.DataFrame:
    nseq = len(records)

    rows = []
    for i in range(aln_len):
        col = [str(r.seq[i]).upper() for r in records]

        nongap = [b for b in col if b not in ("-", ".")]
        coverage = len(nongap) / nseq

        # Consider A/C/G/T only for entropy and SNP count
        acgt = [b for b in nongap if b in ("A", "C", "G", "T")]
        base_counts = Counter(acgt)

        # SNP count: number of distinct A/C/G/T alleles minus 1 (0 means invariant)
        distinct = len(base_counts)
        snp_count = max(0, distinct - 1)

        entropy = shannon_entropy_basefreq(base_counts)

        rows.append(
            {
                "pos": i + 1,  # 1-based
                "coverage": coverage,
                "snp_count": snp_count,
                "entropy": entropy,
            }
        )

    return pd.DataFrame(rows)


def rolling_windows(df_pos: pd.DataFrame, window_sizes=(25, 50)) -> pd.DataFrame:
    aln_len = int(df_pos["pos"].max())

    out = []
    for w in window_sizes:
        for start in range(1, aln_len - w + 2):
            end = start + w - 1
            sl = df_pos[(df_pos["pos"] >= start) & (df_pos["pos"] <= end)]

            out.append(
                {
                    "window": w,
                    "start": start,
                    "end": end,
                    "mean_coverage": sl["coverage"].mean(),
                    "mean_entropy": sl["entropy"].mean(),
                    "sum_snp_count": sl["snp_count"].sum(),
                }
            )
    return pd.DataFrame(out)


def find_core_bounds(df_pos: pd.DataFrame, cov_thr: float, min_run: int = 50):
    """
    Find left/right bounds of the 'core' region where per-position coverage >= cov_thr.
    min_run requires a consecutive run of high-coverage positions to avoid single-column noise.
    Returns (core_left, core_right) as 1-based inclusive positions.
    """
    high = (df_pos["coverage"] >= cov_thr).tolist()
    aln_len = int(df_pos["pos"].max())

    # Find left: first position starting a run of length min_run
    core_left = None
    run = 0
    for i, ok in enumerate(high, start=1):  # 1-based
        run = run + 1 if ok else 0
        if run >= min_run:
            core_left = i - min_run + 1
            break

    # Find right: last position ending a run of length min_run
    core_right = None
    run = 0
    for i, ok in enumerate(reversed(high), start=1):
        run = run + 1 if ok else 0
        if run >= min_run:
            core_right = aln_len - (i - 1)
            break

    if core_left is None or core_right is None or core_left >= core_right:
        raise ValueError(
            f"Could not determine core bounds (cov_thr={cov_thr}, min_run={min_run}). "
            "Try lowering cov_thr or min_run."
        )

    return core_left, core_right


def label_end_region(df_win: pd.DataFrame, core_left: int, core_right: int, end_bp: int) -> pd.DataFrame:
    """
    Label windows near the ends of the *core* region.
    core_left/core_right are 1-based inclusive positions.
    """
    # Clamp so end_bp can't push beyond core bounds
    left_end = min(core_right, core_left + end_bp - 1)
    right_start = max(core_left, core_right - end_bp + 1)

    def region(row):
        if row["start"] >= core_left and row["end"] <= left_end:
            return "left"
        if row["start"] >= right_start and row["end"] <= core_right:
            return "right"
        return "none"

    df_win = df_win.copy()
    df_win["end_region"] = df_win.apply(region, axis=1)
    return df_win


def main():
    ap = argparse.ArgumentParser(
        description="Compute per-position coverage/SNP/entropy and rolling window summaries from an aligned FASTA."
    )
    ap.add_argument("--aln", required=True, help="Aligned FASTA (MAFFT/MUSCLE output).")
    ap.add_argument("--outprefix", required=True, help="Output prefix path, e.g. outputs/mit")
    ap.add_argument(
        "--end_bp",
        type=int,
        default=300,
        help="Define 'near ends' as within this many bp of the CORE region ends.",
    )
    ap.add_argument("--cov_thr", type=float, default=0.90, help="Coverage threshold for primer candidate windows.")
    ap.add_argument("--entropy_thr", type=float, default=0.20, help="Mean entropy threshold for primer candidate windows.")
    ap.add_argument("--windows", default="25,50", help="Comma-separated window sizes, e.g. 25,50")
    ap.add_argument(
        "--min_run",
        type=int,
        default=50,
        help="Consecutive high-coverage positions required to define core bounds.",
    )
    args = ap.parse_args()

    window_sizes = tuple(int(x) for x in args.windows.split(",") if x.strip())
    outprefix = Path(args.outprefix)
    outprefix.parent.mkdir(parents=True, exist_ok=True)

    records, aln_len = load_alignment(args.aln)

    # 1) per-position
    df_pos = per_position_metrics(records, aln_len)
    core_left, core_right = find_core_bounds(df_pos, cov_thr=args.cov_thr, min_run=args.min_run)
    core_len = core_right - core_left + 1
    print(f"Core bounds (1-based): {core_left}-{core_right} (len {core_len})")

    df_pos.to_csv(f"{outprefix}.per_position.tsv", sep="\t", index=False)

    # Write core bounds to a small TSV (handy for logs / pipeline)
    with open(f"{outprefix}.core_bounds.tsv", "w") as out:
        out.write("core_left\tcore_right\tcore_len\tcov_thr\tmin_run\n")
        out.write(f"{core_left}\t{core_right}\t{core_len}\t{args.cov_thr}\t{args.min_run}\n")

    # 2) windows
    df_win = rolling_windows(df_pos, window_sizes=window_sizes)
    df_win = label_end_region(df_win, core_left=core_left, core_right=core_right, end_bp=args.end_bp)
    df_win.to_csv(f"{outprefix}.windows.tsv", sep="\t", index=False)

    # 3) primer candidates — the DESIGN-side output. End-anchored windows only
    # (within end_bp of core_left/core_right). This is what primer3 selection
    # reads when picking long-amplicon flanking primers.
    primer = df_win[
        (df_win["end_region"].isin(["left", "right"])) &
        (df_win["mean_coverage"] >= args.cov_thr) &
        (df_win["mean_entropy"] <= args.entropy_thr)
    ].copy()

    # Sort: best (lowest entropy) then highest coverage
    primer = primer.sort_values(["end_region", "mean_entropy", "mean_coverage"], ascending=[True, True, False])
    primer.to_csv(f"{outprefix}.primer_candidates.tsv", sep="\t", index=False)

    # 4) interior candidates — the BENCHMARK-side output. Same entropy +
    # coverage gates, no end_region filter. Used to ask "is this published
    # primer site a defensible candidate by entropy/coverage standards?"
    # without the design-bias of restricting to amplicon ends. Step 2 of
    # Phase 5 (the Imwong forward-primer position lives at end_region=none).
    interior = df_win[
        (df_win["mean_coverage"] >= args.cov_thr) &
        (df_win["mean_entropy"] <= args.entropy_thr)
    ].copy()
    interior = interior.sort_values(["window", "mean_entropy", "mean_coverage"],
                                     ascending=[True, True, False])
    interior.to_csv(f"{outprefix}.interior_candidates.tsv", sep="\t", index=False)

    print(f"Alignment length: {aln_len} bp; sequences: {len(records)}")
    print(f"Wrote: {outprefix}.per_position.tsv")
    print(f"Wrote: {outprefix}.windows.tsv")
    print(f"Wrote: {outprefix}.core_bounds.tsv")
    print(f"Wrote: {outprefix}.primer_candidates.tsv  ({len(primer)} rows, design-side)")
    print(f"Wrote: {outprefix}.interior_candidates.tsv ({len(interior)} rows, benchmark-side)")


if __name__ == "__main__":
    main()