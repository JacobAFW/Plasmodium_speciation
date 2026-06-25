#!/usr/bin/env python3
"""
interior_threshold_sweep.py — calibrate the interior-side candidate-window
thresholds. Sweeps `(cov_thr, min_run)` over a grid, invokes seq_entropy.py
per cell, and records metrics from the resulting interior_candidates.tsv.

Provenance: Phase 5 Step 2 (HANDOFF). The interior-side output was added to
seq_entropy.py in Step 2; this script picks the default thresholds for the
benchmark-side calls that downstream chapters will reference.

Note: each cell runs seq_entropy.py to a temp directory, so the canonical
outputs/entropy/*.interior_candidates.tsv files are not touched.

Usage:
  interior_threshold_sweep.py \
    --mit-aln outputs/alignment/ma_mit.target.fasta \
    --s18-aln outputs/alignment/ma_18S.target.fasta \
    --cov-thrs 0.70,0.80,0.85,0.90 \
    --min-runs 15,18,21,25 \
    --entropy-thr 0.20 \
    --windows 25,50 \
    --imwong-start 1683 --imwong-end 1740 \
    --out outputs/sensitivity/interior_threshold_sweep.tsv
"""
import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
SEQ_ENTROPY = SCRIPT_DIR / "seq_entropy.py"


def sweep_cell(marker: str, aln_path: str, cov_thr: float, min_run: int,
               entropy_thr: float, windows: str,
               imwong_start: int, imwong_end: int) -> dict:
    """Run seq_entropy.py with given params; read interior_candidates and
    summarise. Returns a dict ready for a DataFrame row."""
    row = {
        "marker": marker,
        "cov_thr": cov_thr,
        "min_run": min_run,
        "entropy_thr": entropy_thr,
        "n_windows": None,
        "n_blocks": None,
        "median_entropy": None,
        "median_coverage": None,
        "imwong_recovered": None,
        "imwong_n_windows": None,
        "core_left": None,
        "core_right": None,
        "core_len": None,
        "error": None,
    }
    with tempfile.TemporaryDirectory() as td:
        prefix = Path(td) / marker
        cmd = [
            sys.executable, str(SEQ_ENTROPY),
            "--aln", aln_path,
            "--outprefix", str(prefix),
            "--cov_thr", str(cov_thr),
            "--min_run", str(min_run),
            "--entropy_thr", str(entropy_thr),
            "--windows", windows,
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            row["error"] = (result.stderr or "")[:300].strip()
            return row

        interior_path = Path(f"{prefix}.interior_candidates.tsv")
        core_path     = Path(f"{prefix}.core_bounds.tsv")
        interior = pd.read_csv(interior_path, sep="\t")
        core     = pd.read_csv(core_path, sep="\t")

        row["n_windows"] = int(len(interior))
        if len(interior):
            row["median_entropy"]  = float(interior["mean_entropy"].median())
            row["median_coverage"] = float(interior["mean_coverage"].median())

            # Count "blocks" = number of distinct (window-size) contiguous
            # runs of admitting windows. Reduces overlap noise when comparing
            # cells — 25 consecutive passing 25-bp windows is one block.
            block_count = 0
            for w_size, grp in interior.groupby("window"):
                starts = sorted(grp["start"].tolist())
                if not starts:
                    continue
                runs = 1
                for prev, cur in zip(starts, starts[1:]):
                    if cur - prev > 1:
                        runs += 1
                block_count += runs
            row["n_blocks"] = block_count

        row["core_left"]  = int(core["core_left"].iloc[0])
        row["core_right"] = int(core["core_right"].iloc[0])
        row["core_len"]   = int(core["core_len"].iloc[0])

        # Imwong recovery applies only to 18S (where the published forward
        # primers bind). For MIT the column is left blank.
        if marker == "18S":
            imwong = interior[(interior["start"] <= imwong_end)
                              & (interior["end"]   >= imwong_start)]
            row["imwong_n_windows"] = int(len(imwong))
            row["imwong_recovered"] = bool(len(imwong) > 0)

    return row


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mit-aln", required=True)
    ap.add_argument("--s18-aln", required=True)
    ap.add_argument("--cov-thrs", default="0.70,0.80,0.85,0.90")
    ap.add_argument("--min-runs", default="15,18,21,25")
    ap.add_argument("--entropy-thr", type=float, default=0.20)
    ap.add_argument("--windows", default="25,50")
    ap.add_argument("--imwong-start", type=int, default=1683)
    ap.add_argument("--imwong-end",   type=int, default=1740)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cov_thrs = [float(x) for x in args.cov_thrs.split(",") if x.strip()]
    min_runs = [int(x)   for x in args.min_runs.split(",") if x.strip()]
    print(f"Grid: {len(cov_thrs)} cov_thrs × {len(min_runs)} min_runs × 2 markers "
          f"= {len(cov_thrs) * len(min_runs) * 2} cells", flush=True)

    rows = []
    for marker, aln in (("mit", args.mit_aln), ("18S", args.s18_aln)):
        for cov in cov_thrs:
            for mr in min_runs:
                row = sweep_cell(
                    marker, aln, cov, mr,
                    entropy_thr=args.entropy_thr,
                    windows=args.windows,
                    imwong_start=args.imwong_start,
                    imwong_end=args.imwong_end,
                )
                rows.append(row)
                tag = f"[{marker}] cov={cov} min_run={mr:>2}"
                if row["error"]:
                    print(f"  {tag}: ERROR — {row['error']}", flush=True)
                else:
                    extra = ""
                    if marker == "18S":
                        extra = f"  imwong={'Y' if row['imwong_recovered'] else 'N'}"
                    print(f"  {tag}: n_windows={row['n_windows']:>6}  "
                          f"n_blocks={row['n_blocks']:>3}  "
                          f"core={row['core_left']}..{row['core_right']}{extra}",
                          flush=True)

    df = pd.DataFrame(rows)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, sep="\t", index=False)
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
