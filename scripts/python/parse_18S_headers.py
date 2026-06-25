#!/usr/bin/env python3
"""
parse_18S_headers.py — extract accession + species from a multi-FASTA's
header lines and write a TSV mapping table.

Provenance: mit_similarity.Rmd lines 198-237 (legacy 18S header parser).

Headers in the 18S DB come in three flavours; the regex ladder below
mirrors the legacy script exactly so the species call is bit-equivalent
to the legacy `outputs/18S_acc_to_header.tsv`.

Usage:
  parse_18S_headers.py --fasta data/reference/18S_ref_db.fasta \\
      --out outputs/cross_species/18S_acc_to_header.tsv
"""
import argparse
import re
from pathlib import Path

import pandas as pd


def extract_species(header: str):
    # 1) Abbrev genus: "P.falciparum"
    m = re.search(r"\bP\.\s*([A-Za-z0-9_-]+)\b", header)
    if m:
        return f"p.{m.group(1).lower()}"
    # 2) Full genus: "Plasmodium vivax"
    m = re.search(r"\bPlasmodium\s+([A-Za-z0-9_-]+)\b", header, flags=re.IGNORECASE)
    if m:
        return f"p.{m.group(1).lower()}"
    # 3) "P. vivax" with a space after the dot
    m = re.search(r"\bP\.\s+([A-Za-z0-9_-]+)\b", header)
    if m:
        return f"p.{m.group(1).lower()}"
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fasta", required=True, help="input multi-FASTA (e.g. 18S_ref_db.fasta)")
    ap.add_argument("--out",   required=True, help="output TSV (acc, species, header)")
    args = ap.parse_args()

    rows = []
    with open(args.fasta) as f:
        for line in f:
            if line.startswith(">"):
                header = line[1:].strip()
                acc = header.split()[0]
                rows.append({
                    "acc":     acc,
                    "species": extract_species(header),
                    "header":  header,
                })

    df = pd.DataFrame(rows)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, sep="\t", index=False)

    n_missing = df["species"].isna().sum()
    print(f"[parse_18S_headers] wrote {len(df)} rows to {args.out}", flush=True)
    if n_missing:
        print(f"[parse_18S_headers] {n_missing} headers had no species call",
              flush=True)


if __name__ == "__main__":
    main()
