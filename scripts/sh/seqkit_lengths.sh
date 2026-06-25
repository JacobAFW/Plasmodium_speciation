#!/usr/bin/env bash
# seqkit_lengths.sh — emit a (id, length) TSV for every sequence in a FASTA.
#
# Provenance: mit_similarity.Rmd lines 414-420.
# MIT used `awk '{print $1,$2}' OFS="\t"` after fx2tab; that's a no-op for the
# default fx2tab output (id, length) so we keep it as a single seqkit call.
#
# Usage:
#   seqkit_lengths.sh <in_fasta> <out_tsv>

set -euo pipefail

IN="${1:?input fasta required}"
OUT="${2:?output tsv required}"

mkdir -p "$(dirname "$OUT")"

seqkit fx2tab -n -l "$IN" > "$OUT"

echo "[seqkit_lengths] $(wc -l < "$OUT") sequences → $OUT"
