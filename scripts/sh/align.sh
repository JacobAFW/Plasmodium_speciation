#!/usr/bin/env bash
# align.sh — run MAFFT --auto on a multi-FASTA.
#
# Provenance: mit_similarity.Rmd lines 422-424 (`mafft --auto`).
# MAFFT is deterministic with --auto on the same input + version, so the
# alignment file is the contract that downstream entropy expects.
#
# Usage:
#   align.sh <in_fasta> <out_fasta> [log_file]

set -euo pipefail

IN="${1:?input fasta required}"
OUT="${2:?output fasta required}"
LOG="${3:-/dev/null}"

mkdir -p "$(dirname "$OUT")"

mafft --auto "$IN" > "$OUT" 2> "$LOG"

echo "[align] $(grep -c '^>' "$OUT") sequences aligned → $OUT"
