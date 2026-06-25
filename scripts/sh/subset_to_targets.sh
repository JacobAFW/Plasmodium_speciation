#!/usr/bin/env bash
# subset_to_targets.sh — subset a multi-FASTA to a comma-separated species
# pattern using seqkit grep.
#
# Provenance: mit_similarity.Rmd lines 405-407 (`seqkit grep -r -p ...`).
# The `-n -i` flags are used for 18S (descriptive headers, case-insensitive
# match anywhere in the header) and dropped for MIT (IDs are self-describing
# and case-stable).
#
# Usage:
#   subset_to_targets.sh <in_fasta> <out_fasta> <species_pat> [seqkit_extra_flags...]

set -euo pipefail

IN="${1:?input fasta required}"
OUT="${2:?output fasta required}"
PAT="${3:?species pattern required}"
shift 3

mkdir -p "$(dirname "$OUT")"

seqkit grep -r -p "$PAT" "$@" "$IN" > "$OUT"

echo "[subset_to_targets] $(grep -c '^>' "$OUT") sequences → $OUT"
