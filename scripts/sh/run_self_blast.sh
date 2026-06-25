#!/usr/bin/env bash
# run_self_blast.sh — build a nucleotide BLAST DB from a multi-FASTA and run
# a self-BLAST against it, emitting the canonical 10-column tabular format.
#
# Provenance: mit_similarity.Rmd lines 8-20 (MIT) and 185-196 (18S).
#
# Usage:
#   run_self_blast.sh <fasta> <db_prefix> <out_tsv>

set -euo pipefail

FASTA="${1:?fasta path required}"
DB_PREFIX="${2:?BLAST db prefix required}"
OUT_TSV="${3:?output tsv path required}"

mkdir -p "$(dirname "$DB_PREFIX")" "$(dirname "$OUT_TSV")"

# Build the DB. -logfile lives next to the DB so logs don't pile up at PWD.
makeblastdb \
  -in "$FASTA" \
  -dbtype nucl \
  -out "$DB_PREFIX" \
  -logfile "${DB_PREFIX}.makeblastdb.log"

# Self-BLAST. The 10 outfmt columns are the contract every downstream rule
# reads — DO NOT reorder without updating cross_species_filter.R too.
blastn \
  -query "$FASTA" \
  -db    "$DB_PREFIX" \
  -out   "$OUT_TSV" \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qlen slen evalue bitscore" \
  -task   blastn

echo "[run_self_blast] $(wc -l < "$OUT_TSV") hits → $OUT_TSV"
