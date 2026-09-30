#!/usr/bin/env bash
# concat_tsv.sh — concatenate TSVs that share a header into one TSV.
#
# Keeps the first file's header, drops the others', and fails loudly if any
# input's header differs — a silent column-order mismatch would corrupt the
# merged table without any downstream rule noticing.
#
# Usage:
#   concat_tsv.sh OUT.tsv IN1.tsv IN2.tsv [...]
set -euo pipefail

if [ "$#" -lt 2 ]; then
  echo "usage: concat_tsv.sh OUT.tsv IN1.tsv [IN2.tsv ...]" >&2
  exit 2
fi

out="$1"; shift
mkdir -p "$(dirname "$out")"

header="$(head -1 "$1")"
for f in "$@"; do
  h="$(head -1 "$f")"
  if [ "$h" != "$header" ]; then
    echo "ERROR: header mismatch in $f" >&2
    echo "  expected: $header" >&2
    echo "  found:    $h" >&2
    exit 1
  fi
done

{
  printf '%s\n' "$header"
  for f in "$@"; do
    tail -n +2 "$f"
  done
} > "$out"

echo "[concat_tsv] wrote $(($(wc -l < "$out") - 1)) rows to $out"
