# 02_cross_species — cross-species filtering and length QC

Five rules:

- `parse_18S_headers` — Python regex ladder over the 18S FASTA's descriptive headers → `outputs/cross_species/18S_acc_to_header.tsv` (acc, species, header). Required by the 18S cross-species filter (the MIT IDs are self-describing so no parser is needed there).
- `cross_species_filter` (called per-marker) — applies the canonical strict filter `pident ≥ 99 AND cov_shorter ≥ 0.99 AND shorter_len ≥ 10`, after dropping self-hits and within-species hits, keeping the top 5 HSPs per query–subject pair. Two outputs:
  - `outputs/cross_species/{marker}_suspicious.tsv` — every HSP that survives the strict filter. Includes a `species_pair` column (alphabetised "psimium vs pvivax") matching the legacy MIT output (the Rmd as checked-in omits this column; see `MEMORY.md` Session 4).
  - `outputs/cross_species/{marker}_closest_per_seq.tsv` — best non-self hit per query, with derived `pct_covered` and `divergence_pct`.
- `length_table` — `seqkit fx2tab -n -l` on the target-species FASTA → `outputs/qc/{marker}_lengths.tsv` (id, len).
- `length_classification` — adds `pct_of_max`, `missing_bases`, `is_partial` (true if `len < 0.90 * max(len)`) → `outputs/qc/{marker}_length_classification.tsv`.

## Provenance

- `parse_18S_headers` ← `mit_similarity.Rmd` lines 198-237.
- `cross_species_filter` ← Rmd 25-139 (MIT) and 242-391 (18S); refactored to a single argv-driven R/tidyverse script that handles both markers via a `--marker` flag.
- `length_table` + `length_classification` ← Rmd 414-460.

## Parameters

Pulled from `config.yaml::suspicious`:

| Parameter | Value | Meaning |
|---|---|---|
| `pident_min` | 99 | percent identity threshold (0–100) |
| `cov_shorter_min` | 0.99 | coverage of the shorter sequence (0–1) |
| `shorter_len_min` | 10 | minimum shorter-sequence length to consider |
| `top_n_hsps` | 5 | how many HSPs to keep per query–subject pair before filtering |

`length_classification` uses a hard-coded 0.90 partial cutoff (matches the legacy).

## Validation

For both markers the strict suspicious + closest tables should match Step 2's replicated outputs. Refactor into R introduces formatting drift (`TRUE`/`FALSE` not `True`/`False`; whole-number floats printed as `1` not `1.0`; trailing whitespace stripped from FASTA IDs by `data.table::fread`) — content-equivalent, validated row-count + content-after-rounding-to-6dp.

## Interpretation

`{marker}_suspicious.tsv` is the canonical cross-species ambiguity flag. The headline result: MIT flags 8 vivax/simium hits; 18S flags zero pairs at the strict threshold. (The looser-threshold 18S re-pass that surfaces partial ambiguities is Phase 4.)

`{marker}_closest_per_seq.tsv` is the per-reference-sequence "what is your closest non-self relative" summary — useful for spotting reference-DB outliers and pairs sitting just below the strict threshold.
