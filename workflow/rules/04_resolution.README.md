# 04_resolution — Phase 4 extension targets

NEW (not in the legacy `.Rmd`). Implements HANDOFF.md Step 4 acceptance criteria 2-5.

Four rules:

- `pair_summary` (per marker) — collapses cross-species BLAST hits to one row per unordered species pair: `pair, n_hits, best_pident, max_aln_len, max_cov_shorter, max_mincov, min_mismatch, min_gapopen, strict_ambiguity, ladder_class`. The strict-ambiguity column re-derives the `pident ≥ 99 AND cov_shorter ≥ 0.99 AND shorter_len ≥ 10` condition; `ladder_class` is the highest-priority class (`high_confidence_possible_ambiguity` > `moderate_possible_ambiguity` > `subject_partial` > `query_partial` > `short_high_identity_local_match` > `lower_concern`) achieved by any HSP in the pair.
- `cross_species_filter_18S_loose` — same script as the strict 18S filter, called with `--pident-min 98 --cov-min 0.85 --cov-metric mincov`. Output `18S_suspicious_loose.tsv`. The strict pipeline's 18S file stays untouched.
- `species_coverage_18S` — counts how many 18S sequences map to each of the 11 canonical target species via `18S_acc_to_header.tsv`. Output `18S_species_coverage.tsv` with `species, n_sequences, is_covered`.
- `cross_marker_resolution` — the v1 deliverable. One row per unordered panel-species pair (`n*(n-1)/2 = 55` for 11 species). Columns: `pair, mit_resolved, 18S_resolved, 18S_resolved_loose, combined_resolved, notes`.

## Resolution semantics

For one marker on one pair (A, B):
- **NA** if either species has zero references for that marker (we can't decide).
- **FALSE** if the pair appears in the marker's strict (or loose, for the loose column) suspicious table.
- **TRUE** otherwise.

`combined_resolved` follows three-valued OR logic: TRUE if either marker is TRUE; FALSE if both are FALSE; NA only if both are NA. Mixed FALSE/NA collapses to FALSE; mixed TRUE/NA collapses to TRUE.

## Notes column

Hard-coded annotations the human reader benefits from at-a-glance:
- "MIT vivax/simium ambiguity (P. simium recently zoonotic from P. vivax)" — flags the canonical example.
- "18S coverage low: A=n, B=m" — fires whenever either species has ≤ 1 reference. Most pairs trigger this on the current panel; it's a confidence signal, not a defect.

## Validation gate

`outputs/cross_species/resolution_table.tsv` must reproduce the vivax/simium MIT ambiguity finding from the legacy outputs. If `simium vs vivax` ever shows `mit_resolved = TRUE`, something is wrong upstream.

## Headline numbers from the current run

| Metric | Value |
|---|---|
| Cross-species pairs (panel) | 55 |
| MIT strict-ambiguity pairs | 1 (`simium vs vivax`) |
| 18S strict-ambiguity pairs | 0 |
| 18S looser-threshold pairs | 1 (`cynomolgi vs fieldi`, partial) |
| 18S coverage covered | 10/11 (missing: `simiovale`) |
| `mit_resolved = TRUE` | 54 / 55 |
| `mit_resolved = FALSE` | 1 / 55 |
| `18S_resolved = TRUE` | 45 / 55 |
| `18S_resolved = NA` | 10 / 55 (any pair with `simiovale`) |
| `combined_resolved = TRUE` | 55 / 55 |

Interpretation: every panel pair is resolvable by at least one marker. The single MIT failure (vivax/simium) is rescued by 18S, supporting the dual-marker rationale. The new looser 18S pair (cynomolgi/fieldi) is a partial-ambiguity finding the strict filter is too aggressive to surface — the report should flag it as "watch this pair" rather than "unresolvable".

## Provenance

`HANDOFF.md` Step 4 + `CONTEXT-analysis.md` v1 acceptance criteria 2-5. No legacy provenance — this is the v1 extension.
