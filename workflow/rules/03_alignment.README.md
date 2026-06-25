# 03_alignment — alignment, entropy, primer-candidate windows

Five rules per marker:

- `subset_to_targets_{marker}` — `seqkit grep -r -p '<species_pat>'` over the full reference FASTA → `data/reference/{marker}.target.fasta`. The 18S call uses `-n -i` (match against the full descriptive header, case-insensitive); MIT uses neither (IDs are self-describing and case-stable). Subsets are *written back* to `data/reference/` because they're long-lived inputs the length-QC and entropy stages both consume.
- `align` — `mafft --auto` on the target subset → `outputs/alignment/ma_{marker}.target.fasta`. Deterministic given the same input + version, so the alignment file is the contract that downstream entropy expects.
- `entropy` — `scripts/python/seq_entropy.py` (copied unchanged from `scripts/legacy/speciation_long/scripts/`) emits four files from one `--outprefix` argument:
  - `{marker}.per_position.tsv` — pos, coverage, snp_count, entropy
  - `{marker}.windows.tsv` — sliding-window scores at multiple window sizes
  - `{marker}.primer_candidates.tsv` — windows that pass low-entropy / high-coverage thresholds, classified by which end of the alignment they sit in (`left` / `right` / `none`)
  - `{marker}.core_bounds.tsv` — alignment-wide left/right core boundaries used by the entropy classifier
- `plot_entropy` — refactored from `scripts/legacy/.../plot_entropy.R` into an argv-driven script. Produces four figures per marker (entropy + coverage, each as `.png` and `.svg`) under `reports/figures/`. Highlights candidate primer windows and internal high-entropy blocks.

## Provenance

`mit_similarity.Rmd` lines 405-407 (subset), 422-424 (mafft), 464 (entropy). The legacy `plot_entropy.R` has been broken out into one parameterised script callable from Snakemake.

## Validation

Alignments and entropy TSVs match the legacy outputs (Step 2). Float-printing drift (`0.9544340029249649` vs `0.954434002924965`) appears in entropy outputs as a numpy-2.x / pandas-3.0 cosmetic difference; values are identical at any reasonable precision.

## Interpretation

`{marker}.per_position.tsv` and `{marker}.windows.tsv` are the canonical inputs for the v1 primer-design conversation. The `primer_candidates.tsv` set is a starting point — actual primer3 selection inside those windows is Phase 5+ work, deferred per `HANDOFF.md`.

The PNGs are for embedding in the Quarto report; the SVGs are the editable handoff format collaborators can tweak in Inkscape / Figma.
