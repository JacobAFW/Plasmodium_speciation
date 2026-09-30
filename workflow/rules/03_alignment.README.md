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

## Reference curation (applied 2026-08-28)

The archive under `data/reference/` is the complete record and is **never
pruned**. Curation happens where the panel ("target") set is built, so every
decision is one config line and reverting it is deleting that line.

### 1. `AB444133` relabelled, not deleted

```
- >Pfieldi_AB444133.1_strain_ATCC_30157
+ >Punverified_ATCC-30157_AB444133
```

Applied in `data/reference/mit_all.fasta` and `mit_all.target.fasta`. The msp1
deposit of the same strain, `AB444066.1`, is relabelled the same way in
`data/reference/msp1/{fieldi_msp1,msp1_ref_db,msp1_ref_db.no_chrom}.fasta`.
The accession is kept in the name so the record stays traceable.

**Why.** ATCC 30157 is not *P. fieldi*. It carries the *P. simiovale* state at
**12 of 12** fixed mitochondrial diagnostic sites (7/7 cox1, 5/5 cytb); it sits
99.883 % from the simiovale type anchor against 99.130 % from its nearest
genuine *fieldi*; and the independent nuclear *msp1* deposit of the same strain
agrees (94.485 % to simiovale, 90.090 % to its nearest genuine fieldi — closer
to simiovale than the two genuine fieldi strains are to *each other*, 91.708 %).
A sequence that genuinely blurred two species would sit *between* them; this one
sits *inside* simiovale, closer to the anchor than the simiovale field isolates
are. Full evidence and the leave-one-out audit that guarded against
cherry-picking:
`simiovale_fieldi_aside/outputs/simiovale_fieldi_report_corrected.md` §4.

**Why `unverified` and not `simiovale`.** Whether ATCC 30157 was misidentified
at source or the deposited sequences were swapped cannot be settled from
sequence data — that needs the stock. Renaming it `Psimiovale_*` would be a
second taxonomic call on the same thin evidence, in the opposite direction.
`unverified` records what is known: the *fieldi* label is wrong, and the right
one is not established.

**How the pipeline treats it.** `unverified` is in `SPECIES_PAT`, so the record
survives into the panel FASTA and the curation decision stays visible in the
data. But it is **not a species**: every consumer that maps a header onto one of
the 11 panel species returns `None` for it —

| consumer | behaviour |
|---|---|
| `insilico_pcr.py::normalise_species` | `None` → recorded `species=NA`, not dropped; excluded from `*_target_species_coverage.tsv` and so from `panel_pass_fail.tsv`'s `n_target_species` |
| `resolution_margin.py` | label `unverified` ∉ `species_targets`; held **out of the BLAST database** by default (`--unverified-in-db exclude`), matching the "cleaned DB" of the corrected analysis |
| `dbs_diagnostic_coverage.py` | excluded — it carries the simiovale allele at every site, so including it would corrupt exactly the margin being measured |
| R `canon()` (`cross_species_filter.R`, `pair_summary.R`, `cross_marker_resolution.R`) | yields `unverified`, which is not in `config.species_targets` and so drops out of every pair grid |

**Consequence to keep in view:** *P. fieldi* now has **one** MIT reference
(`Pfieldi_AB354574_1`). Its leave-one-out call is undefined by construction, and
every fieldi pair in `resolution_table_margin_mit.tsv` is measured from the
other species' references only — reported as `resolved_one_side_only`, never as
plain `resolved`. Acquiring a second genuine *fieldi* mitogenome is the single
highest-value addition to this reference set.

`scripts/R/msp1_amplicon_specificity.R` is the one consumer that does **not**
see the relabel: it derives species from the *filename*, and `AB444066.1` still
sits in `fieldi_msp1.fasta`. Moving records between per-species files is a
structural change, deliberately not made here. Treat that script's `fieldi` rows
as containing one unverified record until it is.

### 2. `AY800109` dropped from the panel set as a duplicate

`Psimiovale_AY800109_1` and `Psimiovale_AB434920_1` are the same sequence: **0
mismatches over all 5809 mutually covered columns** after rotation (these are
fragments of a circular molecule deposited from different start points).
`AB434920` is kept — it is 89 bp longer and is the annotated type anchor.

Two copies of one haplotype double-weight it in every BLAST margin and overstate
within-simiovale tightness. Dropping it changes **no** diagnostic-site call at
cox1, cytb, cox3 or the whole genome — only *n*. Configured at
`config.reference_curation.mit.exclude_ids`; the record remains in
`mit_all.fasta`.

### Why the alignment was edited in place rather than re-derived

`outputs/alignment/ma_mit.target.fasta` was curated by row surgery — one header
renamed, one row deleted — not by re-running MAFFT. Verified afterwards:

- all **6747** alignment columns preserved; **0** columns became all-gap;
- every retained sequence **byte-identical** to before;
- `tpl_to_aln` **identical** with and without the duplicate, so **not one primer
  coordinate moved**;
- the consensus template changes at exactly **1 base of 5932** (template
  position 1409, a plurality tie-break), and that position falls inside **no**
  designed primer footprint in any of the three MIT pair tables.

So the Phase-5 primer tables remain valid against the curated alignment, which
is what makes the corrected resolution and DBS-coverage numbers comparable to
the ones they replace. Re-running MAFFT on the 47-record set is the cleaner
route and belongs with the panel down-select, when primer3 is re-run anyway.
See the re-run warning in `12_resolution_margin.README.md`.
