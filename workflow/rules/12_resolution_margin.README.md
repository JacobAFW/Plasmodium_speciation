# 12_resolution_margin — resolution as a margin, and the DBS diagnostic-window check

Two rules, one purpose: replace a resolution criterion that could not fail
anything with one that measures what the assay actually does, and re-answer the
`fieldi`/`simiovale` question on the corrected reference set.

| Rule | Question | Output |
|---|---|---|
| `resolution_margin` | for every species pair, by how much does own-species identity beat other-species identity? | `outputs/cross_species/resolution_table_margin_mit.tsv`, `loo_assignment_mit.tsv` |
| `dbs_diagnostic_coverage` | does each candidate amplicon still carry the cox1 and/or cytb sites that separate `fieldi` from `simiovale` — and does it still separate them? | `outputs/primer_design/dbs/dbs_diagnostic_coverage.tsv` |

Run them with:

```bash
source envs/activate.sh
snakemake -s workflow/Snakefile --configfile workflow/config.yaml \
  --cores 4 resolution_margin_all
```

Deliberately **not** in `rule all` — this is a correction pass over Phase 5's
existing primer tables, not something every run must redo.

> **Before you run it via snakemake, read the re-run warning at the bottom of
> this file.** The rules themselves are cheap; what is not cheap is the upstream
> cascade snakemake may decide to run first.

## Why the old criterion had to go

The project scored a pair resolvable if it carried **≥ 1 fixed diagnostic site**
inside the amplicon. On this panel that returns **55/55 under every amplicon
option** — Tier A, DBS single, DBS tiled, and the worst-case candidate interval.
A criterion that never fails and never discriminates between options cannot
inform a decision.

It was also wrong where it mattered. `fieldi vs simiovale` had 5–7 diagnostic
sites and was recorded `TRUE` in both `resolution_table.tsv` and
`mit_amplicon_resolution_tradeoff.tsv` — while a real best-hit read-out called
**both** *P. fieldi* references *P. simiovale*. A handful of fixed sites
scattered across ~5 kb does not outvote thousands of shared positions.

## The two numbers, and what "resolved" is allowed to mean

**`loo_margin`** (per reference) — best own-species %id minus best
other-species %id, over the whole panel. Negative means that isolate would be
called as another species.

**`pair_margin_pct`** (per species pair) — the *smallest* own-minus-other margin
over every reference of both species, considering only those two. Negative or
zero means no best-hit rule can separate the pair.

Neither is a threshold. Both are reported to four decimals and the sign is the
verdict.

Verdicts are deliberately four-valued, not two:

| Verdict | Meaning |
|---|---|
| `resolved` | margin > 0, measured from **both** species' references |
| `resolved_one_side_only` | margin > 0, but one species has a single reference so it contributes no margin of its own |
| `CONFUSABLE` | margin ≤ 0 — at least one reference is pulled across the boundary |
| `untestable_sole_reference` | both species are single-reference; nothing can be measured |

`resolved_one_side_only` is not pedantry. After the 2026-08-28 curation
*P. fieldi* has **one** MIT reference, so every fieldi pair is measured from the
other species' references alone. "No simiovale reference is pulled to fieldi" is
a real result; it is silent on whether the single fieldi reference is pulled to
simiovale. Flattening that to `resolved` would repeat exactly the optimism this
metric exists to remove.

A leave-one-out call is a **conservative lower bound** for a species with few
references — removing one of two removes half that species' representation. It
is the right test for "can the panel self-validate", and the `simium`/`vivax`
result does not depend on it: that follows directly from between-species
identity exceeding within-species identity, which is reference-count
independent.

## Result on the corrected reference set (2026-08-28)

55 pairs × 3 MIT amplicon options. Identical verdict counts under all three:

| Verdict | n |
|---|---|
| `resolved` | 35 |
| `resolved_one_side_only` | 18 |
| `CONFUSABLE` | 1 |
| `untestable_sole_reference` | 1 |

| pair | Tier A | DBS single | DBS tiled | verdict |
|---|---|---|---|---|
| `fieldi` vs `simiovale` | **+0.889** | **+1.140** | **+0.874** | `resolved_one_side_only` — **was `CONFUSABLE` (−0.925 / −1.176 / −0.925)** |
| `cynomolgi` vs `fieldi` | +0.964 | +0.930 | +0.927 | `resolved_one_side_only` — was `CONFUSABLE` under DBS single (−0.281) |
| `simium` vs `vivax` | −0.037 | 0.000 | −0.035 | `CONFUSABLE` under every option — **confirmed genuinely unresolvable** |
| `falciparum` vs `fieldi` | — | — | — | `untestable_sole_reference` (both species n = 1) |

`fieldi`/`simiovale` flipping sign is the whole point of the reference fix. The
+0.889 % Tier A margin independently reproduces the corrected analysis's
0.867 % minimum between-species mitochondrial divergence.

`falciparum` vs `fieldi` being untestable is an artefact of the metric, not a
risk: the two sit ~11 % apart, so a real read cannot plausibly cross. It is
reported as untestable rather than resolved because the panel cannot
*self-validate* it — both species have one MIT reference.

Leave-one-out accuracy on the curated set, over 46 references of which 2 are
sole-reference and therefore not evaluable: **43/44 Tier A (97.7 %)**, **43/44
DBS single**, **41/44 DBS tiled**. Every remaining misassignment is at the
`simium`/`vivax` boundary — the `fieldi` crossings are gone.

## The DBS constraint that changed

The original simiovale/fieldi investigation concluded that all discriminating
power for the pair sat in **cytb** (`AB434920` 4780–5928), and that became the
constraint on the DBS ≤ 4 kb tier: reach cytb or lose the pair.

**That premise was an artefact.** The cox1 locus set had pooled cytb-only
fragments with cox1-only deposits, leaving **8 mutually covered columns out of
6676** — a scan that could only ever return zero. Rebuilt properly, cox1 carries
**7** fixed sites (3378–4394) on *better* sampling (7 simiovale / 18 fieldi,
across 3 and 5 independent depositing groups) than cytb's **5** (4881–5835, 6/7).

So the constraint is **cox1 OR cytb**, and it is materially looser. The windows
and their sites live in `config.diagnostic_windows`, sourced from
`simiovale_fieldi_report_corrected.md` §6a/6b.

### Result: the DBS ≤ 4 kb tier is not the binding constraint

**All 60 filtered DBS-single pairs retain fieldi/simiovale separation**, margins
**+0.749 % to +1.157 %** (median +1.084 %). Same for all 60 tiled pairs
(+0.759 % to +1.095 %) and all 40 Tier A pairs (+0.877 % to +0.898 %).

| tier | complete cox1 only | complete cox1 + cytb | partial (no complete window) | all separate the pair |
|---|---|---|---|---|
| MIT_DBS_single | **33** | 20 | 7 | yes (60/60) |
| MIT_DBS_tiled | 2 | 40 | 18 | yes (60/60) |
| MIT_TierA | 20 | 20 | 0 | yes (40/40) |

**The make-or-break case.** `mit_dbs_0021` — the lowest-penalty DBS-single pair,
the one you would actually order — covers the **complete cox1 window (7/7
sites)** and only **2 of 5 cytb sites**. Under the old "must reach cytb"
constraint it would have been rejected. Under the corrected constraint it
separates the pair with a **+1.140 % margin — the best of any option, better
than Tier A's +0.889 %.**

The 7 DBS-single pairs with no complete window all sit at the 5′ end
(`AB434920` ~9–3478) and catch only site 3378. They still separate the pair
(+0.749 % to +0.774 %), because at ~0.87 % whole-molecule divergence the margin
is driven by overall divergence, not by site count. Site coverage and margin
agree in direction here; where they ever disagree, **the margin is what a real
read-out does.**

## Method

- **Amplicon interval** = the lowest-`pair_penalty` filtered pair per option
  (the pair you would order); tiled = the union of the best pair per tile, with
  overlapping tiles **merged before extraction** so the shared stretch is not
  double-weighted in every identity computed.
- **Spanning test**: a reference must carry real sequence at *every* primer
  footprint (± 30 alignment columns). Non-spanning references stay in the
  database — a real read can still hit them — but their own call is marked.
- **Coordinates**: primer tables are in ungapped consensus-template
  coordinates; the template and its bijective alignment↔template map are rebuilt
  with `run_primer3.build_consensus_template`, the same construction the design
  used, so the projection is drift-free. Verified: `mit_dbs_0021` → alignment
  [2794, 6163], reproducing the speciation-power audit exactly.
- **BLAST**: `blastn -task blastn -dust no -evalue 1e-5 -max_target_seqs 500`,
  self-hits discarded, best hit by bitscore. A hit must cover ≥ 80 % of the
  query (union of that pair's HSP query intervals, so a stack of short spurious
  local matches cannot manufacture coverage).
- **Species labels** come from the self-describing MIT id prefix
  (`Psimiovale_AB434920_1` → `simiovale`), the same prefix-strip route the R
  scripts use. **MIT only** — 18S headers are descriptive and need the accession
  map; 18S is untouched by the reference curation, and its resolution is already
  characterised in `outputs/primer_design/speciation_power_report.md` (6/13
  leave-one-out, driven by A-type/S-type paralog divergence exceeding
  between-species divergence — 18S is a genus-level confirmatory marker here,
  not a speciation marker).

## Relationship to `resolution_table.tsv`

`outputs/cross_species/resolution_table.tsv` is the older, threshold-based
product (`cross_marker_resolution.R`, driven by the `pident ≥ 99` **and**
`cov_shorter ≥ 0.99` strict-suspicious gate). It is **not** regenerated here and
it still records `fieldi vs simiovale` as `mit_resolved = TRUE`, which is
optimistic — `AB444133`'s cross-species hit cleared the identity gate at
99.848 % but failed the coverage gate at 88.55 %, so no row ever reached
`mit_suspicious.tsv` and the pair was scored `TRUE` by default.

`resolution_table_margin_mit.tsv` supersedes it for MIT. The two should be
reconciled when the R resolution chain is next re-run; until then, **the margin
table is the one to quote.** Relaxing that coverage gate, or adding a plain
nearest-neighbour check, is the standing fix — a pair whose members' closest
match is another species should be flagged regardless of HSP coverage.

## Re-run warning

Both rules read `outputs/alignment/ma_mit.target.fasta`. That file was curated
in place on 2026-08-28 (one row renamed, one duplicate row deleted; **all 6747
alignment columns preserved**, every retained sequence byte-identical), rather
than re-derived with MAFFT, precisely so that the existing primer tables stay
valid — `tpl_to_aln` is identical with and without the duplicate, so not one
primer coordinate moved, and the single consensus-template base that changes
(template position 1409) is inside no designed primer footprint.

Because the alignment is newer than the entropy outputs, a `snakemake` run
re-runs `entropy` and then `primer3_mit` / `primer3_dbs_mit_*`, regenerating the
primer tables. **That cascade was run on 2026-08-28 and verified to change
nothing:**

- all 160 MIT pair ids, amplicon lengths and template coordinates identical
  (`dbs_diagnostic_coverage.tsv` is content-identical across every pair;
  `mit_dbs_0021` still f_start 2032 / r_end 5355 / 3324 bp);
- `resolution_table_margin_mit.tsv` **byte-identical**;
- `panel_pass_fail.tsv` **byte-identical** after re-deriving the target half
  from the regenerated tables.

That is what the coordinate invariance predicts: the one consensus-template base
that changes (position 1409) lies inside no primer footprint, so primer3 sees an
identical design problem. The primer tables on disk are now genuinely derived
from the curated alignment, so the pipeline is internally consistent.

**Still stale by timestamp only:** `offtarget_amplification.tsv` and the
per-genome amplicon tables under `outputs/specificity/offtarget/`. `snakemake
specificity_all` will want to redo them (~16 min, plus a ~1.5 GB blastdb) because
the primer tables are newer. The primer *sequences* are provably unchanged, so
that run would reproduce the same numbers — it was deliberately not run here.
Check with `-n` before any full run.
