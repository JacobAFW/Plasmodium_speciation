# 11_specificity — Phase 5 Step 5: in-silico PCR specificity

Simulates PCR for every filtered primer pair and records what amplifies.
Step 5 has two halves; **both are implemented.**

| Half | Question | Subject | Status |
|---|---|---|---|
| (a) target amplification | does the pair amplify every target species? | `data/reference/*.target.fasta` | **done** — `target_amplification*` rules |
| (b) off-target amplification | does it amplify anything it shouldn't? | `data/staged/*.fna` per `MANIFEST.tsv` | **done** — `offtarget_amplification*` rules, 12 of 13 genomes staged |

Half (b) is a rule, not a new engine: `scripts/python/insilico_pcr.py` already
takes an arbitrary subject FASTA and a `--subject-kind` label. All four primer
sets — Tier A MIT (40), DBS single (60), DBS tiled (60), 18S (20) = 180 pairs —
go against every staged genome in **one pass per genome**: the pairs are cheap
and a 3 Gbp subject is not, so splitting by marker would pay for the blastn
search twice and buy nothing.

*P. adleri* is **not staged** — `taxon:Plasmodium adleri` is not a valid
GenBank taxon name, so there is no assembly to fetch. It is treated as
optional; the screen runs on the 12 that are there and every table records
which genomes were tested (`n_offtarget_tested`).

## The engine is blastn, not MFEprimer

`HANDOFF.md` Step 5 and `CONTEXT-software.md` both say to install
`bioconda::mfeprimer`. **That package does not exist.** A `micromamba
repoquery search -c bioconda mfeprimer` returns `No entries matching` for
osx-arm64, osx-64, noarch *and* linux-64 (checked 2026-08-17), so the
instruction was never runnable — and the corresponding line in
`envs/install.sh` would have aborted the script under `set -e`.

Jacob chose (2026-08-17) to build Step 5 on the already-installed
blastn 2.16.0 rather than fetch an MFEprimer release binary from GitHub. The
substitution is confined to *how a binding site is scored*; the table shapes
and the pass/fail threshold in `HANDOFF.md` are unchanged.

## Method

1. **Locate.** Every primer of every pair is blastn'd (`-task blastn-short`,
   `-word_size 7`, `-dust no`, `-evalue 1000`) against the subject. blast is a
   fast, sensitive *locator only* — none of its scores reach the verdict.
2. **Re-derive exactly.** Each hit's subject coordinates are extended to the
   primer's full footprint, that window is pulled from the subject,
   reverse-complemented if the primer binds the minus strand, and compared
   base-by-base with the primer. This is what produces the mismatch counts,
   including at the 3' terminus. Gapped binding is not modelled: an indel
   inside a ~20-mer duplex does not prime in practice.
3. **Pair into amplicons.** A plus-strand binder upstream of a minus-strand
   binder, product length measured 5'-end to 5'-end. Both role assignments are
   tried, so a reference stored in either orientation is found.
4. **Score.** See below.

## The efficiency model (substitute for MFEprimer's `max_efficiency`)

Per primer:

- A mismatch anywhere in the **3'-terminal `three_prime_window` bases (5)**
  scores **0**. Taq cannot extend from an unpaired 3' end, so such a site does
  not prime however well the rest of the duplex anneals.
- Otherwise, primer3 computes the thermodynamics of the *actual* duplex
  (`calc_heterodimer`, so interior mismatches depress Tm rather than being
  counted heuristically) and the score is a logistic in that Tm:

  ```
  eff = 1 / (1 + exp(-(Tm_duplex - anneal_temp_c) / logistic_k))
  ```

Per amplicon, `efficiency = eff_forward × eff_reverse`.

`anneal_temp_c` defaults to **55 °C** — roughly 5 °C below the 58–62 °C Tm
window the pairs were designed to, i.e. a realistic annealing step. That makes
`eff = 0.5` mean *"this duplex melts at the annealing temperature"*, which is
why `HANDOFF.md`'s `max_efficiency > 0.5` threshold carries over unchanged and
still means something physical. Calibration reference: a perfectly-matched
21-mer pair scores ≈ 0.71; a pair with 1–2 interior mismatches on one primer
drops to ≈ 0.06–0.12.

**This is not MFEprimer's model and the numbers are not interchangeable with
published MFEprimer efficiencies.** It is a documented, monotone stand-in over
the same physical quantity. Treat the *ordering* and the threshold crossings as
the result, not the third decimal place.

## Word size: 7 on the reference set, 10 on the genomes

The locator's word size is the one parameter that differs between the halves,
and it is a scale problem, not a scientific one. Measured on E. coli with all
180 pairs: `word_size 7` returns **619,744** raw hits — a random 7-mer every
~16 kb — and scaling that to the 3.1 Gbp human genome gives ~4 × 10⁸ hit rows,
which exhausts memory long before it finishes.

`offtarget_word_size: 10` is not a guess:

- Every primer here is **21–24 bp**, so a site with **at most one mismatch
  anywhere** still contains an exact run of ≥ 10 and is still located.
- One mismatch is already far past the point of mattering. By the calibration
  above, a single interior mismatch drops pair efficiency to ≈ 0.1, so **any
  hit that can cross the 0.5 flag threshold needs both primers binding
  essentially perfectly** — a full-length exact match, which word size 10 finds
  trivially.

So the flag decision is untouched; what is given up is some zero-efficiency
junk in the `borderline` bucket. Measured on E. coli: word 7 → 3 amplicons,
word 10 → 1, and **both dropped amplicons scored efficiency 0.0**.

Two other changes were needed to make a genome-scale subject tractable, neither
of which touches a verdict: `insilico_pcr.py` now counts mismatches by string
comparison *before* paying for primer3's thermodynamic alignment (the sites
this skips were being discarded by `--max-mismatch` anyway — 5.4× faster on
E. coli, identical output), and `--db-dir` keeps the blastn database between
runs so a re-run does not rebuild a 3 Gbp index.

The off-target product cap is `offtarget_max_product_bp: 8000` — the larger of
the two markers' caps, because one pass covers both. Erring long is the
conservative direction for a screen that is looking *for* off-target products.

## Inputs

| Marker | Primer pairs | Subject |
|---|---|---|
| mit | `mit_pairs_filtered.tsv` (Tier A), `dbs/mit_dbs_single_pairs_filtered.tsv`, `dbs/mit_dbs_tiled_pairs_filtered.tsv` | `mit_all.target.fasta` (48 seqs, 11 species) |
| 18S | `18S_pairs_filtered.tsv` | `18S_ref_db.target.fasta` (17 seqs, 10 species) |
| off-target (all 180 pairs at once) | all four tables above | every `staged = TRUE` row of `data/staged/MANIFEST.tsv` (12 genomes) |

Reference → species uses the canonical `case_when` ladder from
`CONTEXT-data.md`, applied to the **FASTA header** rather than a pre-parsed
genus token — that is what makes `Plasmodium cf. malariae type 2`
(`AF488000.1`) resolve to `malariae` instead of the `p.cf` a naive genus regex
yields. `simiovale` is tested before `ovale` so it cannot be swallowed — as of
2026-08-28 that is also the order `CONTEXT-data.md` documents, so this is no
longer a local deviation. Headers that map to no panel species return `None`
and are recorded as `species=NA` rather than dropped; that is how
`Punverified_ATCC-30157_AB444133` is kept visible in the evidence tables
without ever being scored as a species.

## Outputs

`outputs/specificity/`

| File | Grain | Notes |
|---|---|---|
| `{marker}_target_amplicons.tsv` | one row per amplicon | the evidence: per-primer mismatches, 3'-end mismatches, duplex Tm/ΔG, per-primer efficiency |
| `{marker}_target_amplification.tsv` | one row per (pair, reference) | HANDOFF Step 5 shape. **Full grid — zero-amplicon cells are present, not dropped**, which is the whole point for a "must amplify" check |
| `{marker}_target_species_coverage.tsv` | one row per (pair, species) | `n_refs`, `n_amplified`, `frac_refs_amplified`, `max_efficiency` |
| `target_amplification.tsv`, `target_species_coverage.tsv` | merged across markers | what downstream reads |
| `offtarget/{genome}.amplicons.tsv` | one row per off-target amplicon | the evidence, per staged genome |
| `offtarget_amplification.tsv` | one row per (pair, off_target_species, subject) | HANDOFF Step 5's off-target shape, zero-filled once per (pair, genome). `flag = amplicon_count > 0 AND max_efficiency > 0.5`. `flag` is evidence, **not** a verdict — only the non-*Plasmodium* subset of flags fails a pair |
| `panel_pass_fail.tsv` | one row per pair | the verdict — see below |

The off-target rollup is **not** `insilico_pcr.py --summary-out`. That is a
full (pair × subject *sequence*) grid, which is right for the target half
("did every target amplify?") and wrong here: a staged genome is thousands of
contigs, so the grid would be 180 pairs × 8,090 A. gambiae contigs = 1.5M rows
of "no, this pair does not amplify contig 7,412 either". The off-target
question is asked per **genome**, and zero is the expected answer, so
`scripts/python/offtarget_summary.py` zero-fills once per (pair, genome). The
off-target species label comes from `MANIFEST.tsv`, not from the panel-species
ladder — the ladder exists to map a header onto one of the 11 *panel* species
and by construction returns `None` for a host, a vector or a non-panel
*Plasmodium*; any subject sequence that *does* hit the ladder is reported as a
warning rather than silently labelled.

## The verdict (`panel_pass_fail.tsv`)

`scripts/python/panel_verdict.py`, one row per pair:

### Off-target means NON-Plasmodium (corrected 2026-08-28)

**The panel is genus-conserved by design.** The whole concept is: amplify any
*Plasmodium* with one conserved pair, then speciate by *reading the amplicon* —
not by which primer fired. A pair that amplifies *P. berghei* or
*P. gallinaceum* is therefore doing exactly what it was designed to do.

Until 2026-08-28 this rubric scored that as a specificity failure, which
inverts the design intent. It marked **34 pairs `fail_offtarget` whose only
flagged off-targets were congeners** — every one of them clean on host, vector
and contaminant. Under the corrected rubric those 34 read as they should: 27
`pass`, 7 `borderline`.

The failure that actually matters is priming on something that is **not
*Plasmodium*** and **is** in the tube: host (human, macaque), vector
(*Anopheles*), contaminant (*E. coli*, *M. tuberculosis*). That is the pass/fail
axis. Congener amplification is still measured and reported in full —
`congener_amplification`, `n_congener_flagged`, `congener_flagged_species` — it
just annotates instead of failing.

Which manifest `purpose` counts as a congener is
`config.specificity.congener_purpose`, not a literal, so a manifest that grows a
new category cannot silently move a verdict.

| Verdict | Meaning |
|---|---|
| `pass` | amplifies every target species AND no **non-*Plasmodium*** off-target above the threshold |
| `fail_no_target` | misses ≥ 1 target species |
| `fail_offtarget` | ≥ 1 **non-*Plasmodium*** genome (host / vector / contaminant) with `amplicon_count > 0 AND max_efficiency > 0.5` |
| `borderline` | a **non-*Plasmodium*** product exists but stays below the threshold — review, do not auto-reject |

Precedence is `fail_no_target` > `fail_offtarget` > `borderline` > `pass`, and
it is only a display order: `target_ok`, `n_offtarget_flagged`,
`n_offtarget_flagged_nonplasmodium`, `n_congener_flagged` and `n_offtarget_any`
are all kept as their own columns, so nothing a verdict outranks is hidden.

**Result as of 2026-08-28** (180 pairs): `fail_offtarget` **0**, `pass` 27
(21 Tier A MIT + 6 DBS single), `borderline` 7, `fail_no_target` 146.
`host_cross_amplification` is **0 pairs** — no pair puts an above-threshold
product on human or macaque. The 7 borderline pairs do produce on-size host /
*E. coli* products, but every one scores efficiency **0.0**: a duplex that
cannot prime. They are flagged for review rather than passed silently, and they
are the pairs to scrutinise first at down-select.

Note that `fail_no_target` is a separate axis and the rubric change does not
touch it — the 146 pairs that miss a target species still miss it.

Two judgement calls worth knowing about:

- **"Amplifies a target"** has two readings and both are reported, because
  they disagree by a lot. `target_ok` (**strict**, drives the verdict) needs
  ≥ 1 on-size amplicon *and* `max_efficiency > 0.5` for every target species —
  this is the definition the Session-14 target-half table was built on (Tier A
  26/40, DBS single 8/60, DBS tiled 0/60, 18S 0/20), so using anything else
  would silently restate that result. `target_ok_amplicon_only` is HANDOFF's
  literal wording, ≥ 1 on-size amplicon, efficiency ignored (40/55/59/8). The
  gap between them *is* the finding: an on-size product made by a duplex that
  melts below the annealing temperature has not amplified anything in a real
  tube. The species involved are carried as `n_species_low_eff` /
  `low_eff_species`.
- **The expected species set is per marker**, taken from what the pair was
  actually tested against rather than from `config["species_targets"]`. The 18S
  reference DB has no simiovale, so scoring 18S pairs out of 11 would fail
  every one of them for a gap in the reference set rather than a defect in the
  primer.

`host_cross_amplification` is broken out as its own column because it is the
failure that matters most: human and macaque DNA dominate a real sample, so a
pair that primes on `homo_sapiens` or `macaca_fascicularis` is the one result
that would sink a pair outright. The flagged off-targets are also split by
manifest `purpose` (`n_offtarget_flagged_nonplasmodium`,
`n_congener_flagged`, `offtarget_flagged_purposes`) so the down-selection can be
made on either reading without re-deriving it.

## Parameters

All in `workflow/config.yaml` under `specificity:` — `anneal_temp_c`,
`logistic_k`, `three_prime_window`, `max_mismatch`, `size_tolerance`,
`blast_word_size`, `efficiency_threshold`, `max_product_bp` (per marker),
`offtarget_word_size`, `offtarget_max_product_bp`.

## How to run

```bash
source envs/activate.sh
# both halves
snakemake -s workflow/Snakefile --configfile workflow/config.yaml \
  --cores 10 specificity_all
# target half only (no staged genomes needed)
snakemake -s workflow/Snakefile --configfile workflow/config.yaml \
  --cores 4 specificity_target_all
```

Deliberately **not** in `rule all`: the off-target half depends on hand-staged
genomes, and a default-rule run must not imply Step 5 is complete.

Cost, measured 2026-08-17 on a 12-core / 32 GB Mac: ~4 min wall for the ten
small genomes plus ~6 min each for the two ~3 Gbp hosts, and ~1.5 GB of
persistent blastn index under `outputs/specificity/offtarget/blastdb/`
(gitignored). Deleting that directory only costs the rebuild.

## How to read the result

A cell fails for one of three distinguishable reasons, and they mean different
things:

- **`no_binding`** — no site pair at all. The primer sites are absent or too
  divergent to locate. Hard miss.
- **`off_size_only`** — sites pair up but the product is outside
  ±`size_tolerance` of the designed length. Usually a spurious pairing.
- **`low_efficiency`** — the amplicon is there and on-size, but
  `max_efficiency ≤ 0.5`: the primers bind with mismatches that put the duplex
  Tm at or below the annealing temperature. **This is the interesting failure**
  — it is a design problem at a specific species, not a coordinate artefact,
  and it is invisible if you only count amplicons.

Specificity is not resolution. A pair can amplify every target species cleanly
and still be unable to *distinguish* two of them — the vivax/simium MIT
ambiguity is the standing example. Read this table next to
`outputs/cross_species/resolution_table.tsv`, never instead of it.
