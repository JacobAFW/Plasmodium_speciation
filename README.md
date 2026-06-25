# Plasmodium long-amplicon speciation pipeline

A reproducible pipeline for designing and validating a long-amplicon Nanopore
(MinION) speciation assay for *Plasmodium*, using mitochondrial (MIT) and 18S
rRNA markers — plus an *msp1* analysis — to distinguish closely related species.

## What this repo is

A Snakemake workflow with its R, Python, and shell scripts that builds BLAST
databases from MIT/18S/*msp1* references, quantifies cross-species ambiguity,
computes per-position alignment entropy for primer-window scoring, assembles a
combined cross-marker resolution table, and renders the whole thing as a Quarto
report. The stack is Snakemake + BLAST+ / seqkit / MAFFT + R (tidyverse) +
Python, with a pinned environment.

## What's included — and what's deliberately not

This repository contains **code only**. By design it does **not** include:

- raw or processed **data** (sequence, genotype, phenotype, tabular records) —
  here: FASTA references, BLAST databases, `data/`, `share/`, and all `.tsv`/`.csv`
- **sample sheets, manifests, or metadata** that link samples to individuals
- any **identifying or sensitive** information — here: the legacy HPC scripts
  (`scripts/legacy/`, `legacy/`) that carry hardcoded usernames and internal
  server paths, and the internal handoff notes (`CONTEXT-*.md`, `HANDOFF.md`,
  `INSTRUCTIONS.md`, `MEMORY.md`, `drafts/`)
- **generated outputs** (`outputs/`, `logs/`, `reports/_book/`,
  `reports/figures/`) — reproducible by running the pipeline
- the **vendored tool environment** (`vvg-box/`, `.snakemake/`,
  `envs/activate.sh`) — recreate it with `envs/install.sh`
- credentials, tokens, or environment files (none were found in the scan)
- third-party manuscript documents (`*.docx`, e.g. a published supplementary
  table)

Data lives outside version control (institutional storage / the references you
supply). The scripts expect it at the paths described in `workflow/config.yaml`.

## Reproducing the analysis

1. **Environment:** `bash envs/install.sh` (versions pinned in
   `envs/environment.lock.yaml`). On macOS, `brew install coreutils` first.
2. **Input data:** supply MIT, 18S, and *msp1* references at the locations named
   in `workflow/config.yaml` (the `paths:` block; `data/reference/` by default).
   `scripts/python/fetch_msp1_refs.py` can fetch *msp1* references from NCBI
   Entrez — set its `EMAIL` value to your own address first (NCBI's terms).
3. **Run:** `snakemake --cores 4`
4. **Report:** `quarto render reports/`

## Structure

```
workflow/
  Snakefile              # top-level workflow
  config.yaml            # parameters, thresholds, species targets, paths
  rules/                 # 01–08 .smk rule files (+ per-rule READMEs)
scripts/
  python/                # reference fetching, entropy, header parsing, sweeps
  R/                     # filtering, summaries, plots, trees, resolution tables
  sh/                    # seqkit / BLAST / alignment helpers
  legacy_replicate.sh    # convenience wrapper
reports/
  _quarto.yml            # Quarto book config
  *.qmd                  # per-stage report chapters (sources only)
envs/
  install.sh             # one-shot environment bootstrap
  environment.lock.yaml  # pinned tool versions
.lintr                   # R lint configuration
```

## License

<Choose a license — e.g. MIT — before publishing. Until a LICENSE file is added,
default "all rights reserved" applies.>
