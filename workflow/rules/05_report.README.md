# 05_report — Quarto book render

One rule:

- `render_report` — runs `quarto render reports/` after every upstream output is in place. Output: `reports/_book/index.html` plus the chapter HTMLs and figures Quarto bundles.

## Inputs

The rule lists every TSV / FASTA / figure produced by Phases 1-4 as an explicit input. Snakemake therefore re-renders whenever any analytical output changes; a stale render against an updated `outputs/` is impossible by construction.

The `jiy519_suppl_supplementary-table1-5 (2).docx` at the project root is also an input so that Chapter 6 can re-parse it at render time.

## Chapters

```
reports/
├── _quarto.yml
├── index.qmd                          # preface, project structure, caveats
├── 01_overview.qmd                    # reference panels, target species, panel sizes
├── 02_cross_species_mit.qmd           # MIT suspicious + closest + per-pair summary
├── 03_cross_species_18S.qmd           # strict empty-set, looser re-pass, coverage
├── 04_alignment_entropy.qmd           # entropy + coverage panels, primer-candidate windows
├── 05_resolution.qmd                  # cross-marker resolution table (v1 deliverable)
├── 06_published_panel_comparison.qmd  # docx parsing + position overlap with our windows
└── 07_methods.qmd                     # command-by-command provenance
```

## Conventions used in chapters

- **Headline numbers come from `outputs/`, not from prior knowledge.** Every count, percentage, threshold appears in a code chunk that reads a TSV; prose then references the chunk output.
- R chunks are hidden by default (`echo: false` in `_quarto.yml`). Unfold them only if you need column types or filter thresholds.
- Each chapter ends with a `## Methods` section linking to its corresponding `workflow/rules/0X_*.README.md`.
- Figures live in `reports/figures/` and are produced by `03_alignment.smk::plot_entropy`. Both `.png` (for embedding) and `.svg` (for handoff editing) are saved.

## Running the render manually

```bash
source envs/activate.sh
quarto render reports/
# or via Snakemake:
snakemake --snakefile workflow/Snakefile reports/_book/index.html --cores 4
```

## Validation

`reports/_book/index.html` opens cleanly in a browser, every chapter renders, every figure resolves, every `outputs/...` reference is live. Chapter 6 surfaces concrete primer sequences from the docx — not boilerplate.
