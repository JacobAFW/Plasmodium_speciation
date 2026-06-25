# 05_report.smk — Phase 5: render the Quarto book.
#
# The book reads every TSV / FASTA / figure in `outputs/` and `reports/figures/`
# at render time, so the rule's `input:` lists every upstream target. That way
# Snakemake re-renders when any output changes, and a stale render against an
# updated `outputs/` becomes impossible.

OUT     = config["paths"]["outputs"]
REPORTS = config["paths"]["reports"]
MARKERS = config["markers"]


def _report_inputs():
    items = [
        # Phase 1
        *[f"{OUT}/blast/self_blast_{m}.tsv" for m in MARKERS],
        # Phase 2
        f"{OUT}/cross_species/18S_acc_to_header.tsv",
        *[f"{OUT}/cross_species/{m}_suspicious.tsv"      for m in MARKERS],
        *[f"{OUT}/cross_species/{m}_closest_per_seq.tsv" for m in MARKERS],
        *[f"{OUT}/qc/{m}_lengths.tsv"                    for m in MARKERS],
        *[f"{OUT}/qc/{m}_length_classification.tsv"      for m in MARKERS],
        # Phase 3
        *[f"{OUT}/alignment/ma_{m}.target.fasta"         for m in MARKERS],
        *[f"{OUT}/entropy/{m}.{kind}.tsv"
          for m in MARKERS
          for kind in ("per_position", "windows", "primer_candidates", "core_bounds")],
        # Figures
        *[f"{REPORTS}/figures/{m}_{plot}.{ext}"
          for m in MARKERS
          for plot in ("entropy_vs_position", "coverage_vs_position")
          for ext in ("png", "svg")],
        # Phase 4
        *[f"{OUT}/cross_species/{m}_pair_summary.tsv" for m in MARKERS],
        f"{OUT}/cross_species/18S_suspicious_loose.tsv",
        f"{OUT}/cross_species/18S_species_coverage.tsv",
        f"{OUT}/cross_species/resolution_table.tsv",
    ]
    return items


REPORT_QMDS = [
    f"{REPORTS}/_quarto.yml",
    f"{REPORTS}/index.qmd",
    f"{REPORTS}/01_overview.qmd",
    f"{REPORTS}/02_cross_species_mit.qmd",
    f"{REPORTS}/03_cross_species_18S.qmd",
    f"{REPORTS}/04_alignment_entropy.qmd",
    f"{REPORTS}/05_resolution.qmd",
    f"{REPORTS}/06_published_panel_comparison.qmd",
    f"{REPORTS}/07_methods.qmd",
]


rule render_report:
    input:
        qmds   = REPORT_QMDS,
        upstream = _report_inputs(),
        docx     = "jiy519_suppl_supplementary-table1-5 (2).docx",
    output:
        index = f"{REPORTS}/_book/index.html",
    log:
        "logs/05_report/render_report.log",
    message:
        "[05_report] quarto render reports/"
    shell:
        r"""
        mkdir -p $(dirname {log})
        quarto render reports/ > {log} 2>&1
        """
