# 08_msp1.smk — MSP1 as a candidate third marker for the long-read panel.
#
# Mirrors the MIT and 18S rule shapes from 03_alignment.smk so MSP1 slots
# into the same alignment → entropy → primer-candidate → plot pipeline
# without modifying the existing rules. The underlying scripts
# (align.sh, seq_entropy.py, plot_entropy.R, subset_to_targets.sh) are the
# same — only the marker tag, the input reference FASTA, and the output
# file stems differ.
#
# Prerequisites:
#   1. data/reference/msp1/msp1_ref_db.fasta exists (output of
#      scripts/python/fetch_msp1_refs.py, run from a machine with NCBI
#      access — the sandbox proxy blocks eutils).
#   2. The existing pipeline (03_alignment, seq_entropy.py) is in place.
#
# Outputs (targets a downstream user would request):
#   data/reference/msp1/msp1_ref_db.target.fasta   — panel-species filter
#   outputs/alignment/ma_msp1.target.fasta         — MAFFT alignment
#   outputs/entropy/msp1.per_position.tsv          — Shannon entropy + coverage
#   outputs/entropy/msp1.windows.tsv               — sliding-window scores
#   outputs/entropy/msp1.primer_candidates.tsv     — candidate primer windows
#   outputs/entropy/msp1.core_bounds.tsv           — alignment core boundaries
#   reports/figures/msp1_entropy_vs_position.{png,svg}
#   reports/figures/msp1_coverage_vs_position.{png,svg}
#
# After this rule file runs, Chapter 5 (cross-marker resolution) and the new
# Chapter 7 (methodological comparison) can be extended to include MSP1 in
# the per-pair resolution table — i.e. answer "does MSP1 resolve pairs that
# MIT + 18S can't?".

DATA_REF = config["paths"]["data_ref"]
OUT      = config["paths"]["outputs"]
REPORTS  = config["paths"]["reports"]

# Same regex as 03_alignment.smk — matches the 11 target species on either
# encoded prefix (Pfalciparum_...) or descriptive header (P.falciparum ...).
SPECIES_PAT_MSP1 = "vivax|falciparum|knowlesi|malariae|ovale|coat|inui|fieldi|cyno|simiovale|simium"


rule subset_to_targets_msp1:
    input:
        fasta = f"{DATA_REF}/msp1/msp1_ref_db.fasta",
    output:
        fasta = f"{DATA_REF}/msp1/msp1_ref_db.target.fasta",
    log:
        "logs/08_msp1/subset_to_targets_msp1.log",
    message:
        "[08_msp1] seqkit grep msp1 (panel species)"
    params:
        pat = SPECIES_PAT_MSP1,
    shell:
        r"""
        mkdir -p $(dirname {log})
        # MSP1 FASTAs come out of fetch_msp1_refs.py with descriptive headers
        # ('>ACC Pspecies msp1 | ...') — same convention as 18S, so -n -i.
        bash scripts/sh/subset_to_targets.sh \
            {input.fasta} {output.fasta} '{params.pat}' -n -i > {log} 2>&1
        """


rule align_msp1:
    input:
        fasta = f"{DATA_REF}/msp1/msp1_ref_db.target.fasta",
    output:
        fasta = f"{OUT}/alignment/ma_msp1.target.fasta",
    log:
        "logs/08_msp1/align_msp1.log",
    message:
        "[08_msp1] mafft --auto msp1"
    shell:
        r"""
        mkdir -p $(dirname {log}) $(dirname {output.fasta})
        bash scripts/sh/align.sh {input.fasta} {output.fasta} {log}
        """


rule entropy_msp1:
    input:
        aln = f"{OUT}/alignment/ma_msp1.target.fasta",
    output:
        per_position      = f"{OUT}/entropy/msp1.per_position.tsv",
        windows           = f"{OUT}/entropy/msp1.windows.tsv",
        primer_candidates = f"{OUT}/entropy/msp1.primer_candidates.tsv",
        core_bounds       = f"{OUT}/entropy/msp1.core_bounds.tsv",
    log:
        "logs/08_msp1/entropy_msp1.log",
    message:
        "[08_msp1] seq_entropy msp1"
    params:
        prefix = f"{OUT}/entropy/msp1",
    shell:
        r"""
        mkdir -p $(dirname {log}) $(dirname {output.per_position})
        python scripts/python/seq_entropy.py \
            --aln {input.aln} --outprefix {params.prefix} > {log} 2>&1
        """


rule plot_entropy_msp1:
    input:
        per_position      = f"{OUT}/entropy/msp1.per_position.tsv",
        windows           = f"{OUT}/entropy/msp1.windows.tsv",
        primer_candidates = f"{OUT}/entropy/msp1.primer_candidates.tsv",
    output:
        entropy_png  = f"{REPORTS}/figures/msp1_entropy_vs_position.png",
        entropy_svg  = f"{REPORTS}/figures/msp1_entropy_vs_position.svg",
        coverage_png = f"{REPORTS}/figures/msp1_coverage_vs_position.png",
        coverage_svg = f"{REPORTS}/figures/msp1_coverage_vs_position.svg",
    log:
        "logs/08_msp1/plot_entropy_msp1.log",
    message:
        "[08_msp1] plot_entropy msp1"
    shell:
        r"""
        mkdir -p $(dirname {log}) $(dirname {output.entropy_png})
        Rscript scripts/R/plot_entropy.R \
            --marker msp1 \
            --per-position {input.per_position} \
            --windows      {input.windows} \
            --primers      {input.primer_candidates} \
            --entropy-png  {output.entropy_png} \
            --entropy-svg  {output.entropy_svg} \
            --coverage-png {output.coverage_png} \
            --coverage-svg {output.coverage_svg} \
            > {log} 2>&1
        """


# Convenience target — request this to materialise all MSP1 artefacts.
rule msp1_all:
    input:
        f"{DATA_REF}/msp1/msp1_ref_db.target.fasta",
        f"{OUT}/alignment/ma_msp1.target.fasta",
        f"{OUT}/entropy/msp1.per_position.tsv",
        f"{OUT}/entropy/msp1.windows.tsv",
        f"{OUT}/entropy/msp1.primer_candidates.tsv",
        f"{OUT}/entropy/msp1.core_bounds.tsv",
        f"{REPORTS}/figures/msp1_entropy_vs_position.png",
        f"{REPORTS}/figures/msp1_coverage_vs_position.png",
