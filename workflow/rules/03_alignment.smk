# 03_alignment.smk — Phase 3: target-species subset → MAFFT alignment →
# per-position entropy + windowed scoring + primer-candidate windows →
# entropy/coverage figures.
#
# Provenance: mit_similarity.Rmd lines 405-407 (seqkit grep), 422-424
# (mafft --auto), 464 (seq_entropy.py), and the legacy plot_entropy.R.

DATA_REF = config["paths"]["data_ref"]
OUT      = config["paths"]["outputs"]
REPORTS  = config["paths"]["reports"]

# Same regex the legacy chunks use — substring-matching against the panel.
SPECIES_PAT = "vivax|falciparum|knowlesi|malariae|ovale|coat|inui|fieldi|cyno|simiovale|simium"


rule subset_to_targets_mit:
    input:
        fasta = f"{DATA_REF}/mit_all.fasta",
    output:
        fasta = f"{DATA_REF}/mit_all.target.fasta",
    log:
        "logs/03_alignment/subset_to_targets_mit.log",
    message:
        "[03_alignment] seqkit grep mit (panel species)"
    params:
        pat = SPECIES_PAT,
    shell:
        r"""
        mkdir -p $(dirname {log})
        # MIT IDs are self-describing & case-stable; no -n -i.
        bash scripts/sh/subset_to_targets.sh \
            {input.fasta} {output.fasta} '{params.pat}' > {log} 2>&1
        """


rule subset_to_targets_18S:
    input:
        fasta = f"{DATA_REF}/18S_ref_db.fasta",
    output:
        fasta = f"{DATA_REF}/18S_ref_db.target.fasta",
    log:
        "logs/03_alignment/subset_to_targets_18S.log",
    message:
        "[03_alignment] seqkit grep 18S (panel species)"
    params:
        pat = SPECIES_PAT,
    shell:
        r"""
        mkdir -p $(dirname {log})
        # 18S has descriptive headers — match against the full header,
        # case-insensitive.
        bash scripts/sh/subset_to_targets.sh \
            {input.fasta} {output.fasta} '{params.pat}' -n -i > {log} 2>&1
        """


TARGET_FASTA = {
    "mit": f"{DATA_REF}/mit_all.target.fasta",
    "18S": f"{DATA_REF}/18S_ref_db.target.fasta",
}


rule align:
    input:
        fasta = lambda w: TARGET_FASTA[w.marker],
    output:
        fasta = f"{OUT}/alignment/ma_{{marker}}.target.fasta",
    log:
        "logs/03_alignment/align_{marker}.log",
    message:
        "[03_alignment] mafft --auto {wildcards.marker}"
    shell:
        r"""
        mkdir -p $(dirname {log}) $(dirname {output.fasta})
        bash scripts/sh/align.sh {input.fasta} {output.fasta} {log}
        """


# seq_entropy.py emits all five files from one --outprefix. We declare them
# all as outputs so Snakemake knows about each artefact.
#
# Two candidate-window outputs share the same upstream computation:
#   primer_candidates.tsv   — design side, end-anchored only.
#   interior_candidates.tsv — benchmark side, no end_region gate. Step 2
#                              of Phase 5 added this for the Imwong-forward
#                              site (sits at end_region=none).
rule entropy:
    input:
        aln = f"{OUT}/alignment/ma_{{marker}}.target.fasta",
    output:
        per_position        = f"{OUT}/entropy/{{marker}}.per_position.tsv",
        windows             = f"{OUT}/entropy/{{marker}}.windows.tsv",
        primer_candidates   = f"{OUT}/entropy/{{marker}}.primer_candidates.tsv",
        interior_candidates = f"{OUT}/entropy/{{marker}}.interior_candidates.tsv",
        core_bounds         = f"{OUT}/entropy/{{marker}}.core_bounds.tsv",
    log:
        "logs/03_alignment/entropy_{marker}.log",
    message:
        "[03_alignment] seq_entropy {wildcards.marker}"
    params:
        prefix = f"{OUT}/entropy/{{marker}}",
    shell:
        r"""
        mkdir -p $(dirname {log}) $(dirname {output.per_position})
        python scripts/python/seq_entropy.py \
            --aln {input.aln} --outprefix {params.prefix} > {log} 2>&1
        """


rule plot_entropy:
    input:
        per_position      = f"{OUT}/entropy/{{marker}}.per_position.tsv",
        windows           = f"{OUT}/entropy/{{marker}}.windows.tsv",
        primer_candidates = f"{OUT}/entropy/{{marker}}.primer_candidates.tsv",
    output:
        entropy_png  = f"{REPORTS}/figures/{{marker}}_entropy_vs_position.png",
        entropy_svg  = f"{REPORTS}/figures/{{marker}}_entropy_vs_position.svg",
        coverage_png = f"{REPORTS}/figures/{{marker}}_coverage_vs_position.png",
        coverage_svg = f"{REPORTS}/figures/{{marker}}_coverage_vs_position.svg",
    log:
        "logs/03_alignment/plot_entropy_{marker}.log",
    message:
        "[03_alignment] plot_entropy {wildcards.marker}"
    shell:
        r"""
        mkdir -p $(dirname {log}) $(dirname {output.entropy_png})
        Rscript scripts/R/plot_entropy.R \
            --marker {wildcards.marker} \
            --per-position {input.per_position} \
            --windows      {input.windows} \
            --primers      {input.primer_candidates} \
            --entropy-png  {output.entropy_png} \
            --entropy-svg  {output.entropy_svg} \
            --coverage-png {output.coverage_png} \
            --coverage-svg {output.coverage_svg} \
            > {log} 2>&1
        """
