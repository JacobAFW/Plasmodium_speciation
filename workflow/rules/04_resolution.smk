# 04_resolution.smk — Phase 4 (extension): per-pair summaries, looser-
# threshold 18S re-pass, 18S species coverage, cross-marker resolution.
#
# Provenance: NEW (HANDOFF.md Step 4 v1 acceptance criteria 2-5).

OUT  = config["paths"]["outputs"]
SUSP = config["suspicious"]


def _pair_summary_inputs(wc):
    base = {"blast": f"{OUT}/blast/self_blast_{wc.marker}.tsv"}
    if wc.marker == "18S":
        base["acc_map"] = f"{OUT}/cross_species/18S_acc_to_header.tsv"
    return base


rule pair_summary:
    input:
        unpack(_pair_summary_inputs),
    output:
        tsv = f"{OUT}/cross_species/{{marker}}_pair_summary.tsv",
    log:
        "logs/04_resolution/pair_summary_{marker}.log",
    message:
        "[04_resolution] pair summary {wildcards.marker}"
    params:
        pident_min       = SUSP["pident_min"],
        cov_shorter_min  = SUSP["cov_shorter_min"],
        shorter_len_min  = SUSP["shorter_len_min"],
        top_n_hsps       = SUSP["top_n_hsps"],
        acc_map_arg = lambda wc, input: (
            f"--acc-map {input.acc_map}" if wc.marker == "18S" else ""
        ),
    shell:
        r"""
        mkdir -p $(dirname {log}) $(dirname {output.tsv})
        Rscript scripts/R/pair_summary.R \
            --marker {wildcards.marker} \
            --blast {input.blast} \
            {params.acc_map_arg} \
            --pident-min {params.pident_min} \
            --cov-shorter-min {params.cov_shorter_min} \
            --shorter-len-min {params.shorter_len_min} \
            --top-n-hsps {params.top_n_hsps} \
            --out {output.tsv} \
            > {log} 2>&1
        """


# Looser-threshold 18S re-pass. Same input as the strict 18S filter; only
# the thresholds change (pident ≥ 98 AND mincov ≥ 0.85). The strict pipeline
# stays untouched.
LADDER_HC = config["ladder"]["high_conf"]


rule cross_species_filter_18S_loose:
    input:
        blast   = f"{OUT}/blast/self_blast_18S.tsv",
        acc_map = f"{OUT}/cross_species/18S_acc_to_header.tsv",
    output:
        suspicious = f"{OUT}/cross_species/18S_suspicious_loose.tsv",
        closest    = f"{OUT}/cross_species/18S_closest_per_seq_loose.tsv",
    log:
        "logs/04_resolution/cross_species_filter_18S_loose.log",
    message:
        "[04_resolution] 18S looser-threshold re-pass (pident≥98 AND mincov≥0.85)"
    params:
        pident_min       = LADDER_HC["pident_min"],
        cov_min          = LADDER_HC["mincov_min"],
        shorter_len_min  = SUSP["shorter_len_min"],
        top_n_hsps       = SUSP["top_n_hsps"],
    shell:
        r"""
        mkdir -p $(dirname {log}) $(dirname {output.suspicious})
        Rscript scripts/R/cross_species_filter.R \
            --marker 18S \
            --blast {input.blast} \
            --acc-map {input.acc_map} \
            --pident-min {params.pident_min} \
            --cov-min {params.cov_min} \
            --cov-metric mincov \
            --shorter-len-min {params.shorter_len_min} \
            --top-n-hsps {params.top_n_hsps} \
            --suspicious-out {output.suspicious} \
            --closest-out    {output.closest} \
            > {log} 2>&1
        """


rule species_coverage_18S:
    input:
        acc_map = f"{OUT}/cross_species/18S_acc_to_header.tsv",
    output:
        tsv = f"{OUT}/cross_species/18S_species_coverage.tsv",
    log:
        "logs/04_resolution/species_coverage_18S.log",
    message:
        "[04_resolution] 18S species coverage check"
    params:
        targets = ",".join(config["species_targets"]),
    shell:
        r"""
        mkdir -p $(dirname {log}) $(dirname {output.tsv})
        Rscript scripts/R/species_coverage_18S.R \
            --acc-map {input.acc_map} \
            --species-targets {params.targets} \
            --out {output.tsv} \
            > {log} 2>&1
        """


rule cross_marker_resolution:
    input:
        mit_suspicious        = f"{OUT}/cross_species/mit_suspicious.tsv",
        s18_suspicious        = f"{OUT}/cross_species/18S_suspicious.tsv",
        s18_suspicious_loose  = f"{OUT}/cross_species/18S_suspicious_loose.tsv",
        mit_lengths           = f"{OUT}/qc/mit_lengths.tsv",
        s18_coverage          = f"{OUT}/cross_species/18S_species_coverage.tsv",
    output:
        tsv = f"{OUT}/cross_species/resolution_table.tsv",
    log:
        "logs/04_resolution/cross_marker_resolution.log",
    message:
        "[04_resolution] cross-marker resolution table (v1 deliverable)"
    params:
        targets = ",".join(config["species_targets"]),
    shell:
        r"""
        mkdir -p $(dirname {log}) $(dirname {output.tsv})
        Rscript scripts/R/cross_marker_resolution.R \
            --species-targets {params.targets} \
            --mit-suspicious  {input.mit_suspicious} \
            --18s-suspicious  {input.s18_suspicious} \
            --18s-suspicious-loose {input.s18_suspicious_loose} \
            --mit-lengths     {input.mit_lengths} \
            --18s-coverage    {input.s18_coverage} \
            --out             {output.tsv} \
            > {log} 2>&1
        """
