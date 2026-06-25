# 02_cross_species.smk — Phase 2: BLAST tabular → cross-species suspicious +
# closest-per-sequence tables, plus the 18S header parser and per-marker
# length QC.
#
# Provenance: mit_similarity.Rmd lines 25-139 (MIT cross-species),
# 198-237 (18S header parser), 242-391 (18S cross-species),
# 148-179 (MIT length QC), 429-460 (18S length QC).

OUT = config["paths"]["outputs"]
SUSP = config["suspicious"]


rule parse_18S_headers:
    input:
        fasta = f"{config['paths']['data_ref']}/18S_ref_db.fasta",
    output:
        tsv = f"{OUT}/cross_species/18S_acc_to_header.tsv",
    log:
        "logs/02_cross_species/parse_18S_headers.log",
    message:
        "[02_cross_species] parse 18S FASTA headers → species map"
    shell:
        r"""
        mkdir -p $(dirname {log}) $(dirname {output.tsv})
        python scripts/python/parse_18S_headers.py \
            --fasta {input.fasta} --out {output.tsv} > {log} 2>&1
        """


# MIT and 18S take different inputs (the 18S filter needs the acc->species
# map) so we keep them as two rules sharing the same script. Snakemake's
# wildcard-driven function inputs let us parameterise the script call.
def _filter_inputs(wc):
    base = {"blast": f"{OUT}/blast/self_blast_{wc.marker}.tsv"}
    if wc.marker == "18S":
        base["acc_map"] = f"{OUT}/cross_species/18S_acc_to_header.tsv"
    return base


rule cross_species_filter:
    input:
        unpack(_filter_inputs),
    output:
        suspicious = f"{OUT}/cross_species/{{marker}}_suspicious.tsv",
        closest    = f"{OUT}/cross_species/{{marker}}_closest_per_seq.tsv",
    log:
        "logs/02_cross_species/cross_species_filter_{marker}.log",
    message:
        "[02_cross_species] cross-species filter {wildcards.marker}"
    params:
        pident_min       = SUSP["pident_min"],
        cov_min          = SUSP["cov_shorter_min"],
        shorter_len_min  = SUSP["shorter_len_min"],
        top_n_hsps       = SUSP["top_n_hsps"],
        # Build the optional --acc-map flag only for 18S.
        acc_map_arg = lambda wc, input: (
            f"--acc-map {input.acc_map}" if wc.marker == "18S" else ""
        ),
    shell:
        r"""
        mkdir -p $(dirname {log}) $(dirname {output.suspicious})
        Rscript scripts/R/cross_species_filter.R \
            --marker {wildcards.marker} \
            --blast {input.blast} \
            {params.acc_map_arg} \
            --pident-min {params.pident_min} \
            --cov-min {params.cov_min} \
            --cov-metric cov_shorter \
            --shorter-len-min {params.shorter_len_min} \
            --top-n-hsps {params.top_n_hsps} \
            --suspicious-out {output.suspicious} \
            --closest-out    {output.closest} \
            > {log} 2>&1
        """


# Length QC reads the target-subset FASTAs produced by Phase 3's seqkit grep.
# Snakemake will resolve the dependency the right way around because the
# subset rule outputs are declared in 03_alignment.smk.

TARGET_FASTA = {
    "mit": f"{config['paths']['data_ref']}/mit_all.target.fasta",
    "18S": f"{config['paths']['data_ref']}/18S_ref_db.target.fasta",
}


rule length_table:
    input:
        fasta = lambda w: TARGET_FASTA[w.marker],
    output:
        tsv = f"{OUT}/qc/{{marker}}_lengths.tsv",
    log:
        "logs/02_cross_species/length_table_{marker}.log",
    message:
        "[02_cross_species] seqkit fx2tab lengths {wildcards.marker}"
    shell:
        r"""
        mkdir -p $(dirname {log}) $(dirname {output.tsv})
        bash scripts/sh/seqkit_lengths.sh {input.fasta} {output.tsv} > {log} 2>&1
        """


rule length_classification:
    input:
        tsv = f"{OUT}/qc/{{marker}}_lengths.tsv",
    output:
        tsv = f"{OUT}/qc/{{marker}}_length_classification.tsv",
    log:
        "logs/02_cross_species/length_classification_{marker}.log",
    message:
        "[02_cross_species] length classification {wildcards.marker}"
    shell:
        r"""
        mkdir -p $(dirname {log}) $(dirname {output.tsv})
        Rscript scripts/R/length_classification.R \
            --in {input.tsv} --out {output.tsv} > {log} 2>&1
        """
