# 11_specificity.smk — Phase 5 Step 5: in-silico PCR specificity.
#
# Step 5 has two halves, both implemented here:
#
#   (a) TARGET amplification — the panel must amplify every target species.
#       Reads only `data/reference/*.target.fasta`.
#   (b) OFF-TARGET amplification — the panel must NOT amplify host, vector,
#       contaminant or non-panel Plasmodium genomes. Reads the genomes marked
#       `staged = TRUE` in `data/staged/MANIFEST.tsv`, which Jacob stages by
#       hand (embargo: the pipeline never downloads sequence data). 12 of the
#       13 manifest rows are staged as of 2026-08-17; *P. adleri* is absent
#       (invalid GenBank taxon name) and is treated as optional, so the screen
#       runs on what is there and the tables record which genomes were tested.
#       Half (b) is a rule, not a new engine — `insilico_pcr.py` already takes
#       an arbitrary subject FASTA and a `--subject-kind` label.
#
# ENGINE — blastn, not MFEprimer. `HANDOFF.md` Step 5 and `CONTEXT-software.md`
# both specify `bioconda::mfeprimer`; that package does not exist on bioconda
# for osx-arm64, osx-64, noarch or linux-64. Jacob chose the blastn route on
# 2026-08-17 over fetching an MFEprimer release binary from GitHub. See the
# README for the efficiency model that substitutes for MFEprimer's
# `max_efficiency`, and why HANDOFF's 0.5 threshold survives the swap.
#
# Both primer tiers are tested against the same target references:
#   Tier A     — outputs/primer_design/{mit,18S}_pairs_filtered.tsv
#   DBS <=4 kb — outputs/primer_design/dbs/mit_dbs_{single,tiled}_pairs_filtered.tsv
# MIT pairs of both tiers go in one blastn pass against the MIT references;
# the 18S pairs are shared across tiers and go against the 18S references.

import csv
import os
import sys

OUT      = config["paths"]["outputs"]
DATA_REF = config["paths"]["data_ref"]
SPEC     = config["specificity"]

_SPEC_DIR = f"{OUT}/specificity"

# Which primer-pair tables feed each marker's target-amplification pass.
_TARGET_PAIRS = {
    "mit": [
        f"{OUT}/primer_design/mit_pairs_filtered.tsv",
        f"{OUT}/primer_design/dbs/mit_dbs_single_pairs_filtered.tsv",
        f"{OUT}/primer_design/dbs/mit_dbs_tiled_pairs_filtered.tsv",
    ],
    "18S": [
        f"{OUT}/primer_design/18S_pairs_filtered.tsv",
    ],
}

# Target reference per marker. These are the target-species-filtered FASTAs
# staged by envs/install.sh, i.e. the same references every earlier phase used.
_TARGET_SUBJECT = {
    "mit": f"{DATA_REF}/mit_all.target.fasta",
    "18S": f"{DATA_REF}/18S_ref_db.target.fasta",
}


rule target_amplification_marker:
    """In-silico PCR of one marker's primer pairs against its target references."""
    input:
        pairs   = lambda w: _TARGET_PAIRS[w.marker],
        subject = lambda w: _TARGET_SUBJECT[w.marker],
    output:
        amplicons = f"{_SPEC_DIR}/{{marker}}_target_amplicons.tsv",
        summary   = f"{_SPEC_DIR}/{{marker}}_target_amplification.tsv",
        coverage  = f"{_SPEC_DIR}/{{marker}}_target_species_coverage.tsv",
    log:
        "logs/11_specificity/target_amplification_{marker}.log",
    message:
        "[11_specificity] target amplification {wildcards.marker} (blastn + primer3)"
    params:
        max_product   = lambda w: SPEC["max_product_bp"][w.marker],
        anneal_temp   = SPEC["anneal_temp_c"],
        logistic_k    = SPEC["logistic_k"],
        tp_window     = SPEC["three_prime_window"],
        max_mismatch  = SPEC["max_mismatch"],
        size_tol      = SPEC["size_tolerance"],
        word_size     = SPEC["blast_word_size"],
    threads: 4
    shell:
        r"""
        mkdir -p $(dirname {log}) $(dirname {output.summary})
        python scripts/python/insilico_pcr.py \
            --pairs {input.pairs} \
            --subject {input.subject} \
            --subject-kind reference \
            --amplicons-out {output.amplicons} \
            --summary-out   {output.summary} \
            --coverage-out  {output.coverage} \
            --max-product        {params.max_product} \
            --size-tolerance     {params.size_tol} \
            --three-prime-window {params.tp_window} \
            --max-mismatch       {params.max_mismatch} \
            --anneal-temp        {params.anneal_temp} \
            --logistic-k         {params.logistic_k} \
            --word-size          {params.word_size} \
            --threads            {threads} \
            > {log} 2>&1
        """


rule target_amplification:
    """Merge the per-marker passes into HANDOFF Step 5's named artefacts."""
    input:
        summaries = expand(f"{_SPEC_DIR}/{{marker}}_target_amplification.tsv",
                           marker=config["markers"]),
        coverages = expand(f"{_SPEC_DIR}/{{marker}}_target_species_coverage.tsv",
                           marker=config["markers"]),
    output:
        summary  = f"{_SPEC_DIR}/target_amplification.tsv",
        coverage = f"{_SPEC_DIR}/target_species_coverage.tsv",
    log:
        "logs/11_specificity/target_amplification_merge.log",
    message:
        "[11_specificity] merge per-marker target amplification tables"
    shell:
        r"""
        mkdir -p $(dirname {log})
        bash scripts/sh/concat_tsv.sh {output.summary}  {input.summaries} >  {log} 2>&1
        bash scripts/sh/concat_tsv.sh {output.coverage} {input.coverages} >> {log} 2>&1
        """


# Convenience target — request this to materialise Step 5's target half.
# Deliberately NOT in `rule all`: the off-target half needs hand-staged
# genomes, and a default-rule run must not imply Step 5 is complete.
rule specificity_target_all:
    input:
        f"{_SPEC_DIR}/target_amplification.tsv",
        f"{_SPEC_DIR}/target_species_coverage.tsv",


# ==========================================================================
# (b) OFF-TARGET amplification
# ==========================================================================
# Every primer set goes against every staged genome in one pass per genome:
# the pairs are cheap and the subject is not, so splitting by marker would buy
# nothing and pay for the blastn search twice.

_STAGED_DIR   = "data/staged"
_MANIFEST     = f"{_STAGED_DIR}/MANIFEST.tsv"
_OFF_DIR      = f"{_SPEC_DIR}/offtarget"
_ALL_PAIRS    = _TARGET_PAIRS["mit"] + _TARGET_PAIRS["18S"]


def _staged_subjects():
    """Manifest rows with `staged = TRUE` AND a file actually on disk.

    Read at DAG-build time. A row flagged TRUE whose FASTA has since been
    deleted is dropped with a warning rather than becoming a missing-input
    error, because the manifest is a hand-maintained ledger and the disk is
    the ground truth — the same rule `offtarget_manifest.py` heals by.
    """
    subjects = {}
    if not os.path.exists(_MANIFEST):
        return subjects
    with open(_MANIFEST) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["staged"].strip().upper() != "TRUE":
                continue
            path = os.path.join(_STAGED_DIR, row["expected_filename"])
            if not os.path.exists(path):
                print(f"[11_specificity] WARNING: {row['genus_species']} is "
                      f"staged=TRUE in the manifest but {path} is missing — "
                      f"skipped", file=sys.stderr)
                continue
            subjects[row["expected_filename"].rsplit(".", 1)[0]] = path
    return subjects


OFFTARGET_SUBJECTS = _staged_subjects()


rule offtarget_amplification_genome:
    """In-silico PCR of every primer pair against one staged off-target genome.

    No `--summary-out`: that is a full (pair x subject SEQUENCE) grid, and a
    staged genome is thousands of contigs. The amplicons ARE the evidence
    here — zero is the expected answer — and the zero-fill happens once per
    (pair, genome) in `offtarget_summary.py`.
    """
    input:
        pairs   = _ALL_PAIRS,
        subject = lambda w: OFFTARGET_SUBJECTS[w.subject],
    output:
        amplicons = f"{_OFF_DIR}/{{subject}}.amplicons.tsv",
    log:
        "logs/11_specificity/offtarget_{subject}.log",
    message:
        "[11_specificity] off-target amplification vs {wildcards.subject}"
    params:
        max_product  = SPEC["offtarget_max_product_bp"],
        anneal_temp  = SPEC["anneal_temp_c"],
        logistic_k   = SPEC["logistic_k"],
        tp_window    = SPEC["three_prime_window"],
        max_mismatch = SPEC["max_mismatch"],
        size_tol     = SPEC["size_tolerance"],
        word_size    = SPEC["offtarget_word_size"],
        db_dir       = f"{_OFF_DIR}/blastdb",
    threads: 8
    shell:
        r"""
        mkdir -p $(dirname {log}) $(dirname {output.amplicons})
        python scripts/python/insilico_pcr.py \
            --pairs {input.pairs} \
            --subject {input.subject} \
            --subject-kind offtarget \
            --amplicons-out {output.amplicons} \
            --max-product        {params.max_product} \
            --size-tolerance     {params.size_tol} \
            --three-prime-window {params.tp_window} \
            --max-mismatch       {params.max_mismatch} \
            --anneal-temp        {params.anneal_temp} \
            --logistic-k         {params.logistic_k} \
            --word-size          {params.word_size} \
            --db-dir             {params.db_dir} \
            --threads            {threads} \
            > {log} 2>&1
        """


rule offtarget_amplification:
    """Roll the per-genome evidence up to (pair, off_target_species, subject)."""
    input:
        amplicons = expand(f"{_OFF_DIR}/{{subject}}.amplicons.tsv",
                           subject=sorted(OFFTARGET_SUBJECTS)),
        pairs     = _ALL_PAIRS,
        manifest  = _MANIFEST,
    output:
        summary = f"{_SPEC_DIR}/offtarget_amplification.tsv",
    log:
        "logs/11_specificity/offtarget_amplification.log",
    message:
        "[11_specificity] roll up off-target amplification across staged genomes"
    params:
        threshold = SPEC["efficiency_threshold"],
    shell:
        r"""
        mkdir -p $(dirname {log})
        python scripts/python/offtarget_summary.py \
            --amplicons {input.amplicons} \
            --pairs     {input.pairs} \
            --manifest  {input.manifest} \
            --efficiency-threshold {params.threshold} \
            --out       {output.summary} \
            > {log} 2>&1
        """


rule panel_pass_fail:
    """Combine both halves into the per-pair verdict."""
    input:
        coverage  = f"{_SPEC_DIR}/target_species_coverage.tsv",
        offtarget = f"{_SPEC_DIR}/offtarget_amplification.tsv",
        manifest  = _MANIFEST,
    output:
        verdict = f"{_SPEC_DIR}/panel_pass_fail.tsv",
    log:
        "logs/11_specificity/panel_pass_fail.log",
    message:
        "[11_specificity] panel pass/fail verdict (target + off-target)"
    params:
        threshold = SPEC["efficiency_threshold"],
        congener  = SPEC["congener_purpose"],
    shell:
        r"""
        mkdir -p $(dirname {log})
        python scripts/python/panel_verdict.py \
            --target-coverage {input.coverage} \
            --offtarget       {input.offtarget} \
            --manifest        {input.manifest} \
            --efficiency-threshold {params.threshold} \
            --congener-purpose {params.congener} \
            --out             {output.verdict} \
            > {log} 2>&1
        """


# Convenience target — the whole of Step 5.
rule specificity_all:
    input:
        rules.specificity_target_all.input,
        f"{_SPEC_DIR}/offtarget_amplification.tsv",
        f"{_SPEC_DIR}/panel_pass_fail.tsv",
