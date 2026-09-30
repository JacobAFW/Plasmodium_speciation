# 12_resolution_margin.smk — Phase 5 correction: resolution measured as an
# ASSIGNMENT MARGIN, and the fieldi/simiovale diagnostic-window check that the
# DBS <=4 kb tier turns on.
#
# WHY THIS RULE FILE EXISTS
# -------------------------
# The project's older resolution criterion — "does this pair carry >= 1 fixed
# diagnostic site inside the amplicon" — returns 55/55 under every amplicon
# option. It cannot discriminate between options, cannot fail a pair, and is
# therefore not decision-relevant. It was also optimistic in the one place it
# mattered: `fieldi vs simiovale` was scored resolvable while a real best-hit
# read-out called BOTH *P. fieldi* references *P. simiovale*.
#
# These two rules replace the threshold with a margin:
#
#   resolution_margin        per reference, and per species pair, how much
#                            own-species identity beats other-species identity.
#                            Negative => a real read from that isolate would be
#                            called wrong. Emits margins, never a pass mark.
#   dbs_diagnostic_coverage  per primer pair, does the amplicon still carry the
#                            cox1 and/or cytb sites that separate
#                            fieldi/simiovale — and does it still separate them
#                            when measured directly.
#
# Both read the CURATED reference set (see the reference-curation note in
# 03_alignment.README.md): AB444133 relabelled `Punverified_*` and held out of
# the assignment database, AY800109 dropped as an exact duplicate.
#
# Deliberately NOT in `rule all`. These are a correction pass over Phase 5's
# existing primer tables, not a step every run has to redo.

DATA_REF = config["paths"]["data_ref"]
OUT      = config["paths"]["outputs"]

_MARGIN_PAIRS = {
    "MIT_TierA":      f"{OUT}/primer_design/mit_pairs_filtered.tsv",
    "MIT_DBS_single": f"{OUT}/primer_design/dbs/mit_dbs_single_pairs_filtered.tsv",
}
_MARGIN_TILED = {
    "MIT_DBS_tiled":  f"{OUT}/primer_design/dbs/mit_dbs_tiled_pairs_filtered.tsv",
}
_ALL_MIT_PAIRS = {**_MARGIN_PAIRS, **_MARGIN_TILED}


rule resolution_margin:
    """Leave-one-out assignment + per-pair margin on the curated MIT set."""
    input:
        alignment = f"{OUT}/alignment/ma_mit.target.fasta",
        pairs     = list(_ALL_MIT_PAIRS.values()),
    output:
        loo  = f"{OUT}/cross_species/loo_assignment_mit.tsv",
        pair = f"{OUT}/cross_species/resolution_table_margin_mit.tsv",
    log:
        "logs/12_resolution_margin/resolution_margin_mit.log",
    message:
        "[12_resolution_margin] leave-one-out assignment margins (MIT)"
    params:
        targets = ",".join(config["species_targets"]),
        pairs   = " ".join(f"--pairs {k}={v}" for k, v in _MARGIN_PAIRS.items()),
        tiled   = " ".join(f"--tiled {k}={v}" for k, v in _MARGIN_TILED.items()),
    shell:
        r"""
        mkdir -p $(dirname {log})
        python scripts/python/resolution_margin.py \
            --alignment {input.alignment} \
            {params.pairs} {params.tiled} \
            --species-targets '{params.targets}' \
            --loo-out  {output.loo} \
            --pair-out {output.pair} \
            > {log} 2>&1
        """


rule dbs_diagnostic_coverage:
    """Does each candidate amplicon keep the cox1 and/or cytb sites that
    separate fieldi from simiovale — and does it still separate them?"""
    input:
        alignment = f"{OUT}/alignment/ma_mit.target.fasta",
        pairs     = list(_ALL_MIT_PAIRS.values()),
        config    = "workflow/config.yaml",
    output:
        tsv = f"{OUT}/primer_design/dbs/dbs_diagnostic_coverage.tsv",
    log:
        "logs/12_resolution_margin/dbs_diagnostic_coverage.log",
    message:
        "[12_resolution_margin] fieldi/simiovale diagnostic-window coverage per pair"
    params:
        pairs = " ".join(f"--pairs {k}={v}" for k, v in _ALL_MIT_PAIRS.items()),
    shell:
        r"""
        mkdir -p $(dirname {log})
        python scripts/python/dbs_diagnostic_coverage.py \
            --config {input.config} \
            --alignment {input.alignment} \
            {params.pairs} \
            --out {output.tsv} \
            > {log} 2>&1
        """


rule resolution_margin_all:
    input:
        rules.resolution_margin.output,
        rules.dbs_diagnostic_coverage.output,
