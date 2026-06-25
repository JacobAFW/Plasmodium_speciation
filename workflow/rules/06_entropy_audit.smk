# 06_entropy_audit.smk — Phase 5 Step 1: audit the entropy/coverage pipeline
# before trusting any downstream threshold sensitivity or primer3 selection.
#
# Three rules:
#   - entropy_audit_mit  : independent recomputation + diff + spot-check +
#                          top/bottom window eyeball.
#   - entropy_audit_18S  : same, plus a 100-bootstrap CI band (small n=17
#                          reference panel).
#   - inspect_imwong_position : drill-down at the published PlasmoM forward
#                          primer binding region (18S aln 1683-1740).
#
# Provenance: HANDOFF.md Phase 5 Step 1.

OUT     = config["paths"]["outputs"]
REPORTS = config["paths"]["reports"]


def _audit_outputs(marker):
    base = {
        "recompute":   f"{OUT}/audit/{marker}_entropy_recompute.tsv",
        "diff":        f"{OUT}/audit/{marker}_entropy_diff.tsv",
        "spotcheck":   f"{OUT}/audit/{marker}_coverage_spotcheck.tsv",
        "eyeball_tsv": f"{OUT}/audit/{marker}_window_eyeball.tsv",
        "eyeball_png": f"{REPORTS}/figures/{marker}_window_eyeball.png",
        "eyeball_svg": f"{REPORTS}/figures/{marker}_window_eyeball.svg",
    }
    if marker == "18S":
        base["bootstrap_tsv"] = f"{OUT}/audit/18S_entropy_bootstrap.tsv"
        base["bootstrap_png"] = f"{REPORTS}/figures/18S_entropy_bootstrap.png"
        base["bootstrap_svg"] = f"{REPORTS}/figures/18S_entropy_bootstrap.svg"
    return base


rule entropy_audit_mit:
    input:
        aln        = f"{OUT}/alignment/ma_mit.target.fasta",
        perpos     = f"{OUT}/entropy/mit.per_position.tsv",
        windows    = f"{OUT}/entropy/mit.windows.tsv",
        candidates = f"{OUT}/entropy/mit.primer_candidates.tsv",
    output:
        recompute   = _audit_outputs("mit")["recompute"],
        diff        = _audit_outputs("mit")["diff"],
        spotcheck   = _audit_outputs("mit")["spotcheck"],
        eyeball_tsv = _audit_outputs("mit")["eyeball_tsv"],
        eyeball_png = _audit_outputs("mit")["eyeball_png"],
        eyeball_svg = _audit_outputs("mit")["eyeball_svg"],
    log:
        "logs/06_entropy_audit/entropy_audit_mit.log",
    message:
        "[06_entropy_audit] independent recomputation + audit checks (mit)"
    shell:
        r"""
        mkdir -p $(dirname {log}) $(dirname {output.recompute})
        Rscript scripts/R/entropy_audit.R \
            --marker mit \
            --alignment {input.aln} \
            --pipeline-perpos     {input.perpos} \
            --pipeline-windows    {input.windows} \
            --pipeline-candidates {input.candidates} \
            --recompute-out {output.recompute} \
            --diff-out      {output.diff} \
            --spotcheck-out {output.spotcheck} \
            --eyeball-tsv   {output.eyeball_tsv} \
            --eyeball-png   {output.eyeball_png} \
            --eyeball-svg   {output.eyeball_svg} \
            > {log} 2>&1
        """


rule entropy_audit_18S:
    input:
        aln        = f"{OUT}/alignment/ma_18S.target.fasta",
        perpos     = f"{OUT}/entropy/18S.per_position.tsv",
        windows    = f"{OUT}/entropy/18S.windows.tsv",
        candidates = f"{OUT}/entropy/18S.primer_candidates.tsv",
    output:
        recompute     = _audit_outputs("18S")["recompute"],
        diff          = _audit_outputs("18S")["diff"],
        spotcheck     = _audit_outputs("18S")["spotcheck"],
        eyeball_tsv   = _audit_outputs("18S")["eyeball_tsv"],
        eyeball_png   = _audit_outputs("18S")["eyeball_png"],
        eyeball_svg   = _audit_outputs("18S")["eyeball_svg"],
        bootstrap_tsv = _audit_outputs("18S")["bootstrap_tsv"],
        bootstrap_png = _audit_outputs("18S")["bootstrap_png"],
        bootstrap_svg = _audit_outputs("18S")["bootstrap_svg"],
    log:
        "logs/06_entropy_audit/entropy_audit_18S.log",
    message:
        "[06_entropy_audit] independent recomputation + audit + bootstrap (18S)"
    shell:
        r"""
        mkdir -p $(dirname {log}) $(dirname {output.recompute})
        Rscript scripts/R/entropy_audit.R \
            --marker 18S \
            --alignment {input.aln} \
            --pipeline-perpos     {input.perpos} \
            --pipeline-windows    {input.windows} \
            --pipeline-candidates {input.candidates} \
            --recompute-out {output.recompute} \
            --diff-out      {output.diff} \
            --spotcheck-out {output.spotcheck} \
            --eyeball-tsv   {output.eyeball_tsv} \
            --eyeball-png   {output.eyeball_png} \
            --eyeball-svg   {output.eyeball_svg} \
            --bootstrap-tsv {output.bootstrap_tsv} \
            --bootstrap-png {output.bootstrap_png} \
            --bootstrap-svg {output.bootstrap_svg} \
            --bootstrap-n   100 \
            --seed 42 \
            > {log} 2>&1
        """


# Drill-down at the 18S alignment positions where the published Imwong
# forward primers bind. Per Chapter 6: PlasmoM_N1F at aln 1683-1702 and
# PlasmoM_N2F at aln 1718-1740. The script re-derives these positions to
# avoid hard-coding bounds.
rule inspect_imwong_position:
    input:
        aln     = f"{OUT}/alignment/ma_18S.target.fasta",
        perpos  = f"{OUT}/entropy/18S.per_position.tsv",
        windows = f"{OUT}/entropy/18S.windows.tsv",
        candidates = f"{OUT}/entropy/18S.primer_candidates.tsv",
    output:
        tsv = f"{OUT}/audit/imwong_forward_position.tsv",
    log:
        "logs/06_entropy_audit/inspect_imwong_position.log",
    message:
        "[06_entropy_audit] Imwong forward-primer position drill-down (18S aln 1683-1740)"
    params:
        primers = "PlasmoM_N1F=ATGGCCGTTTTTAGTTCGTG,PlasmoM_N2F=GTTAATTCCGATAACGAACGAGA",
        # Threshold pairs: (cov_thr, ent_thr). v1 default is (0.90, 0.20).
        # Each successive pair relaxes one or both knobs.
        thresholds = "0.90:0.20,0.85:0.25,0.85:0.30,0.80:0.30,0.80:0.40,0.70:0.50",
    shell:
        r"""
        mkdir -p $(dirname {log}) $(dirname {output.tsv})
        Rscript scripts/R/inspect_imwong_position.R \
            --alignment {input.aln} \
            --perpos    {input.perpos} \
            --windows   {input.windows} \
            --primers   '{params.primers}' \
            --thresholds '{params.thresholds}' \
            --out       {output.tsv} \
            > {log} 2>&1
        """
