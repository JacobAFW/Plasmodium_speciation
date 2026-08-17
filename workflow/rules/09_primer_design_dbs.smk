# 09_primer_design_dbs.smk — Phase 5 Step 3b: DBS-compatible MIT amplicon tier
# (≤4 kb). Parallel to `08_primer_design.smk`; that Tier-A ~5.5 kb rule stays
# authoritative for whole-blood samples.
#
# Why: the assay is scaling to in-country partners who work from dried
# blood spots (DBS), which constrain ONT amplicons to <=4 kb. MIT's Tier-A ~5.5 kb amplicon is DBS-unsafe; 18S (~1.3–1.4 kb)
# is already DBS-ready and is REUSED here unchanged — no new 18S rule.
#
# Two candidate MIT geometries (Jacob picks from the trade-off table):
#   single  — one ≤4 kb amplicon,      product-size 3300–4000 bp
#   tiled   — two amplicons that jointly tile the MIT core with a modest
#             overlap, each ≤4 kb,      product-size 3000–3500 bp per tile
#
# The trade-off is materialised in the resolution table (Task C): per
# species pair × option → {resolvable, n_diagnostic_sites}. The Tier-A
# column is reconciled against `outputs/cross_species/resolution_table.tsv`.
#
# Downstream (Task D): multiplex compatibility over the FULL DBS panel
# (DBS-MIT singles + DBS-MIT tiles + 18S filtered) is produced by
# `multiplex_compat.py`, reused without edit.

OUT     = config["paths"]["outputs"]

_DBS_DIR = f"{OUT}/primer_design/dbs"

# Locked per-option product-size ranges (see README § "Bands used per option").
_SINGLE_MIN, _SINGLE_MAX = 3300, 4000
_TILED_MIN,  _TILED_MAX  = 3000, 3500


rule primer3_dbs_mit_single:
    input:
        aln        = f"{OUT}/alignment/ma_mit.target.fasta",
        candidates = f"{OUT}/entropy/mit.interior_candidates.tsv",
    output:
        pairs    = f"{_DBS_DIR}/mit_dbs_single_pairs.tsv",
        filtered = f"{_DBS_DIR}/mit_dbs_single_pairs_filtered.tsv",
        unfilled = f"{_DBS_DIR}/mit_dbs_single_pairs_unfilled.tsv",
    log:
        "logs/09_primer_design_dbs/primer3_dbs_mit_single.log",
    message:
        "[09_primer_design_dbs] primer3 mit single ≤4 kb (product-size 3300-4000)"
    params:
        product_min = _SINGLE_MIN,
        product_max = _SINGLE_MAX,
    shell:
        r"""
        mkdir -p $(dirname {log}) $(dirname {output.pairs})
        python scripts/python/run_primer3_dbs.py single \
            --aln          {input.aln} \
            --candidates   {input.candidates} \
            --pairs-out    {output.pairs} \
            --filtered-out {output.filtered} \
            --unfilled-out {output.unfilled} \
            --product-min  {params.product_min} \
            --product-max  {params.product_max} \
            --num-return   20 \
            > {log} 2>&1
        """


rule primer3_dbs_mit_tiled:
    input:
        aln        = f"{OUT}/alignment/ma_mit.target.fasta",
        candidates = f"{OUT}/entropy/mit.interior_candidates.tsv",
        core       = f"{OUT}/entropy/mit.core_bounds.tsv",
    output:
        pairs    = f"{_DBS_DIR}/mit_dbs_tiled_pairs.tsv",
        filtered = f"{_DBS_DIR}/mit_dbs_tiled_pairs_filtered.tsv",
        unfilled = f"{_DBS_DIR}/mit_dbs_tiled_pairs_unfilled.tsv",
    log:
        "logs/09_primer_design_dbs/primer3_dbs_mit_tiled.log",
    message:
        "[09_primer_design_dbs] primer3 mit tiled ~3–3.5 kb (tile1 start×mid, tile2 mid×end)"
    params:
        product_min = _TILED_MIN,
        product_max = _TILED_MAX,
    shell:
        r"""
        mkdir -p $(dirname {log}) $(dirname {output.pairs})
        python scripts/python/run_primer3_dbs.py tiled \
            --aln          {input.aln} \
            --candidates   {input.candidates} \
            --core-bounds  {input.core} \
            --pairs-out    {output.pairs} \
            --filtered-out {output.filtered} \
            --unfilled-out {output.unfilled} \
            --product-min  {params.product_min} \
            --product-max  {params.product_max} \
            --num-return   20 \
            > {log} 2>&1
        """


rule mit_amplicon_resolution_tradeoff:
    # Task C. Species-pair diagnostic-position counts per option, reconciled
    # against v1's mit_resolved column. The load-bearing decision artifact.
    input:
        aln            = f"{OUT}/alignment/ma_mit.target.fasta",
        tier_a         = f"{OUT}/primer_design/mit_pairs_filtered.tsv",
        single         = f"{_DBS_DIR}/mit_dbs_single_pairs_filtered.tsv",
        tiled          = f"{_DBS_DIR}/mit_dbs_tiled_pairs_filtered.tsv",
        v1_resolution  = f"{OUT}/cross_species/resolution_table.tsv",
    output:
        tsv = f"{_DBS_DIR}/mit_amplicon_resolution_tradeoff.tsv",
    log:
        "logs/09_primer_design_dbs/resolution_tradeoff.log",
    message:
        "[09_primer_design_dbs] mit amplicon resolution trade-off (tier_a vs single vs tiled)"
    shell:
        r"""
        mkdir -p $(dirname {log}) $(dirname {output.tsv})
        python scripts/python/run_primer3_dbs.py resolution \
            --aln              {input.aln} \
            --tier-a-pairs     {input.tier_a} \
            --single-pairs     {input.single} \
            --tiled-pairs      {input.tiled} \
            --resolution-table {input.v1_resolution} \
            --out              {output.tsv} \
            > {log} 2>&1
        """


rule dbs_multiplex_compatibility:
    # Task D. DBS panel = DBS-MIT singles + DBS-MIT tiles + existing 18S
    # filtered pairs. Reuses `multiplex_compat.py` (import from 08's script
    # collection, do NOT edit — same file is being used by the Tier-A rule).
    input:
        mit_single = f"{_DBS_DIR}/mit_dbs_single_pairs_filtered.tsv",
        mit_tiled  = f"{_DBS_DIR}/mit_dbs_tiled_pairs_filtered.tsv",
        s18        = f"{OUT}/primer_design/18S_pairs_filtered.tsv",
    output:
        tsv = f"{_DBS_DIR}/dbs_multiplex_compatibility.tsv",
    log:
        "logs/09_primer_design_dbs/dbs_multiplex_compatibility.log",
    message:
        "[09_primer_design_dbs] dbs panel multiplex compatibility (mit-single + mit-tiled + 18S)"
    shell:
        r"""
        mkdir -p $(dirname {log}) $(dirname {output.tsv})
        python scripts/python/multiplex_compat.py \
            --pairs {input.mit_single} {input.mit_tiled} {input.s18} \
            --out   {output.tsv} \
            > {log} 2>&1
        """


# Convenience target — request this to materialise all Phase 5 Step 3b outputs.
rule primer_design_dbs_all:
    input:
        f"{_DBS_DIR}/mit_dbs_single_pairs.tsv",
        f"{_DBS_DIR}/mit_dbs_single_pairs_filtered.tsv",
        f"{_DBS_DIR}/mit_dbs_single_pairs_unfilled.tsv",
        f"{_DBS_DIR}/mit_dbs_tiled_pairs.tsv",
        f"{_DBS_DIR}/mit_dbs_tiled_pairs_filtered.tsv",
        f"{_DBS_DIR}/mit_dbs_tiled_pairs_unfilled.tsv",
        f"{_DBS_DIR}/mit_amplicon_resolution_tradeoff.tsv",
        f"{_DBS_DIR}/dbs_multiplex_compatibility.tsv",
