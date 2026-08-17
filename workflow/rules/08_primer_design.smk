# 08_primer_design.smk — Phase 5 Step 3: primer3-driven primer-pair selection.
#
# Consumes outputs/entropy/{marker}.primer_candidates.tsv (end-anchored only —
# do NOT read interior_candidates.tsv end-anchored by design), builds a
# per-marker consensus template with a bijective aln↔template map, projects
# the left/right end-anchored blocks onto the template, and drives primer3
# per (left_block × right_block) combination.
#
# Per-marker product-size range (do NOT hardcode 4-5 kb globally; product-size is derived per-marker below).
#
# Informed caveat: I derived the achievable envelope from the actual left/
# right block separations on the majority-consensus template (see
# 08_primer_design.README.md § "Per-marker product-size range").
#
#   - MIT: HANDOFF nominates 4000-5000 bp. Achievable envelope from the
#     end-anchored blocks is 5230-5824 bp — the left block sits at aln
#     749-1048 and the two right blocks at aln 6357-6531 / 6605-6656, so
#     the amplicon effectively spans (almost) the full 5908 bp core.
#     Setting 4000-5000 would return zero pairs. **Widened to 5000-5900
#     bp** to cover the envelope; the "long-amplicon flanking primers"
#     intent is preserved.
#
#   - 18S: HANDOFF says derive from core_bounds. Core is aln 317-2218
#     (core_len 1902); achievable envelope from the two left blocks
#     (aln 317-368, 487-616) and two right blocks (aln 1988-2128,
#     2194-2218) is 1209-1682 bp. **Set to 1200-1700 bp.**
#
# Other primer3 params per HANDOFF Step 3: Tm 58-62 °C, GC 40-60 %, length
# 18-25 bp. Overridable on the primer3 driver's CLI.

OUT     = config["paths"]["outputs"]
REPORTS = config["paths"]["reports"]

# Locked per-marker product-size ranges. See rule note above.
_PRODUCT_MIN = {"mit": 5000, "18S": 1200}
_PRODUCT_MAX = {"mit": 5900, "18S": 1700}


rule primer3_mit:
    input:
        aln        = f"{OUT}/alignment/ma_mit.target.fasta",
        candidates = f"{OUT}/entropy/mit.primer_candidates.tsv",
    output:
        pairs    = f"{OUT}/primer_design/mit_pairs.tsv",
        filtered = f"{OUT}/primer_design/mit_pairs_filtered.tsv",
        unfilled = f"{OUT}/primer_design/mit_pairs_unfilled.tsv",
        template = f"{OUT}/primer_design/mit_template.fasta",
    log:
        "logs/08_primer_design/primer3_mit.log",
    message:
        "[08_primer_design] primer3 mit (product-size 5000-5900; widened from HANDOFF 4-5 kb — see README)"
    params:
        product_min = _PRODUCT_MIN["mit"],
        product_max = _PRODUCT_MAX["mit"],
        tm_min      = 58.0,
        tm_max      = 62.0,
        gc_min      = 40.0,
        gc_max      = 60.0,
        len_min     = 18,
        len_max     = 25,
        num_return  = 20,
    shell:
        r"""
        mkdir -p $(dirname {log}) $(dirname {output.pairs})
        python scripts/python/run_primer3.py \
            --marker mit \
            --aln {input.aln} \
            --candidates {input.candidates} \
            --pairs-out    {output.pairs} \
            --filtered-out {output.filtered} \
            --unfilled-out {output.unfilled} \
            --template-out {output.template} \
            --product-min {params.product_min} \
            --product-max {params.product_max} \
            --tm-min {params.tm_min} --tm-max {params.tm_max} \
            --gc-min {params.gc_min} --gc-max {params.gc_max} \
            --len-min {params.len_min} --len-max {params.len_max} \
            --num-return {params.num_return} \
            > {log} 2>&1
        """


rule primer3_18S:
    input:
        aln        = f"{OUT}/alignment/ma_18S.target.fasta",
        candidates = f"{OUT}/entropy/18S.primer_candidates.tsv",
    output:
        pairs    = f"{OUT}/primer_design/18S_pairs.tsv",
        filtered = f"{OUT}/primer_design/18S_pairs_filtered.tsv",
        unfilled = f"{OUT}/primer_design/18S_pairs_unfilled.tsv",
        template = f"{OUT}/primer_design/18S_template.fasta",
    log:
        "logs/08_primer_design/primer3_18S.log",
    message:
        "[08_primer_design] primer3 18S (product-size 1200-1700; derived from actual block separations)"
    params:
        product_min = _PRODUCT_MIN["18S"],
        product_max = _PRODUCT_MAX["18S"],
        tm_min      = 58.0,
        tm_max      = 62.0,
        gc_min      = 40.0,
        gc_max      = 60.0,
        len_min     = 18,
        len_max     = 25,
        num_return  = 20,
    shell:
        r"""
        mkdir -p $(dirname {log}) $(dirname {output.pairs})
        python scripts/python/run_primer3.py \
            --marker 18S \
            --aln {input.aln} \
            --candidates {input.candidates} \
            --pairs-out    {output.pairs} \
            --filtered-out {output.filtered} \
            --unfilled-out {output.unfilled} \
            --template-out {output.template} \
            --product-min {params.product_min} \
            --product-max {params.product_max} \
            --tm-min {params.tm_min} --tm-max {params.tm_max} \
            --gc-min {params.gc_min} --gc-max {params.gc_max} \
            --len-min {params.len_min} --len-max {params.len_max} \
            --num-return {params.num_return} \
            > {log} 2>&1
        """


rule multiplex_compatibility:
    input:
        mit_filtered = f"{OUT}/primer_design/mit_pairs_filtered.tsv",
        s18_filtered = f"{OUT}/primer_design/18S_pairs_filtered.tsv",
    output:
        tsv = f"{OUT}/primer_design/multiplex_compatibility.tsv",
    log:
        "logs/08_primer_design/multiplex_compatibility.log",
    message:
        "[08_primer_design] pairwise multiplex compatibility (Tm spread, cross-dimer ΔG, len diff)"
    shell:
        r"""
        mkdir -p $(dirname {log}) $(dirname {output.tsv})
        python scripts/python/multiplex_compat.py \
            --pairs {input.mit_filtered} {input.s18_filtered} \
            --out   {output.tsv} \
            > {log} 2>&1
        """


# Convenience target — request this to materialise all Phase 5 Step 3 outputs.
rule primer_design_all:
    input:
        f"{OUT}/primer_design/mit_pairs.tsv",
        f"{OUT}/primer_design/mit_pairs_filtered.tsv",
        f"{OUT}/primer_design/mit_pairs_unfilled.tsv",
        f"{OUT}/primer_design/18S_pairs.tsv",
        f"{OUT}/primer_design/18S_pairs_filtered.tsv",
        f"{OUT}/primer_design/18S_pairs_unfilled.tsv",
        f"{OUT}/primer_design/multiplex_compatibility.tsv",
