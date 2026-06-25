# 07_threshold_sensitivity.smk — Phase 5 Step 2: calibrate the interior-side
# candidate-window thresholds.
#
# Single rule. Calls scripts/python/interior_threshold_sweep.py over a grid
# of (cov_thr, min_run) values for both markers and reports per-cell metrics
# from the resulting interior_candidates.tsv emissions. Each cell runs
# seq_entropy.py to a temp directory, so the canonical outputs/entropy/
# files are untouched.
#
# Provenance: HANDOFF.md Phase 5 Step 2. Step 1 audit established that the
# Imwong forward-primer region (18S aln 1683-1740) is rejected by the
# end_region gate, not by entropy/coverage thresholds. Step 2 confirms which
# (cov_thr, min_run) combinations surface the region in the new interior
# output.

OUT = config["paths"]["outputs"]


rule interior_threshold_sweep:
    input:
        mit_aln = f"{OUT}/alignment/ma_mit.target.fasta",
        s18_aln = f"{OUT}/alignment/ma_18S.target.fasta",
    output:
        tsv = f"{OUT}/sensitivity/interior_threshold_sweep.tsv",
    log:
        "logs/07_threshold_sensitivity/interior_threshold_sweep.log",
    message:
        "[07_threshold_sensitivity] sweep (cov_thr × min_run) over interior candidates"
    params:
        cov_thrs    = "0.70,0.80,0.85,0.90",
        min_runs    = "15,18,21,25",
        entropy_thr = 0.20,
        windows     = "25,50",
        imwong_start = 1683,
        imwong_end   = 1740,
    shell:
        r"""
        mkdir -p $(dirname {log}) $(dirname {output.tsv})
        python scripts/python/interior_threshold_sweep.py \
            --mit-aln {input.mit_aln} \
            --s18-aln {input.s18_aln} \
            --cov-thrs {params.cov_thrs} \
            --min-runs {params.min_runs} \
            --entropy-thr {params.entropy_thr} \
            --windows {params.windows} \
            --imwong-start {params.imwong_start} \
            --imwong-end   {params.imwong_end} \
            --out {output.tsv} \
            > {log} 2>&1
        """
