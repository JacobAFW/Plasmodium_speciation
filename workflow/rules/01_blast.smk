# 01_blast.smk — Phase 1: build BLAST DBs from staged references and run a
# self-BLAST per marker. Outputs are the contract every Phase-2 rule reads.
#
# Provenance: mit_similarity.Rmd lines 8-20 (MIT) and 185-196 (18S).
#
# Naming: the DB is a derivative artefact, so it lives under outputs/blast/
# rather than overwriting data/reference/. Canonical pattern:
#   outputs/blast/db_{marker}.{ndb,nhr,nin,njs,not,nsq,ntf,nto}
# We declare .nsq as the rule's primary output; makeblastdb writes the
# whole set atomically.

DATA_REF = config["paths"]["data_ref"]
OUT      = config["paths"]["outputs"]

REF_FASTA = {
    "mit": f"{DATA_REF}/mit_all.fasta",
    "18S": f"{DATA_REF}/18S_ref_db.fasta",
}


rule make_blast_db:
    input:
        fasta = lambda w: REF_FASTA[w.marker],
    output:
        primary = f"{OUT}/blast/db_{{marker}}.nsq",
    log:
        "logs/01_blast/make_blast_db_{marker}.log",
    message:
        "[01_blast] makeblastdb {wildcards.marker}"
    params:
        prefix = f"{OUT}/blast/db_{{marker}}",
    shell:
        r"""
        mkdir -p $(dirname {log}) $(dirname {output.primary})
        makeblastdb -in {input.fasta} -dbtype nucl -out {params.prefix} \
            -logfile {log}
        """


rule self_blast:
    input:
        fasta   = lambda w: REF_FASTA[w.marker],
        db      = f"{OUT}/blast/db_{{marker}}.nsq",
    output:
        tsv = f"{OUT}/blast/self_blast_{{marker}}.tsv",
    log:
        "logs/01_blast/self_blast_{marker}.log",
    message:
        "[01_blast] self-BLAST {wildcards.marker}"
    params:
        db_prefix = f"{OUT}/blast/db_{{marker}}",
    shell:
        r"""
        mkdir -p $(dirname {log}) $(dirname {output.tsv})
        bash scripts/sh/run_self_blast.sh \
            {input.fasta} {params.db_prefix} {output.tsv} > {log} 2>&1
        """
