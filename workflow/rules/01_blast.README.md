# 01_blast — self-BLAST

Two rules per marker, run twice (once for `mit`, once for `18S`):

- `make_blast_db` — `makeblastdb -dbtype nucl` over `data/reference/<fasta>` → `outputs/blast/db_{marker}.{ndb,nhr,nin,njs,not,nsq,ntf,nto}`. The DB lives in `outputs/` rather than overwriting `data/reference/` so derivative artefacts stay out of the staged-source area.
- `self_blast` — `blastn -task blastn` of the FASTA against its own DB → `outputs/blast/self_blast_{marker}.tsv`. Tabular outfmt with the canonical 10 columns: `qseqid sseqid pident length mismatch gapopen qlen slen evalue bitscore`. The column order is the contract every Phase-2 reader assumes.

## Provenance

`mit_similarity.Rmd` lines 8-20 (MIT) and 185-196 (18S). Both lifted unchanged.

## Validation

`outputs/blast/self_blast_{marker}.tsv` should match `scripts/legacy/speciation_long/outputs/self_blast{,_18S}.tsv`. Step 2 already verified this: 23 506 MIT rows match the legacy row count exactly; 4 of those rows differ at the gapopen↔mismatch boundary by 1 (BLAST 2.16.x patch-level drift; aln_len and bitscore unchanged). 18S is byte-identical.

## Interpretation

Each row is one HSP (high-scoring pair). Multiple HSPs per query–subject pair are normal. `pident`, `length`, `qcov`/`scov`, `cov_shorter` — never use any of these alone. The downstream cross-species filter combines `pident ≥ 99` and `cov_shorter ≥ 0.99` to flag near-identical pairs.
