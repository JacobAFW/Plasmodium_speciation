#!/usr/bin/env bash
# legacy_replicate.sh — single-file driver that replays mit_similarity.Rmd
# end-to-end against the project-staged refs in data/reference/, writing
# outputs to outputs/legacy_replicate/.
#
# Step 2 of HANDOFF.md (Approach A — verify the install supports every
# tool the Rmd needs, before refactoring into Snakemake).
#
# Run from the project root:
#   bash scripts/legacy_replicate.sh

set -euo pipefail

PROJECT_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )/.." && pwd )"
cd "$PROJECT_ROOT"

# shellcheck disable=SC1091
set +eu; source envs/activate.sh; set -eu

REF=data/reference
OUT=outputs/legacy_replicate
mkdir -p "$OUT/tmp_db"

LEGACY_SCRIPTS="$PROJECT_ROOT/scripts/legacy/speciation_long/scripts"

echo "==> [1/12] makeblastdb (MIT)"
makeblastdb -in "$REF/mit_all.fasta" -dbtype nucl -out "$OUT/tmp_db/homology_db" \
  -logfile "$OUT/tmp_db/makeblastdb_mit.log"

echo "==> [2/12] blastn self-BLAST (MIT)"
blastn -query "$REF/mit_all.fasta" -db "$OUT/tmp_db/homology_db" \
  -out "$OUT/self_blast.tsv" \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qlen slen evalue bitscore" \
  -task blastn

echo "==> [3/12] MIT cross-species filter (suspicious + closest)"
python <<'PY'
import pandas as pd, re

df = pd.read_csv(
    "outputs/legacy_replicate/self_blast.tsv", sep="\t",
    names=["query","subject","pident","aln_len","mismatch","gapopen",
           "qlen","slen","evalue","bitscore"],
)
species_pat = re.compile(
    r"(vivax|falciparum|knowlesi|malariae|ovale|coat|inui|fieldi|cyno|simiovale|simium)",
    flags=re.IGNORECASE,
)
df = df[df["query"].str.contains(species_pat, na=False)
        & df["subject"].str.contains(species_pat, na=False)].copy()
df["query_species"] = (df["query"].str.extract(r"^([^_]+)")[0].str.lower()
                       .str.replace("pcynomologi","pcynomolgi", regex=False))
df["subject_species"] = (df["subject"].str.extract(r"^([^_]+)")[0].str.lower()
                         .str.replace("pcynomologi","pcynomolgi", regex=False))
df = df[(df["query"]!=df["subject"]) & (df["query_species"]!=df["subject_species"])].copy()
df["qcov"] = df["aln_len"]/df["qlen"]
df["scov"] = df["aln_len"]/df["slen"]
df["shorter_len"] = df[["qlen","slen"]].min(axis=1)
df["cov_shorter"] = df["aln_len"]/df["shorter_len"]
df_top5 = (df.sort_values(["query","subject","bitscore","aln_len"],
                          ascending=[True,True,False,False])
             .groupby(["query","subject"], as_index=False).head(5))
PID_THR=99; COV_SHORTER_THR=0.99; MIN_SHORTER_LEN=10
suspicious = df_top5[(df_top5["pident"]>=PID_THR)
                     & (df_top5["cov_shorter"]>=COV_SHORTER_THR)
                     & (df_top5["shorter_len"]>=MIN_SHORTER_LEN)].copy()
# species_pair: alphabetically-sorted unordered pair. Present in the legacy
# output but missing from the Rmd as checked in — known legacy/Rmd drift.
suspicious["species_pair"] = suspicious.apply(
    lambda r: " vs ".join(sorted([r["query_species"], r["subject_species"]])), axis=1
)
print("MIT suspicious:", len(suspicious))
suspicious.to_csv("outputs/legacy_replicate/suspicious_cross_species_near_identity.tsv",
                  sep="\t", index=False)

best = (df.sort_values(["query","bitscore","aln_len"], ascending=[True,False,False])
          .groupby("query", as_index=False).first())
best["pct_covered"] = best["cov_shorter"]*100
best["divergence_pct"] = 100 - best["pident"]
summary = best[["query","query_species","subject","subject_species","pident",
                "divergence_pct","aln_len","pct_covered","mismatch","gapopen",
                "bitscore"]].sort_values("pident", ascending=False)
summary.to_csv("outputs/legacy_replicate/closest_cross_species_match_per_sequence.tsv",
               sep="\t", index=False)
PY

echo "==> [4/12] makeblastdb (18S)"
makeblastdb -in "$REF/18S_ref_db.fasta" -dbtype nucl -out "$OUT/tmp_db/homology_db_18S" \
  -logfile "$OUT/tmp_db/makeblastdb_18S.log"

echo "==> [5/12] blastn self-BLAST (18S)"
blastn -query "$REF/18S_ref_db.fasta" -db "$OUT/tmp_db/homology_db_18S" \
  -out "$OUT/self_blast_18S.tsv" \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qlen slen evalue bitscore" \
  -task blastn

echo "==> [6/12] 18S header → species map"
python <<'PY'
import pandas as pd, re

def extract_species(header: str):
    m = re.search(r"\bP\.\s*([A-Za-z0-9_-]+)\b", header)
    if m: return f"p.{m.group(1).lower()}"
    m = re.search(r"\bPlasmodium\s+([A-Za-z0-9_-]+)\b", header, flags=re.IGNORECASE)
    if m: return f"p.{m.group(1).lower()}"
    m = re.search(r"\bP\.\s+([A-Za-z0-9_-]+)\b", header)
    if m: return f"p.{m.group(1).lower()}"
    return None

rows = []
with open("data/reference/18S_ref_db.fasta") as f:
    for line in f:
        if line.startswith(">"):
            header = line[1:].strip()
            acc = header.split()[0]
            rows.append({"acc": acc, "species": extract_species(header), "header": header})
pd.DataFrame(rows).to_csv("outputs/legacy_replicate/18S_acc_to_header.tsv",
                           sep="\t", index=False)
PY

echo "==> [7/12] 18S cross-species filter (suspicious + closest)"
python <<'PY'
import pandas as pd, re

df = pd.read_csv("outputs/legacy_replicate/self_blast_18S.tsv", sep="\t",
                 names=["query","subject","pident","aln_len","mismatch","gapopen",
                        "qlen","slen","evalue","bitscore"])
m = pd.read_csv("outputs/legacy_replicate/18S_acc_to_header.tsv", sep="\t")
df = (df.merge(m[["acc","species"]].rename(columns={"acc":"query","species":"query_species"}),
               on="query", how="left")
        .merge(m[["acc","species"]].rename(columns={"acc":"subject","species":"subject_species"}),
               on="subject", how="left"))
df["query"]   = df["query"]   + "_" + df["query_species"].fillna("unknown")
df["subject"] = df["subject"] + "_" + df["subject_species"].fillna("unknown")
df = df.drop(columns=["query_species","subject_species"])

species_pat = re.compile(
    r"(vivax|falciparum|knowlesi|malariae|ovale|coat|inui|fieldi|cyno|simiovale|simium)",
    flags=re.IGNORECASE,
)
df = df[df["query"].str.contains(species_pat, na=False)
        & df["subject"].str.contains(species_pat, na=False)].copy()
df["query_species"]   = (df["query"].str.split("_", n=1).str[1].fillna("unknown")
                          .str.lower().str.replace("pcynomologi","pcynomolgi", regex=False))
df["subject_species"] = (df["subject"].str.split("_", n=1).str[1].fillna("unknown")
                          .str.lower().str.replace("pcynomologi","pcynomolgi", regex=False))
df = df[(df["query"]!=df["subject"]) & (df["query_species"]!=df["subject_species"])].copy()
df["qcov"] = df["aln_len"]/df["qlen"]
df["scov"] = df["aln_len"]/df["slen"]
df["shorter_len"] = df[["qlen","slen"]].min(axis=1)
df["cov_shorter"] = df["aln_len"]/df["shorter_len"]
df_top5 = (df.sort_values(["query","subject","bitscore","aln_len"],
                          ascending=[True,True,False,False])
             .groupby(["query","subject"], as_index=False).head(5))
PID_THR=99; COV_SHORTER_THR=0.99; MIN_SHORTER_LEN=10
suspicious = df_top5[(df_top5["pident"]>=PID_THR)
                     & (df_top5["cov_shorter"]>=COV_SHORTER_THR)
                     & (df_top5["shorter_len"]>=MIN_SHORTER_LEN)].copy()
print("18S suspicious:", len(suspicious))
suspicious.to_csv("outputs/legacy_replicate/suspicious_cross_species_near_identity_18S.tsv",
                  sep="\t", index=False)

best = (df.sort_values(["query","bitscore","aln_len"], ascending=[True,False,False])
          .groupby("query", as_index=False).first())
best["pct_covered"] = best["cov_shorter"]*100
best["divergence_pct"] = 100 - best["pident"]
summary = best[["query","query_species","subject","subject_species","pident",
                "divergence_pct","aln_len","pct_covered","mismatch","gapopen",
                "bitscore"]].sort_values("pident", ascending=False)
summary.to_csv("outputs/legacy_replicate/closest_cross_species_match_per_sequence_18S.tsv",
               sep="\t", index=False)
PY

echo "==> [8/12] seqkit grep — target-species subsets"
SPP_PAT="vivax|falciparum|knowlesi|malariae|ovale|coat|inui|fieldi|cyno|simiovale|simium"
seqkit grep -n -r -i -p "$SPP_PAT" "$REF/18S_ref_db.fasta" \
  > "$OUT/18S_ref_db.target.fasta"
seqkit grep -r -p "$SPP_PAT" "$REF/mit_all.fasta" \
  > "$OUT/mit_all.target.fasta"

echo "==> [9/12] seqkit fx2tab — sequence lengths"
seqkit fx2tab -n -l "$OUT/mit_all.target.fasta" \
  | awk 'BEGIN{FS="\t"; OFS="\t"} {print $1,$2}' > "$OUT/mit_lengths.tsv"
seqkit fx2tab -n -l "$OUT/18S_ref_db.target.fasta" > "$OUT/18S_lengths.tsv"

echo "==> [10/12] mafft --auto"
mafft --auto "$OUT/mit_all.target.fasta"     > "$OUT/ma_mit.target.fasta"  2>"$OUT/tmp_db/mafft_mit.log"
mafft --auto "$OUT/18S_ref_db.target.fasta"  > "$OUT/ma_18S.target.fasta"  2>"$OUT/tmp_db/mafft_18S.log"

echo "==> [11/12] length classification"
python <<'PY'
import pandas as pd
for marker, infile, outfile in [
    ("MIT", "outputs/legacy_replicate/mit_lengths.tsv",
            "outputs/legacy_replicate/mit_length_classification.tsv"),
    ("18S", "outputs/legacy_replicate/18S_lengths.tsv",
            "outputs/legacy_replicate/18S_length_classification.tsv"),
]:
    df = pd.read_csv(infile, sep="\t", names=["id","len"])
    max_len = df["len"].max()
    df["pct_of_max"] = df["len"]/max_len
    df["missing_bases"] = max_len - df["len"]
    df["is_partial"] = df["pct_of_max"] < 0.90
    print(f"{marker} partial:", df["is_partial"].sum(), "of", len(df))
    df.to_csv(outfile, sep="\t", index=False)
PY

echo "==> [12/12] seq_entropy.py — per_position / windows / primer_candidates / core_bounds"
python "$LEGACY_SCRIPTS/seq_entropy.py" --aln "$OUT/ma_mit.target.fasta" --outprefix "$OUT/mit"
python "$LEGACY_SCRIPTS/seq_entropy.py" --aln "$OUT/ma_18S.target.fasta" --outprefix "$OUT/18S"

echo
echo "==> Done. Replicated outputs in $OUT/"
ls -la "$OUT/"
