#!/usr/bin/env bash
# FastaLake demo — the whole idea in about a minute, on a laptop.
#
# Builds a sequence lake, filters it to peptide-supported sequences, then curates
# a SEPARATE database for each of three samples from their own de novo evidence.
#
#   ./demo/run_demo.sh
#
# No cluster, no SLURM, no network, no Python, no MMseqs2, no SAGE.
# Set MMSEQS=/path/to/mmseqs to additionally run the optional clustering stage.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DEMO="$ROOT/demo"
OUT="${1:-$DEMO/output_$(date -u +%Y%m%dT%H%M%SZ)}"
BIN_FE="${FASTALAKE_BIN_DIR:-$ROOT/rust/fasta_extractor_v2/target/release}"
BIN_PE="${FASTALAKE_BIN_DIR:-$ROOT/rust/parsimony_engine/target/release}"

rule() { printf '\n\033[1m%s\033[0m\n' "$*"; }

rule "FastaLake demo"
echo "repo: $ROOT"

# ---------------------------------------------------------------- build
if [ ! -x "$BIN_FE/lake_builder" ] || [ ! -x "$BIN_PE/parsimony_engine" ]; then
  rule "Building the two required crates (first run only, ~1-2 min)"
  command -v cargo >/dev/null || { echo "ERROR: cargo not found. Install Rust: https://rustup.rs"; exit 1; }
  cargo build --release --locked --quiet --manifest-path "$ROOT/rust/fasta_extractor_v2/Cargo.toml"
  cargo build --release --locked --quiet --manifest-path "$ROOT/rust/parsimony_engine/Cargo.toml"
fi

[ ! -e "$OUT" ] || { echo "ERROR: output exists: $OUT; choose a new directory" >&2; exit 2; }
mkdir -p "$OUT"

# ---------------------------------------------------------------- stage 0
rule "Stage 0 - build the sequence lake"
echo "Pools sources into one FASTA and gives every unique sequence a SHA-256 name,"
echo "so the same sequence gets the same identifier no matter which source it came from."
"$BIN_FE/lake_builder" \
  --source "DEMO:1:$DEMO/lake/reference_lake.fasta" \
  --uniform-tag DEMO --track-all-sources \
  --output "$OUT/lake.fasta" \
  --output-headers "$OUT/lake.joint_headers.tsv" \
  --output-stats "$OUT/lake.stats.txt" >/dev/null
echo "  lake: $(grep -c '^>' "$OUT/lake.fasta") unique sequences"
# The lake keys every sequence by its hash, so an accession no longer says which
# catalogue it came from. The joint-header table still does: derive the
# accession -> catalogue map from it, so the per-sample breakdown below is a real
# composition instead of an honest but useless UNCLASSIFIED.
awk -F'\t' 'NR > 1 {
  # all_source_headers = "TAG header|||TAG header|||..."; classify the first
  # (highest-priority) source by the original accession after its tag.
  # The separator is a regex in every awk (mawk in the container, gawk on
  # workstations); "|||" as a plain string is not a valid regex in mawk.
  split($4, entries, /[|][|][|]/); orig = entries[1]; sub(/^[^ ]+ /, "", orig)
  src = "UNCLASSIFIED"
  if (orig ~ /^(sp|tr)\|/)         src = "Human"
  else if (orig ~ /^MGYG/)         src = "UHGG"
  else if (orig ~ /^GMGC/)         src = "GMGC"
  else if (orig ~ /^GMSC/)         src = "smORF"
  else if (orig ~ /^SMORF/)        src = "smORF"
  else if (orig ~ /Assembly_/)     src = "Assembly"
  print $3 "\t" src
}' "$OUT/lake.joint_headers.tsv" > "$OUT/source_map.tsv"
echo "  source map: $(wc -l < "$OUT/source_map.tsv") accessions -> catalogue"

# ---------------------------------------------------------------- stage 2
rule "Stage 2 - evidence filter (cohort level)"
echo "Keeps only lake sequences that actually contain an observed de novo peptide."
"$BIN_FE/fasta_extractor" \
  --database "$OUT/lake.fasta" \
  --peptides-dir "$DEMO/predictions" \
  --output "$OUT/evidence_lake.fasta" \
  --output-hits "$OUT/hits.tsv" \
  --min-length 8 >/dev/null 2>&1
echo "  evidence lake: $(grep -c '^>' "$OUT/evidence_lake.fasta") of $(grep -c '^>' "$OUT/lake.fasta") sequences survived"

# ---------------------------------------------------------------- stage 3 (optional)
DB_FOR_SAMPLES="$OUT/evidence_lake.fasta"
if [ -n "${MMSEQS:-}" ] && command -v "$MMSEQS" >/dev/null; then
  rule "Stage 3 - cluster97 (optional, MMseqs2 detected)"
  "$MMSEQS" easy-cluster "$OUT/evidence_lake.fasta" "$OUT/clu97" "$OUT/tmp" \
    --min-seq-id 0.97 -c 0.8 --cov-mode 1 >/dev/null 2>&1
  DB_FOR_SAMPLES="$OUT/clu97_rep_seq.fasta"
  echo "  clustered to $(grep -c '^>' "$DB_FOR_SAMPLES") representatives"
else
  rule "Stage 3 - cluster97 (SKIPPED)"
  echo "  MMseqs2 not requested. Set MMSEQS=/path/to/mmseqs to enable."
fi

# ---------------------------------------------------------------- stages 4-5
rule "Stages 4-5 - per-sample curation + parsimony"
echo "This is the point of FastaLake: each sample gets its OWN database,"
echo "built from its OWN evidence, not one shared database for everybody."
printf '\n  %-10s %10s %10s %10s\n' SAMPLE PEPTIDES "DB(seqs)" "PARSIMONY"
printf '  %-10s %10s %10s %10s\n' ---------- -------- -------- ---------
for p in "$DEMO"/predictions/*_predictions.csv; do
  s=$(basename "$p" _predictions.csv)
  mkdir -p "$OUT/$s/peps"; cp "$p" "$OUT/$s/peps/"
  "$BIN_FE/fasta_extractor" --database "$DB_FOR_SAMPLES" --peptides-dir "$OUT/$s/peps" \
      --output "$OUT/$s/${s}_db.fasta" --min-length 8 >/dev/null 2>&1
  "$BIN_PE/parsimony_engine" --fasta "$OUT/$s/${s}_db.fasta" --peptides "$p" \
      --sample "$s" --output-dir "$OUT/$s/inference" --strategy razor --tiebreak hash-acc \
      --min-length 8 --source-map "$OUT/source_map.tsv" >/dev/null
  npep=$(( $(wc -l < "$p") - 1 )); ndb=$(grep -c '^>' "$OUT/$s/${s}_db.fasta" 2>/dev/null || echo 0)
  npars=$(grep -c '^>' "$OUT/$s/inference/${s}_razor.fasta")
  printf '  %-10s %10s %10s %10s\n' "$s" "$npep" "$ndb" "$npars"
done

rule "Done"
echo "Outputs in $OUT. Each SAMPLE_*/ holds a database curated for that sample alone."
echo "Compare them - they are not the same file, and that is the whole idea:"
echo "the database becomes a per-sample measurement instead of a fixed assumption."
