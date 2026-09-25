#!/bin/bash
# Smoke: weighted parsimony shifts protein selection.
# 2 proteins: protA has unique peptide PEPA, protB has unique peptide PEPB.
# info_score min_info default = 0.1.
# Without weights: both peptides score 1.0 -> both proteins selected.
# With weights:    protB's peptide weighted 0.05 -> info=0.05 < 0.1 -> dropped.
set -euo pipefail
SCRATCH=$(mktemp -d)
echo "Retained fixture directory: $SCRATCH"

cat > "$SCRATCH/db.fasta" <<'EOF'
>protA
ACGFKAAAAAAAAAAAAAAPEPTIDEAAAAAAA
>protB
LMNQRTYAAAAAAAAAAAAAAGOODSEQEAAAAAAA
EOF

# AlphaNovo-style CSV (header: peptide, score, sequence_with_bags optional)
cat > "$SCRATCH/peptides.csv" <<'EOF'
peptide_prediction_detokenized_unmodified,score
PEPTIDEAAAAAAA,0.9
GOODSEQEAAAAAAA,0.9
EOF

cat > "$SCRATCH/weights.tsv" <<'EOF'
peptide	weight
PEPTIDEAAAAAAA	1.0
GOODSEQEAAAAAAA	0.05
EOF

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIN="${FASTALAKE_BIN_DIR:-$SCRIPT_DIR/../target/release}/parsimony_engine"

mkdir -p "$SCRATCH/out_no_w" "$SCRATCH/out_with_w"

echo "=== Run 1: info_score WITHOUT weights ==="
"$BIN" --strategy info-score --tiebreak hash-acc \
  -f "$SCRATCH/db.fasta" \
  -p "$SCRATCH/peptides.csv" \
  -o "$SCRATCH/out_no_w" \
  -s smoke_unw \
  --min-length 9 \
  --max-length 50 \
  -t 1 2>&1 | tail -3

NO_W=$(grep -c "^>" "$SCRATCH/out_no_w"/smoke_unw_info_score.fasta 2>/dev/null || echo 0)
echo "  selected (no weights): $NO_W"

echo
echo "=== Run 2: info_score WITH weights (protB peptide downweighted 0.05) ==="
"$BIN" --strategy info-score --tiebreak hash-acc \
  -f "$SCRATCH/db.fasta" \
  -p "$SCRATCH/peptides.csv" \
  -o "$SCRATCH/out_with_w" \
  -s smoke_w \
  --min-length 9 \
  --max-length 50 \
  --weights "$SCRATCH/weights.tsv" \
  -t 1 2>&1 | tail -3

WITH_W=$(grep -c "^>" "$SCRATCH/out_with_w"/smoke_w_info_score.fasta 2>/dev/null || echo 0)
echo "  selected (with weights): $WITH_W"

echo
echo "=== ASSERTIONS ==="
fail=0
[[ $NO_W -eq 2 ]] || { echo "FAIL: expected 2 selected without weights, got $NO_W"; fail=1; }
[[ $WITH_W -eq 1 ]] || { echo "FAIL: expected 1 selected with weights (protA only), got $WITH_W"; fail=1; }

if grep -q "^>protA" "$SCRATCH/out_with_w"/smoke_w_info_score.fasta 2>/dev/null \
   && ! grep -q "^>protB" "$SCRATCH/out_with_w"/smoke_w_info_score.fasta 2>/dev/null; then
  :
else
  echo "FAIL: with-weights output should contain protA but not protB"
  cat "$SCRATCH/out_with_w"/smoke_w_info_score.fasta 2>/dev/null
  fail=1
fi

[[ $fail -eq 0 ]] && echo "ALL CHECKS PASSED" || { echo "TEST FAILED"; exit 1; }
