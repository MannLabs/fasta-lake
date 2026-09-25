#!/bin/bash
# Smoke test for shared_peptide_classifier on a hand-crafted fixture.
#
# Cluster A has 3 members. Peptide P1 is in the conserved region (all 3 hit) -> SHARED.
# Peptide P2 is in member1's divergent region only (1 hit) -> UNIQUE.
# Cluster B has 2 members. Peptide P3 is in member4 only (1/2 = 50% < 80%) -> UNIQUE.
# Cluster C is referenced but missing from cluster.tsv -> peptide P4 dropped silently.

set -euo pipefail
SCRATCH=$(mktemp -d)
echo "Retained fixture directory: $SCRATCH"

cat > "$SCRATCH/lake.fasta" <<'EOF'
>member1
AAACCCDDDEEEPEPTIDEAAA
GGGFFFFKKKKLLLLM
>member2
AAACCCDDDEEEXXXAAA
GGGFFFFKKKKLLLLM
>member3
AAACCCDDDEEEYYYAAA
GGGFFFFKKKKLLLLM
>member4
NNNNHHHHWWWWQQQQRARERAREPEPSEQK
>member5
ZZZZZZZZZZ
EOF

cat > "$SCRATCH/cluster.tsv" <<'EOF'
clusterA	member1
clusterA	member2
clusterA	member3
clusterB	member4
clusterB	member5
EOF

cat > "$SCRATCH/rescued.tsv" <<'EOF'
peptide	cluster_rep	score
GGGFFFFK	clusterA	0.9
PEPTIDE	clusterA	0.85
PEPSEQK	clusterB	0.8
GHOSTPEP	clusterC	0.7
EOF

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIN="${FASTALAKE_BIN_DIR:-$SCRIPT_DIR/../target/release}/shared_peptide_classifier"
"$BIN" \
  --rescued "$SCRATCH/rescued.tsv" \
  --cluster-tsv "$SCRATCH/cluster.tsv" \
  --lake "$SCRATCH/lake.fasta" \
  --output "$SCRATCH/lifted.tsv" \
  --shared-threshold 0.8 \
  --threads 2

echo
echo "=== EXPECTED ==="
echo "GGGFFFFK     clusterA  3  3  member1,member2,member3  true   (shared, conserved region)"
echo "PEPTIDE      clusterA  3  1  member1                  false  (unique to member1)"
echo "PEPSEQK      clusterB  2  1  member4                  false  (1/2 = 50% < 80%)"
echo "GHOSTPEP     clusterC  -- dropped (cluster not in cluster.tsv)"
echo
echo "=== ACTUAL ==="
cat "$SCRATCH/lifted.tsv"

# Assertions
fail=0
grep -P '^GGGFFFFK\tclusterA\t3\t3\t' "$SCRATCH/lifted.tsv" >/dev/null \
  && grep -P 'true$' "$SCRATCH/lifted.tsv" | grep -P '^GGGFFFFK\t' >/dev/null \
  || { echo "FAIL: GGGFFFFK expected SHARED on all 3 members"; fail=1; }
grep -P '^PEPTIDE\tclusterA\t3\t1\tmember1\tfalse$' "$SCRATCH/lifted.tsv" >/dev/null \
  || { echo "FAIL: PEPTIDE expected UNIQUE on member1 only"; fail=1; }
grep -P '^PEPSEQK\tclusterB\t2\t1\tmember4\tfalse$' "$SCRATCH/lifted.tsv" >/dev/null \
  || { echo "FAIL: PEPSEQK expected UNIQUE on member4 only"; fail=1; }
if grep -q '^GHOSTPEP' "$SCRATCH/lifted.tsv"; then
  echo "FAIL: GHOSTPEP should have been dropped (no cluster)"; fail=1
fi

[[ $fail -eq 0 ]] && echo "ALL CHECKS PASSED" || { echo "TEST FAILED"; exit 1; }
