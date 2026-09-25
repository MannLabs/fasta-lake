# shared_peptide_classifier

Projects k-mer-rescued peptides from cluster representatives onto member proteins and classifies each UNIQUE or SHARED (>=80% of members) — the conservation call `parsimony_engine --weights` consumes.

FastaLake v5.4 helper. Given peptides anchored to MMseqs2 cluster representatives by `mass_kmer_anchor` (or any rescue tool), lift each peptide back to the unclustered member proteins that actually contain it. Tags each `(peptide, rep)` pair as **SHARED** or **UNIQUE** based on the fraction of cluster members in which the peptide appears. Downstream, `parsimony_engine --weights` consumes a per-peptide weight derived from this output to penalise non-conserved evidence.

## Inputs / Outputs

| Flag | Schema |
|---|---|
| `--rescued` | TSV with header. Must contain columns named by `--peptide-col` (default `peptide`) and `--rep-col` (default `cluster_rep`). Extra columns ignored. |
| `--cluster-tsv` | 2-column TSV `representative<TAB>member` (one row per member; reps repeat). MMseqs2 `createtsv` output is the canonical format. |
| `--lake` | Unclustered FASTA. Only the first whitespace-delimited token of the header is used as the protein id. |
| `--output` | TSV with columns: `peptide, cluster_rep, n_cluster_members, n_with_hit, members_with_hit, is_shared`. |

## Build

```bash
cargo build --release --manifest-path rust/shared_peptide_classifier/Cargo.toml
```

## Usage

```bash
./target/release/shared_peptide_classifier \
    --rescued     results/rescued_peptides.tsv \
    --cluster-tsv results/lake_cluster.tsv \
    --lake        results/lake_unclustered.fasta \
    --output      results/lifted_peptides.tsv \
    --shared-threshold 0.8 \
    --threads 8
```

`--shared-threshold` (default 0.8) is the fraction-of-members hit at or above which a peptide is tagged `is_shared=true`. Below threshold the peptide is `is_shared=false` (UNIQUE), and the matching members are listed in `members_with_hit`.

## Tests

```bash
cargo build --release --manifest-path rust/shared_peptide_classifier/Cargo.toml
bash rust/shared_peptide_classifier/tests/fixture_smoke.sh
```

The fixture builds a 3-protein cluster (1 conserved peptide → SHARED, 1 divergent peptide → UNIQUE), a 2-protein cluster at 50% hit rate (UNIQUE), and a missing cluster (peptide silently dropped), then asserts the expected output rows.

## Notes

- The lake is streamed once; only sequences belonging to clusters referenced by `--rescued` are materialised in memory.
- Sequences are uppercased and stripped to ASCII alphabetic on load — peptide search is exact substring via `memchr::memmem`.
- Peptides anchored to a cluster_rep that is missing from `--cluster-tsv` are silently dropped (warning lines to stderr show the count).
- See FastaLake CLAUDE.md for the v5.4 weighted-parsimony pipeline this feeds into; the per-peptide weight TSV is derived from `n_with_hit / n_cluster_members` and consumed by `parsimony_engine --weights`.
