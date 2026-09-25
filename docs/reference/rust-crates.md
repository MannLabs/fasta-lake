# Rust crates and executables

The compute-heavy stages run in six Rust crates under `rust/`. Each crate has its own README with the exact command lines, and `--help` on every executable is the authoritative option list. The Python package finds the executables through an explicit path, `FASTALAKE_BIN_DIR`, the repository build directories, and then `PATH` (see [Environment](../get-started/environment.md)).

| Crate | Executables | Stage | README |
|---|---|---|---|
| `fasta_extractor_v2` | `fasta_extractor`, `lake_builder`, `hash_fasta`, `fasta_sort_by_id`, `mass_kmer_anchor`, `mass_corasick`, `mass_corasick_v2`, `mass_corasick_v3`, `parsimony` | Reference lake (stage 0), pooled and per-acquisition evidence extraction (stages 1 and 2), mass-tolerant rescue | [fasta_extractor_v2](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/README.md), [parameters](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/README_parameters.md) |
| `parsimony_engine` | `parsimony_engine` | Protein inference (stage 3): razor and the alternative strategies | [parsimony_engine](https://github.com/MannLabs/fasta-lake/blob/main/rust/parsimony_engine/README.md) |
| `aggregation_engine` | `aggregation_engine` | Cross-acquisition matrices from Sage output (stage 5) | [aggregation_engine](https://github.com/MannLabs/fasta-lake/blob/main/rust/aggregation_engine/README.md) |
| `fasta_export` | `fasta_export` | Engine-aware FASTA export and entrapment FDP from Sage results | [fasta_export](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_export/README.md) |
| `shared_peptide_classifier` | `shared_peptide_classifier` | UNIQUE or SHARED calls for rescued peptides, consumed by `parsimony_engine --weights` | [shared_peptide_classifier](https://github.com/MannLabs/fasta-lake/blob/main/rust/shared_peptide_classifier/README.md) |
| `lake_stats` | `lake_stats` | Tryptic peptide-space saturation of a FASTA | [lake_stats](https://github.com/MannLabs/fasta-lake/blob/main/rust/lake_stats/README.md) |

## Build

The toolchain is pinned in `rust-toolchain.toml`; every crate has a `Cargo.lock`. The three crates the study runner needs:

```bash
for crate in fasta_extractor_v2 parsimony_engine aggregation_engine; do
  cargo build --release --locked --manifest-path "rust/$crate/Cargo.toml"
done
```

Continuous integration runs `cargo fmt --check`, `cargo clippy -D warnings`, `cargo test` and `cargo deny` for all six crates on every push. The Docker image contains release builds of all of them.

## Determinism

Every executable writes its outputs in a fixed order, independent of thread count and directory listing order, so identical inputs give byte-identical files. The [algorithms](../concepts/algorithms.md) page states the tie rules that make this possible.
