# parsimony_engine

Stage 3 protein inference for FastaLake. Six strategies — `razor`, `species_budget`, `uniform`, `uniform_2pep`, `info_score`, `bayesian` — across the same peptide-to-protein mapping. Replaces the Python parsimony for cohort-scale runs: razor on a 3.16 M-protein Stage 2 input finishes in ≈39 s on a single core. Optionally consumes per-peptide weights from `shared_peptide_classifier` to penalise non-conserved cluster evidence.

## Inputs / Outputs

**Inputs**
- `-f / --fasta` — Stage 2 sample-specific FASTA.
- `-p / --peptides` — de novo predictions CSV (AlphaNovo / DIANovo / Casanovo / plain text, auto-detected).
- `--weights` *(optional)* — TSV with header `peptide<TAB>weight`. Only `info_score` and `bayesian` use it; `razor` / `uniform*` ignore it. Missing peptides default to weight 1.0.

**Outputs** (`-o / --output-dir`):
- `{sample}_{strategy}.fasta` — the inferred protein subset.
- `{sample}_{strategy}_stats.json` — per-strategy diagnostics (`stage2_proteins`, `proteins_with_evidence`, `selected_proteins`, `reduction_factor`, `unique_peptides`, `shared_peptides`, `db_breakdown` per source DB, plus species-budget extras when applicable).
- `{sample}_{strategy}_provenance.jsonl` — one JSON record per selected protein listing its supporting peptides, scores, exact-vs-rescued admission, and contaminant flags.

## Build

```bash
cargo build --release --manifest-path rust/parsimony_engine/Cargo.toml
```

## Usage

```bash
# Razor inference, top-30% peptides, sample S1
./target/release/parsimony_engine \
    --strategy razor \
    -f results/stage2/S1.fasta \
    -p data/predictions/S1.csv \
    -o results/parsimony/S1/ \
    -s S1 \
    --min-length 9 --max-length 50 \
    --top-percent 0.30 \
    --threads 1

# v5.4 weighted info_score: shared_peptide_classifier-derived weights drop non-conserved evidence
./target/release/parsimony_engine \
    --strategy info-score \
    -f results/stage2/S1.fasta \
    -p data/predictions/S1.csv \
    --weights results/peptide_weights.tsv \
    -o results/parsimony_weighted/S1/ \
    -s S1 \
    --threads 4
```

## Tests

```bash
# Rust unit/integration tests
cargo test --manifest-path rust/parsimony_engine/Cargo.toml --release

# Smoke test: --weights demonstrably shifts protein selection
cargo build --release --manifest-path rust/parsimony_engine/Cargo.toml
bash rust/parsimony_engine/tests/fixture_weights.sh
```

The `fixture_weights.sh` smoke test builds two proteins each carrying a unique peptide. Without weights, `info_score` selects both. With one peptide downweighted to 0.05 (below the default `min_info=0.1`), only the un-downweighted protein is selected.

## Notes

- `--weights` is gated on the upstream `shared_peptide_classifier` output: derive a per-peptide weight from `n_with_hit / n_cluster_members` and pass the resulting TSV here. This is the v5.4 cluster-aware inference path described in the FastaLake CLAUDE.md session memories.
- The `--min-length` default is **9** (matches the Python wrapper and the v5.3 entrapment-FDR fix). Earlier versions defaulted to 7; do not lower it without re-tuning entrapment.
- I/L normalisation is on by default (`--normalize-il true`); peptides and protein sequences are folded I → L before Aho–Corasick matching.
- Strategy semantics, tie-breaking rules, and the Bayesian EM posterior are formalised in §5.2 of the FastaLake Methods paper. The `db_breakdown` block in `stats.json` reports unique / shared-only / total counts per source DB (FASTALAKE, sp, tr, MGYG, GMGC, etc.).
- For peptide–protein mapping the engine uses the same Aho–Corasick automaton as `fasta_extractor_v2`, so identical input gives identical peptide → protein assignments across the two binaries.
