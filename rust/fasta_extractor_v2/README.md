# fasta_extractor_v2

The production peptide-to-protein toolkit for FastaLake. One crate, six binaries: the main `fasta_extractor` (Stage 1 + Stage 2 Aho–Corasick extraction), the `lake_builder` (Stage 0 SHA-256 deduplicated reference lake), the cross-language `hash_fasta` checker, mass-tolerant rescue (`mass_kmer_anchor`, `mass_corasick_v3`), and a standalone Bayesian/info-score `parsimony` binary. All binaries memory-map their FASTAs and parallelise with Rayon.

## Binaries

| Binary | Stage | Purpose |
|---|---|---|
| `fasta_extractor` | 1 / 2 | Aho–Corasick peptide → protein extraction with I/L normalization, top-X% per-sample score filter, hits/evidence TSVs, optional include FASTA, optional entrapment spike-in (shuffled or noble). Auto-detects AlphaNovo / DIANovo / Casanovo / plain-text input. |
| `lake_builder` | 0 | Two-phase reference-lake builder. Streams sources, SHA-256 hashes, dedups, writes one record per unique sequence keeping the highest-priority source's header. Joint-header TSV records all collisions. |
| `hash_fasta` | — | Emits `sha256_hex<TAB>accession<TAB>normalized_length` per sequence. Byte-identical to `fasta_lake.hashing.sequence_hash`; pinned by `tests/test_hashing_crosslang.py`. |
| `mass_kmer_anchor` | rescue | 5-mer mass-window matching for peptides that fail exact match (de novo errors). Indexes 0.001 Da bins; min-vote thresholded; every candidate verified at ±20 ppm against the measured precursor mass (`prec_mz`/`prec_charge`) or, for the two-column input FastaLake uses, the theoretical mass of the de novo sequence (`mass_source` column says which). Opt-in. |
| `mass_corasick_v3` | rescue | Per-sample parallel rescue with `prep` / `rescue` / `merge` subcommands. Built for SLURM arrays — 306 samples × 273 M lake in ~45 min wall. |
| `parsimony` | 3 | Standalone EM-based protein parsimony (Bayesian / info_score / razor). Superseded by the dedicated `parsimony_engine` crate for cohort-scale runs but kept here for dependency parity. |

## Build

```bash
cargo build --release --manifest-path rust/fasta_extractor_v2/Cargo.toml
```

All six binaries land in `target/release/`.

## Usage

```bash
# Stage 0: build a deduplicated reference lake from 5 sources
./target/release/lake_builder \
    -s "PRODIGAL:1:data/assemblies:*_proteins.faa" \
    -s "UNIPROT_HUMAN:2:data/uniprot_human.fasta" \
    -s "UHGP:3:data/uhgp-100.faa" \
    -o results/lake.fasta \
    --output-headers results/lake.joint_headers.tsv \
    --output-stats results/lake.stats.txt \
    -t 32

# Stage 1: cohort evidence lake (pooled predictions vs. reference lake)
./target/release/fasta_extractor \
    -p data/predictions/ \
    -d results/lake.fasta \
    -o results/evidence_lake.fasta \
    -s FASTALAKE \
    --min-length 9 --max-length 50 \
    --top-percent 0.30 \
    --output-hits results/hits.tsv \
    -t 24

# Rescue: mass-tolerant 5-mer matching for peptides that failed exact match
./target/release/mass_kmer_anchor \
    --database results/evidence_lake.fasta \
    --predictions data/failed_peptides.csv \
    --output results/rescued.tsv \
    --min-votes 5 --skip-exact -t 4

# Cross-language hash check (used by Python<->Rust dedup tests)
./target/release/hash_fasta -i results/lake.fasta -o results/lake.hashes.tsv
```

## Tests

```bash
cargo test --manifest-path rust/fasta_extractor_v2/Cargo.toml --release
```

Integration tests in `tests/integration_test.rs` cover I/L normalization, top-percent filter, entrapment, hits/evidence output, and round-tripping. The cross-language hash test lives at `fasta_lake/tests/test_hashing_crosslang.py` and invokes `hash_fasta` against the Python reference.

## Notes

- `lake_builder` exits with status **141 (SIGPIPE)** on a successful run when downstream consumers close pipes early; SLURM chains use `afterany` (not `afterok`) per the FastaLake CLAUDE.md memory `feedback_afterany_not_afterok.md`.
- v5.3 and later use **SHA-256** (not MD5) for sequence dedup. Pre-v5.3 lakes must be rebuilt before mixing with current outputs.
- The `--min-length` default was raised from 8 to 9 in v5; SAGE-tryptic digest produces peptides ≥ 9 and shorter de novo peptides admit low-evidence proteins that inflate entrapment FDP. See `DESIGN_V5.md`.
- `mass_kmer_anchor` and `mass_corasick_v3` are **opt-in**. K-mer rescue is DVP-specific and can inflate FDP in other contexts; the FastaLake paper benchmarks the WITH/WITHOUT k-mer pipeline pair on every dataset (memory `project_kmer_rescue_comparison.md`).
- Hashing contract (uppercased A–Z only; drops whitespace, digits, stop codons, punctuation; preserves B/Z/X/U/O) is shared byte-for-byte with `fasta_lake/hashing.py`.
