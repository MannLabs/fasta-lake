# fasta_export

Engine-aware FASTA writer for FastaLake. Heals headers from heterogeneous source databases (UniProt sp/tr, UHGP MGYG, GMGC, GENCODE, Ensembl, NCBI, Prodigal, gnomAD variants) into the form each search engine expects, and computes entrapment FDP from SAGE results. Two subcommands: `export` and `entrapment`.

## Inputs / Outputs

**`export` subcommand**
- Input: any FASTA. Headers are parsed; `GN=` is honoured; gene name is fallback-extracted from the entry name (sp/tr), gnomAD parent accession, or the GMGC tail token.
- Output: FASTA with healed headers. For `--engine msfragger`, target sequences are followed by reversed decoys with prefix `rev_` (one decoy per target).

**`entrapment` subcommand**
- Input: directory of per-sample subdirectories each containing `results.sage.tsv`.
- Output: TSV to stdout with `sample, n_target, n_entrapment, fdp_lower, fdp_combined`. `fdp_combined` uses Wen, Keich & Noble 2025 Eq 1: `FDP = N_E × (1 + 1/r) / (N_T + N_E)` with `r = entrap_db_size / target_db_size`.

## Build

```bash
cargo build --release --manifest-path rust/fasta_export/Cargo.toml
```

## Usage

```bash
# Heal headers for SAGE: ensures GN= present; tolerates non-sp/tr accessions
./target/release/fasta_export export \
    -i results/lake.fasta \
    -o results/lake.sage.fasta \
    --engine sage

# Heal for DIA-NN: for non-sp/tr proteins, gene must be the FIRST WORD of the description
./target/release/fasta_export export \
    -i results/lake.fasta \
    -o results/lake.diann.fasta \
    --engine diann

# Heal for MSFragger and append reversed rev_ decoys
./target/release/fasta_export export \
    -i results/lake.fasta \
    -o results/lake.msfragger.fasta \
    --engine msfragger

# Entrapment FDP from per-sample SAGE results
./target/release/fasta_export entrapment \
    -i results/sage_per_sample/ \
    --prefix Entrap_ \
    --target-db-size 8123456 \
    --entrap-db-size 5000 \
    -q 0.01 \
    > results/entrapment_fdp.tsv
```

## Tests

```bash
cargo test --manifest-path rust/fasta_export/Cargo.toml --release
```

The unit tests in `src/main.rs` cover header parsing across all 9 supported formats and the per-engine healing rules.

## Notes

- DIA-NN ignores `GN=` for non-sp/tr accessions (UHGP/GMGC/Prodigal/etc.) — the gene name must be the first whitespace-delimited token of the description. The `diann` engine path enforces this; SAGE/MSFragger paths only ensure `GN=` is present. See FastaLake CLAUDE.md memory `project_header_fix.md` for why this matters: every gene-level result computed before this fix is wrong.
- Decoy generation reverses sequence characters in place and prefixes the accession with `rev_`. SAGE and DIA-NN generate their own decoys internally, so the `msfragger` path is the only one that emits them.
- Entrapment counting reads `proteins, spectrum_q, label` from `results.sage.tsv`; rows with `label = -1` (decoys) are skipped. A peptide is entrapment if any group member starts with `--prefix`.
- `fdp_lower = N_E / (N_T + N_E)` is reported alongside `fdp_combined` but is only useful as a failure check (it cannot prove FDR control).
