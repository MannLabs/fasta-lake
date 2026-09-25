# lake_stats

Compute the **D/K saturation** of tryptic peptides in a FASTA database — the fraction of all possible length-N peptides (K = 20^N) that actually appear (D) after in-silico tryptic digestion. Used to generate Extended Data Figure 1 in the FastaLake paper. Mirrors `sequence_space_saturation.py` in parallel Rust (~50–100× speedup over the Python reference).

## Inputs / Outputs

**Input** (`-f / --fasta`): any FASTA. Headers are ignored; sequences are uppercased and stripped to ASCII alphabetic.

**Tryptic rule:** cleave after K/R unless followed by P. The last fragment is always emitted.

**Filter:** drop peptides containing any of `B J O U X Z * -` (ambiguous / non-standard).

**Output** (`-o / --output`): TSV with one row per requested length:

```
database	min_length	K_theoretical_20pow_N	D_unique_peptides_observed	saturation_percent	n_proteins	elapsed_seconds
lake	6	64000000	5821334	9.095834	273000000	412.3
```

## Build

```bash
cargo build --release --manifest-path rust/lake_stats/Cargo.toml
```

## Usage

```bash
# Default: lengths 6,7,8,9,10,11,12
./target/release/lake_stats \
    -f results/lake.fasta \
    -o results/lake_stats.tsv \
    -n MicrobPredict_lake \
    -t 32

# Just the lengths SAGE actually digests by default (≥9)
./target/release/lake_stats \
    -f results/lake.fasta \
    -o results/lake_stats_sage.tsv \
    -l 9,10,11,12,13,14 \
    -n MicrobPredict_lake \
    -t 32
```

## Tests

```bash
cargo test --manifest-path rust/lake_stats/Cargo.toml --release
```

`tests/golden.rs` builds a 2-protein FASTA with hand-enumerated tryptic cuts and asserts the resulting D values exactly: e.g. `MAKPEPTIDERWAYR` digests to `MAKPEPTIDER` (the K@3 isn't cut because it's followed by P) and `WAYR`.

## Notes

- D values are reported per length, not cumulatively — `D` for length 9 does not include length-9 substrings of longer peptides.
- The Python reference (`sequence_space_saturation.py`) is the ground truth; the Rust binary's golden test is verified against hand-enumerated truth, and full-FASTA results match the Python reference bit-for-bit on the small fixtures the paper uses.
- For typical reference-lake-scale inputs (200 M+ proteins) memory grows with D — set `-l` conservatively. Length 14 on a 273 M-protein lake is the practical ceiling on a 760 GB node.
- `K = 20^N` undercounts the true theoretical space slightly (ignores non-canonical residues), but this is the convention used in the saturation literature and in the paper figures.
