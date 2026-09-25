# aggregation_engine

Stage 5 cross-sample aggregation for FastaLake metaproteomics. Reads per-sample SAGE/DIA-NN output and emits dual-index matrices (peptide × sample, protein-group × sample) plus optional functional roll-ups (eggNOG: KEGG, COG, CAZy, Pfam, GOs, BRITE, EC). Protein-group keys are canonicalised — sorted, deduplicated, source-tag-stripped, semicolon-joined — so the same group is identical across samples and reruns. Each group is classified as Human / Microbial / Mixed.

## Inputs / Outputs

**Input directory layout** (`-i`): one subdirectory per sample, each containing either
- `lfq.tsv` (DIA-NN-style LFQ: needs `peptide`, `proteins`, `q_value`, `*.mzML` intensity column), or
- `results.sage.tsv` (SAGE PSM: needs `peptide`, `proteins`, `spectrum_q`, `peptide_q`, `protein_q`, `ms2_intensity`, `label`).

**Outputs** (`-o`): `{dataset}_*` matrices at q ≤ 0.01, q ≤ 0.05, and unfiltered. PSM mode repeats the matrix set for `spectrum_q`, `peptide_q`, `protein_q`. Also written: `{dataset}_protein_group_metadata.tsv` (group key, classification, representative member) and a per-sample summary TSV. With `--functional-lookup`, additional `{dataset}_functional_{level}_{plurality|splitintensity}.tsv` matrices appear for each requested level.

**Functional lookup TSV** (`--functional-lookup`): tab-separated, must include a column named `accession` plus any of the requested level columns. Sentinel values (`-`, `--`, empty, `NA`, `None`) are dropped. Multi-value cells split on `;` or `,`.

## Build

```bash
cargo build --release --manifest-path rust/aggregation_engine/Cargo.toml
```

## Usage

```bash
# LFQ aggregation with eggNOG functional roll-up
./target/release/aggregation_engine \
    -i  results/sage_per_sample/ \
    -o  results/aggregated/ \
    -d  MicrobPredict \
    -m  FastaLake \
    --data-type lfq \
    --functional-lookup data/eggnog_emapper.tsv \
    -t 16

# PSM aggregation, all q levels, no functional
./target/release/aggregation_engine \
    -i  results/sage_per_sample/ \
    -o  results/aggregated/ \
    -d  CAPSCAN \
    --data-type psm \
    --no-functional \
    -t 16
```

## Tests

```bash
cargo test --manifest-path rust/aggregation_engine/Cargo.toml --release
```

## Notes

- Default functional levels = `OG,KEGG_ko,KEGG_Pathway,KEGG_Module,KEGG_Reaction,COG_category,CAZy,Pfam,GOs,BRITE,EC`. Pass `--functional-levels ""` to suppress, or `--no-functional` to disable entirely.
- Two functional aggregation modes always run side by side: `plurality` (assign full group intensity to plurality-winner term — overcounts by design, conservative) and `splitintensity` (split equally across N unique values members point to — mass-balanced).
- `--bh-correction` is **retired** and is a hard error. SAGE's q-values are already FDR-controlled (target–decoy), `lfq.tsv` is written only after that filtering, and aggregation computes no new p-value — so there is no multiple-testing family to correct. Use the entrapment FDP (`fasta_export entrapment-fdp`) for cohort-level error control.
- Source-tag stripping: a leading `TAG|PRIORITY|...` form (uppercase letters, digits, underscores) is removed so that `FASTALAKE|3|sp|P02768|ALBU_HUMAN` and `sp|P02768|ALBU_HUMAN` collapse to the same group member.
