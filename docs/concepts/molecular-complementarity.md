> For the current additive metaG/metaT workflow, portable v2 sources and explicit available-reference mode, start with [Molecular inputs](../how-to/molecular-inputs.md). The original v1 complete-reference contract described below remains supported.

# Molecular evidence and de novo complementarity

These optional tools add specimen-matched molecular evidence after FastaLake's
final de novo protein selection and compare the resulting Sage searches.
The normal de novo workflow is unchanged. They support DNA-only inputs and
paired DNA/RNA inputs represented against one verified gene reference.

This is a supported optional workflow. Start with [how to supply your own
specimen FASTA and TPM files](../how-to/molecular-inputs.md). The [completed DNA benchmark](molecular-benchmark.md)
documents measured gains and losses; cutoff choice and whole-workflow error
calibration remain separate scientific questions.

## What to compare

Keep sequence availability, sequence selection, identification and quantification
as separate measurements. For the same acquisition, compare a complete matched
metagenome search, de novo-guided FastaLake, and FastaLake with TPM additions.
The complete metagenome comparator and a TPM-filtered metagenome are different
conditions. Preserve the same host/contaminant policy, spectra, Sage executable,
search settings and quantitative settings between the arms.

An experimental 2x2 crosses metagenome sequence input with FastaLake
curation. An outcome-overlap 2x2 instead asks whether a protein sequence has
accepted evidence in each search. An unrun experimental cell is not an empty
outcome category. The outcome category “neither” is defined only within a declared
sequence universe.

Report metagenome-only and FastaLake-only accepted peptides in both directions.
Then measure which previously missed observations are recovered by TPM addition,
including any losses among the original FastaLake observations. Adding a FASTA
record is not evidence that the protein was identified.

## Add exact molecular sequences

Install FastaLake as described in the main README, then run:

```bash
fasta-lake molecular add \
  --baseline final_razor.fasta \
  --source specimen_source.json --specimen specimen_A \
  --dna-percent 1 --rna-absolute 1 --out union_top1_rna1
```

Omit `--rna-absolute` to disable RNA selection, including for DNA-only datasets.
The DNA rule uses a percentage of positive original genes, before exact
sequence deduplication. The historical default rank is
`max(1, floor(n_positive * percent / 100))`; all boundary ties pass.
`--rank-rounding ceil` is an explicit alternative. Zero percent disables the
DNA percentile rule. For an absolute DNA threshold use, for example,
`--dna-percent 0 --dna-absolute 3`. Absolute DNA and RNA thresholds both use
strict `>`; percentile selection uses `>=` at the selected positive cutoff.

The source manifest has this shape; replace each digest with the actual SHA-256:

```json
{
  "schema": "fastalake.molecular-source.v1",
  "specimen": "specimen_A",
  "reference_id": "assembly_A_version_1",
  "reference_status": "verified_original_quantification_reference",
  "tpm": {"path": "genes_tpm.tsv", "sha256": "ACTUAL_SHA256"},
  "proteins": {"path": "original_genes.faa.gz", "sha256": "ACTUAL_SHA256"},
  "provenance": {"path": "reference_validation.json", "sha256": "ACTUAL_SHA256"}
}
```

Paths are relative to the manifest unless absolute. The TPM table must have
unique `prodigal_protein`, `metaG_tpm`, and `metaT_tpm` columns; extra columns are
retained in the original input. Empty/NA values mean unmeasured, whereas numeric
zero is an observed zero. Negative, nonfinite and malformed numeric values fail.
All TPM gene IDs must equal all original protein FASTA IDs, not cluster IDs.
Duplicate sequences with different original gene IDs are allowed. Duplicate gene
IDs, changed checksums, foreign specimens, and missing genes fail before output.

The provenance record is required and hashed. Its scientific interpretation
is the input provider's responsibility: establish that the quantification
reference really produced these genes/proteins, with the exact specimen linkage,
software settings, genetic code and any translation/index checks. Setting the
status string or matching gene names alone does not establish that provenance.
Do not mark an uncertain candidate reference as verified. A source with known
unresolved genes cannot be used for a complete sweep.

Proteins pass if any encoding gene passes. Gene TPM values are not summed for
selection. Full protein I and L remain different; exact identity follows the
FastaLake 5.3 normalization contract. Existing baseline accessions and headers
are retained. New sequences receive `FMO_SHA256_<full digest>` identifiers.
Every original gene/header remains in `gene_evidence.tsv`, including missing
RNA and genes that do not pass. A no-addition result preserves the baseline
FASTA content exactly. The output is reparsed to verify the complete union.

Use `search.fasta` directly for Sage after this step.

!!! warning

    Do not apply de novo parsimony or clustering again, because it can remove
    TPM-only additions.

`COMPLETE.json` records parameters, gene counts, hashes and the actual union
size. An existing destination is refused; failed attempts remain inspectable.

If a new assembly adds sequences absent from the original reservoir, expose that
same expanded reservoir to the fresh de novo comparator before interpreting the
effect as TPM selection alone. Report reference-refresh effects separately.

## Compare two searches quantitatively

Each search directory must contain `config.json`, `results.sage.tsv` and
`lfq.tsv`. Its configuration must reference the FASTA actually searched and one
acquisition's mzML. Relative paths are resolved against the search directory.

```bash
fasta-lake molecular compare \
  --left searches/metagenome --right searches/denovo \
  --left-label metaG --right-label de_novo --out metag_vs_denovo

fasta-lake molecular compare \
  --left searches/denovo --right searches/tpm_union \
  --left-label de_novo --right-label tpm_union --out addition_effect
```

Add `--compress-tables` for gzip TSV outputs; `COMPLETE.json` records their exact
filenames and checksums.

The comparison rejects different non-file Sage settings and different spectra.
For the same resolved spectrum path it records that identity check; for different
paths it requires equal content hashes. Freeze and record the Sage executable
and source data independently when constructing a new experiment.

| Output | Interpretation |
| --- | --- |
| `protein_overlap.tsv` | Both/left-only/right-only/neither by complete protein SHA-256 over the union of searched target sequences. Includes full sequences, original aliases, headers, accepted peptides, database availability, shared evidence and the actual support peptides not accepted by the other arm. |
| `peptide_overlap.tsv` | Target canonical peptide overlap at peptide q <= 0.01 by default; modifications removed and I/L-equivalent sequences merged for identification counts. |
| `lfq_overlap.tsv` | Positive finite reported peptidoform/charge rows whose canonical peptide passed the PSM whitelist. Modifications and original I/L variants are retained. Missing measurements remain blank. |
| `lfq_excluded.tsv` | Accepted-peptide rows with missing or nonpositive intensity; retains the raw value and exclusion reason. |
| `COMPLETE.json` | Exact definitions, counts, input/output hashes, code hash, intensity agreement and the number of ambiguous LFQ rows excluded from agreement. |

Protein evidence is sequence-compatible accepted PSM support, not unique protein
identification or a new protein-group FDR. Peptides are checked for containment in
every reported protein. The single-sequence-candidate count describes the
reported PSM candidate set only; it does not establish biological uniqueness.
Do not compare arbitrary representative IDs from independently grouped searches.
Use a matched acquisition roster and rerun the supported study grouping per arm
for quantified group counts.

Sage 0.14.6 uses charge `-1` for a feature integrated across charge states; this
code is preserved as a combined feature, not interpreted as a physical ion charge.
The receipt reports the observed charge mode. Missing and nonpositive intensities
are excluded from quantification and retained in a separate raw-value ledger.
Negative values are never log-transformed, clipped or interpreted as abundance. See the versioned
[Sage output implementation](https://github.com/lazear/sage/blob/v0.14.6/crates/sage-cli/src/output.rs).

For LFQ, multiple reported I/L variants can fall in the same isobaric family and
carry different intensities. The tool preserves them without summing or choosing
one. Intensity agreement uses exact shared peptidoform/charge rows only when that
I/L family has one retained row in each arm. Excluded families remain in the
ledger. LFQ's own `q_value` is not used as an additional gate. The acceptance rule
is the declared canonical-peptide whitelist plus positive finite intensity.
Pearson agreement is calculated on log2 intensity; Spearman agreement and raw
median log2 differences are reported separately. These are agreement measures,
not accuracy estimates or evidence of quantitative gain on undetected features.

## Sweep and validation

Keep all raw gene values. Vary one DNA threshold axis and, where measured, the
RNA axis; retain the de novo-only control. Deduplicate identical FASTAs across
grid settings only with an explicit mapping to their shared measured search.
Keep donors and their repeat acquisitions together when separating tuning and
evaluation. Show recovery and losses alongside database size and quantitative
coverage. An optimum chosen and evaluated on the same samples is exploratory.

Run the planted-defect checks with:

```bash
python -m pytest tests/test_molecular_exact.py
```

The suite covers rank edges, ties, disabled/missing layers, strict RNA thresholds,
duplicate genes/sequences, wrong specimens, changed/rehashed incomplete sources,
accession collisions, no-addition identity, shared peptide ambiguity, changed
search settings, false peptide membership, duplicate LFQ rows and I/L families.

A [bundled synthetic example](https://github.com/MannLabs/fasta-lake/blob/main/examples/molecular/README.md) exercises the source
manifest and exact union without private inputs or a search engine.

See the dated software validation record (in the private review record, not part of this repository) for
the tested environment and real-data preflights. The linked benchmark record
supersedes its historical account of unfinished DNA comparisons.

## Earlier clustered workflow

The historical `fasta-lake fl-mo` commands remain available for reproduction
of the archived clustered FL_MO analysis. They are omitted from the main help.
For new exact-sequence additions use `fasta-lake molecular` as described above.
An AND evidence rule is a selection rule, not an independently calibrated FDR.
