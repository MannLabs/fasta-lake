# Quantification details

Start with [Quantify a study](../how-to/quantification.md) for method selection and commands.

## Which key is shared between samples?

The merge key is `(study_group, ion)`, where `ion` is the modified peptide with
I/L equivalence plus charge. The assignment lookup uses the unmodified canonical
peptide. Charge and modification states remain different quantitative traces;
they do not count as different canonical peptide sequences.

Sage can already combine charge states during LFQ extraction. Its `charge=-1`
sentinel is preserved as one combined-charge trace; those values cannot be split
back into individual charges. The audited CAMPI input uses this encoding. Mixing
combined-charge and charge-resolved quantities for the same modified peptide
within a study is rejected. [Sage release notes](https://github.com/lazear/sage/blob/v0.14.6/CHANGELOG.md)

| Study group | Ion | Sample 1 | Sample 2 |
|---|---|---:|---:|
| A | peptide-1, charge 2 | 100 | 200 |
| A | peptide-1, charge 3 | 30 | missing |
| A | peptide-2, charge 2 | missing | 20 |

These illustrative values are not experimental results. A sum gives 130 and
220. directLFQ needs the individual ion rows to estimate their profiles. Sending
an already summed group matrix loses that information.

The same canonical peptide has one assigned reporting group throughout a fixed
study, even when its local reported proteins differ between samples. A join on
local Sage membership strings would split such evidence. Copying a shared
peptide's full intensity into every compatible protein would count it repeatedly.
The group label is a reporting convention and does not establish unique molecular
origin. Full and matched rosters, or different database arms, need separately
identified dictionaries; their group labels are not automatically interchangeable.

## Inspect the peptides, including those without LFQ measurements

Grouping runs write:

| File | Key and contents |
|---|---|
| `peptide_evidence.tsv.gz` | One row per acquisition, study group and accepted canonical peptide; original target accession union, local accession count and study accession count |
| `peptide_ambiguity.tsv` | Histogram of study-wide reported accession counts over distinct accepted canonical peptides |
| `feature_assignments.tsv.gz` | Every accepted positive LFQ source feature, charge, modification, line number, local memberships and assigned group |

`peptide_evidence.tsv.gz` includes accepted identifications without an accepted
positive LFQ feature. The feature ledger alone cannot represent those IDs.
Evidence and quantities join on `(sample, study_group)`; retain the peptide column
to compare supporting sequences between samples. Reduce peptide rows to a set
or a count before attaching one group quantity, to avoid replicating and then
summing that quantity once per peptide.

Example: count peptides reported against more than x proteins:

```python
import pandas as pd
h = pd.read_csv('study_groups/peptide_ambiguity.tsv', sep='\t')
x = 5
n = h.loc[h.n_reported_proteins_study > x, 'n_canonical_peptides'].sum()
print(n, n / h.n_canonical_peptides.sum())
```

These counts refer to distinct target accessions reported in accepted PSMs.
They are neither an exhaustive remapping against every searched FASTA nor a count
of distinct protein sequences, species or reporting groups. Database redundancy
can affect them. A peptide can have one locally reported protein but several
across the study. Identification ambiguity and quantitative consistency are
separate properties.

## Inspect what directLFQ retained

The directLFQ destination additionally contains:

| File | Interpretation |
|---|---|
| `ions.aq_reformat.tsv` | The exact input matrix, indexed by `protein` (study group) and `ion`; 0 means missing |
| `ion_mapping.tsv` | Ion-to-canonical-peptide mapping, numbers of input and retained acquisitions, and presence in the output ion roster |
| `peptide_support.tsv.gz` | Each measured canonical peptide by acquisition/group; numbers of input and retained ions, input intensity sum and whether the group retains a quantity |
| `ions.aq_reformat.tsv.ion_intensities.tsv` | Backend-aligned ion values; they are not raw intensities or additive protein contributions |
| `ions.aq_reformat.tsv.protein_intensities.tsv` | Unmodified backend protein output |
| `ions.aq_reformat.tsv.run_config.tsv` | Backend settings from the temporary execution path; the stable input is the adjacent `ions.aq_reformat.tsv` |

Left-join identification evidence to peptide support on
`(sample, study_group, canonical_peptide)` to distinguish an identification with
no accepted LFQ measurement from measured peptide evidence lost during alignment.
The source evidence files remain in the input study directory, recorded in
`summary.json`.

The exporter checks unique source rows, fixed peptide assignments, positive
finite values, valid charges and exact acquisition rosters. I/L-equivalent
repeated ion measurements are summed before import; the actual importer must
preserve every key and intensity. Exported ions must reproduce every source
study-group sum. Outputs cannot introduce positive ion/group cells where input
measurements were missing. Source files are hashed before and after processing;
backend version, implementation hashes, settings and dependency versions are saved.
These checks validate the handoff, not independent identification confidence.

## What does “more stringent” mean here?

The backend uses `min_nonan=1`, sample normalization enabled, 50 sample anchors
and 10 ion anchors, matching the audited settings (thread count is configurable).
It retains at most 100 input ions per group. The cap favors completeness, then
the sum of normalized log2 intensities. This is not a filter on peptide q-values.

During quadratic alignment, traces with fewer than two observed values are set
to missing. For ion alignment, this means observations across acquisitions. The
linear alignment path for ions outside the anchor subset has different behavior.
Thus `min_nonan=1` is not a guarantee that every originally measured ion or
single-acquisition group retains a quantity. Disconnected traces also need care:
small disconnected components can retain values without an estimated between-
component offset, whereas a linear trace with no overlap to its reference can
become missing. Presence in an output alone does not establish a supported ratio.

A lost quantity is not evidence that its peptide identification was false.
Inspect coverage and supporting ions alongside precision. Do not fill losses
with raw sums and call the hybrid matrix directLFQ: the scales and estimators differ.

In the previous 14-run manuscript comparison, directLFQ improved HSTOOL median
CV on fixed shared complete groups: 41.46% to 31.13% for technical replicates and
59.37% to 37.21% for workflow replicates. Improvements remained after matching
sample normalization and retained ion evidence. Median paired coverage losses
ranged from 2.29% to 11.86% across runs. These are descriptive precision and
coverage results; concentration accuracy and preservation of biological effects
were not established. The sum remains the default.

## QuantSelect compatibility

[QuantSelect](https://github.com/MannLabs/QuantSelect) learns ion quality using
additional features and directLFQ alignment. The reviewed implementation expects
fragment intensity, correlation, mass-error and related feature tables from a
DIA workflow. Current Sage DDA `lfq.tsv` files provide precursor MS1 quantities
and do not supply those fragment-level inputs. QuantSelect is consequently not
exposed as a selectable backend here. A DDA adaptation would need real additional
features and separate validation; duplicating MS1 values into fragment columns
would not implement the published method.

Method sources: [directLFQ publication](https://pmc.ncbi.nlm.nih.gov/articles/PMC10315922/),
[directLFQ code](https://github.com/MannLabs/directlfq),
[QuantSelect feature loader at the reviewed revision](https://github.com/MannLabs/QuantSelect/blob/05a62d00de9001e2ab4e940bc6be5159c1ab644c/quantselect/loader.py),
[QuantSelect preprint](https://doi.org/10.64898/2026.01.08.698341).

Automatic [AlphaPeptTools downstream analysis](../how-to/downstream-analysis.md) exports AnnData, coverage/completeness QC, explicit log2 PCA and declared replicate-set CV after either quantification method.

## Automatic plots

The grouping tool, the manifest runner and `fasta-lake quantify-study` all
write AlphaPeptTools QC/PCA under `qc/` in their result directory. For the
runner, this is `out/study_groups/qc/`. Plots follow the selected quantities;
the original sum matrix and missing observations remain traceable.
`--skip-qc` explicitly opts out for a core-only environment (including Python 3.10).
The Python `quantify_study` function remains a low-level quantification operation.
