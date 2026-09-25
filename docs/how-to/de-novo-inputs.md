# Reuse existing DIANovo predictions

DIANovo inference can be computationally expensive. Reuse a completed, checked
prediction file when possible. FastaLake reads its CSV and curates the candidate
database; it does not run DIANovo or reproduce its model inference.

Take the `pred_seq` and `pred_prob` columns from
[DIANovo](https://github.com/hearthewind/dianovo), developed for *Disentangling the
Complex Multiplexed DIA Spectra in De Novo Peptide Sequencing*, and convert them
to a two-column `sequence,score` CSV (`pred_seq` becomes `sequence`, `pred_prob`
becomes `score`). The study runner accepts only that generic CSV or AlphaNovo's
`peptide_prediction_detokenized_unmodified,score` CSV and rejects other headers.
Give the converted CSV to the prediction input in the
[study tutorial](../get-started/first-study.md). Preserve the original
prediction file, model/preprocessing settings and acquisition identity.

## Two formats, explicit confidence rules

| Input | Example | FastaLake interpretation |
|---|---|---|
| Sequence generation | `AcD` and `[0, 0.9, 0.9]` | Three residue probabilities; mean 0.6 |
| Complete path output | `A c D` and `0 0 0.9 0.9` | One extra leading path marker, then three residue probabilities; mean 0.6 |
| Incomplete path | `A 113.084 D` | Excluded: unresolved mass remains between the residues |

Only canonical amino-acid letters are accepted after converting case and removing
ASCII whitespace. Lowercase modified residues retain their underlying amino acid.
Numeric gaps and modification strings are not erased or inferred. For a path
vector, the extra leading zero is removed only when it is space-delimited and
has exactly one more entry than the complete peptide. Bracket/comma vectors need
exactly one entry per residue. Values must be finite and between zero and one.
Real zero probabilities participate in the mean; an all-zero mean is excluded by
the normal positive-confidence selection rule.

Extraction and razor selection log excluded unresolved/invalid prediction rows.
Rust reads the complete file during selection.
The same parsing contract is used by the experimental Python helper.

Zero probabilities are kept and numeric gaps are preserved.

## DIA scope

A smaller database can feed spectral-library prediction and DIA searching.
Reuse of one acquisition for de novo selection and searching does not establish
whole-workflow FDR. Validate acquisition windows, modifications and calibration
before interpreting results. See [development scope](../project/development.md).

## Find the code

- [Python parser](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/dianovo.py): complete sequence and probability checks.
- [Shared Rust parser](https://github.com/MannLabs/fasta-lake/blob/main/rust/common/dianovo.rs): extraction and parsimony input.
- [Regression tests](https://github.com/MannLabs/fasta-lake/blob/main/tests/test_dianovo_inputs.py): zero-aware target ranking,
  path markers and rejection of unresolved or malformed input.

The format contract was checked against upstream commit
`cc413f25140e47a5f6a4234bae15da0ef525ac3a`, specifically
[`sequence_generation_inference_utils.py`](https://github.com/hearthewind/dianovo/blob/cc413f25140e47a5f6a4234bae15da0ef525ac3a/denovo_from_feature_detection/rnova/task/sequence_generation_inference_utils.py)
and
[`optimal_path_inference_classic.py`](https://github.com/hearthewind/dianovo/blob/cc413f25140e47a5f6a4234bae15da0ef525ac3a/denovo_from_feature_detection/rnova/task/optimal_path_inference_classic.py).
The parser is a FastaLake adapter; DIANovo's model implementation is not copied.
