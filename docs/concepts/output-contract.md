# One study-group table from direct peptide assignment

The current release implements the direct assignment agreed on 10 September 2026.
Each retained LFQ feature follows its canonical peptide directly to the first
selected study representative that explains it. There is one quantitative table,
one matrix and one per-acquisition summary under this definition. Full local
memberships are retained separately for interpretation.

The pre-search razor, search databases, Sage searches, acceptance threshold and
greedy study cover are unchanged. The former LFQ detour through a whole-protein
plurality map is removed. In this route, the lexical local representative is no longer a
quantitative partition or a second protein-group counter.

## HSTOOL population

HSTOOL is exactly 58 acquisitions. Differential centrifugation is excluded before
building either dictionary or quantifying either approach. Use the frozen
`publication/analyses/simplified_grouping_2026-09-10/hstool_samples.txt`; the
paired design contains 12 technical, 15 workflow, 10 spatial, 11 intra-individual
and 10 inter-individual acquisitions. Filtering display rows after inferring a
broader dictionary does not meet this definition.

## Run on completed searches

```bash
python tools/group_study.py --search /path/to/completed/acquisitions \
    --out /path/to/new/grouped-output --config-name results.json \
    --samples-file publication/analyses/simplified_grouping_2026-09-10/hstool_samples.txt
```

Use `config.json` (the default) for the supported release runner's searches. Use
`results.json` explicitly for the retained publication runs. Each configuration
must identify the actual searched FASTA by absolute path. The output directory
must be new and outside the search-input tree. This command does not run a search.

## Output contract

| File | Meaning |
|---|---|
| `study_group_long.tsv` | One row per acquisition and quantified study group; distinct canonical peptide count, LFQ feature count and summed intensity |
| `study_group_matrix.tsv` | Exactly the same groups and intensities as a wide matrix; unquantified cells are empty |
| `acquisition_summary.tsv` | Group presence and evidence counts derived from that same quantitative table |
| `feature_assignments.tsv.gz` | Every eligible positive feature, original file line, peptide/charge, reported proteins, full member set and direct assigned group |
| `peptide_to_study_group.tsv` | The direct pooled peptide assignment |
| `study_groups.tsv` | Selected representative, selection rank and peptide-support counts |
| `input_manifest.json` | Hashes of configurations, searched FASTAs, search/LFQ inputs and implementation |
| `peptide_evidence.tsv.gz` | Accepted canonical peptides by acquisition/group, including unquantified IDs, local target memberships and local/study accession counts |
| `peptide_ambiguity.tsv` | Histogram of study-wide reported accession counts per canonical peptide |
| `summary.json` | Population, rule, confidence scope and conservation checks; written last |

`study_group` is the selected accession naming a group. It does not claim that
the accession is the uniquely identified biological source of all group signal.
The greedy cover is a heuristic, and equal-gain ties still use accession ordering.
This first change does not make the inference independent of naming or study scope.

Local memberships no longer create separate quantitative rows. They remain in
the feature ledger, alongside flags showing whether the representative was in
that acquisition's FASTA and in that feature's reported member set. A peptide
can have an explaining representative elsewhere in the study without establishing
that representative's identification in the contributing acquisition.

The implementation checks accession-to-sequence consistency across all FASTAs,
each observation's membership in its own searched database, and containment of
every assigned peptide in its representative sequence, with I/L equivalent.
These are evidence-integrity checks, not FDR calibration or unique-source proof.

## Populations and biological interpretation

Input target rows pass finite peptide q in [0, threshold]. Modification brackets
are removed and I/L collapsed before pooling and the local LFQ peptide gate.
Positive finite LFQ rows contribute once each; charge/modification features are
counted separately from canonical peptides. Peptide acceptance is not replaced
with the LFQ row's own q-value. Zero-output acquisitions remain in the acquisition
summary and matrix columns. Inferred groups and quantified groups are reported
separately because some accepted evidence has no positive quantified feature.

This route's per-acquisition group count means quantified presence under the
declared study dictionary. It is not the historical independently inferred
`pg_parsimony` endpoint. The v3_simplified manuscript (kept privately)
uses the completed four-cohort reanalysis, with updated figure inputs and
explicit endpoint definitions. The current working draft is
6.1 (manuscript drafts are kept privately).

Representative annotations still require interpretation against the contributing
evidence. Direct peptide assignment repairs the sequence-support detour. It does
not resolve all shared taxonomic origins, missing annotations or group confidence.
Functional contributor and correlation sensitivity have now been measured in the
completed 58-acquisition comparison below. Biological source accuracy remains a
separate question.

## Completed HSTOOL comparison

The paired post-search rerun is complete. Both the original uploaded rule and this
rule used the exact 58-acquisition roster before inference. The original protein
map matched the independent cohort audit; all direct assignments passed sequence
compatibility and the group table was independently recounted from its features.

- Full report: in the private review record (not part of this repository)
- Workbook with definitions and all principal endpoints: in the private review record (not part of this repository)
- Independent delivery checks: in the private review record (not part of this repository)

The package includes selection/search inventories, observed peptides per accession,
assigned peptides per quantified group, coverage, group completeness, correlations
with both method-specific and identical row populations, and KO/species sensitivity.
Search results were reused. No raw de novo prediction or Sage search was rerun.
The checks support simpler sequence-consistent reporting. They do not establish
unique biological origins or turn higher group counts into additional identifications.

Result source: `comparison/SUMMARY.json` in the linked run directory, SHA256
`5c1b7eac6f9e426a5a886a24d3f2687cd125582efadf21830ad08f199f0d88d4`.
Generating code commits and input hashes accompany each run.

## Provenance

Direct peptide assignment replaced the dual-key rule on 11 September 2026.
`tools/group_study.py` writes the new schema and the example consumers follow it.
Outputs produced before that date record the previous method.

## Optional quantification

The sum remains the default. [Selectable quantification](../how-to/quantification.md) adds
directLFQ over the same assignments, with the full ion input, peptide retention
and coverage diagnostics in a new destination.

The grouping command also writes basic QC and eligible PCA plots under `qc/`.
Use Python 3.11+ with the analysis extra, or the Docker image. A core-only
installation can explicitly request `--skip-qc`. See [QC/PCA](../how-to/downstream-analysis.md).
