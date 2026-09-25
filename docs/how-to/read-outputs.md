# Read the outputs

This page lists the files a completed run writes, what each one means and the
first check to make on it. You need a completed run directory from the manifest
runner.

## Open these first

Inside your completed run directory:

| Start with | What to look for |
|---|---|
| `study_groups/qc/QC_overview.png` | Do acquisitions have similar coverage? Which groups are missing? Was PCA eligible? |
| `study_groups/study_group_matrix.tsv` | One reporting group per row, acquisitions in columns; the default quantities are sums. |
| `study_groups/peptide_evidence.tsv.gz` | Which accepted peptide sequences support each group and acquisition? |
| `COMPLETE.json` | Did the requested workflow finish, and did it include searching? |

For a selected quantification method, follow the matrix path in the run's
quantification summary; directLFQ writes a separate matrix for the same dictionary.
A blank cell means no positive accepted quantity was recorded under that rule.
It does not show that the protein was biologically absent.

To see real filenames before your first study, open
[the demo outputs](../get-started/demo.md#open-these-first).
An assistant can walk through them using the
[read-results task](assistants.md#ready-to-copy-task-files).

## Run outputs

| Output | Meaning | First check |
|---|---|---|
| `PROVENANCE.json` | Parameters, input identities and executable identities | Correct acquisition pairing and versions |
| `commands.json` and stage logs | Actual execution and failures | Every required stage returned zero |
| `evidence.fasta` and sidecars | Cohort evidence-supported reference subset | Nonempty output and plausible retained size |
| `samples/<sample>/sample_lake.fasta` | Acquisition-specific drawn sequence candidates | Reduction relative to the evidence lake |
| `samples/<sample>/inference/<sample>_razor.fasta` | Exact selected search candidates | Preserve together with the SAGE configuration |
| `search/<sample>/results.sage.tsv` | Spectrum matches and confidence fields | Target label, appropriate q-value and counting unit |
| `search/<sample>/lfq.tsv` | Quantified peptide features | Acquisition column, finite intensities and matching PSM table |
| `aggregate/` | Peptide-feature and protein-member-set matrices | Threshold and grouping mode in filenames/metadata |
| `COMPLETE.json` | Successful termination of the requested runner scope | Whether databases-only was enabled |

The core LFQ aggregator requires the corresponding results table and exactly one acquisition intensity column per search directory. Multiple intensity columns fail clearly instead of silently selecting the first. The peptide-feature matrix is gated on search-level peptide_q, and the protein-member-set matrix on search-level protein_q. The core aggregator does not gate these matrices on the MS1-specific LFQ q_value column. Their total signal can differ because their identification gates differ. Aggregate observations are filtered before summation. For example, intensities 10 and 1,000 from observations with q-values 0.001 and 0.9 must contribute 10 at a 0.01 threshold, even when they share a reporting key.

!!! warning
    A missing matrix entry or zero output intensity indicates that no accepted positive contribution was recorded under that rule. It is not proof of biological absence. Different spectra, ionization, sampling depth, database inclusion and filtering can all create missingness. Do not replace missing values with small constants without a separate, justified downstream analysis plan.

## Direct peptide assignment to study groups

The manifest runner performs this step automatically, followed by selected
quantification and QC/PCA under `study_run_01/study_groups/qc/`.
The default method is sum; use `--quantification directlfq` to choose directLFQ.
`QC_overview.png` and `.pdf` show coverage, completeness and PCA eligibility.
Eligible studies also receive `PCA.png` and `.pdf`. With fewer than three
nonempty acquisitions or two complete varying groups, the QC summary explains
why PCA was skipped. No missing quantities are imputed.

For existing searches, grouping and plots can also be generated directly:

```bash
python tools/group_study.py --search study_run_01/search \
  --out study_groups_01 --peptide-q 0.01
```

This post-search operation reads the exact searched FASTAs from each configuration and builds an observed target peptide-to-protein graph at the declared peptide-q threshold. Greedy set cover selects study representatives by maximum additional peptide coverage, with lexical ties. The full rule and ambiguity handling are in the algorithms document.

The output has one row per acquisition and quantified study group in `study_group_long.tsv`, the matching `study_group_matrix.tsv`, and full local memberships in `feature_assignments.tsv.gz`. Each feature follows its canonical peptide directly to the first selected representative containing it. `n_peptides` counts distinct I/L-canonical unmodified sequences, whereas `n_features` counts included LFQ rows. Different charge states or modifications can therefore contribute several features to one peptide count. Intensity is summed once per included feature and checked for conservation.

The study utility uses a target peptide-identification gate and positive finite LFQ intensities. Its default does not impose an additional threshold on the MS1-specific `q_value` column. This definition differs from a claim that every quantification feature has MS1 q ≤ 0.01. State both identification and quantification filters when reporting a quantitative result.

!!! warning
    A study group is a pooled reporting construct. A representative can change if the cohort, accepted peptide graph or tie rule changes. Retain the peptide-to-study map, the declared acquisition list, and the full feature ledger with local-member and local-FASTA flags. Never combine study keys from independently built dictionaries as if they were a universal identifier system.
