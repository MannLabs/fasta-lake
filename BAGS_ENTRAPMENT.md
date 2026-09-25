# End-to-end absent-species diagnostic

This experiment addresses the reuse of spectra for selection and search. It is
separate from the running yield comparison and requires independently documented
absence. No absent-species calibration result has been obtained by adding this code.

The full target reservoir and full absent-species proteins are combined **before
retrieval**. Original accessions and the existing deterministic selection rule are
retained. Accession collisions abort; control labels are kept outside selection.
The exact, bag, mass and combined arms share the same augmented reservoir, original
predictions, selection rules, spectra and Sage settings. The two shuffled-query
arms remain retrieval diagnostics. Controls are ordinary targets for Sage; Sage
also generates its usual decoys from each selected database. No control protein
is forced to survive selection or added only to a final selected FASTA.

The runner retains input and retrieved/selected control counts, accepted target
PSMs, peptide identities, shared sequences and reciprocal gains/losses. Thresholds
are fixed at 0.001, 0.005, 0.01, 0.02 and 0.05. Canonical peptide sets use peptide_q;
PSM assignments use spectrum_q. Modifications are removed and I/L is collapsed
for peptide identity, consistently with the yield experiment. Neither endpoint
is a protein-group error estimate.

Control-only identification requires exact peptide containment in the full control
FASTA and no containment in either the full target reservoir or the complete
possible-present reference, including contaminants. Membership is assessed against
preselection sequences, so discarding a real protein cannot make its shared peptide
a false control hit. Shared peptides remain in a separate category. A peptide not
found in any declared source aborts the diagnostic. Target compatibility alone is
not proof of a correct spectrum assignment.

The reported control-only/accepted fraction is a descriptive observed fraction,
not an estimate of target-only FDR, an upper bound or a calibration pass. Zero
surviving controls may indicate insufficient exposure. Interpretation must consider
control size, sequence similarity, searchable peptide opportunities, depletion,
replication and uncertainty. No universal multiplier is applied to convert species
counts or selected protein counts into FDR. Entrapment estimators depend on their
construction and assumptions; see [Wen et al., Nature Methods (2025)](https://www.nature.com/articles/s41592-025-02719-x).
All reports retain `workflow_fdr: UNVERIFIED` pending independent scientific review.

## Inputs and execution

Use a fresh manifest following `tools/run_bags.py`, with full SHA-256 identities
for every input (including the complete unselected reservoir), binaries and code.
Freeze this checkout's code, including both runners, rather than editing a manifest
belonging to active jobs. The archived exact FASTA is not the comparator here:
exact selection is rerun from the same augmented reservoir as the relaxed arms.
The source reservoir must precede *all* spectrum-driven filtering, including any
cohort-wide prefilter. Do not supply an acquisition-selected or pooled-evidence lake.

The control-specification JSON has exactly three file metadata objects (each
containing `path`, `bytes`, `sha256`):

- `entrapment_fasta`: full proteins of independently absent species.
- `possible_present_fasta`: complete plausible present proteomes and contaminants.
- `absence_evidence`: JSON with `absence_confirmed: true`, `samples` (exact acquisition
  IDs), `species`, `source`, `reviewed_by` and `contaminant_policy`.

The reviewer must establish absence from the experiment's design and provenance;
the software only checks that the assertion is present and covers the acquisition.
Missing detection in this workflow is not evidence of absence. Keep specimen and
control identities private unless their distribution is authorized.

```bash
python tools/run_bag_entrapment.py frozen_manifest.json 0 controls.json new_entrapment_output
```

Submit a representative pilot under the study's scheduler limits before a broader
run. This writes a full augmented reservoir and requires space for that copy, hit
tables and all searches. It streams sequence membership with an Aho–Corasick query
index; it does not load the full target reservoir into a Python dictionary. Control
accessions are held in memory. A completed mechanical fixture or pilot does not
establish biological absence, statistical power or workflow FDR. Retain failed and
unfavourable results. Review the whole population and the added-identification
subset separately; a favorable result in one does not validate the other. Cohort
pooling, reporting groups, quantification and their error rates are separate
endpoints; this runner evaluates acquisition-level peptide and PSM identifications.

## Separate shuffled-query sensitivity check

The existing native scan contains route tags for each seed's exact, bag and mass
matches. `tools/run_bag_nulls.py` applies the same fixed-count selection separately
to the six previously missing null arms and runs matched Sage searches. Its
`candidate`, `search` and `report` stages use a separate frozen manifest with
`source_manifest` metadata pointing to the original study. Reports reuse completed
original combined-control searches. Retain all comparisons, including increased
null hits and reduced authentic yield. This is a diagnostic of matching/selection
promiscuity, not independent final-identification truth.
