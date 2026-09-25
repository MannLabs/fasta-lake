# Measured evidence for molecular additions

FastaLake supports metaG and metaT additions on top of the de novo result.
The completed CapScan evidence tests metaG (DNA); the complete PREDICT metaT
(RNA) cohort analysis is not measured here. The input and selection code supports
both layers through the same specimen gene mapping.

## Shared-reference, held-out CapScan comparison

The completed benchmark covers 504 conditions in 28 acquisitions from 15 donors.
Ten donors selected cutoffs; five held-out donors provide the values below.
De novo and molecular selection start from the same available specimen reference,
and every arm receives the same human background.

| Condition | Positive identified LFQ features relative to de novo |
| --- | ---: |
| De novo baseline | 100.0% |
| De novo + all-positive DNA TPM | 107.1% |
| De novo + metabolite-linked candidates | 101.3% |
| DNA TPM alone, all positive genes | 106.7% |
| DNA TPM alone, top 25% of positive genes | 57.6% |

All-positive DNA was selected on tuning donors only for both the DNA-only and
de novo-plus-DNA families. It is a tested-panel result, not a universal cutoff
recommendation. De novo uses about 2,345 non-host candidate sequences versus
65,056 for all-positive DNA, as donor-weighted mean counts. Both include the
same 20,353 human sequences.

Adding DNA to de novo gains about 1,601 canonical peptide identities and loses
654 after donor weighting; metabolite additions gain 230 and lose 89. A FASTA
union preserves candidate sequences. It does not guarantee retention of every
previous search identification.

The reviewed result, its independent validation and the donor values are kept with the private review record (not part of this repository).

## Earlier broader-reference comparison

The earlier CapScan comparison allowed broader public/cohort candidates in the
de novo arm. Across 28 acquisitions and 15 equally weighted donors, all-positive
DNA alone reached 63.95% of de novo feature counts and adding DNA reached 103.75%.
These are different candidate-reference and evaluation populations, not a change
in performance on an unchanged baseline. The
original evidence bundle (in the private review record, not part of this repository)
preserves that curve and its reciprocal gains/losses.

## Interpretation

Counts use the original single-acquisition, accepted-peptide-whitelisted positive
LFQ features. They are distinct from directLFQ quantities, corrected joint LFQ or
validated match-between-runs. Neither software execution nor nominal peptide
q values establish whole-selection-procedure FDR.

The compound-to-KO results describe compatible evidence. Metabolites used to
select candidates are not independent validation of those candidates. No flux,
producing species or unique enzyme identity follows from a mapping.

PREDICT uses a fixed mapped compound panel across its saved patients. Its current
metabolite arm therefore does not establish a benefit from personalized compound
presence. Complete held-out metaT cohort results are not measured here.

See the [input guide](../how-to/molecular-inputs.md) for the installed additive workflow
and the v7 manuscript proposal (in the private review record, not part of this repository)
for wording linked to these source records.
