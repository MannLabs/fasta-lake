# What has been tested?

Use this page to choose a check and understand what its result means.
For an explanation of the method, start with the
[visual walkthrough](walkthrough.md). For a working installation, run the
[measured CAMPI example](../get-started/first-study.md).

## Four kinds of evidence

| Question | Evidence you can inspect | What it supports | What it leaves open |
|---|---|---|---|
| Does the code follow its declared rules? | [Python/Rust tests and CI](https://github.com/MannLabs/fasta-lake/actions/workflows/ci.yml); worked [grouping example](walkthrough.md) | Checked input handling, deterministic rules, sequence compatibility and accounting | A correct implementation can still have scientific limitations |
| Can I run the workflow? | [CAMPI example](https://github.com/MannLabs/fasta-lake/blob/main/examples/campi/README.md), its expected results and your own `VALIDATION.json`; [Docker checks](https://github.com/MannLabs/fasta-lake/actions/workflows/docker.yml) | Installation and execution on the supplied small measured inputs | Full-acquisition resource needs and performance on your study |
| What happened on measured study data? | [Dated acceptance experiment](validation-acceptance-2026-09.md): 20 FastaLake acquisitions plus 10 matched-metagenome searches | Execution, descriptive recovery, overlap and quantitative variation in that declared subset | A general advantage across cohorts or an isolated causal effect of curation |
| How reliable are the biological identifications? | Search configuration and confidence fields; [error-control framework](error-control.md); study-specific control analyses | Explicit search acceptance and the evidence needed to evaluate stronger claims | Whole-workflow calibration, unique biological origins and study-group FDR are not established by the checks above |

These are different questions. Passing a test or reproducing a matrix is useful
evidence for its stated purpose; neither is an independent calibration experiment.

## A reproducible starting point

The distributed CAMPI teaching example uses two measured acquisition windows.
It includes the spectra, predictions, reference, runner and expected results.
You can inspect every stage and retain the exact commands and input identities.

1. [Run the example](../get-started/first-study.md).
2. Check `VALIDATION.json` and the process exit status.
3. Open the selected FASTA, peptide evidence and group matrix.
4. Follow a contribution using the [aggregation walkthrough](walkthrough.md#3-how-do-measured-features-become-quantities).
5. Pilot one representative full acquisition before choosing resources for a study.

The cropped example teaches execution and output interpretation. Its numbers
are not full-run abundance estimates or a benchmark of de novo model accuracy.
[Resource measurements](performance.md) describe specific workloads, not guaranteed
maximum requirements.

## Questions to carry into a scientific study

| Question | Why it matters | Evidence that would address it |
|---|---|---|
| Does selecting candidates using the same acquisition affect error calibration? | The subsequent search is conditional on an adaptively selected database | Controls with justified absence/exchangeability assumptions and designs separating construction from evaluation |
| How does reliability vary with peptide length? | Shorter strings are less specific, and search populations influence confidence estimation | Length-resolved controls; clearly distinguish a subset of an existing search from a new search with a different minimum length |
| Does a group have a unique protein, species or functional origin? | Several reference sequences can explain the same peptide evidence | Full local memberships, sequence support, annotation provenance and appropriate ambiguity/sensitivity analyses |
| Will results transfer to a new cohort or predictor? | Reference content, acquisition quality and prediction errors can change | Held-out studies with declared inputs and paired gains as well as losses |
| Are missing quantities evidence of biological absence? | Sampling, selection, search acceptance and quantification can each produce missingness | Acquisition-level evidence and a separately justified downstream missing-data model |

The [error-control page](error-control.md) explains the statistical scope.

### Additional discoveries need a separate check

A peptide q-value threshold applies to the population for which the search
estimated it. The identities added relative to another search form a selected
subset; they do not automatically inherit that population's error guarantee.
Report additional identities, reciprocal losses and shared identities separately,
using a matched comparator and the same declared acceptance rule. Identification
counts at a nominal threshold are not verified gains at matched true error rates.

Database selection can remove a sequence that would otherwise compete to explain
a spectrum. Retaining more candidates, using a stricter q threshold, or increasing
the minimum peptide length does not by itself establish calibration. A proposed
change needs evaluation of both error and retained true yield, followed by fresh
evidence that was not used to choose the change. An empty accepted set is
abstention and supplies no empirical precision information about recovered
identities. Simulation truth can diagnose a failure mechanism; representative
measured controls and justified control comparability remain separate requirements.

The [output contract](output-contract.md) explains the reporting scope.
The [manuscript reproduction guide](../reference/manuscript-reproduction.md)
keeps paper settings and producing versions distinct from a current software run.

## Keep the historical record interpretable

The earlier validation page is retained as the
[September 2026 acceptance record](validation-acceptance-2026-09.md), including
its measurements, adverse observations, failed attempts, retries and limitations.
Its ten-acquisition HSTOOL subset is not the later 58-acquisition reporting study.
Its test counts and dependency audit describe that dated build, not every later
commit. Use the CI result for the exact revision you run.

For reproduction, record the revision, input identities, searched FASTAs,
search settings, acquisition roster and grouping/quantification policy.
[Report a study](../how-to/reporting.md) lists the files to preserve.

## Help improve the method

A useful collaborator report names the claim or behavior, the exact inputs and
revision, what was observed, and a small shareable reproducer when possible.
Failures and unexpected results are useful. Keep participant data out of public
issues; provide an approved minimal example. See
[contributing](../project/contributing.md) for the development workflow.
