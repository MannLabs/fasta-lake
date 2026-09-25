# What deserves further development

FastaLake is optimised for a cluster (see [laptop or cluster](../how-to/compute.md#laptop-or-cluster)). Improvements that make substantial studies easier to
run, inspect and reproduce are welcome. Laptop support benefits from the same careful use of
memory and storage, with its practical scope guided by measurements.

FastaLake should make its assumptions and unfinished work easy to inspect.
Passing integration tests establishes that the recorded workflow runs and
preserves its contracts. Scientific calibration and general performance claims
need their own evidence.

## Questions guiding the next steps

| Question | What is known, and where help is useful |
|---|---|
| Where can memory be reduced most usefully? | Chunked extraction is small; razor's protein–peptide relationships dominate the measured full run. Compact representations are the first target |
| Can large runs avoid repeated setup? | Each reference piece rebuilds the prediction matcher. Reusing one matcher across pieces could save work while keeping the same global evidence rules |
| Can an interrupted study restart with confidence? | Stage receipts provide a foundation. Reuse should require matching inputs, settings, executables and outputs before an explicit resume feature is offered |
| When do molecular additions improve a study? | Paired searches reveal gains and losses. Both must be retained, together with reference coverage and held-out evaluation, as the evidence grows |
| How can plots help readers assess their own data? | Missingness, annotation coverage and the eligible PCA feature core belong beside the plots. Shared-peptide and functional ambiguity should remain visible |
| Which small improvements could benefit the surrounding tools? | Reproduced directLFQ dependency and AlphaPeptTools styling suggestions are ready for discussion; the existing Sage correction deserves attribution |

These questions are an invitation to improve the workflow together. The tables
below connect each proposal to an observed need and a checkable result.

## A dedicated review of the implementation

A systematic, line-by-line review deserves its own time and benchmark
record. The current measurements identify promising targets; the size of any
further improvement remains to be measured.

Start with the [execution-order function map](../reference/function-map.md), then follow each
stage through allocations, copies, input scans, algorithmic complexity and
failure/restart behavior. Profile representative full inputs before changing a
stage. Review Rust graph storage and matcher reuse first, followed by Python
parsing, reporting and optional analysis where their measured costs justify it.

Keep the accepted source revision and outputs as the comparison baseline. For
each proposed change, record the input, bottleneck, implementation, correctness
checks and before/after time, peak memory and disk use. Preserve exact scientific
contracts or document an intentional method change separately. Small, reviewable
changes make it easier to retain improvements and trace published results.

## Priorities in FastaLake

| Priority | Observation | Useful next change | Acceptance condition |
|---|---|---|---|
| Razor memory | A full-reference pilot exceeded 8 GiB; extraction itself used much less memory | Intern protein/peptide IDs and use compact adjacency arrays | Same selected sequences, tie handling and output order on adversarial and full inputs; lower measured peak RAM |
| Repeated extraction setup | Each reference piece repeats prediction parsing/ranking and matcher construction | A long-lived Rust matcher over successive pieces | One immutable global prediction roster; identical merged evidence across piece sizes |
| Restarting long runs | Completed stages have receipts, but the runner deliberately requires a new destination | Explicit resume from verified stage identities | Reject changed inputs, settings or binaries; interruption/restart must reproduce a fresh run |
| Portable resource guidance | Full-data measurements are from Linux Slurm; Apple Silicon uses the x86_64 image through emulation | Measure complete acquisitions on representative laptops | Record actual hardware, wall time, RAM and disk; retain failed runs |
| Reference-wide molecular support | Full-reference sequence matching can dominate a large comparison | Profile scans and retained compatibility sets; compact only measured bottlenecks | Preserve exhaustive peptide specificity, denominator and gained/lost support accounting |
| Optional analysis | Complete-case PCA can retain a small feature core; functional attribution remains ambiguous | Better diagnostics and explicit, separately validated analysis choices | Export eligibility/coverage and missingness; do not turn imputed values or representative labels into evidence |

The [performance map](../concepts/performance.md) links these priorities to current code and
measurements. They are development proposals, not unmeasured speedup claims.

## Scientific questions remain explicit

- **Adaptive selection and error rates.** Search-level q values do not calibrate
  the whole data-dependent selection/search workflow. Independent entrapment and
  held-out tests remain separate from software acceptance.
- **Molecular additions.** DNA/RNA/compound support selects candidates. Matched
  search comparisons must retain reciprocal losses, reference coverage and the
  difference between selection support and measured peptide support.
- **Biological attribution.** Shared peptides, incomplete annotation and repeated
  donors limit species/function claims. Experimental downstream helpers expose
  these choices; QC/PCA does not establish a biological mechanism.

The scientific status and the downstream function review are kept with the private review record (not part of this repository).

## Contributions to the surrounding tools

Checked against current upstream source on 15 September 2026. Reproducers use
the stated released versions; source inspection of a development branch is
recorded separately. These findings are focused proposals for review.

| Project | Finding | Disposition |
|---|---|---|
| Sage | The full target/decoy LFQ return key is already corrected upstream | Credit [the existing fix](https://github.com/lazear/sage/commit/0e2c7a6d44d935dfab4230e747f9e12e1a993ce2). The older FastaLake benchmark runtime and audited backport (in the private review record, not part of this repository) retain explicit identities; no duplicate issue is needed |
| directLFQ | A clean 0.3.3 core installation passes dependency checking but `import directlfq.lfq_manager` fails because PyYAML is absent. [Current core requirements](https://github.com/MannLabs/directlfq/blob/c1b5b650b61557c00fafc1f436eeeff6552b25de/requirements/requirements_loose.txt) still omit it | Propose declaring the imported dependency and testing a core-only install. FastaLake already declares PyYAML explicitly |
| AlphaPeptTools | `create_figure` in 0.4.0 changes global Matplotlib settings; [current source](https://github.com/MannLabs/alphapepttools/blob/df8714f3b9410d7c80a4b328b1ce4bf24e8d5c05/src/alphapepttools/pl/figure.py) retains this behavior | Propose scoped figure styling or an explicit global-style option. FastaLake scopes its own plotting settings; this is an interoperability improvement, not a numerical error |

The source records and minimal upstream requests (in the private review record, not part of this repository)
include duplicate-search scope and reproducer results. No upstream submission is
implied by a draft. Follow each project's contribution process and preserve its
attribution.
