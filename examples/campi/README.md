# Run two real CAMPI acquisitions

This example runs database construction, acquisition-specific razor selection, SAGE searching, label-free quantification and study grouping on two measured acquisition windows. The spectra, matching de novo predictions, reference FASTA and expected results are included. No laboratory storage paths or additional dataset downloads are required.

The inputs are five-minute windows from S01 and S02 in the public CAMPI SIHUMIx mock-community dataset, [PXD023217](https://www.ebi.ac.uk/pride/archive/projects/PXD023217). They are teaching inputs, not the complete CAMPI experiment. See [data provenance and attribution](DATA_PROVENANCE.md).

## One command

From the top level of the extracted FastaLake source bundle:

```bash
python examples/campi/run_example.py --out campi_run_01 --threads 2
```

First use requires **Linux x86_64, Python 3.12 with venv/pip, Cargo and a C compiler**, plus internet access to obtain Python and Rust dependencies. Python packages are installed from the hash-locked `requirements/container.txt` (compiled for Python 3.12 on Linux x86_64, the same set the container and CI use), so the example runs against known dependency versions. The pinned Rust toolchain is 1.91.1. SAGE 0.14.6 for Linux x86_64 is already included, with its MIT licence. Other operating systems are not validated by this example. The Python wheel alone does not contain these datasets or Rust executables: use the source bundle.

The command makes a private Python environment and compiles FastaLake inside `campi_run_01/_runtime`. It leaves the shared source tree and earlier outputs intact. Installation and compilation take longer than the analysis, which took approximately five seconds for these inputs after setup in the release environment. A 6 GiB allocation was sufficient for the tested setup and workflow; this is not a measured minimum or a full-reference resource estimate. Allow several gigabytes of writable space for the retained build directory and environment.

On an HPC system, submit the command through the site's scheduler. On a workstation, run it in a terminal. A completed installation can be reused:

```bash
python examples/campi/run_example.py --out campi_run_02 --threads 4 \
  --runtime campi_run_01/_runtime
```

Reuse verifies the source fingerprint and executable checksums. It does not repeat dependency installation. Keep the original runtime at its existing path because Python virtual environments are not relocatable. Choose a fresh output name for every run. An existing folder is rejected, including a failed or incomplete run.

## What is actually included

| Input | S01 | S02 |
|---|---:|---:|
| Retention-time interval, minutes | 45 ≤ RT < 50 | 45 ≤ RT < 50 |
| MS1 spectra retained | 381 | 375 |
| MS2 spectra retained | 5,689 | 5,615 |
| Matching prediction rows | 5,689 | 5,615 |
| Additional precursor MS1 outside interval | 1 | 0 |

The mzML files retain measured binary-array payloads and native scan identifiers. Sequential XML indices were rebuilt after cropping. The prediction CSVs preserve the original unmodified peptide prediction, score and scan metadata; they were not replaced by database-search identifications. All prediction scores are retained in the files, including scores that the workflow subsequently excludes.

`inputs/reference.fasta` contains all 29,673 target entries in the deposited `SIHUMI_DB1UNIPROT.faa`. The 29,673 deposited `DECOY_` entries were excluded so SAGE can generate decoys for each selected database. Target entries were not chosen by inspecting these searches. The reference retains the deposited target and contaminant content; it is not asserted that every entry belongs to a strain actually present. Exact sequence normalization and deduplication yield 29,635 reference sequences.

## Follow the analysis

1. **Reference identity.** The lake builder normalizes sequences, assigns content-derived accessions and preserves the original source headers in a separate mapping table. Identical normalized sequences share an accession. I/L-equivalent sequences can match the same peptide without becoming the same exact reference identity.
2. **Evidence selection.** Eligible prediction rows have positive finite scores, length 9–50 and pass the low-complexity filters. The top 30% rank cutoff is computed within each cropped acquisition; cutoff ties are retained. Pooled containment produces an evidence lake of 538 sequences.
3. **Acquisition-specific selection.** Each acquisition draws sequences using its own retained evidence. Rust razor inference then assigns all eligible mapped inference peptides without the 30% rank restriction, using the documented unique-support, total-support and hash-accession tie rules. A protein survives when it receives at least one assignment. This is not a guaranteed minimum set, and a protein with only shared peptide support can survive.
4. **Search and quantification.** SAGE 0.14.6 searches each mzML against its own selected FASTA, generates decoys and performs LFQ. The exact DDA settings are written to each search configuration and are described in the [Methods](../../docs/concepts/methods.md). The aggregator filters accepted target observations before adding intensity.
5. **Study reporting.** A greedy cover over pooled accepted canonical peptides creates study reporting groups within the searched sequence universe. The output retains both the original sample member-set key and the assigned study key. This does not calculate a new study-wide FDR.

For a worked example of shared peptides and the distinction between razor assignment and greedy set cover, read [Algorithms](../../docs/concepts/algorithms.md). The small synthetic demo is easier for tracing every individual assignment; this example tests the same workflow with measured spectra.

## Expected results and their units

| Measurement | S01 | S02 |
|---|---:|---:|
| Drawn protein sequences | 370 | 316 |
| Razor-selected search sequences | 348 | 287 |
| Target PSM rows, spectrum_q ≤ 0.01 | 758 | 593 |
| Distinct canonical target peptides, peptide_q ≤ 0.01 | 548 | 426 |

The reference run has 819 distinct canonical peptides across the two acquisitions, 498 sequences in the union of searched databases, and 448 study reporting groups. There are 577 rows in its long table of sample and study keys. The included intensity totals 87,045,942,192.60388 before and after study assignment, in the units returned by SAGE LFQ.

These are different counting units. The 348 and 287 selected sequences are search hypotheses, not identified protein counts. PSMs use the spectrum q-value, whereas canonical peptides use the peptide q-value. Canonical peptide counting removes bracketed modifications and equates I with L. A group count is meaningful only with its grouping rule and evidence universe. LFQ totals from cropped windows do not estimate the full acquisitions' protein abundance.

## Read the outputs

| Path inside the chosen output directory | Purpose |
|---|---|
| `VALIDATION.json`, `COMPLETE.json` | Measured results, comparison with expected results, and successful completion |
| `commands.json`, `*.log` | Exact commands, timings and retained diagnostic logs |
| `lake.fasta`, `source_headers.tsv` | Exact reference identities and original-header mappings |
| `workflow/PROVENANCE.json` | Input and executable checksums, parameters and sample paths |
| `workflow/evidence.fasta` | Pooled construction-evidence lake |
| `workflow/samples/S01/` and `S02/` | Drawn databases, selected FASTAs and inference sidecars |
| `workflow/search/S01/` and `S02/` | Search configuration, SAGE identification and LFQ results |
| `workflow/aggregate/` | Filtered acquisition-level aggregation outputs |
| `study/summary.json`, `study/study_group_long.tsv` | Study accounting and one quantitative study-group table |

Open `study/study_group_long.tsv` in a spreadsheet or TSV reader. Each row reports one acquisition/study-group combination, its assigned peptides and features, and intensity. Original local memberships are retained in `feature_assignments.tsv.gz`. The accompanying `expected/reference_dual_key_long.tsv` preserves the original-rule reference execution. It is historical evidence, not the output schema of this branch. The example still checks the unchanged inferred-cover count against that reference; quantified row counts can change.

To check a completed analysis again:

```bash
python examples/campi/check_results.py --run campi_run_01
```

## What PASS means

Input checksums must agree before execution. Reference, pooled evidence, drawn and selected FASTAs must match the expected SHA-256 hashes and sequence counts exactly. Search PSM/peptide and study-group counts may differ by at most 2%, with a minimum allowance of two counts. Each acquisition's accepted canonical peptide set must have Jaccard similarity of at least 0.98 to its reference set. Jaccard is the size of the intersection divided by the size of the union.

Both acquisitions must contribute to the study table, counting keys must be present, and included intensity must be conserved to relative tolerance 10⁻¹². These checks permit small downstream floating-point or threshold differences while requiring identical construction databases. A PASS is a software acceptance result, not a statement that allowed numerical differences are biologically negligible. `expected/results.json` records the exact baseline and tolerances.

If a check fails, read `VALIDATION.json` and the last stage log. Do not edit expected results merely to make the run pass. Verify the bundled source and input checksums, executable versions and parameters first. If setup failed before the workflow began, retain its logs and start again in a new directory after correcting the missing prerequisite. Do not reuse a runtime without its completed `RUNTIME.json`.

## Scientific limits and extension to your data

The retention-time interval was specified before inspecting the new search outcomes. Cropping changes prediction rank distributions, search score calibration, retention-time modelling and available chromatographic peaks. The added precursor MS1 preserves a scan reference; it does not restore full peak boundaries. Therefore these results need not equal the corresponding interval extracted from a full-acquisition analysis. Two windows cannot establish repeatability across laboratories, biological representativeness or independent FDR calibration.

The supplied predictions were generated earlier with AlphaNovo. Running this example does not download a neural-network checkpoint or reproduce raw conversion or model inference. The recorded checkpoint and decoding configuration are in the provenance note; their relationship to model training data is not established here. This example is consequently unsuitable as an independent de novo accuracy benchmark.

For your own complete acquisitions, use [the user guide](../../docs/how-to/acquisition-manifest.md) and `tools/run_manifest.py` with your own manifest and appropriate reference lake. Do not repurpose these expected counts as thresholds for unrelated samples. To recreate the bundled crops from the original full inputs, the optional `prepare_inputs.py` accepts a TSV containing exactly S01 and S02, original prediction CSV and mzML paths, plus `--reference`, `--out`, `--start-minutes 45` and `--end-minutes 50`. Paths resolve relative to the manifest. The bundled predictions and mzML files are already prepared; recreating them is not required to run the example.

## Basic QC output

The normal workflow automatically writes `study/qc/QC_overview.png` and
`.pdf`, sample/group QC tables and `summary.json`.
`study/` is a relative alias of `workflow/study_groups/`.
These two-acquisition examples explicitly skip PCA; studies with at least three
nonempty acquisitions and two complete varying groups also produce
`PCA.png` and `PCA.pdf`. The synthetic three-acquisition example exercises
that path for both sum and directLFQ.
