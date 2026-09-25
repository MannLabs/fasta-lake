# HSTOOL technical-acquisition example

The measured inputs and full-acquisition manifests are not distributed in this
public source package. The runner and expected regression records are retained
for authorized users who already hold the matching inputs. Use the bundled
[public CAMPI example](../campi/README.md) for installation checks.

Run from the software directory:

~~~bash
python3 examples/hstool/run_example.py --out hstool_run_01 --threads 2
~~~

Or reuse the runtime created by the CAMPI example:

~~~bash
python3 examples/hstool/run_example.py --out hstool_run_02 --threads 2 --runtime campi_run_01/_runtime
~~~

Linux x86_64, Python 3.12, Cargo and a C compiler are required for first setup. Python dependencies install from the hash-locked `requirements/container.txt`; Rust dependencies from the committed lockfiles. SAGE 0.14.6 is bundled under ../campi/runtime with its licence. The source and binary fingerprints must match when reusing a runtime. Choose a fresh output directory; failed output is preserved.

## Included inputs and scope

Two measured 45 <= RT < 50 minute windows from TechnicalReplicates 1 and 2 retain 6,429 and 6,432 MS2 spectra, 430 MS1 spectra each (including one preceding precursor spectrum), and one original prediction per MS2. Native scan IDs and measured binary-array payloads are preserved. See inputs/input_provenance.json and the per-spectrum checksum tables. Raw-file SHA-1 metadata was absent from these mzML headers; original mzML and prediction SHA-256 values are recorded instead.

The compact reference contains 706,838 sequences extracted with the reviewed exact matcher from the frozen ten-acquisition HSTOOL evidence lake. The two cropped prediction sets supplied the top-30%, 9–50-residue matching evidence. This reference is therefore evidence-conditioned. Repeating evidence filtering on it retains all its sequences, as expected. It is not the original 61 GB reservoir and cannot test original-reservoir breadth, independent spectral accuracy, full-acquisition yield or whole-procedure FDR.

Reference bytes are stored in six independently compressed parts so every tracked file is below GitHub's per-file limit. The runner reconstructs them in its output directory and verifies the uncompressed SHA-256 before lake construction. All inputs have SHA-256 entries in inputs/SHA256SUMS.

The example runs lake construction, evidence filtering, per-acquisition extraction, Rust razor, Sage, aggregation and post-search study grouping. It does not perform de novo model inference. Technical replicate naming follows the archived acquisition identity; definitive donor/preparation/injection relationships require the source records.

## Expected measurements

| Acquisition | Selected sequences | Target PSMs at spectrum q <= 0.01 | Canonical target peptides at peptide q <= 0.01 |
|---|---:|---:|---:|
| TR01 | 607 | 925 | 544 |
| TR02 | 596 | 858 | 514 |

The baseline has 550 study groups. Expected database hashes and peptide identity sets are included. The verifier requires exact database identities and checks search counts within 2% (at least two counts), peptide Jaccard >=0.98, nonempty output, finite intensity and conservation of included study intensity. The first measured run defined regression expectations; subsequent executions are checked against them. Expectations are not independent biological truth.

Outputs include workflow/PROVENANCE.json, commands and stage logs, searched FASTAs, Sage results and LFQ, aggregate matrices, study/study_group_long.tsv, and VALIDATION.json. COMPLETE.json is written only after the expected-result checks pass.

## Full acquisition records

full_run_records/ contains the portable ten-acquisition input manifest and frozen review settings. These larger inputs are not part of this teaching bundle. The September 6 full-acquisition validation tables are in ../../docs/guide/data. Running cropped windows does not reproduce their counts or CVs.

HSTOOL input redistribution is outside this source package. FastaLake's MIT licence does not relicense source data or sequence catalogues; retain the original catalogue attribution and source provenance.

## Basic QC output

The normal workflow automatically writes `study/qc/QC_overview.png` and
`.pdf`, sample/group QC tables and `summary.json`.
`study/` is a relative alias of `workflow/study_groups/`.
These two-acquisition examples explicitly skip PCA; studies with at least three
nonempty acquisitions and two complete varying groups also produce
`PCA.png` and `PCA.pdf`. The synthetic three-acquisition example exercises
that path for both sum and directLFQ.
