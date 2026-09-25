# CAMPI example: provenance, modifications and attribution

## Public source and attribution

The measured spectra and deposited reference originate from the Critical Assessment of MetaProteome Investigation project, [PRIDE PXD023217](https://www.ebi.ac.uk/pride/archive/projects/PXD023217). Cite:

Van Den Bossche, T. et al. Critical Assessment of MetaProteome Investigation (CAMPI): a multi-laboratory comparison of established workflows. *Nature Communications* **12**, 7305 (2021). [doi:10.1038/s41467-021-27542-8](https://doi.org/10.1038/s41467-021-27542-8).

The PRIDE project API records the deposit under [CC0](https://creativecommons.org/publicdomain/zero/1.0/). The selected file records and licence field were retrieved on 6 September 2026 and are retained in `public_source_records.json`. The reference includes UniProt-derived content; credit the [UniProt Consortium](https://www.uniprot.org/) and preserve its [CC BY 4.0 attribution notice](https://www.uniprot.org/help/license/). The deposited FASTA and original sequence headers identify the source material. No fresh UniProt download or current database version is implied.

The FastaLake code licence does not replace the dataset or third-party licences. SAGE is distributed under the MIT licence: see `runtime/SAGE_LICENSE.txt` and the licence inside its unmodified official release archive. SAGE source and release attribution are at [lazear/sage v0.14.6](https://github.com/lazear/sage/releases/tag/v0.14.6).

## Original files

| File in PXD023217 | SHA-1 recorded by PRIDE |
|---|---|
| `S01.raw` | `810578a8620e323432464bb1c35d0338b84d9c30` |
| `S02.raw` | `5a53f1679052874d89bd061e51c57bc8b32edcbf` |
| `SIHUMI_DB1UNIPROT.faa` | `65b2440a81c56b9f210f77183e482ffc17e7a277` |

The original deposited files are available through the PRIDE project and its file records. These SHA-1 values are used to connect the inputs to the repository metadata, not as a modern security guarantee. Full SHA-256 digests of the original local mzML, predictions, reference and bundled derivatives are in `inputs/input_provenance.json`.

The two original mzML headers record raw-file SHA-1 values matching the corresponding PRIDE raw-file entries. The original reference FASTA itself was hashed and matched to PRIDE's deposited file checksum. The example derives from S01/S02 mock-community acquisitions; it does not contain the CAMPI F-series faecal acquisitions. The original full raw files are not included because the runnable mzML windows are sufficient for this example.

## Spectrum preparation

The interval 45 ≤ retention time < 50 minutes was selected before the new searches were inspected. Preparation retained all MS1 and MS2 spectra in that interval and any additional MS1 spectrum referenced as a precursor. This adds one preceding MS1 to S01. Measured binary-array payloads and native Thermo scan IDs are unchanged. Spectrum indices were rebuilt sequentially, and the original chromatogram and random-access index sections were omitted to produce non-indexed mzML. Instrument and conversion metadata from the original mzML are retained in each header.

After writing, the files were reparsed and every retained spectrum's binary-array payload was checked against its source. Per-spectrum checksums, native IDs, retention times and precursor references are in `inputs/S01_spectra.tsv` and `inputs/S02_spectra.tsv`. This verifies preservation of the recorded arrays; it is not a fresh validation of the earlier raw conversion.

## Prediction origin and pairing

The original predictions were generated on 31 July 2026 using the AlphaNovo isolation-window implementation recorded as commit `af3ac89b`. The recorded checkpoint filename is `multislope_fragann_comp02_peptdeep_retrained_2pp/model-epoch12-step058020-peprecall0.5156.ckpt`. The saved execution logs identify constrained decoding, ten beams, variable-length decoding disabled, precursor matching in `window` mode, batch size 64 and execution on an NVIDIA L40S. These are prior model outputs supplied as inputs to FastaLake; neither the model weights nor an independently reproducible training pipeline are included.

The AlphaNovo converter uses the integer in the native `scan=` identifier for `spec_idx`. The preparation code selected prediction rows by that native scan number and required exactly one matching row per retained MS2 spectrum. Precursor m/z agreed within 0.059 ppm for both files, consistent with the recorded numeric precision. The six-column bundled CSV preserves `sequence`, `score`, `spec_idx`, `prec_charge`, `prec_mz` and `beam_rank`; `sequence` copies the original unmodified peptide prediction field. Original scores and peptides were not replaced by search results or filtered to successful identifications.

The relationship between these acquisitions and the model's training data has not been established in this review. The example tests the downstream software from fixed predictions and spectra, not independent model accuracy. Full raw conversion, inference and model-training reproducibility require additional upstream assets and metadata.

## Reference preparation and runtime

The original FASTA contains 59,346 records. Preparation retained all 29,673 original target records and excluded only the 29,673 entries whose headers begin `DECOY_`. Original target sequence text and descriptive headers were retained. The runner subsequently normalizes exact sequence identity and produces 29,635 distinct reference sequences with a source-header sidecar. SAGE generates its own `rev_` decoys for each selected search database.

No search-outcome filter was used to construct this target reference. Its deposited organism and contaminant coverage is retained, including reference breadth beyond a claim of one sequence per known mock-community strain. It should not be treated as an exhaustive ground-truth list of proteins actually present.

The bundled SAGE archive is the official Linux x86_64 GNU release for version 0.14.6. Its SHA-256 is `3492d04b9c922d47ee44fd7fd4248b946d8afb7e88893d1d97c403bb3760be30`. FastaLake's Rust executables are compiled from the supplied source with locked dependencies; Python runtime dependencies are installed from the hash-locked `requirements/container.txt` (Python 3.12, Linux x86_64). The runner records the source fingerprint and executable checksums inside the new output directory.
