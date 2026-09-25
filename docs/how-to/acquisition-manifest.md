# Prepare an acquisition manifest

FastaLake starts from one de novo prediction file and one mzML per acquisition, listed in a tab-separated manifest. This page is the contract for that manifest and for the prediction files, followed by the preflight and the run.

## The manifest

Create a tab-separated file with these columns:

```text
sample	predictions	mzml
Sample01	data/Sample01.csv	data/Sample01.mzML
Sample02	data/Sample02.csv	data/Sample02.mzML
```

The separators above must be actual tabs when saved. Paths resolve relative to the manifest's directory. Absolute paths also work. Each sample identifier must be unique and use letters, digits, underscores, dots or hyphens, starting with a letter or digit. Do not use a subject identifier alone when that subject has several acquisitions.

## Prediction files

The supported prediction headers are `sequence,score` or `peptide_prediction_detokenized_unmodified,score`. Supply uppercase, unmodified amino-acid sequences. The Rust extractor accepts some numeric modification decorations, but the inference parser does not share that cleaning contract. Use unmodified sequences throughout this workflow. Named modifications and general ProForma syntax are not a supported de novo interchange format. Convert such predictions explicitly and retain the conversion record. In particular, simply dropping punctuation from a named modification can change the apparent peptide sequence.

Scores must be numeric and higher must mean better. This is essential because the extraction rule ranks scores in descending order. A score whose meaning is a loss or negative log likelihood cannot be supplied without a documented transformation. Rows with missing, nonfinite or nonpositive scores are not eligible for extraction. A top fraction is a within-acquisition rank filter, not an error probability and not a calibrated confidence threshold.

## Pairing predictions with spectra

Check the prediction-to-mzML pairing using acquisition metadata. The runner verifies file existence and identity by checksums, but cannot determine from a plausible filename whether a CSV belongs to the correct injection. For matched metagenome comparisons, verify the additional link from acquisition to biological sample and DNA assembly. A valid file path does not establish a correct biological pairing.

## Preflight

```bash
python tools/run_manifest.py --manifest samples.tsv --lake lake_build_01/lake.fasta \
  --out study_run_01 --sage /path/to/sage --threads 4 --validate-only
```

Preflight checks the manifest structure, unique identifiers, nonempty input files, prediction headers, required binaries, SAGE version and a fresh output destination. It does not parse every spectrum or guarantee that a large job fits memory. It also does not establish that a reference catalogue covers the sample.

## Execute

```bash
python tools/run_manifest.py --manifest samples.tsv --lake lake_build_01/lake.fasta \
  --out study_run_01 --sage /path/to/sage --threads 4 \
  --top-fraction 0.30 --min-length 9 --max-length 50
```

The runner hashes all inputs, including the full lake, before starting. Hashing a large FASTA adds I/O time and provides an exact input identity. The command records parameter values and executable hashes in `PROVENANCE.json`, and stage commands, return codes and elapsed times in `commands.json`.

For database construction without searching, add `--databases-only`. The current manifest format still requires the mzML path so the intended acquisition pairing is recorded. `COMPLETE.json` states whether searching was omitted. A databases-only completion must not be described as an end-to-end search validation.
