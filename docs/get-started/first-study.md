# Your first FastaLake result

This tutorial takes two measured CAMPI acquisition windows through database
selection, Sage searching and a study-level result. The small bundled example
lets you learn the workflow before processing full acquisitions.
For the scientific rationale and terminology, see
[the principle behind FastaLake](../concepts/principle.md).

## Prepare the environment

Build the image once as described in [Install](install.md), then create a results directory:

```bash
mkdir -p fastalake-results
```

## Run measured spectra

```bash
docker run --rm --platform linux/amd64 --network none \
  --memory 4g --memory-swap 4g --cpus 2 \
  --user "$(id -u):$(id -g)" \
  --mount "type=bind,source=$PWD/fastalake-results,target=/work" \
  fastalake:docker fastalake-example campi --out /work/campi_01 --threads 2
```

The inputs are real spectra and matching predictions, cropped to two acquisition
windows. The command checks its outputs against the bundled expectations and
records `VALIDATION.json`. Reuse the image for your next run; choose a new output
name to preserve the first result.

## Open the result

Under `fastalake-results/campi_01/`:

| Output | What to look for |
|---|---|
| `study/qc/QC_overview.png` | Intensity, measured groups and completeness across acquisitions |
| `study/study_group_matrix.tsv` | Group quantities; a blank cell is unquantified |
| `study/peptide_evidence.tsv.gz` | Peptide identities supporting each reporting group |
| `workflow/samples/S01/inference/S01_razor.fasta` | The compact target database selected for S01 |
| `workflow/search/S01/config.json` | The exact search settings |
| `workflow/search/S01/results.sage.tsv` | PSM-level search results |
| `VALIDATION.json` | Checks, counts and scope of the example |

PCA needs at least three nonempty acquisitions and two complete, varying groups.
This two-acquisition example gets QC and an explicit PCA exclusion. Your eligible
study also gets `PCA.png`, `PCA.pdf` and coordinate tables automatically.
The separately labelled [synthetic analysis example](https://github.com/MannLabs/fasta-lake/blob/main/examples/analysis/run_example.py)
exercises that three-acquisition route.

## Use your own acquisitions

Create `inputs/samples.tsv` with real tab separators and one row per
acquisition. Repeated measurements of a specimen still have distinct sample IDs.
The [bundled acquisition manifest](https://github.com/MannLabs/fasta-lake/blob/main/examples/campi/inputs/samples.tsv) is a
working TSV you can copy and adapt.

```text
sample	predictions	mzml
sample_01	sample_01_predictions.csv	sample_01.mzML
sample_02	sample_02_predictions.csv	sample_02.mzML
```

Place the files and prepared `lake.fasta` in `inputs/`. Paths are relative to the
manifest. Predictions are an AlphaNovo CSV or a generic `sequence,score` CSV, with peptide
sequence and confidence fields described in [Prepare an acquisition manifest](../how-to/acquisition-manifest.md);
DIANovo output is converted first, as described in [Reuse existing DIANovo predictions](../how-to/de-novo-inputs.md).
FastaLake starts from these predictions; it does not run the de novo model.

Check inputs, executables and settings before starting:

```bash
docker run --rm --platform linux/amd64 --network none \
  --memory 11g --memory-swap 11g --cpus 2 \
  --user "$(id -u):$(id -g)" \
  --mount "type=bind,source=$PWD/inputs,target=/data,readonly" \
  --mount "type=bind,source=$PWD/fastalake-results,target=/work" \
  fastalake:docker python /opt/fastalake/tools/run_manifest.py \
  --manifest /data/samples.tsv --lake /data/lake.fasta \
  --out /work/study_01 --laptop --validate-only
```

Remove `--validate-only` to execute. Docker enforces the memory limit; `--laptop`
chooses conservative threads and reference pieces using the visible capacity.
The [compute guide](../how-to/compute.md) explains what the measured limits cover.

!!! warning
    An 11 GiB allocation is not a guarantee that every dataset fits. Keep the failed logs if a run is killed for memory and size the next attempt from the measured peaks.

The runner writes the matrix under `study_01/study_groups/` and QC under
`study_01/study_groups/qc/`. Add `--quantification directlfq` to choose directLFQ;
the default is sum. Matched DNA/RNA additions require the extra specimen/source
columns described in [Molecular inputs](../how-to/molecular-inputs.md).

### Three terms used in the results

- **Acquisition:** one measured mzML file. A specimen can have several acquisitions.
- **Razor database:** a compact set of complete proteins selected to explain the
  retained de novo peptide evidence. This is the target database Sage searches.
- **Study group:** one reporting representative with its assigned peptide evidence
  across acquisitions. Shared peptides can leave the biological source ambiguous.

## When something needs attention

The [troubleshooting page](../how-to/troubleshooting.md) lists the first messages a run can produce and what to do about each.

