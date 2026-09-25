# Run FastaLake with Docker

The image is built or downloaded as described in [Install](install.md). This page shows how to run the bundled examples, your own study and the optional stages inside it.

## Run the real-data example

From the repository root on macOS or Linux:

```bash
mkdir -p fastalake-results
docker run --rm --platform linux/amd64 --network none \
  --user "$(id -u):$(id -g)" \
  --mount "type=bind,source=$PWD/fastalake-results,target=/work" \
  fastalake:docker fastalake-example campi --out /work/campi_01 --threads 2
```

A successful run prints `PASS` after comparing its results with the bundled expected results. The output directory on your computer contains `VALIDATION.json`, stage logs, selected FASTAs, search results and study groups. Each run needs a new output directory.

Run the second example using the same image:

```bash
docker run --rm --platform linux/amd64 --network none \
  --user "$(id -u):$(id -g)" \
  --mount "type=bind,source=$PWD/fastalake-results,target=/work" \
  fastalake:docker fastalake-example campi --out /work/campi_01 --threads 2
```

For the small synthetic example, replace the command after `fastalake:docker` with `bash /opt/fastalake/demo/run_demo.sh /work/demo_01`. The expected evidence lake has 14 sequences and the three selected databases contain 4, 4 and 3 sequences.

## Use your own data

Mount your input directory at `/data` with `readonly`, and your output directory at `/work`. Paths inside the manifest must resolve inside the container, either relative to the manifest or using `/data/...` paths. Host paths outside these mounts are unavailable.

```bash
docker run --rm --platform linux/amd64 --network none \
  --user "$(id -u):$(id -g)" \
  --mount "type=bind,source=$PWD/inputs,target=/data,readonly" \
  --mount "type=bind,source=$PWD/fastalake-results,target=/work" \
  fastalake:docker python /opt/fastalake/tools/run_manifest.py \
  --manifest /data/samples.tsv --lake /data/lake.fasta \
  --out /work/study_01 --sage /usr/local/bin/sage --threads 4 --validate-only
```

Remove `--validate-only` to execute. The runner automatically writes grouping, quantification and QC/PCA under `/work/study_01/study_groups/`. Inputs are mzML and de novo prediction files; raw conversion and de novo prediction are upstream steps. A full reference lake still needs the memory recorded in the [compute guide](../how-to/compute.md).

!!! tip "Large references"
    Add `--reference-chunk-mib 256` to the runner command to build the evidence lake in bounded pieces. The pieces merge before per-acquisition inference and molecular additions; the selected databases are identical.

## Quantification, QC and TPM additions

The image includes the quantification and QC dependencies. Optional eggNOG and
ESM-C require their separately configured local tools/resources; see
[annotation](../how-to/annotation.md). After the CAMPI example,
replace the command after `fastalake:docker` in the run pattern above with:

```bash
fasta-lake quantify-study --groups /work/campi_01/study \
  --method directlfq --out /work/campi_lfq --threads 2
```

The command also writes plots under `/work/campi_lfq/qc/`.

For your own molecular evidence, mount the directory containing its source
manifest and referenced files at `/data`, then run:

```bash
fasta-lake molecular add --baseline /data/final_razor.fasta \
  --source /data/specimen_A/source.json --specimen specimen_A \
  --dna-percent 1 --out /work/specimen_A_union
```

The [molecular input walkthrough](../how-to/molecular-inputs.md) explains the source files
and how these inputs were linked and imported in the benchmark.
