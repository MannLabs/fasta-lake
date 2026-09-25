# Add matched DNA and RNA with Docker

This page adds matched DNA and RNA measurements to a study run in Docker. It
follows the [first-study tutorial](../get-started/first-study.md), using its
`fastalake:docker` image. De novo selection remains the baseline; matched
molecular candidates are added before Sage searches the final FASTA.

## Prepare the local inputs

Use this project layout, with your actual files:

```text
inputs/
  lake.fasta
  samples.tsv
  raw/
    run_A_predictions.csv
    run_A.mzML
    specimen_A_genes.faa
    specimen_A_tpm.tsv
    specimen_A_reference.json
```

The [input specification](molecular-inputs.md#1-prepare-one-source-per-specimen)
defines the shared gene IDs, DNA/RNA columns and provider provenance record.
Unavailable measurements are missing, not zero. The commands below require
matching TPM/protein gene populations; declared incomplete references have a
separate [coverage policy](molecular-inputs.md).

Create destinations and prepare a checked, portable source:

```bash
mkdir -p prepared results
docker run --rm --platform linux/amd64 --network none \
  --user "$(id -u):$(id -g)" \
  --mount "type=bind,source=$PWD/inputs,target=/data,readonly" \
  --mount "type=bind,source=$PWD/prepared,target=/prepared" \
  fastalake:docker fasta-lake molecular prepare-source \
  --specimen specimen_A --reference-id assembly_A_v1 \
  --tpm /data/raw/specimen_A_tpm.tsv --proteins /data/raw/specimen_A_genes.faa \
  --provenance /data/raw/specimen_A_reference.json --out /prepared/specimen_A
```

## Link the acquisition and check the run

Save this header and row in `inputs/samples.tsv` with tab separators:

```text
sample	predictions	mzml	molecular_source	specimen
run_A	raw/run_A_predictions.csv	raw/run_A.mzML	../prepared/specimen_A/source.json	specimen_A
```

Paths resolve beside the manifest. Repeat the source/specimen pair for additional
acquisitions of that specimen; use distinct acquisition IDs.

```bash
docker run --rm --platform linux/amd64 --network none \
  --memory 11g --memory-swap 11g --cpus 2 \
  --user "$(id -u):$(id -g)" \
  --mount "type=bind,source=$PWD/inputs,target=/data,readonly" \
  --mount "type=bind,source=$PWD/prepared,target=/prepared,readonly" \
  --mount "type=bind,source=$PWD/results,target=/work" \
  fastalake:docker python /opt/fastalake/tools/run_manifest.py \
  --manifest /data/samples.tsv --lake /data/lake.fasta --out /work/study_01 \
  --laptop --quantification directlfq \
  --molecular-dna-percent 1 --molecular-rna-percent 1 --validate-only
```

The two 1% rules are illustrative choices, applied independently to positive
DNA and RNA values. Choose and record your rules before evaluating outcomes.
Set a layer's percentage to zero to disable it. An 11 GiB container limit does
not guarantee that every reference or molecular union fits; consult the
[resource guide](compute.md).

## Execute and inspect

Remove `--validate-only` from the same command to execute. Use a new output name
for each run. Inputs and prepared sources stay mounted read-only.

Under `results/study_01/`:

| Output | Meaning |
|---|---|
| `samples/run_A/molecular/search.fasta` | Preserved de novo targets plus exact additions |
| `samples/run_A/molecular/gene_evidence.tsv` | Each gene's measurements, selection and sequence identity |
| `search/run_A/` | Sage results and the actual search configuration |
| `study_groups/quantification/study_group_matrix.tsv` | directLFQ quantities for the fixed study groups |
| `study_groups/qc/` | Coverage, missingness and PCA when the study is eligible |

One acquisition receives QC but cannot support PCA. The included
[multi-omics demo](https://github.com/MannLabs/fasta-lake/blob/main/examples/multiomics/README.md) provides a tested teaching
workflow with real spectra and explicitly synthetic molecular values.
