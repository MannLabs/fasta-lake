# Annotate your selected proteins

FastaLake can run local eggNOG annotation and then generate ESM-C embeddings for
proteins without your requested functional terms. These are optional steps after
searching and grouping. Normal construction does not download or load either backend.

## Inputs you already have

A completed study includes `study_groups/representatives.fasta`: exactly the
reporting representatives selected by grouping, with original sequence identity.
The summary records its hash. The peptide and protein-member ledgers remain the
source for ambiguity; an annotation on the representative does not establish the
origin of every shared peptide or annotate every compatible member.

You can also supply your own target-only selected FASTA with `annotate eggnog`.
Exact sequence aliases are collapsed for annotation and preserved in `aliases.tsv`.
I/L variants remain different proteins. The default guard allows 100,000 input
records; this step is for selected proteins, not automatic annotation of an entire
worldwide reference catalogue.

## Set up local tools and data

Use eggNOG-mapper 2.1.12 with database 5.0.2 and its DIAMOND backend, following
[the upstream installation instructions](https://github.com/eggnogdb/eggnog-mapper/wiki/Installation).
The mapper can use its own Python environment through `emapper_python`.
FastaLake checks its reported version and supplies the database directory explicitly.

Download the three required databases and, optionally, model weights with
[the resource helper](references.md). The database remains on disk;
`--dbmem` is not used. The default DIAMOND block size is 0.5, with two threads,
`sensitive` search, automatic taxonomic scope and all target orthologues. Automatic
scope accommodates mixed microbial/host inputs; it does not silently restrict
every sample to bacteria. A specific taxonomic scope is an explicit analysis choice.

For ESM-C, install the optional Python dependencies in the environment running
FastaLake (Python as required by the [installation page](../get-started/install.md#source-installation)):

```bash
python -m pip install '.[analysis,directlfq,embeddings]'
```

This pins `esm==3.2.3` and `torch==2.13.0`. The standard Docker image contains
construction, Sage, directLFQ and AlphaPeptTools; ESM/eggNOG software and their
large external resources are optional additions. The commands are available in
the package and explain missing dependencies; weights are never fetched during analysis.

## Write the settings file

Save `annotation.json` beside your local resource folder. Paths are resolved
relative to this JSON file, so the same configuration works from another directory.

```json
{
  "emapper": "/path/to/eggnog-mapper/emapper.py",
  "emapper_python": "/path/to/eggnog-environment/bin/python",
  "data_dir": "resources/eggnog-data",
  "threads": 2,
  "sensitivity": "sensitive",
  "tax_scope": "auto",
  "unknown_fields": ["KEGG_ko", "EC"],
  "generate_embeddings": true,
  "model": "esmc_300m",
  "checkpoint": "resources/esmc-300m/esmc_300m_2024_12_v0.pth",
  "device": "cpu",
  "max_residues": 1024
}
```

Set `generate_embeddings` to `false` to run annotation alone. An **unknown** means
none of the selected fields has an assignment. The defaults are KO and EC;
GO, PFAM, COG category and Description are other explicit choices. A missing
annotation is not evidence of a novel sequence or function.

## Run annotation on a study

Run on an existing study:

```bash
fasta-lake annotate study --groups results/study_groups \
  --settings annotation.json --out results/annotation
```

Or add `--annotation-settings annotation.json` to your normal
`python tools/run_manifest.py ...` command. After searching, grouping and QC,
the runner annotates its verified representatives and records the annotation
completion receipt. Invalid settings fail before the workflow starts.

## Standalone FASTA and embedding routes

```bash
fasta-lake annotate eggnog --fasta selected.fasta --out annotation \
  --emapper /path/to/emapper.py --emapper-python /path/to/python \
  --data-dir resources/eggnog-data --threads 2

fasta-lake annotate embed-unknown --annotation annotation --out embeddings \
  --checkpoint resources/esmc-300m/esmc_300m_2024_12_v0.pth \
  --model esmc_300m --device cpu --threads 2
```

`annotate eggnog --embed-unknown --checkpoint ...` combines both steps directly.
Outputs are new directories; failed runs retain their logs and lack the applicable
completion receipt.

## Understand the outputs

| File | Meaning |
|---|---|
| `queries.fasta`, `aliases.tsv` | Exact unique sequences and original accession mapping |
| `annotation.emapper.annotations` | Full original eggNOG annotations |
| `annotations.tsv` | Every input query, including no hits and hits without requested terms |
| `unknown.fasta` | Exact query roster lacking the declared terms |
| `PROVENANCE.json`, `RESOURCES.json`, `COMPLETE.json` | Tool/database/input hashes, settings, measured usage and completion |
| `embeddings/vectors.tsv`, `embeddings/vectors/*.npy` | Sequence-to-vector map and numeric vectors |
| `STUDY_COMPLETE.json` | Completed study annotation, tied to its dictionary and settings |

ESM-C uses all residues through nonoverlapping windows. It excludes BOS/EOS
special tokens and combines windows by residue count. This avoids silently
truncating long proteins, but changes long-range context: comparison vectors must
use the same model, window size and pooling. Old prefix-only/all-token embeddings
are not interchangeable. CPU is the default; CUDA must be chosen explicitly.

Embeddings are representations for subsequent analysis, not inferred functional
labels. The historical selected ESM-C manuscript panel remains a separate,
closed scientific analysis; adding reusable software does not reopen its claims.
