# Choose search settings

This page shows how to choose the digestion rule and search mode for the
manifest runner, and how to compare the results. Cleavage specificity and
precursor mass window answer different questions. FastaLake exposes them
independently in the manifest runner; all combinations use the same database
construction and optional molecular additions.

## Digestion and search mode

| Option | Setting | Intended use |
|---|---|---|
| `--digestion tryptic` | Two tryptic termini; cleavage after K/R, restricted before P | Standard trypsin DDA; default |
| `--digestion semi-tryptic` | At least one tryptic terminus | Inspect peptides with one non-tryptic end |
| `--digestion unspecific` | No enzyme constraint | Peptide processing and endogenous peptide exploration |
| `--search-mode closed` | Precursor ±20 ppm, fragment ±20 ppm | Narrow search; default |
| `--search-mode open` | Precursor −500 to +100 Da, fragment ±20 ppm | Explore mass differences outside declared variable modifications |

Protein termini are valid enzymatic boundaries. Search peptide lengths default
to 5–50; use `--search-min-length` and `--search-max-length` to change them.
These are separate from the de novo evidence length filters.

## Run a declared comparison

From the repository root, after installation and preparing the ordinary manifest:

```bash
python tools/run_manifest.py --manifest samples.tsv --lake lake.fasta \
  --out semitryptic_01 --digestion semi-tryptic --threads 2

python tools/run_manifest.py --manifest samples.tsv --lake lake.fasta \
  --out open_01 --search-mode open --threads 2

python tools/run_manifest.py --manifest samples.tsv --lake lake.fasta \
  --out unspecific_01 --digestion unspecific \
  --search-min-length 7 --search-max-length 30 --threads 2
```

Use `--validate-only` first. Every search writes a complete `config.json`;
the run records selected modes and the Sage binary hash. Open searches retain
fixed C carbamidomethylation, disable variable modifications and isotope-error
expansion, and disable retention-time prediction. LFQ and generated decoys
remain enabled. The mode is a concrete starting configuration, not an
automatically optimized instrument or modification model.

These settings target the bundled Sage 0.14.6. The installed
[`sage_config`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/search.py) function returns a fresh complete JSON
dictionary for callers who need to inspect settings. See the
[Sage configuration documentation](https://sage-docs.vercel.app/docs/configuration)
and [versioned source](https://github.com/lazear/sage/tree/v0.14.6).

## Read the biology carefully

A non-tryptic end can reflect endogenous processing, sample handling or search
error. It does not identify a particular active protease. Report terminal
positions, peptide support, acquisition pairing and the database context.
An open-search mass difference is not automatically a localized modification.

Compare identification gains and losses on the same acquisitions. Broader
settings increase search space, time and often memory; start with a bounded
pilot. A narrow-search memory measurement does not establish the cost of a
semi-tryptic, unspecific or open search. Independent whole-workflow error
calibration remains separate from execution tests.

## Executable teaching check

After the [CAMPI example](https://github.com/MannLabs/fasta-lake/blob/main/examples/campi/README.md), run:

```bash
python examples/searches/run_example.py --campi campi_run_01 \
  --out search_choices_01 --threads 2
```

The check runs all six combinations on the measured windows with 7–30-aa search
peptides. It verifies the target FASTA bytes remain unchanged, acquisition
identity is retained, and each route produces study QC. It records outputs and
does not assume that a broader setting must yield more identifications.
