---
title: FastaLake
hide:
  - navigation
  - toc
---

# FastaLake

**Compact protein databases, guided by your spectra.**

FastaLake uses de novo peptide predictions to select complete proteins from large
sequence reservoirs. Search the selected databases with Sage, then group and
quantify the peptide evidence across your study.

[Install FastaLake](get-started/install.md){ .md-button .md-button--primary }
[Run your first example](get-started/first-study.md){ .md-button }

<figure markdown="1">
![From peptide evidence to a study result](assets/fastalake-funnel.svg){ width="560" .no-lightbox }
</figure>

## What you can do

- Build an acquisition-specific protein database from public sequence catalogues.
- Run selection, database searching and study-level quantification from a manifest.
- Add matched DNA and RNA evidence when available.
- Use chunked lookup for large reservoirs on smaller machines.
- Explore protein groups, peptide evidence and QC plots.

## Start here

<div class="grid cards" markdown>

-   :material-download:{ .lg .middle } **Installation**

    ---

    Use the bundled Docker image, install from source, or run on a cluster.

    [Choose your installation :octicons-arrow-right-24:](get-started/install.md)

-   :material-play-circle:{ .lg .middle } **Quickstart**

    ---

    Run a small CAMPI example and open your first protein-group matrix.

    [Get a first result :octicons-arrow-right-24:](get-started/first-study.md)

</div>

## Bring your own study

Start with mzML spectra, de novo prediction CSVs and a protein FASTA reservoir.
An acquisition manifest connects the files for each sample. RAW conversion and
de novo inference happen upstream.

[Prepare a manifest](how-to/acquisition-manifest.md) ·
[Choose references](how-to/references.md) ·
[Read the results](how-to/read-outputs.md)

Learn the method through the [visual walkthrough](concepts/walkthrough.md),
follow the [notebook tutorials](tutorials/index.md), or use a
[task file with your coding assistant](how-to/assistants.md).

## Citation and credits

Cite the exact software version used; [CITATION.cff](https://github.com/MannLabs/fasta-lake/blob/main/CITATION.cff)
contains the software citation. The manuscript analyses are in
[fastalake-paper](https://github.com/treitpeter/fastalake-paper).

FastaLake uses [Sage](https://github.com/lazear/sage),
[directLFQ](https://github.com/MannLabs/directlfq) and
[AlphaPeptTools](https://github.com/MannLabs/alphapepttools).
Code is MIT licensed; bundled data and dependencies retain their own terms.
