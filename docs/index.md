---
title: FastaLake
hide:
  - navigation
---

# FastaLake

FastaLake uses de novo peptide evidence to construct compact protein databases
for proteomics and metaproteomics. Complete proteins are selected from large
sequence repositories for each mass spectrometry acquisition. The resulting
databases are searched with [Sage](https://github.com/lazear/sage), followed by
protein grouping and quantification within each study. Matched metagenomic and
metatranscriptomic sequences can be included when available.

<figure markdown="1">
![FastaLake selection funnel](assets/fastalake-funnel.svg){ width="560" .no-lightbox }
</figure>

## Getting started { #run-it-in-30-seconds }

The Docker image bundles Python, the compiled Rust tools, Sage, directLFQ, AlphaPeptTools and the example inputs.

```bash
git clone https://github.com/MannLabs/fasta-lake.git
cd fasta-lake
bash tools/check_local_docker.sh docker_check_01
```

The first build downloads dependencies. The examples then run offline within 4 GiB of RAM and two CPUs. Success ends with `PASS`, and every table, plot and log stays in `docker_check_01/`.

## Set up with an assistant

Working with an agent or LLM? Give it the repository's `AGENTS.md` and a
[task file](how-to/assistants.md#ready-to-copy-task-files) for what you want to do.
Each task names the inputs to supply and the outputs to check, with example
files to adapt.

## Documentation guide

<div class="grid cards" markdown>

-   :material-play-circle:{ .lg .middle } **Example data**

    ---

    Every bundled example, what it produces and how to read the result.

    [:octicons-arrow-right-24: The demo](get-started/demo.md)

-   :material-flask:{ .lg .middle } **Study analysis**

    ---

    Acquisition manifest, preflight, result folders, first study.

    [:octicons-arrow-right-24: Your first study](get-started/first-study.md)

-   :material-dna:{ .lg .middle } **Matched metagenome or metatranscriptome**

    ---

    One FASTA and one TPM table per specimen, added after de novo selection.

    [:octicons-arrow-right-24: Molecular inputs](how-to/molecular-inputs.md)

-   :material-memory:{ .lg .middle } **A large reference and limited memory**

    ---

    Chunked extraction, measured limits, laptop or cluster.

    [:octicons-arrow-right-24: Compute choices](how-to/compute.md)

-   :material-robot:{ .lg .middle } **A cluster, with an assistant**

    ---

    Scheduler-first prompts for ChatGPT, Claude or a coding agent.

    [:octicons-arrow-right-24: Working with an assistant](how-to/assistants.md)

-   :material-help-circle:{ .lg .middle } **Methods and interpretation**

    ---

    The principle, the algorithms and the output contract.

    [:octicons-arrow-right-24: Concepts](concepts/principle.md)

</div>

## How it is tested

Continuous integration runs on every push: the six Rust crates (fmt, clippy, tests), the Python suite under Python 3.10 and 3.12 with warnings treated as errors, the analysis and directLFQ suites, dependency and licence checks, the source manifest, and the public CAMPI SIHUMIx acquisition windows end to end (PRIDE data). The Docker workflow builds the image and runs every bundled example offline. [Validation](concepts/validation.md) records what the examples check; [Development](project/development.md) records what remains open.

## Credits

FastaLake uses [Sage](https://github.com/lazear/sage), [directLFQ](https://github.com/MannLabs/directlfq) and [AlphaPeptTools](https://github.com/MannLabs/alphapepttools). Please credit them alongside the FastaLake revision and settings. Code is under the MIT licence. Bundled data and third-party dependencies keep their own terms. A manuscript is in preparation. Cite the repository revision until a DOI exists.

## AI-assisted development

Much of FastaLake's code and documentation was developed with Codex Astra
(reasoning effort: xhigh), Opus 4.8, Opus 5 and Fable. The authors remain
responsible for the software, design decisions, validation and scientific
claims. AI assistance is not evidence that an implementation is correct.
The [validation documentation](concepts/validation.md) describes the checks
and their limits.
