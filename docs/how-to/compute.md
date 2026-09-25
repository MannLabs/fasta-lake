# What computer do I need?

Try the bundled examples on a laptop. For your own study, check memory and disk
before choosing a machine, then measure one representative acquisition.

## Measured memory, time and disk

Small references can run locally. Large lakes also need memory for peptide
matching and protein selection. Reading the reference in pieces reduces the
extraction footprint, but later stages still need their own allocation.

| Measured workload | Memory | Time and disk |
|---|---|---|
| Pooled extraction from 61.3 / 71.7 GB references, without chunking | About 67 / 81 GiB peak RSS | Two datasets; not a general sizing formula |
| Chunked full workflow, two acquisitions, 61.26 GB reference | 8.395 GiB largest-process peak; 11 GiB allocation | 38 min 52 s; 9.64 GB outputs |
| Chunked full workflow, 24 acquisitions, same reference | 10.233 GiB largest-process peak; 11 GiB allocation | 3 h 35 min; 105.12 GB inputs + 84.31 GB outputs |

These are Linux cluster measurements. They exclude RAW conversion, de novo
prediction and initial reservoir construction. Sequence length, peptide
patterns, match density and search settings affect cost; the numbers do not
predict laptop speed.

A separate 20-acquisition series took about 1.3–4.3 minutes per acquisition for
sample-specific extraction, selection and search. One matched-metagenome search
exceeded 32 GiB and completed with 96 GiB. Keep failed logs when sizing a retry.

Every run needs a new directory. Budget for inputs, retained candidate FASTAs
and outputs, plus container images and build files. The runner preserves
previous results and does not resume a partial directory.

## Laptop or cluster?

| Your task | Starting point |
|---|---|
| Try the examples | Docker on a laptop; example containers use 4 GiB and two CPUs. The initial build needs additional resources. |
| Use a modest reference | Profile one full acquisition locally. |
| Use a large reference with limited RAM | Try chunked extraction, then measure selection and Sage. |
| Build large catalogs or process many acquisitions | Use a cluster; size one job before increasing concurrency. |
| Inspect finished results | Open the study matrix and QC locally. |

For a pilot comparable to the chunked run above, a 32 GiB machine with about
11 GiB available to FastaLake is a plausible starting point. A 16 GiB machine's
default planning share falls below the measured selection requirement.

## Plan and validate a run

After [installation](../get-started/install.md), from the repository root:

```bash
fasta-lake plan-resources
python tools/run_manifest.py --manifest samples.tsv --lake reference_lake.fasta \
  --laptop --out pilot_01 --validate-only
```

Remove `--validate-only` to execute. `--laptop` plans threads and reference
pieces from visible capacity; `--memory-budget-gib 11` sets a planning budget.
It does **not** enforce a memory cap. Docker or the scheduler supplies the hard
limit.

The runner processes acquisitions sequentially. Threads share process data;
separate jobs can duplicate it. See
[resource controls](../reference/resource-controls.md) for overrides and
memory-accounting details, or use the
[assistant's cluster task](assistants.md#ready-to-copy-task-files).

## Build only the evidence lake

```bash
fasta-lake extract-chunked --database reference_lake.fasta \
  --predictions predictions/ --out evidence_01 --chunk-mib 256 --threads 2
```

The result is `evidence_01/evidence.fasta`. Every reference piece uses the same
prediction roster and ranking rules; results merge before selection parsimony.
Chunk size limits reference pieces, not the whole process.

## Full-reference measurements

The two- and 24-acquisition runs used 233,497,702 proteins, two Rust threads and
256 MiB reference pieces. Both completed through sum, directLFQ and QC/PCA
eligibility checks. In the 24-acquisition run, pooled extraction took 80.8 minutes
and peaked at 394.5 MiB; selection parsimony was the memory bottleneck.
An earlier 8 GiB pilot failed at that step.

Inputs and retained outputs totaled 189.42 GB (176.41 GiB), before images,
build caches or original catalog downloads. These runs demonstrate modest-RAM
cluster execution; laptop wall time and Apple Silicon emulation were not measured.
The [detailed resource record](../reference/resource-controls.md) preserves
settings, accounting scope and validation boundaries.

<span id="what-fits-on-a-laptop"></span>
