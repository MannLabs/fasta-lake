# Resource controls and measured scope

For choosing a machine, start with [What computer do I need?](../how-to/compute.md).

## Run your study

Use the Docker image or the source installation (Python as required by the
[installation page](../get-started/install.md#source-installation)) with analysis
dependencies and the Rust/Sage executables:

```bash
fasta-lake plan-resources
python tools/run_manifest.py --manifest samples.tsv --lake reference_lake.fasta \
  --laptop --out laptop_run_01
```

Laptop mode detects visible memory, CPU affinity, container ancestor limits and
Slurm allocations. It starts with 35% of visible RAM, at most 80% of currently
available memory, and half the available CPU slots. It enables reference pieces
of at most 256 MiB, smaller when the planning budget is small. One Rust process
uses the chosen threads; acquisitions and workflow steps run sequentially.
BLAS and NumPy helper threads are limited to one to avoid nested CPU pools.

On an otherwise available eight-core, 32 GiB machine, this chooses four Rust
threads and an 11.2 GiB planning budget. To request 11 GiB explicitly:

```bash
python tools/run_manifest.py --manifest samples.tsv --lake reference_lake.fasta \
  --laptop --memory-budget-gib 11 --out laptop_run_02
```

Use `--memory-fraction 0.5` or `--cpu-fraction 0.6` to change the shares.
`--threads` and `--reference-chunk-mib` override their automatic choices within
the visible CPU slots and reference-piece budget. `RESOURCE_PLAN.json` records
the observations and final settings; `--validate-only` prints them without
creating an output directory. Explicit thread/chunk options remain available
without laptop mode.

!!! warning
    The planning budget is not an enforced memory limit. A Rust thread shares
    much of its process's data; a second independent job duplicates those structures.
    Dividing 11 GiB by one matcher's small-input RSS does not establish how many
    razor or Sage jobs fit. Those stages need full-input measurements and headroom.
    For a hard total limit use Docker `--memory 11g --memory-swap 11g` or Slurm
    `--mem=11G`. In that constrained environment, automatic planning sees the smaller
    allocation. Capacity and current availability can change after preflight.

This runs chunked pooled extraction, chunked per-acquisition extraction, razor,
Sage, study grouping, sum quantification and automatic QC/PCA. Add
`--quantification directlfq` for directLFQ. Specimen-qualified metaG/metaT
columns and molecular options work through the same runner: additions still
follow the completed de novo razor baseline and precede Sage.

The reference must already be a prepared lake with its intended sequence/header
identities. This option does not redesign initial catalogue normalization and
global deduplication, run a de novo prediction model, or convert RAW files.

## Build only the evidence lake

```bash
fasta-lake extract-chunked --database reference_lake.fasta \
  --predictions predictions/ --out evidence_01 --chunk-mib 256 --threads 2
```

Use `evidence_01/evidence.fasta` as the merged product. The full acquisition
manifest, source identities and selection settings remain necessary for
interpreting or reusing a study-specific evidence lake.

- Each reference piece contains complete FASTA records. Its size can exceed the
  MiB target by its final record; a separate 500,000-record limit bounds the
  per-record offset table. `--chunk-records` adjusts that standalone limit.
- Every piece uses the same complete prediction roster. Prediction ranking
  remains per sample with the same cutoff ties. Predictions are not ranked
  separately based on which proteins happen to be in a piece.
- The existing Rust matcher produces each selected piece. Their FASTA and
  evidence rows merge in original reference order. Exact sequences, duplicate
  records and header bytes behave as in unchunked extraction; no extra
  deduplication, clustering or chunk-local razor is introduced.
- I/L equivalence applies to peptide matching. Original protein sequences and
  headers are preserved.
- A piece with no matches is valid. A full study stops if its merged candidate
  lake is empty.
- Input/output hashes, commands, piece statistics and logs accompany the merged
  FASTA and protein evidence table. `--output-hits` also writes peptide hit rows.
  No successful completion receipt is written after an input change or failure.

## Memory, time and disk are separate limits

The MiB option bounds the reference pieces fed to the matcher; it is not a
total-memory cap. The peptide matching index, all retained prediction evidence,
dense hits within one piece, per-acquisition razor, Sage and quantification need
their own memory. Work proceeds sequentially across pieces and acquisitions;
two threads provide a conservative starting point.

Re-reading/ranking the prediction roster for each piece and hashing complete
inputs adds work. Smaller pieces can save memory while taking longer. An
overnight workload must be measured with representative full acquisitions;
cropped teaching examples do not establish full-cohort runtime.

Successfully consumed scratch reference pieces are removed. Selected piece
copies are also reclaimed after being appended to the merged product. Piece
receipts/logs remain; `--keep-pieces` retains selected FASTAs/tables for inspection
at the cost of duplicate storage. The input reference, predictions, spectra and
earlier output directories are preserved. Failed outputs remain diagnosable.
Per-acquisition drawn lakes and normal search outputs also consume disk.

On Linux/macOS, `commands.json` records stage wall time, CPU time and the largest
process RSS reported by the operating system. Chunk receipts record the same
measurements for each matcher invocation. This RSS value is not simultaneous
memory summed across all processes. Scheduler memory allocations and enforcement
records are separate evidence.

## Full-reference measurements

These low-priority Linux Slurm runs used two Rust threads, 256 MiB reference
piece targets and the unchanged reference of 233,497,702 protein records
(61.26 GB). The settings were explicit; automatic laptop planning was tested
separately. Both completed through sum, directLFQ and QC/PCA eligibility checks.

| Full workflow | Requested RAM | Elapsed time | Largest-process peak | Retained outputs |
|---|---|---|---|---|
| Two full acquisitions | 11 GiB | 38 min 52 s | 8.395 GiB | 9.64 GB |
| 24 full acquisitions: 12 technical + 12 workflow replicates | 11 GiB | 3 h 35 min 03 s | 10.233 GiB | 84.31 GB |

The 24-acquisition pooled extraction took 80.8 minutes across 469 pieces and
peaked at 394.5 MiB. Razor was the memory bottleneck. An earlier 8 GiB pilot
failed during razor; its failed receipt remains part of the resource record.
Chunking alone does not make that stage fit a smaller allocation.

The full batch's prepared inputs occupied 105.12 GB, including the reference.
Inputs plus retained outputs therefore occupied 189.42 GB (176.41 GiB),
before Docker images, build caches, original catalogue downloads or other work.
Plan disk space as carefully as RAM. Consumed reference pieces are reclaimed,
but the per-acquisition candidate lakes are retained for inspection.

This demonstrates a few-dozen-acquisition run with modest RAM on a cluster.
It does not measure laptop wall time or Apple Silicon emulation. A 32 GiB
machine with an approximately 11 GiB working budget is a plausible starting
point for a comparable pilot; a 16 GiB machine's default 35% budget is smaller
than the measured razor requirement. Smaller references or cluster selection
may be more practical there. Other inputs and wider searches can need more RAM.

The resource acceptance record, in the private review record that is not part of
this repository, retains source identities, all 24 acquisitions, timings, memory-accounting scope,
PCA eligibility and the failed 8 GiB pilot. Slurm's requested allocation, its
actual cgroup limit and largest-process RSS are recorded separately. RAW
conversion, de novo model inference and initial global catalogue construction
are outside this profile.

## What has been checked

The chunked and unchunked routes produce identical selected FASTA bytes, peptide
hit tables and per-protein evidence under changes of piece boundaries, including
ties, I/L behavior, multiline/CRLF input, repeated sequence headers and no-hit
pieces. Both measured teaching examples additionally reproduce the per-acquisition
razor FASTAs and resulting study matrices, with automatic QC.

The recorded test counts are on the [validation page](../concepts/validation.md). Separate installed-package and
Docker checks execute measured teaching examples, additive metaG/metaT searches,
all six search choices, directLFQ, QC and eligible PCA.
