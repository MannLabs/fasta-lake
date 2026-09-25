# Where performance work will help

Measure the whole stage before choosing Rust, NumPy or Numba. Faster numerical
loops do not make disk scans, annotation database lookups or a model's attention
matrix smaller. The current pipeline already puts its main sequence matching
and razor inference in Rust.

## Measurements so far

These measurements are from low-priority Linux Slurm runs with two Rust threads.
GiB means \(2^{30}\) bytes. Largest-process RSS is not total simultaneous process memory.
They describe the stated inputs, not every reference or instrument.

| Work | Measured input | Result |
|---|---|---|
| Chunked pooled extraction | Original 61.26 GB reference; 233,497,702 records, two full HSTOOL acquisitions | 469 pieces; 1,716 s; 333.285 MiB largest-process RSS |
| First acquisition extraction | Same evidence lake and one full prediction file | 108.5 s; about 446 MiB; 5,803,867 candidate proteins |
| Razor inference alone | Those 5,803,867 candidates, unchanged Rust algorithm | 88.80 s; 8.358 GiB; 11,172 selected proteins in a 16 GiB allocation |
| Full workflow at 8 GiB | Same full reference/input pilot | Failed at razor; the successful extraction does not make this a full-workflow pass |
| Molecular construction | Real PREDICT specimen, DNA+RNA+metabolite selection | 19.89 s; about 0.738 GiB; exact prior FASTA parity |
| Combined small-data acceptance | Cropped real examples, molecular additions, directLFQ and QC | Passed under a 4 GiB job allocation; not a full-acquisition throughput benchmark |

The fresh two-acquisition 11 GiB full workflow passed in 38:52, with an
8.395 GiB largest-process peak, sum/directLFQ and QC. The complete 24-acquisition batch
also passed: 3:35:03, 10.233 GiB largest-process peak, 84.31 GB retained
outputs. Its pooled extraction took 4,849.6 s and 394.5 MiB across 469 pieces.
See [laptop mode](../how-to/compute.md) for the full disk and interpretation scope. RAW conversion, de novo model inference and initial
worldwide catalogue assembly are outside these measurements.

## Optimization map

| Priority | Stage and current implementation | Candidate improvement | Correctness gate |
|---|---|---|---|
| 1 | Rust razor keeps protein strings and both peptide↔protein graphs | Intern IDs once; use compact integer adjacency lists and avoid repeated sequence copies | Identical selected accessions, exact output order/bytes, I/L matching and hash tie behavior across adversarial and full inputs |
| 2 | Chunk orchestration rereads/ranks predictions and builds an automaton for each piece | Long-lived Rust extraction over successive pieces, with one immutable ranked peptide roster | Global ranking and ties identical; no per-piece filtering or inference; input hashes preserved |
| 3 | Full-reference IO and checksums | Profile shared-storage throughput and pass count; combine only scans with equivalent integrity guarantees | Preserve complete-file hashing and before/after mutation detection |
| 4 | Reference-wide molecular support uses the C Aho–Corasick matcher plus Python sets | Stream sequence records; intern digest/peptide IDs; retain compatibility sets only where needed | Identical reference-wide specificity, exhaustive support rows, candidates and gains/losses |
| 5 | Python study graph and repeated TSV parsing | Profile `_rows`, `study_dictionary` and `_quantify`; use integer IDs or columnar batches if dominant | Same greedy ordering, peptide ownership, missingness and deterministic sums |
| 6 | Experimental donor-blocked permutation statistics | Numba can compile the dense within-group distance loop and avoid per-permutation submatrix allocations | Same seeded permutations and blocked label multisets; statistic/p-value parity; no `fastmath` reassociation |
| 7 | QC/PCA and directLFQ | Keep existing compiled NumPy/BLAS and backend operations; limit nested thread pools | Same measured/missing mask, sample roster, PCA method and numerical tolerances |
| 8 | eggNOG/DIAMOND and ESM-C | Batch exact unique selected sequences; reuse only content/version/settings-matched results | Complete annotation/unknown ledgers; fixed database/model identity and pooling/length policy |

Priorities 1–5 remain an engineering map. Priority 6 has an optional tested Numba
kernel. No other new optimization is claimed from this list.
Numba is most promising for repeated numeric loops over fixed-dtype arrays.
It is a poor first choice for CSV dictionaries, string sets, subprocess control
or code already executing in Rust/C/BLAS. Model inference belongs in PyTorch's
device kernels. Test startup/compilation cost as well as steady-state throughput.

## Threads versus independent jobs

Rust threads share a process and its graph. Running four separate inference
processes can duplicate most of the graph memory. On a 32 GiB laptop with an
11 GiB working budget, the measured 8.358 GiB razor case supports one heavy
process at a time. Four available CPU threads do not justify four heavy jobs. `--laptop` therefore uses sequential acquisitions, a bounded
reference piece size and limited helper threads. Docker or Slurm enforces RAM.

## Reproduce a profile before editing

Use representative complete acquisitions, fixed prediction rosters, a new output
directory and the same selected settings. Record source revision, input/binary
hashes, wall/CPU time, peak RSS, thread limits and output disk usage.

On Linux, for a normal run:

```bash
/usr/bin/time -v python tools/run_manifest.py \
  --manifest samples.tsv --lake lake.fasta --out profile_01 \
  --laptop --threads 2 > profile_01.log 2>&1
```

For a Python-only operation, use `python -m cProfile -o stage.prof ...` and inspect
`python -m pstats stage.prof`. Profile Python and native stages separately, because Python
profiling sees the wait for a Rust process, not its internal graph operations.
Include cold-start and warm-start timings when testing Numba. Keep the original
and optimized outputs for direct comparison, including failed cases.

## Measured optional Numba kernel

For 999 within-block permutations of synthetic 80-feature profiles, on one CPU
process: 100 samples took 0.084 s in NumPy, 0.472 s on the first Numba call
(including compilation), and 0.059 s on a warm call. At 400 samples, NumPy took
0.346 s and warm Numba 0.294 s. All permuted statistics agreed to 1e-12 relative
or absolute tolerance; the resulting p-values were identical. NumPy remains the
default because the improvement is modest and first-use compilation costs matter.

## Search and annotation measurements

All six combinations of tryptic/semi-tryptic/unspecific and closed/open searches
passed on both real CAMPI acquisition windows. For 7–30-aa search peptides,
unspecific closed searches peaked at 5.373/4.583 GiB and unspecific open at
3.269/2.810 GiB. The tryptic and semi-tryptic cases stayed below 0.82 GiB in this
example. A separate 7–15-aa teaching variant passed all six choices under 4 GiB.
These are input-specific measurements, not general upper bounds for Sage.

An actual three-sequence annotation-gap pilot used the full eggNOG 5.0.2 database
on disk, two threads and DIAMOND block size 0.5: 1,176 s and a 1.709 GiB
largest-process peak for eggNOG. ESM-C 300M inference for the same three sequences,
including a 1,200-aa synthetic control split into 256-residue windows, took 14.1 s
on CPU after checkpoint validation. These timings do not predict a full cohort's
annotation throughput. Database hashing and model loading have separate costs.
