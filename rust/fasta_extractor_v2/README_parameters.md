# FastaLake pipeline parameters — what every knob does

> **Canonical, exhaustive reference: `--help` on each binary, and [`publication/chain/STAGES.md`](../../publication/chain/STAGES.md)
> for which stage each parameter belongs to.** (A `TUNABLE_PARAMETERS.md` was referenced here but
> never existed in any commit; the link is corrected rather than left dangling.)
> That file lists EVERY flag and config option with its Rust-binary / `config.py` / as-run default and a
> `file:line` citation for each. This page is the per-stage prose summary; if the two ever disagree, the
> code-cited defaults in the source govern.

Every step that turns the sequence lake into a per-sample database exposes parameters. They trade depth
against cross-sample comparability and false-discovery control. Defaults are the published-reproducing
settings. Pairs with `CANONICAL_FASTALAKE_PIPELINE.md` (the fixed recipe).

## Stage 0 — lake build (`lake_builder`)

| Flag | Default | What it does | Effect |
|---|---|---|---|
| `-s TAG:PRIO:PATH[:GLOB]` | (required) | a source. Lower priority keeps the header on a sequence collision. | composition. |
| `--reassign-header-tags TAGS` | none | sources whose headers are replaced by `TAG_<sha256>` ids. Use for metaG/assembly (contig ids collide across samples); catalogue accessions are already unique. | gives identical sequences one id across samples and studies (cross-study join). |
| `--header-tracking-priority N` | 3 | sources with priority <= N get full joint-header tracking in the sidecar. | provenance only. |
| HASH_VERSION | 5.3 (sha256) | the sequence-hashing contract, byte-identical to `fasta_lake/hashing.py`, pinned by `tests/test_hashing_crosslang.py`. | the id; do not change without re-pinning. |

**INPUT ORDER (the D23 lesson, important):** `lake_builder` writes the deduplicated lake in INPUT order
(Phase 2 re-reads sources sequentially; for a single source the output is that file's order, verified).
mmseqs cluster-representative selection is input-order-dependent, so the lake's order changes which
representative wins, and therefore per-sample depth, by up to ~18% with everything else held identical.
**Practical rule:** to reproduce the published depth, build the SHA-256 lake from the native COMPLETE lake
as a single source (preserves its order) rather than by re-unioning component sources (which emit in
source-priority order). No order flag is needed on `lake_builder` itself; the order it is given is the
order it keeps.

## Stage 1 — evidence filter (`fasta_extractor`)

| Flag | Default | What it does | Effect |
|---|---|---|---|
| `--top-percent` | 0.30 | keep the top fraction of each sample's de novo peptides by confidence score. | larger = broader evidence, bigger DB, more de-novo error; smaller = fewer, most-confident calls. |
| `--min-length` | 9 | minimum de novo peptide length admitted. | lower (7) = more ids and higher FDP (short peptides match by chance); higher (11) = lower FDP. |
| `--max-length` | 50 | maximum peptide length. | rarely tuned. |
| `--normalize-il` | true | treat I and L as equivalent (they are isobaric). | matching. |
| `--entrapment-method {shuffled,noble}` `--entrapment-n` | off | inject entrapment proteins for FDR validation. **Inject BEFORE selection (into the lake) for a fair FDP**; appended after selection it cannot be nominated and the FDP is optimistic. | FDP estimate. |

## Stage 1.5 — clustering (`mmseqs easy-cluster`)

| Flag | Default | What it does | Effect |
|---|---|---|---|
| `--min-seq-id` | 0.97 | identity threshold to cluster (cluster97). | higher (0.99) = more, smaller clusters, more depth, less reproducibility; lower (0.90) = orthologs collapse, more reproducible, less strain resolution. |
| `-c` | 0.8 | coverage threshold (fraction of residues aligned). | mild; higher = stricter merging. |
| `--cov-mode` | **0** | what `-c` is measured against: 0 = query AND target, 1 = target only, 2 = query only. **Canonical FL = 0** (the published-reproducing build; the earlier "1" was a different complete-lake sweep). | cov-mode 0 keeps fragments separate (more reps, slightly more depth); 1 collapses them. |
| `--shuffle` (mmseqs flag, **NOT a FastaLake binary flag** — wired only in `projects/predict_orderpreserve_real/scripts/03_cluster.sh` via `MMSEQS_SHUFFLE`; absent from the e2e harness and from `config.py`) | mmseqs default 1; that script defaults it to **0** | mmseqs shuffles the input DB by default, which is what makes the result input-order-sensitive. `--shuffle 0` makes mmseqs respect the input order deterministically. | **`--shuffle 0` = the order-preserving / reproducible setting** in `03_cluster.sh` (`MMSEQS_SHUFFLE=1` restores mmseqs default). |

## Stage 3 — inference (`parsimony_engine --strategy`)

| Strategy | What it does |
|---|---|
| `razor` (default) | shared peptide → the protein with the most unique peptides. |
| `uniform` | keep every protein with >= 1 assigned peptide (inclusive). |
| `species_budget` | cap shared-peptide spread by the number of unique-validated species (UHGG/MGYG taxonomy); preserves strain specificity, lower depth. |
| `uniform_2pep` | require >= 2 peptides (conservative, FDP floor). |
| `info_score` | score = sum(denovo x weight)/n_proteins per peptide; degenerate peptides count less. |
| `bayesian` | EM posterior over iterations. |

Plus `--min-length 9 --max-length 50 --normalize-il` (as Stage 1).

## Stage 5 — aggregation (`aggregation_engine`)

| Flag | Default | What it does | Effect |
|---|---|---|---|
| `--group-key-mode {member-set,representative}` | member-set | how a protein group is keyed: full member set (max count) vs a single canonical anchor (stable across samples). | member-set = per-sample granularity; representative = cross-sample comparability. |
| `--q-thresholds` | 0.01,0.05,unfiltered | protein/peptide q-values at which matrices are emitted (free; one run). | depth vs false discovery, read from the same run. |
| `--bh-correction` | **RETIRED** | hard error; SAGE q-values are already FDR-controlled and aggregation tests no new hypothesis. | use entrapment FDP instead. |

## Evidence gate (`fasta_lake.multi_omics`, Methods M5.5; PREDICT only)

`de-novo-only` (default, no omics) vs `OR` (de novo or matched metaG/metaT) vs `AND` (de novo and omics).
Toward depth: de-novo-only / OR. Toward lowest FDP: AND.
