# Explore functional profiles

These optional helpers turn a fixed study matrix and your explicit annotations
into interpretable summaries. Basic QC/PCA already runs with normal FastaLake
quantification; this page covers additional experimental functional analyses.

The design follows AlphaPeptTools' readable, functional style. Its verified PCA
backend is used directly. FastaLake supplies the metaproteomics-specific accounting
for ambiguous term mappings, unannotated signal and repeated-measure designs.

## Choose a matrix and one annotation level

Use either the sum matrix or the selected directLFQ matrix. The workflow does not
requantify, impute, batch-correct or change peptide assignments. Supply a long TSV:

```text
study_group	term
representative_A	K00001
representative_A	K00002
representative_B	K00001
```

Use real reporting-group accessions from your matrix. Each row declares one mapping
at one level, such as KO, EC or genus; do not mix levels in a profile. Duplicate
pairs count once. A group with two terms contributes half its observed signal to
each. This is a visible attribution convention, not evidence for equal catalytic
contributions. A representative's annotation does not resolve shared-protein origin.

```bash
fasta-lake explore-study \
  --matrix results/study_groups/quantification/study_group_matrix.tsv \
  --annotations group_terms.tsv --level KO --out functional_exploration
```

The annotation table can leave groups unmapped. Unknown group IDs fail so a stale
or foreign mapping cannot silently disappear.

## What each function does

| Function | Output and assumption |
|---|---|
| `aggregate_terms` | Sum signal by distinct terms using equal splitting; retain unannotated signal and verify conservation |
| `relative_signal` | Fractions of measured signal, including the unannotated bucket when present |
| `diversity_summary` | Observed richness, Shannon/effective Shannon diversity and Simpson concentration at the supplied feature level |
| `bray_curtis` | Distances between annotated signal compositions; empty acquisitions excluded explicitly |
| `clr_complete` | Natural-log ratios on features positive in every included sample; no pseudocount or imputation |
| `permanova` | Explicit one-way group comparison with reproducible permutations and optional donor blocks |
| `explore_study` | Portable tables, AnnData, coverage/diversity plots and eligible CLR PCA using AlphaPeptTools |

The high-level workflow excludes the unannotated bucket from diversity, distances
and CLR PCA. Coverage always uses all measured signal. Thus differences in
functional profiles describe the annotated portion, whose size is reported beside
the plot. Diversity of KOs is not species diversity; measured intensity is not
biomass, enzyme activity or pathway flux.

Outputs include `term_signal.tsv`, `relative_signal.tsv`, `annotation_coverage.tsv`,
`diversity.tsv`, `bray_curtis.tsv`, `functional_profile.h5ad`, and
`Functional_QC.png/.pdf`. CLR/PCA tables and `PCA_clr.png/.pdf` appear when eligible;
`COMPLETE.json` records exclusions. PCA needs at least three annotated acquisitions
and two varying complete CLR features. Always inspect the retained feature core.

## Optional repeated-measure comparison

Metadata must label exactly the matrix acquisitions:

```text
sample	condition	donor
D1_before	before	D1
D1_after	after	D1
D2_before	before	D2
D2_after	after	D2
```

```bash
fasta-lake explore-study --matrix matrix.tsv --annotations group_terms.tsv \
  --level KO --out paired_exploration --metadata metadata.tsv \
  --group-column condition --block-column donor --permutations 999 --seed 23
```

Only requesting `--group-column` runs PERMANOVA. Within-donor permutations require
conditions that can actually move within donors. Donor-constant groups cannot be
tested with that blocking scheme. Labels must be exchangeable under your design;
this software does not establish that assumption. PERMANOVA can reflect dispersion
as well as location differences and does not establish a causal effect. Its p-value
is a single exploratory comparison, not a multiple-testing correction.

`PERMANOVA.json` and every permuted statistic are exported. NumPy is the default;
`--backend numba` uses the same seeded permutations and a tested numeric kernel.
Its measured gain was modest (about 1.2–1.4× after compilation in the tested
100/400-sample examples), so it is optional. See [performance](../concepts/performance.md).

## Small runnable example

```bash
python examples/functional/run_example.py --out functional_example
```

The example is explicitly synthetic and checks conservation, annotation coverage,
PCA and within-donor permutations. It is a software demonstration, not a biological
finding. The [function map](../reference/function-map.md) links to the implementation and tests.
