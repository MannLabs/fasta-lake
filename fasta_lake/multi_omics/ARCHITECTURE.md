# FastaLake Multi-Omics Module — Architecture

**Status:** IMPLEMENTED (2026-06-17). The FL_MO evidence gate is now committed
runnable code in this package:
- `evidence_gate.py` — OR/AND gate logic + per-sample top-1% metaG threshold
  (pure functions, frozen `fasta_lake.hashing` SHA-256 contract).
- `runner.py` — file-level Stage-3 runner (`run_fl_mo_sample`): reads the
  sample's de-novo FASTA + TPM TSV + MEGAHIT FASTA + hash_to_rep + cluster rep
  FASTA, writes the augmented FASTA and `evidence_lookup.tsv`.
- `reproduce.py` — regenerates the published razor per-sample median (13,776)
  from the committed verification TSVs.
- CLI: `fasta-lake fl-mo gate ...` and `fasta-lake fl-mo reproduce`
  (also `python -m fasta_lake.multi_omics`). Opt-in: default de-novo-only runs
  are unaffected.
- Tests: `tests/test_multi_omics.py`.

The logic was lifted from the verified project scripts
`projects/groundup_2026-05-12/MP/scripts/17_FL_MO_persample_augment.py` (OR-gate),
`manuscript/v2.2_analyses/scripts/A2_AND_gate_build_sample_fasta.py` (AND-gate),
and `.../c86_sweep_build_fasta.py` (threshold sweep).

**Date:** 2026-05-14 (spec) / 2026-06-17 (implementation)
**Author:** P. Treit + Claude (planning session)

## Why

Cheng et al. (MIM, Cell 2025) anchor metaproteomic database curation on metagenomic abundance. FastaLake currently anchors on de novo peptide evidence. The orthogonality of these information sources is empirically testable: a protein with no de novo evidence but high metaG TPM is a candidate that de novo missed for measurement reasons, not biological-absence reasons.

The OLD `FINAL_PUBLICATION_HEIDELBERG` pipeline (Nov 2025) already proved this works on 280/292 PREDICT samples via a pre-clustering OR-gate (denovo OR metaG_TPM≥top1% OR metaT_TPM>1.0). The current `groundup_2026-05-12` pipeline lost this feature. The multi-omics module restores it as a first-class FastaLake feature with SHA-256-keyed cross-sample identity preservation.

## Architecture

```
┌────────────────────────────────────────────────────────────────────┐
│  Stage 0 — Lake builder (existing)                                 │
│  Sources: UHGG + GMGC + GMSC + UniProt + cohort-MetaG              │
│  Output:  273M-record SHA-256-keyed lake                           │
└────────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌────────────────────────────────────────────────────────────────────┐
│  Stage 1 — Evidence filter (existing)                              │
│  Cohort-pooled de novo peptides → 29.3M-record evidence lake       │
└────────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌────────────────────────────────────────────────────────────────────┐
│  Stage 1.5 — Cluster97 (existing)                                  │
│  MMseqs2 easy-cluster → 7.8M cluster representatives               │
│  Sidecar: cluster97_cluster.tsv (rep → members)                    │
└────────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌────────────────────────────────────────────────────────────────────┐
│  Stage 2 — Per-sample extraction (existing + NEW multi-omics flag) │
│  Inputs (existing):                                                │
│    - cluster97 rep FASTA                                           │
│    - sample's de novo peptide CSV                                  │
│  Inputs (NEW):                                                     │
│    - sample's per-MEGAHIT predictedGenesFromAssembly.<S>.fasta     │
│    - sample's TPM TSV (metaG_tpm, metaT_tpm per k141_X_Y)          │
│    - flag --multi-omics-tpm <DIR>                                  │
│  Logic:                                                            │
│    A.  Run de novo extraction as today → set DENOVO_PROTEINS       │
│    B.  Build (k141_X_Y → SHA-256) map from sample's own FASTA     │
│    C.  Bridge: SHA-256 → cluster97 rep (via lake's joint_headers   │
│        + cluster mapping)                                          │
│    D.  Filter sample's TPM: pass = (metaG_tpm ≥ top1) OR           │
│        (metaT_tpm > 1.0)                                           │
│    E.  Translate passing k141_X_Y → cluster97 reps                 │
│    F.  Output FASTA = DENOVO_PROTEINS ∪ TPM_RESCUED_REPS           │
│  Output sidecar: evidence_lookup.tsv                               │
│    Columns: cluster_rep_id, has_denovo, denovo_max_score,          │
│             denovo_n_peptides, has_metaG_top1, metaG_tpm,          │
│             has_metaT_1, metaT_tpm, evidence_sources, included     │
└────────────────────────────────────────────────────────────────────┘
                              │
                              ▼
                        [Stage 3-7 unchanged]
```

## The lookup-table is the core artifact

Per-sample `evidence_lookup.tsv` is **the contract** between Stage 2 and downstream analysis. It is the answer to "where did this protein come from?" and enables:

1. Per-source attribution of detected protein groups (Fig 3-style breakdown)
2. Per-source entrap FDP audit ("are metaG-rescued proteins disproportionately entrapped?")
3. Cross-sample evidence-source agreement statistics
4. Manuscript Methods §M3 reproducibility — every protein has a documented provenance

## Reproducibility contract

- SHA-256 sequence hashing uses the FastaLake `lake_builder::normalize_sequence` contract verbatim
- top1 threshold is computed per-sample from the sample's own TPM distribution (excluding zeros)
- The `included` boolean in evidence_lookup.tsv is a deterministic function of `(denovo_evidence, metaG_tpm vs top1, metaT_tpm)`
- Test fixture: `tests/multi_omics/test_or_gate.py` covers
   - de novo-only protein (no TPM): included via denovo
   - high-metaG only protein (no de novo): included via metaG rescue
   - low-TPM, low-de-novo: excluded
   - both: included with `evidence_sources=denovo+metaG`
   - top1 threshold sensitivity (top0.1% vs top1% vs top5%)

## Source files identified for implementation

| Purpose | Path |
|---|---|
| Per-sample MEGAHIT assemblies | `Peter_Microbpred/FINAL_PUBLICATION_HEIDELBERG/0_DATA/<run_id>/predictedGenesFromAssembly.<sample>.fasta` |
| Per-sample TPM | `Peter_Microbpred/FINAL_PUBLICATION_HEIDELBERG/0_DATA/<run_id>/<sample>_TPM.tsv` |
| Augmented TPM (with cluster reps, OLD lake) — reference only | `Peter_Microbpred/Clean/cleaned_databases/nMuPr_predictedGenesfromAssembly_Study/augmented_tpm_files_with_reps/<sample>_TPM_augmented.tsv` (280 files) |
| Master pool (OLD lake; for nMuPr_id → seq lookup if needed) | `Peter_Microbpred/Multi-Omics/fasta4_final/nMuPr_lake.fasta` |
| Hash → header map (MD5, OLD lake) | `Peter_Microbpred/FINAL_PUBLICATION_HEIDELBERG/all_study_protein_sequences_MuPr_Lake/hash_to_headers.tsv` |
| Assembly_protein_map.tsv | `Peter_Microbpred/Clean/cleaned_databases/nMuPr_predictedGenesfromAssembly_Study/assembly_protein_map.tsv` |
| Current MP cluster97 representatives | `groundup_2026-05-12/MP/1_5_cluster97/FL/cluster97_rep_seq.fasta` |
| Current MP cluster mapping | `groundup_2026-05-12/MP/1_5_cluster97/FL/cluster97_cluster.tsv` |
| Current MP joint headers | `groundup_2026-05-12/MP/0_lakes/metag_joint_headers.tsv` |

## NOT the way

We are NOT reusing the OLD augmented_tpm_files_with_reps. The OLD files reference OLD lake cluster reps (`GMGC10.298_358_274.UNKNOWN`-style); the NEW lake's cluster reps use different IDs. Bridging through the OLD lake adds a dependency and a potential semantic drift.

We ARE generating fresh augmented TPM files from raw inputs (per-sample MEGAHIT FASTA + per-sample TPM TSV) using the FastaLake SHA-256 hashing contract, producing a `<sample>_evidence_lookup.tsv` that is bit-deterministic given fixed inputs.

## Coverage

- 280/292 samples have metaG_tpm (96%)
- 206/292 samples have metaT_tpm (71%)
- 12 samples have NEITHER → fall back to denovo-only (per-sample fallback path documented)

## Package layout (NEW)

```
fasta_lake/multi_omics/
├── __init__.py
├── ARCHITECTURE.md             (this file)
├── tpm_reader.py               (load per-sample TPM TSVs, compute top1 threshold)
├── evidence_or_gate.py         (the OR-gate logic)
├── lookup_generator.py         (build evidence_lookup.tsv)
└── cli.py                      (CLI wrapper)

rust/fasta_extractor_v2/src/
├── multi_omics.rs              (Rust port; pairs with --multi-omics-tpm flag)

tests/multi_omics/
├── test_or_gate.py
├── test_evidence_lookup.py
└── fixtures/
    └── tiny_sample/            (1 sample, 100 proteins, hand-crafted TPM)
```

## Manuscript hooks

- Results §"Three-arm benchmark": add FL+MO 4th-arm paragraph (stub already drafted at `manuscript/sources/FL_MO_paragraph_stub.md`)
- Methods §M3 Stage-2: document `--multi-omics-tpm` flag
- Methods §M11 (new sub-section): "Multi-omics evidence integration"
- Discussion: explicit MIM (Cheng 2025) head-to-head — "MIM anchors on abundance, FastaLake-MO triangulates abundance + expression + de novo"
- Extended Data: per-sample evidence-source breakdown plot, FL_MO vs FL entrap FDP delta

## Open questions for v6

1. Should the OR-gate become an AND-gate at scale (high-precision mode)? — likely a config option
2. Should we sweep multiple top-1% thresholds (top0.1%, top1%, top5%) and report sensitivity? — yes, in supplementary
3. Should we extend to cohort-level (rather than per-sample) TPM normalization? — methodological choice; can be a future paper
