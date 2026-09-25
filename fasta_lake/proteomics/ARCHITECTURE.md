# FastaLake Proteomics Mode — Architecture

**Status:** scaffolded 2026-06-09 on branch `dev/proteomics-mode`. Not yet wired into CLI.
**Author:** P. Treit + Claude
**Motivation:** single-organism (human/host) proteomics, NOT metaproteomics.

## Why a separate mode

FastaLake was built for metaproteomics: *"out of a vast, unknown sequence space,
which proteins are in this sample?"* De novo evidence subsets an intractable pool.

Single-organism proteomics is nearly the inverse problem:
- The reference (human SwissProt ~20K) is already near-complete and small — the
  completeness gap FastaLake exists to close barely exists.
- What's missing is the **proteoform dimension** (coding variants, isoforms, IG/HLA,
  novel ORFs), which is individual-specific and combinatorial.
- Pooling **population** variants (e.g. gnomAD, 1.37M seqs) into a single-individual
  search is dead weight: ~none beyond the individual's own are in them → pure
  decoy competition → the Böcker "clogging" effect.

**Empirical proof (2026-06-09, this branch's motivating finding):** the naive
transplant (FL lake = SP + TrEMBL isoforms + his WGS + 1.37M gnomAD) gave an apparent
1.66× peptide / 3.66× gene "depth advantage" over canonical at q≤0.01 — but it is
concentrated at the FDR boundary (median FL-only peptide Q = 3.1e-3 vs shared 5.1e-6)
and the advantage shrinks monotonically with FDR (1.66× → 1.37× → and reverses to
0.85× at q≤0.001). NOTE: the strict-q probe is itself confounded by cross-DB q-value
calibration; the correct verdict tool is entrapment FDP (Point 4), not a sliding q.
Either way: a naive depth claim is not defensible. The defensible value is
**specificity** — variant/IG/HLA identities canonical search is blind to *by
construction*.

## Reframe (the design principle)

Swap the prior: **population/catalogue prior → personal-genome prior.**
- Metaproteomics: orthogonal evidence = de novo (no genome).
- Proteomics: primary orthogonal evidence = the individual's own WGS/RNA-seq
  (a far stronger, lower-noise prior than population panels or de novo alone).
  De novo demotes to a secondary discovery channel for what annotation can't
  predict (novel ORFs, editing, IG junctions), keeping its database-free,
  non-FDR-circular property.

## The six architecture changes (priority order)

1. **Personal genome-anchored extension — not population pooling.**
   Variant set from the individual's WGS/WES/RNA-seq. No matched genome → gate a
   population panel by de novo + cap/weight by allele frequency as a separate
   low-prior class. Variant DB drops from millions to thousands → FDR stays clean.
   _First rerun:_ rebuild the personal lake as SP + TrEMBL isoforms + his WGS only
   (drop the gnomAD population glob). Module: `personal_lake.py`.

2. **Cascade / LAYERED TIERS — don't co-search; always keep the canonical backbone.**
   This was the original design (Peter_Probant/sage_comparison March 2026) and is the
   fix to BOTH 2026-regression bugs: every tier is **canonical + an incremental layer**,
   so the canonical proteome is never lost and each layer's gain is attributable.
   - Tier 1: compact canonical reference → standard 1% FDR (clean depth/quant backbone).
   - Tier 2..N: small high-prior layers added ON TOP of canonical (personal WGS variants,
     smart-selected expressed isoforms, IG/HLA alleles), the discovery layer searched with
     **class-specific FDR** so it neither inflates nor borrows strength from canonical.

   **Existing tier ladder (reuse, do NOT reinvent):**
   `benchmark_fastas/` + `50_build_benchmark_fastas.py`:
   | tier | anchored on canonical? | content |
   |---|---|---|
   | 01_canonical (20,420) | — (is canonical) | SwissProt canonical |
   | 02_canonical_plus_wgs (38,081) | YES | + personal WGS variants |
   | 06_canonical_plus_fl_sp (~20.5K) | YES | + FL-discovered SP |
   | 07_canonical_wgs_fl_sp (~38.2K) | YES | + WGS + FL-SP |
   | 09_smart_isoforms / 10_canonical_plus_smart | YES (10) | smart-isoform selection (Point 6) |
   | 11_canonical_wgs_healed (38,081) | YES | + WGS + header-healed gene rollup |
   | 12_plasma_fastalake (~45–62K) | YES | full layered plasma lake |
   | 03_fastalake_sp_only (328–532) / 04_fastalake_full | **NO — TRAP** | bare de-novo extraction, no backbone |

   **2026 REGRESSION:** `personal_proteome_2026/4_diann/fl/cohort_FL.fasta` is the
   `03/04` bare-extraction pattern (only 314 SwissProt survivors) → it lost the canonical
   backbone AND its 92%-accession-GN= TrEMBL entries inflated the gene count to a 94%-fake
   3,759 (real genes only 242, vs canonical 889). The fair comparator is `06/07/10/11/12`.
   Module: `cascade.py` (tier builder, canonical-anchored by construction; ports the
   50_build_benchmark_fastas.py layers to the current 10-probant DIA-NN 2.2.0 run).

3. **Point de novo at the unexplained, not at subsetting.**
   De novo on spectra canonical failed to explain; keep de novo peptides absent from
   canonical even after I/L fold; confirm against personal genome / 3-frame
   translation / variant catalogues → targeted variant detector. Module: `residual_denovo.py`.

4. **Entrapment as a standard output, not a special experiment.**
   Empirical FDP via foreign-proteome spike (P. furiosus / shuffled) baked into every
   run's report. Reuse `fasta_lake/entrapment.py` (estimate_fdp / validate_entrapment_fdr,
   Wen-Keich-Noble 2025) + the proven `Peter_Probant/sage_comparison/2*_entrap*` scripts.
   Module: `entrapment_report.py` (thin wrapper over the existing entrapment.py).

5. **Emit a personalized predicted spectral library (DIA), not just a FASTA.**
   Curated canonical + confident personal variants, predicted RT/IM/fragments.
   Cleaner/faster than DIA-NN library-free against a huge FASTA. Module: `speclib.py`.

6. **Fix the inference/grouping layer — the honest win.**
   TrEMBL's 93%-missing GN= inflates protein groups. Provide complete, gene-resolved,
   isoform-aware sequences so grouping/gene-rollup is correct. Reuse
   `fasta_lake/headers/` + `inference/` + Probant `54/57_*heal*`. Module: `grouping.py`.

## Proven basis to fold in (do NOT reinvent)
`Peter_Probant/sage_comparison/`:
- canonical-vs-FL: `01/02/03`, `36-39_diann*`
- entrapment FDP: `20/21/22/23`  ← Point 4 already implemented here
- header-heal + isoform: `54/57`  ← Point 6
- benchmark FASTA tiers: `50_build_benchmark_fastas.py`

## Package layout (NEW)
```
fasta_lake/proteomics/
├── __init__.py
├── ARCHITECTURE.md          (this file)
├── personal_lake.py         (Point 1)
├── cascade.py               (Point 2 — class-specific FDR)
├── residual_denovo.py       (Point 3)
├── entrapment_report.py     (Point 4 — wraps entrapment.py)
├── speclib.py               (Point 5)
└── grouping.py              (Point 6)
```

## First validation rerun (this branch)
Rebuild personal-anchored lake (Point 1) → DIA-NN FL arm → entrapment FDP (Point 4).
Hypothesis: dropping the 1.37M population gnomAD removes the boundary-clogging; FL's
1% IDs become more confident and entrapment FDP stays ≤1%. Compare against the
existing canonical arm (job 5005509) and the gnomAD-polluted FL arm (job 5005510).
