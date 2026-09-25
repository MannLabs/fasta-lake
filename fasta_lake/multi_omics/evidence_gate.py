"""
FL_MO multi-omics evidence gate (Methods §M5.5).

This is the committed, runnable implementation of the FastaLake multi-omics
(FL_MO) arm reported in the manuscript. It restores, as a first-class library
function, the per-sample OR-gate / AND-gate logic that produced the verified
FL_MO numbers (defect-ledger D3: razor per-sample median 13,776 protein groups,
+11-13% across inference strategies, OR-gate median 13,170, AND-gate 3,544).

The original verified code lived in project-specific scripts:
  - projects/groundup_2026-05-12/MP/scripts/17_FL_MO_persample_augment.py
    (OR-gate: builds the per-sample evidence_lookup.tsv from raw inputs)
  - manuscript/v2.2_analyses/scripts/A2_AND_gate_build_sample_fasta.py
    (AND-gate: subsets an existing evidence_lookup.tsv)
  - manuscript/v2.2_analyses/scripts/c86_sweep_build_fasta.py
    (threshold sweep; generalises both gates)

This module re-expresses that logic against the frozen FastaLake hashing
contract (`fasta_lake.hashing`) with no new dependencies, so the FL_MO numbers
can be regenerated from committed code.

The gate logic (Methods §M5.5)
------------------------------
A FL_MO run augments per-sample extraction (Stage 3). For each cluster
representative r, include r if ANY of:

  (a) a sample de novo peptide matches r exactly (the normal de-novo evidence
      that the Stage-3 extractor already produces); OR
  (b) the sample's MEGAHIT-assembled protein hashing (SHA-256, frozen
      normalize contract) to r's cluster has metaG_tpm >= the per-sample
      top-1% threshold (computed on non-zero metaG_tpm; ties -> smallest tied
      value so all ties are included); OR
  (c) the same hash-matched assembly protein has metaT_tpm > 1.0.

OR-gate  = (a) OR (b) OR (c)            -- default for FL_MO
AND-gate = (a) AND ((b) OR (c))        -- the low-FDR option

Cross-sample identity is by SHA-256 protein hash via the FL_MO hash_to_rep
table (assembly protein sequence hash -> cluster representative id).

FASTA headers are NEVER edited: rescued representative sequences are copied
verbatim from the cluster97 representative FASTA (CLAUDE.md contract).
"""

from __future__ import annotations

from dataclasses import dataclass, field

from fasta_lake.hashing import sequence_hash

__all__ = [
    "Gate",
    "EvidenceRow",
    "top1_threshold",
    "passes_tpm",
    "build_evidence_lookup",
    "gate_members",
    "EVIDENCE_LOOKUP_COLUMNS",
]

# Default thresholds from Methods §M5.5.
DEFAULT_METAG_TOP_PCT = 1.0  # top 1% of non-zero metaG_tpm
DEFAULT_METAT_ABS = 1.0  # metaT_tpm strictly greater than 1.0

EVIDENCE_LOOKUP_COLUMNS = [
    "cluster_rep_id",
    "source",
    "has_denovo",
    "has_metaG_top1",
    "has_metaT_1",
    "metaG_tpm",
    "metaT_tpm",
    "orig_k141",
]


class Gate:
    """Gate variants."""

    OR = "or"
    AND = "and"


@dataclass
class EvidenceRow:
    """One cluster representative's per-sample evidence attribution."""

    cluster_rep_id: str
    source: str  # "DENOVO" or "TPM_RESCUE"
    has_denovo: bool
    has_metaG_top1: bool
    has_metaT_1: bool
    metaG_tpm: float = 0.0
    metaT_tpm: float = 0.0
    orig_k141: list[str] = field(default_factory=list)

    @property
    def has_tpm(self) -> bool:
        """Report whether either historical DNA or RNA support flag is set."""
        return self.has_metaG_top1 or self.has_metaT_1

    def in_gate(self, gate: str) -> bool:
        """Whether this representative is admitted under the given gate."""
        if gate == Gate.OR:
            # (a) de novo, OR (b) metaG-top1, OR (c) metaT>1
            return self.has_denovo or self.has_tpm
        if gate == Gate.AND:
            # (a) de novo AND ((b) OR (c))
            return self.has_denovo and self.has_tpm
        raise ValueError(f"unknown gate {gate!r}; use Gate.OR or Gate.AND")

    def to_tsv(self) -> str:
        """Serialize this legacy gate row, formatting TPM to four decimal places."""
        return "\t".join(
            [
                self.cluster_rep_id,
                self.source,
                "True" if self.has_denovo else "False",
                "True" if self.has_metaG_top1 else "False",
                "True" if self.has_metaT_1 else "False",
                f"{self.metaG_tpm:.4f}" if self.metaG_tpm else "0",
                f"{self.metaT_tpm:.4f}" if self.metaT_tpm else "0",
                ";".join(self.orig_k141),
            ]
        )


def top1_threshold(metag_tpm_values, top_pct: float = DEFAULT_METAG_TOP_PCT) -> float:
    """Per-sample metaG cutoff: smallest value v* such that values >= v* are the
    top `top_pct` percent of the NON-ZERO distribution.

    Ties resolve to the smallest tied value, so every protein tied at the
    boundary is admitted (Methods §M5.5). Returns 0.0 if nothing is non-zero
    (i.e. no metaG rescue possible for the sample).

    Matches the verified contract in 17_FL_MO_persample_augment.py:
        nonzero = sorted(values_desc); top1 = nonzero[ceil(n*pct/100) - 1]
    re-expressed in ascending order with explicit tie-to-smallest handling.
    """
    nonzero = sorted(v for v in metag_tpm_values if v > 0.0)
    if not nonzero:
        return 0.0
    n = len(nonzero)
    # number of entries in the top band (at least 1), matching int(n*0.01) index
    # in the original descending formulation: top1 = desc[max(0, int(n*pct)-1)].
    k = max(1, int(n * top_pct / 100.0))
    # the k-th largest value, i.e. ascending index n - k
    cutoff = nonzero[n - k]
    # extend cutoff down to the smallest value tied with it (include all ties)
    i = n - k
    while i > 0 and nonzero[i - 1] == cutoff:
        i -= 1
    return nonzero[i]


def passes_tpm(
    metag_tpm: float, metat_tpm: float, metag_cutoff: float, metat_abs: float = DEFAULT_METAT_ABS
) -> tuple[bool, bool]:
    """Return (passes_metaG_top1, passes_metaT_1) for a single assembly protein."""
    pass_mg = metag_cutoff > 0.0 and metag_tpm >= metag_cutoff
    pass_mt = metat_tpm > metat_abs
    return pass_mg, pass_mt


def build_evidence_lookup(
    denovo_rep_ids,
    tpm_rows,
    assembly_seqs,
    hash_to_rep,
    metag_top_pct: float = DEFAULT_METAG_TOP_PCT,
    metat_abs: float = DEFAULT_METAT_ABS,
):
    """Build the per-sample evidence attribution table (the FL_MO core).

    Parameters
    ----------
    denovo_rep_ids : iterable of str
        Cluster representative ids the Stage-3 de-novo extractor selected for
        this sample (evidence (a)).
    tpm_rows : iterable of (assembly_protein_id, metaG_tpm, metaT_tpm)
        The sample's per-assembly TPM table (e.g. <sample>_TPM.tsv).
    assembly_seqs : mapping assembly_protein_id -> amino-acid sequence (str)
        The sample's MEGAHIT predicted-genes FASTA, indexed by id.
    hash_to_rep : mapping sha256_hex -> cluster_rep_id
        FL_MO hash_to_rep table (cross-sample identity bridge).
    metag_top_pct, metat_abs : float
        Gate thresholds (defaults are the published §M5.5 values).

    Returns
    -------
    (rows, stats) where rows is a list[EvidenceRow] and stats is a dict.

    A representative may be admitted by both de novo and TPM evidence; in that
    case it appears once with source=DENOVO and the TPM flags set.
    """
    denovo_rep_ids = set(denovo_rep_ids)

    # 1. per-sample metaG top1% cutoff on non-zero metaG_tpm
    cutoff = top1_threshold((mg for _, mg, _ in tpm_rows), metag_top_pct)

    # 2. TPM-passing assembly proteins -> hash -> rep, aggregating per rep
    rescued: dict[str, dict] = {}
    n_tpm_rows = 0
    n_passing = 0
    n_missing_asm = 0
    n_missing_lake = 0
    for prot_id, mg, mt in tpm_rows:
        n_tpm_rows += 1
        pass_mg, pass_mt = passes_tpm(mg, mt, cutoff, metat_abs)
        if not (pass_mg or pass_mt):
            continue
        n_passing += 1
        seq = assembly_seqs.get(prot_id)
        if seq is None:
            n_missing_asm += 1
            continue
        rep = hash_to_rep.get(sequence_hash(seq))
        if rep is None:
            n_missing_lake += 1
            continue
        cur = rescued.get(rep)
        if cur is None:
            cur = {"metaG": False, "metaT": False, "metaG_tpm": 0.0, "metaT_tpm": 0.0, "origs": []}
            rescued[rep] = cur
        cur["metaG"] = cur["metaG"] or pass_mg
        cur["metaT"] = cur["metaT"] or pass_mt
        cur["metaG_tpm"] = max(cur["metaG_tpm"], mg)
        cur["metaT_tpm"] = max(cur["metaT_tpm"], mt)
        cur["origs"].append(prot_id)

    # 3. emit rows. de-novo reps first (source=DENOVO), carrying TPM flags if
    #    the same rep was also TPM-rescued; then TPM-only rescued reps. Both
    #    blocks are emitted in sorted id order: the inputs are a set and a dict
    #    filled in file order, and the lookup TSV must be byte-identical across
    #    runs and PYTHONHASHSEED values.
    rows: list[EvidenceRow] = []
    for rep in sorted(denovo_rep_ids):
        extra = rescued.get(rep)
        rows.append(
            EvidenceRow(
                cluster_rep_id=rep,
                source="DENOVO",
                has_denovo=True,
                has_metaG_top1=bool(extra and extra["metaG"]),
                has_metaT_1=bool(extra and extra["metaT"]),
                metaG_tpm=extra["metaG_tpm"] if extra else 0.0,
                metaT_tpm=extra["metaT_tpm"] if extra else 0.0,
                orig_k141=extra["origs"] if extra else [],
            )
        )
    for rep, info in sorted(rescued.items()):
        if rep in denovo_rep_ids:
            continue
        rows.append(
            EvidenceRow(
                cluster_rep_id=rep,
                source="TPM_RESCUE",
                has_denovo=False,
                has_metaG_top1=info["metaG"],
                has_metaT_1=info["metaT"],
                metaG_tpm=info["metaG_tpm"],
                metaT_tpm=info["metaT_tpm"],
                orig_k141=info["origs"],
            )
        )

    stats = {
        "metaG_top1_cutoff": cutoff,
        "tpm_rows": n_tpm_rows,
        "tpm_passing": n_passing,
        "missing_assembly": n_missing_asm,
        "missing_lake": n_missing_lake,
        "denovo_reps": len(denovo_rep_ids),
        "tpm_rescue_reps_new": sum(1 for r in rescued if r not in denovo_rep_ids),
        "or_gate_members": sum(1 for r in rows if r.in_gate(Gate.OR)),
        "and_gate_members": sum(1 for r in rows if r.in_gate(Gate.AND)),
    }
    return rows, stats


def gate_members(rows, gate: str = Gate.OR):
    """Return the set of cluster_rep_ids admitted under `gate`."""
    return {r.cluster_rep_id for r in rows if r.in_gate(gate)}
