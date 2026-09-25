"""
FastaLake multi-omics (FL_MO) evidence gate — Methods §M5.5.

Opt-in arm that augments per-sample de-novo extraction with metagenomic /
metatranscriptomic abundance evidence. Default (de-novo-only) runs are not
affected: this module only runs when explicitly invoked.

Public API
----------
- evidence_gate : the gate logic (OR / AND), pure functions + EvidenceRow
- runner        : file-level Stage-3 runner (run_fl_mo_sample)
- reproduce     : regenerate the published 13,776 razor median from committed TSVs

See ARCHITECTURE.md for the original design and source-script provenance.
"""

from fasta_lake.multi_omics.evidence_gate import (
    EVIDENCE_LOOKUP_COLUMNS,
    EvidenceRow,
    Gate,
    build_evidence_lookup,
    gate_members,
    passes_tpm,
    top1_threshold,
)
from fasta_lake.multi_omics.runner import run_fl_mo_sample

__all__ = [
    "Gate",
    "EvidenceRow",
    "EVIDENCE_LOOKUP_COLUMNS",
    "top1_threshold",
    "passes_tpm",
    "build_evidence_lookup",
    "gate_members",
    "run_fl_mo_sample",
]
