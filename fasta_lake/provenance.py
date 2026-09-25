"""
Per-protein provenance (v5.1).

Writes a JSONL sidecar alongside each per-sample FASTA. One line per selected
protein, with enough metadata to answer:

- Which de novo peptides admitted this protein?
- What are their scores?
- Was this protein k-mer rescued or exact-matched?
- How many peptides support it?
- Which inference strategy selected it?

Output format: `{sample_id}_{strategy}_provenance.jsonl`

Each line:
    {"protein_id": "sp|P04637|P53_HUMAN",
     "n_peptides": 3,
     "peptides": ["MEEPQSDPSVEPPLSQ", "ETPPLGE...", ...],
     "scores":   [0.97, 0.88, 0.85],
     "mean_score": 0.90,
     "best_score": 0.97,
     "admission": "exact",             # "exact" | "kmer_rescue"
     "strategy": "species_budget",
     "contaminant": false}

Compact, grep-able, line-per-protein.
"""

from __future__ import annotations

import json
import statistics
from pathlib import Path

from fasta_lake.contaminants import is_contaminant_id


def write_provenance(
    out_path: Path,
    selected_proteins: set[str],
    prot2pep: dict[str, set[str]],
    peptide_scores: dict[str, float],
    strategy: str,
    kmer_rescued_proteins: set[str] | None = None,
) -> int:
    """
    Write per-protein provenance sidecar.

    Parameters
    ----------
    out_path : Path
        Output JSONL file.
    selected_proteins : set[str]
        Protein IDs in the final inference output.
    prot2pep : dict[str, set[str]]
        Protein -> peptides mapping (I/L-normalized peptides).
    peptide_scores : dict[str, float]
        Peptide -> score (I/L-normalized keys).
    strategy : str
        Inference strategy name.
    kmer_rescued_proteins : set[str] or None
        Proteins added by k-mer rescue. If None, assumes none.

    Returns
    -------
    int
        Number of provenance records written.
    """
    if kmer_rescued_proteins is None:
        kmer_rescued_proteins = set()

    n_written = 0
    with open(out_path, "w") as f:
        for protein_id in sorted(selected_proteins):
            peps = sorted(prot2pep.get(protein_id, set()))
            if not peps:
                # Protein selected but no peptides in map — defensive
                record = {
                    "protein_id": protein_id,
                    "n_peptides": 0,
                    "peptides": [],
                    "scores": [],
                    "mean_score": None,
                    "best_score": None,
                    "admission": (
                        "kmer_rescue" if protein_id in kmer_rescued_proteins else "exact"
                    ),
                    "strategy": strategy,
                    "contaminant": is_contaminant_id(protein_id),
                }
            else:
                scores = [peptide_scores.get(p, 0.5) for p in peps]
                record = {
                    "protein_id": protein_id,
                    "n_peptides": len(peps),
                    # Cap peptide list to avoid huge JSONL — include all scores though
                    "peptides": peps[:20],
                    "peptides_truncated": len(peps) > 20,
                    "scores": scores[:20],
                    "mean_score": round(statistics.mean(scores), 4),
                    "best_score": round(max(scores), 4),
                    "admission": (
                        "kmer_rescue" if protein_id in kmer_rescued_proteins else "exact"
                    ),
                    "strategy": strategy,
                    "contaminant": is_contaminant_id(protein_id),
                }
            f.write(json.dumps(record) + "\n")
            n_written += 1
    return n_written


def summarise_provenance(jsonl_path: Path) -> dict:
    """Aggregate a provenance JSONL into summary counts.

    Useful for reports.
    """
    n_total = 0
    n_kmer_rescued = 0
    n_contaminant = 0
    peptide_counts = []
    scores = []
    with open(jsonl_path) as f:
        for line in f:
            r = json.loads(line)
            n_total += 1
            if r.get("admission") == "kmer_rescue":
                n_kmer_rescued += 1
            if r.get("contaminant"):
                n_contaminant += 1
            peptide_counts.append(r.get("n_peptides", 0))
            if r.get("best_score") is not None:
                scores.append(r["best_score"])

    return {
        "n_proteins": n_total,
        "n_kmer_rescued": n_kmer_rescued,
        "n_contaminant": n_contaminant,
        "median_peptides_per_protein": (statistics.median(peptide_counts) if peptide_counts else 0),
        "mean_best_score": (round(statistics.mean(scores), 4) if scores else None),
        "proteins_with_1_peptide": sum(1 for c in peptide_counts if c == 1),
        "proteins_with_2plus_peptides": sum(1 for c in peptide_counts if c >= 2),
    }
