"""
ICS-guided variant *proposal* for FastaLake. Research module, proposal-only.

STATUS: this module proposes sequence variants; it does not validate them.
It is not part of the published FastaLake method, is not reachable from the
``fasta-lake`` command line, and its output must not be used as a per-sample
search database. A proposed variant is a hypothesis about where a database
sequence may be wrong; confirming it requires a re-search and an ICS
comparison that this module does not perform.

What it does: after a SAGE search, ``identify_refinement_candidates`` selects
FDR-passing PSMs with low ion coverage (ICS), and ``generate_variant_fasta``
writes, for each candidate protein, the single amino-acid substitution with
the best BLOSUM95 (or de novo confusion) score at the peptide position. The
original sequence is replaced in the output and the accession is suffixed
``__VAR`` so the file cannot be mistaken for the input.

What it does not do: compute ICS for the variant, compare it with the
original, or accept only improving variants. ``run_refinement`` therefore
raises ``NotImplementedError`` until that validation exists.

Reference: Miura et al. (2025) for the ICS definition.
"""

from __future__ import annotations

import logging
from collections import defaultdict
from pathlib import Path

from fasta_lake.blosum import generate_variants
from fasta_lake.helpers import load_fasta_proteins
from fasta_lake.ics import add_ics, clean_sequence

logger = logging.getLogger(__name__)


def identify_refinement_candidates(
    sage_results_path: str | Path,
    ics_threshold: float = 0.4,
    fdr_threshold: float = 0.01,
    max_candidates: int = 500,
) -> list[dict]:
    """
    Find low-ICS PSMs that are candidates for variant refinement.

    Parameters
    ----------
    sage_results_path : path
        SAGE results.sage.tsv
    ics_threshold : float
        PSMs below this ICS are candidates. Default 0.4.
    fdr_threshold : float
        Only consider PSMs passing FDR. Default 0.01.
    max_candidates : int
        Limit candidates to avoid combinatorial explosion. Default 500.

    Returns
    -------
    list of dict
        Each with: peptide, protein_id, position, ics, hyperscore
    """
    import pandas as pd

    df = pd.read_csv(sage_results_path, sep="\t")
    df = add_ics(df, engine="sage")

    # Filter: passing FDR + low ICS
    fdr_col = "spectrum_q" if "spectrum_q" in df.columns else "spectrum_fdr"
    mask = (df[fdr_col] <= fdr_threshold) & (df["ics"] < ics_threshold)
    candidates = df[mask].copy()

    logger.info(
        "Refinement candidates: %d/%d PSMs (ICS < %.2f, FDR <= %.2f)",
        len(candidates),
        len(df),
        ics_threshold,
        fdr_threshold,
    )

    if candidates.empty:
        return []

    # Deduplicate by peptide (keep best hyperscore)
    candidates = candidates.sort_values("hyperscore", ascending=False)
    candidates = candidates.drop_duplicates(subset=["peptide"], keep="first")

    # Limit
    candidates = candidates.head(max_candidates)

    result = []
    for _, row in candidates.iterrows():
        result.append(
            {
                "peptide": clean_sequence(row["peptide"]),
                "peptide_raw": row["peptide"],
                "protein_id": row["proteins"].split(";")[0],
                "ics": row["ics"],
                "hyperscore": row["hyperscore"],
            }
        )

    logger.info("Selected %d unique peptides for refinement", len(result))
    return result


def generate_variant_fasta(
    fasta_path: str | Path,
    candidates: list[dict],
    matrix: str = "blosum95",
    min_score: int = -1,
    max_variants_per_position: int = 3,
    output_path: str | Path | None = None,
) -> dict[str, str]:
    """
    Write a proposal FASTA: one best-scoring substitution per candidate protein.

    For each candidate peptide:
    1. Locate it in the parent protein sequence (exact, then I/L-folded)
    2. Generate BLOSUM/denovo-guided single-AA variants at each peptide position
    3. Keep the single variant with the highest substitution score

    Every candidate protein that maps is modified, unconditionally. No ICS
    or spectral evidence is consulted for the variant; the substitution score
    is the only criterion. The output has the same number of records as the
    input, modified records carry the accession suffix ``__VAR``, and it is
    a proposal for a validating re-search, not a search database.

    Parameters
    ----------
    fasta_path : path
        Per-sample FASTA to refine.
    candidates : list of dict
        From ``identify_refinement_candidates()``.
    matrix : str
        ``'blosum95'`` or ``'denovo'``.
    output_path : path, optional
        Write variant FASTA. If None, returns dict only.

    Returns
    -------
    dict
        Statistics: n_variants, n_proteins_modified, etc.
    """
    proteins = load_fasta_proteins(fasta_path)

    # Map peptides to their protein and position
    variants_by_protein: dict[str, list[tuple[str, str, int]]] = defaultdict(list)
    n_mapped = 0

    for cand in candidates:
        pep = cand["peptide"]
        prot_id = cand["protein_id"]

        if prot_id not in proteins:
            continue

        seq = proteins[prot_id]
        pos = seq.find(pep)
        if pos < 0:
            # Try I/L normalized
            seq_norm = seq.replace("I", "L")
            pep_norm = pep.replace("I", "L")
            pos = seq_norm.find(pep_norm)

        if pos < 0:
            continue

        n_mapped += 1

        # Generate variants for each position in the peptide
        best_variant = None
        best_score = -999

        for pep_pos in range(len(pep)):
            protein_pos = pos + pep_pos
            for var_seq, notation, blosum_score in generate_variants(
                seq,
                protein_pos,
                min_score=min_score,
                matrix=matrix,
                max_per_position=max_variants_per_position,
            ):
                if blosum_score > best_score:
                    best_score = blosum_score
                    best_variant = (var_seq, notation, blosum_score)

        if best_variant:
            variants_by_protein[prot_id].append(best_variant)

    # Apply best variant per protein (highest BLOSUM score)
    n_modified = 0
    variant_proteins = dict(proteins)  # copy

    for prot_id, var_list in variants_by_protein.items():
        best = max(var_list, key=lambda x: x[2])
        variant_proteins[prot_id] = best[0]
        n_modified += 1

    logger.info(
        "Variants: %d candidates mapped, %d proteins modified",
        n_mapped,
        n_modified,
    )

    # Write output
    if output_path:
        output_path = Path(output_path)
        with open(output_path, "w") as f:
            for prot_id in sorted(variant_proteins):
                seq = variant_proteins[prot_id]
                # Mark modified proteins
                suffix = "__VAR" if prot_id in variants_by_protein else ""
                f.write(f">{prot_id}{suffix}\n")
                for i in range(0, len(seq), 60):
                    f.write(seq[i : i + 60] + "\n")

        logger.info("Written variant FASTA: %s (%d proteins)", output_path, len(variant_proteins))

    return {
        "candidates": len(candidates),
        "mapped": n_mapped,
        "proteins_modified": n_modified,
        "total_proteins": len(variant_proteins),
        "matrix": matrix,
    }


def run_refinement(
    sage_results: str | Path,
    sample_fasta: str | Path,
    output_fasta: str | Path,
    ics_threshold: float = 0.4,
    fdr_threshold: float = 0.01,
    matrix: str = "blosum95",
    max_candidates: int = 500,
) -> dict:
    """
    Not implemented: a validated refinement pass does not exist yet.

    The validated pass would (1) propose variants with
    ``generate_variant_fasta``, (2) re-search the proposal FASTA, (3) compute
    ICS per variant PSM, and (4) accept only variants whose ICS improves on the
    original. Steps 2 to 4 are not implemented. Calling this function raises
    ``NotImplementedError`` so that an unvalidated proposal FASTA is never
    produced under the name "refinement". Use
    ``identify_refinement_candidates`` and ``generate_variant_fasta`` directly
    to obtain proposals.
    """
    raise NotImplementedError(
        "run_refinement is proposal-only: variant validation by re-search and "
        "ICS comparison is not implemented. Call identify_refinement_candidates() "
        "and generate_variant_fasta() to obtain unvalidated proposals."
    )
