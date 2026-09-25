"""
Four inference strategy functions — verbatim logic from e2e-tested scripts.

Each function takes an ``InferenceInput`` and returns an ``InferenceResult``.
The algorithms are preserved exactly as they passed 12/12 golden reference
checks against the NM2026 dataset.
"""

from __future__ import annotations

import logging
from collections import Counter, defaultdict

from fasta_lake.helpers import get_db_type, get_species_id
from fasta_lake.inference.base import InferenceInput, InferenceResult

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Species Budget (metaproteomics-optimized)
# ---------------------------------------------------------------------------


def infer_species_budget(inp: InferenceInput) -> InferenceResult:
    """
    Species-budget parsimony: cross-species aware shared peptide resolution.

    Algorithm
    ---------
    1. Classify peptides as unique (1 protein) or shared (>1 protein).
    2. For within-species shared peptides: razor to protein with most unique peptides.
    3. For cross-species shared peptides: keep best representative from EACH species.
    4. Select all proteins with at least 1 assigned peptide.

    When de novo scores are available (``inp.peptide_scores``), tie-breaking
    uses score-weighted unique peptide evidence instead of raw counts.
    A protein supported by 5 peptides at score 0.95 ranks higher than
    one supported by 5 peptides at score 0.70.

    Parameters
    ----------
    inp : InferenceInput
        Pre-computed peptide-protein mappings and sequences.

    Returns
    -------
    InferenceResult
        Selected protein IDs and statistics.
    """
    pep2prot, prot2pep = inp.pep2prot, inp.prot2pep
    scores = inp.peptide_scores
    use_scores = bool(scores)
    tax_map = inp.taxonomy_map
    do_pruning = inp.taxonomy_pruning and bool(tax_map)

    # Step 1: Group proteins by species
    # Use taxonomy_map if provided (resolves assembly proteins to species),
    # otherwise fall back to get_species_id (parses from protein ID format).
    protein_species: dict[str, str] = {}
    species_proteins: dict[str, set[str]] = defaultdict(set)
    n_from_map = 0
    n_from_id = 0
    unresolved: list[str] = []
    for protein_id in prot2pep:
        if tax_map and protein_id in tax_map:
            species = tax_map[protein_id]
            n_from_map += 1
        else:
            species = get_species_id(protein_id)
            if species is None:
                unresolved.append(protein_id)
                continue
            n_from_id += 1
        protein_species[protein_id] = species
        species_proteins[species].add(protein_id)

    if unresolved:
        # A species budget is undefined for a protein whose species is unknown:
        # pooling such proteins into one bucket silently ran razor under the
        # name species_budget on every content-hash lake. Refuse instead.
        unresolved.sort()
        shown = ", ".join(unresolved[:5]) + (" …" if len(unresolved) > 5 else "")
        raise ValueError(
            f"species_budget: {len(unresolved)} of {len(prot2pep)} proteins with evidence "
            f"have no resolvable species (e.g. {shown}). Provide a taxonomy map "
            "(--taxonomy-map, protein_id<TAB>species) covering them, or use "
            "--strategy razor / uniform, which do not need species."
        )

    # Step 2: Count unique peptides per protein (and score-weighted sum)
    unique_count: Counter[str] = Counter()
    unique_score_sum: dict[str, float] = defaultdict(float)
    unique_peptides: set[str] = set()
    for pep, prots in pep2prot.items():
        if len(prots) == 1:
            prot = next(iter(prots))
            unique_count[prot] += 1
            unique_score_sum[prot] += scores.get(pep, 1.0)
            unique_peptides.add(pep)

    # Tie-breaking key: score-weighted if scores available, count-based otherwise
    def _protein_sort_key(p: str) -> tuple:
        """Rank by optional unique-score support, then evidence counts and accession."""
        if use_scores:
            return (
                -unique_score_sum.get(p, 0.0),
                -unique_count[p],
                -len(prot2pep.get(p, set())),
                p,
            )
        return (
            -unique_count[p],
            -len(prot2pep.get(p, set())),
            p,
        )

    # Step 2b: Species abundance profile (for optional taxonomy-based pruning)
    species_unique_count: Counter[str] = Counter()
    for prot, count in unique_count.items():
        species_unique_count[protein_species[prot]] += count

    # Step 3: Assign peptides using species budget logic
    protein_assigned: dict[str, set[str]] = defaultdict(set)
    n_unique = 0
    n_within_species = 0
    n_cross_species = 0
    n_cross_species_pruned = 0

    for pep, prots in pep2prot.items():
        if len(prots) == 1:
            prot = next(iter(prots))
            protein_assigned[prot].add(pep)
            n_unique += 1
        else:
            species_set = {protein_species[p] for p in prots}
            if len(species_set) == 1:
                best = sorted(prots, key=_protein_sort_key)[0]
                protein_assigned[best].add(pep)
                n_within_species += 1
            else:
                # Optional: frequentist taxonomy pruning.
                # Only active when taxonomy_map is provided AND taxonomy_pruning=True.
                # Uses unique peptide counts per species as abundance prior.
                # Species with < 1% of the dominant species' count are excluded.
                if do_pruning:
                    max_count = max(species_unique_count.get(sp, 0) for sp in species_set)
                    min_threshold = max(1, int(max_count * 0.01))
                    qualified = {
                        sp for sp in species_set if species_unique_count.get(sp, 0) >= min_threshold
                    }
                    if not qualified:
                        qualified = species_set
                    if len(qualified) < len(species_set):
                        n_cross_species_pruned += 1
                else:
                    qualified = species_set

                for species in qualified:
                    species_prots = [p for p in prots if protein_species[p] == species]
                    best = sorted(species_prots, key=_protein_sort_key)[0]
                    protein_assigned[best].add(pep)
                n_cross_species += 1

    # Step 4: Select proteins with assigned peptides
    selected = set(protein_assigned.keys())

    # Database breakdown
    db_stats: dict[str, dict[str, int]] = defaultdict(lambda: {"unique": 0, "shared_only": 0})
    for protein_id in selected:
        db_type = get_db_type(protein_id)
        if unique_count[protein_id] > 0:
            db_stats[db_type]["unique"] += 1
        else:
            db_stats[db_type]["shared_only"] += 1

    species_retained = {protein_species[p] for p in selected}

    proteins_with_unique = sum(1 for c in unique_count.values() if c > 0)

    stats = {
        "strategy": "species_budget",
        "score_weighted": use_scores,
        "taxonomy_map_provided": bool(tax_map),
        "taxonomy_pruning": do_pruning,
        "proteins_with_evidence": len(prot2pep),
        "selected_proteins": len(selected),
        "unique_peptides": len(unique_peptides),
        "shared_peptides": len(pep2prot) - len(unique_peptides),
        "within_species_shared": n_within_species,
        "cross_species_shared": n_cross_species,
        "cross_species_pruned": n_cross_species_pruned,
        "proteins_with_unique_peptides": proteins_with_unique,
        "species_detected": len(species_proteins),
        "species_retained": len(species_retained),
        "species_from_taxonomy_map": n_from_map,
        "species_from_protein_id": n_from_id,
        "db_breakdown": {db: dict(s) for db, s in db_stats.items()},
    }

    return InferenceResult(selected_proteins=selected, stats=stats)


# ---------------------------------------------------------------------------
# Razor (MaxQuant-style)
# ---------------------------------------------------------------------------


def infer_razor(inp: InferenceInput) -> InferenceResult:
    """
    Razor parsimony: assign shared peptides to protein with most unique peptides.

    This is the classical MaxQuant approach used in standard proteomics.
    When de novo scores are available, uses score-weighted tie-breaking.

    Parameters
    ----------
    inp : InferenceInput
        Pre-computed peptide-protein mappings and sequences.

    Returns
    -------
    InferenceResult
        Selected protein IDs and statistics.
    """
    pep2prot, prot2pep = inp.pep2prot, inp.prot2pep
    scores = inp.peptide_scores
    use_scores = bool(scores)

    # Count unique peptides per protein (and score-weighted sum)
    unique_count: Counter[str] = Counter()
    unique_score_sum: dict[str, float] = defaultdict(float)
    unique_peptides: set[str] = set()
    for pep, prots in pep2prot.items():
        if len(prots) == 1:
            prot = next(iter(prots))
            unique_count[prot] += 1
            unique_score_sum[prot] += scores.get(pep, 1.0)
            unique_peptides.add(pep)

    def _protein_sort_key(p: str) -> tuple:
        """Rank the Python razor candidates by support, with lexical accession ties."""
        if use_scores:
            return (
                -unique_score_sum.get(p, 0.0),
                -unique_count[p],
                -len(prot2pep.get(p, set())),
                p,
            )
        return (
            -unique_count[p],
            -len(prot2pep.get(p, set())),
            p,
        )

    # Assign peptides using razor logic
    protein_assigned: dict[str, set[str]] = defaultdict(set)
    for pep, prots in pep2prot.items():
        if len(prots) == 1:
            prot = next(iter(prots))
            protein_assigned[prot].add(pep)
        else:
            best = sorted(prots, key=_protein_sort_key)[0]
            protein_assigned[best].add(pep)

    selected = set(protein_assigned.keys())

    # Database breakdown
    db_stats: dict[str, dict[str, int]] = defaultdict(lambda: {"unique": 0, "shared_only": 0})
    for protein_id in selected:
        db_type = get_db_type(protein_id)
        if unique_count[protein_id] > 0:
            db_stats[db_type]["unique"] += 1
        else:
            db_stats[db_type]["shared_only"] += 1

    proteins_with_unique = sum(1 for c in unique_count.values() if c > 0)

    stats = {
        "strategy": "razor",
        "proteins_with_evidence": len(prot2pep),
        "selected_proteins": len(selected),
        "unique_peptides": len(unique_peptides),
        "shared_peptides": len(pep2prot) - len(unique_peptides),
        "proteins_with_unique_peptides": proteins_with_unique,
        "db_breakdown": {db: dict(s) for db, s in db_stats.items()},
    }

    return InferenceResult(selected_proteins=selected, stats=stats)


# ---------------------------------------------------------------------------
# Uniform (baseline — keep everything)
# ---------------------------------------------------------------------------


def infer_uniform(inp: InferenceInput) -> InferenceResult:
    """
    Uniform parsimony: keep ALL proteins with any peptide evidence.

    Baseline approach with no filtering. Useful for comparison and
    for letting downstream tools (e.g., SAGE) handle protein inference.

    Parameters
    ----------
    inp : InferenceInput
        Pre-computed peptide-protein mappings and sequences.

    Returns
    -------
    InferenceResult
        Selected protein IDs and statistics.
    """
    pep2prot, prot2pep = inp.pep2prot, inp.prot2pep

    unique_peptides: set[str] = set()
    shared_peptides: set[str] = set()
    for pep, prots in pep2prot.items():
        if len(prots) == 1:
            unique_peptides.add(pep)
        else:
            shared_peptides.add(pep)

    selected = set(prot2pep.keys())

    # Database breakdown
    db_stats: dict[str, dict[str, int]] = defaultdict(lambda: {"unique": 0, "shared_only": 0})
    for protein_id in selected:
        db_type = get_db_type(protein_id)
        has_unique = any(len(pep2prot[pep]) == 1 for pep in prot2pep[protein_id])
        if has_unique:
            db_stats[db_type]["unique"] += 1
        else:
            db_stats[db_type]["shared_only"] += 1

    stats = {
        "strategy": "uniform",
        "proteins_with_evidence": len(prot2pep),
        "selected_proteins": len(selected),
        "unique_peptides": len(unique_peptides),
        "shared_peptides": len(shared_peptides),
        "db_breakdown": {db: dict(s) for db, s in db_stats.items()},
    }

    return InferenceResult(selected_proteins=selected, stats=stats)


# ---------------------------------------------------------------------------
# Uniform 2-peptide (conservative baseline)
# ---------------------------------------------------------------------------


def infer_uniform_2pep(
    inp: InferenceInput,
    min_peptides: int = 2,
) -> InferenceResult:
    """
    Uniform 2-peptide: require at least ``min_peptides`` per protein.

    More conservative than uniform — reduces false positives by
    requiring multiple peptide evidence.

    Parameters
    ----------
    inp : InferenceInput
        Pre-computed peptide-protein mappings and sequences.
    min_peptides : int
        Minimum number of peptides required (default: 2).

    Returns
    -------
    InferenceResult
        Selected protein IDs and statistics.
    """
    pep2prot, prot2pep = inp.pep2prot, inp.prot2pep

    unique_peptides: set[str] = set()
    shared_peptides: set[str] = set()
    for pep, prots in pep2prot.items():
        if len(prots) == 1:
            unique_peptides.add(pep)
        else:
            shared_peptides.add(pep)

    selected = {prot for prot, peps in prot2pep.items() if len(peps) >= min_peptides}
    excluded = len(prot2pep) - len(selected)

    # Database breakdown
    db_stats: dict[str, dict[str, int]] = defaultdict(lambda: {"selected": 0, "excluded": 0})
    for protein_id, peps in prot2pep.items():
        db_type = get_db_type(protein_id)
        if protein_id in selected:
            db_stats[db_type]["selected"] += 1
        else:
            db_stats[db_type]["excluded"] += 1

    stats = {
        "strategy": "uniform_2pep",
        "min_peptides": min_peptides,
        "proteins_with_evidence": len(prot2pep),
        "selected_proteins": len(selected),
        "excluded_proteins": excluded,
        "unique_peptides": len(unique_peptides),
        "shared_peptides": len(shared_peptides),
        "db_breakdown": {db: dict(s) for db, s in db_stats.items()},
    }

    return InferenceResult(selected_proteins=selected, stats=stats)


# ---------------------------------------------------------------------------
# Info Score (information-theoretic parsimony)
# ---------------------------------------------------------------------------


def infer_info_score(
    inp: InferenceInput,
    min_info: float = 0.1,
) -> InferenceResult:
    """
    Information-content parsimony: keep proteins above an info score threshold.

    Each protein gets an info score that quantifies how much *informative*
    evidence supports it from the de novo predictions:

        info_score(protein) = sum over matching peptides of:
            denovo_score(peptide) / n_proteins_matching(peptide)

    A peptide with score 0.95 matching 1 protein contributes 0.95.
    The same peptide matching 500 proteins contributes 0.0019 each.
    High-confidence unique peptides dominate; shared peptides contribute
    proportionally to their ambiguity.

    This avoids the problems of:
    - razor: arbitrarily picks one winner, loses real species
    - uniform_2pep: keeps millions of proteins with shared-only evidence
    - species_budget: requires species-level taxonomy (not always available)

    The threshold ``min_info`` controls the tradeoff:
    - 0.01: ~200K proteins (broad, includes low-evidence)
    - 0.10: ~50K proteins  (balanced, good for SAGE)
    - 0.50: ~6K proteins   (very clean, may miss real low-abundance species)

    When no scores are available, falls back to: 1/n_proteins per peptide.

    Parameters
    ----------
    inp : InferenceInput
        Pre-computed peptide-protein mappings and sequences.
    min_info : float
        Minimum information score threshold (default: 0.1).

    Returns
    -------
    InferenceResult
        Selected protein IDs and statistics.
    """
    pep2prot, prot2pep = inp.pep2prot, inp.prot2pep
    scores = inp.peptide_scores
    use_scores = bool(scores)

    # Compute info score per protein
    info_scores: dict[str, float] = {}
    for prot_id, peps in prot2pep.items():
        info = 0.0
        # DETERMINISM FIX: `peps` is a set of peptide strings, so the summation
        # order depended on the per-process hash seed. Float addition is not
        # associative, so info_scores differed in the low bits between runs --
        # verified: reported info_score p75 read 7.4578 on 9/11 hash seeds and
        # 7.4577 on the other 2, and every seed produced a different digest of
        # the full info_scores vector. Summing in sorted peptide order makes the
        # value a function of the peptide set alone. (The `>= min_info` cut can
        # also flip for a protein sitting exactly on the threshold.)
        for pep in sorted(peps):
            n_match = len(pep2prot.get(pep, set()))
            if n_match == 0:
                continue
            denovo_score = scores.get(pep, 1.0) if use_scores else 1.0
            info += denovo_score / n_match
        info_scores[prot_id] = info

    # Select proteins above threshold
    selected = {p for p, score in info_scores.items() if score >= min_info}

    # Compute unique vs shared peptide stats
    unique_peptides: set[str] = set()
    for pep, prots in pep2prot.items():
        if len(prots) == 1:
            unique_peptides.add(pep)

    # Proteins with unique evidence in selected set
    selected_with_unique = sum(
        1 for p in selected if any(len(pep2prot[pep]) == 1 for pep in prot2pep[p])
    )

    # Database breakdown
    db_stats: dict[str, dict[str, int]] = defaultdict(lambda: {"selected": 0, "excluded": 0})
    for protein_id in prot2pep:
        db_type = get_db_type(protein_id)
        if protein_id in selected:
            db_stats[db_type]["selected"] += 1
        else:
            db_stats[db_type]["excluded"] += 1

    # Info score distribution of selected
    selected_scores = sorted([info_scores[p] for p in selected], reverse=True)
    percentiles = {}
    if selected_scores:
        import statistics

        percentiles = {
            "median": round(statistics.median(selected_scores), 4),
            "p25": round(selected_scores[int(len(selected_scores) * 0.75)], 4),
            "p75": round(selected_scores[int(len(selected_scores) * 0.25)], 4),
            "max": round(selected_scores[0], 4),
            "min": round(selected_scores[-1], 4),
        }

    stats = {
        "strategy": "info_score",
        "min_info": min_info,
        "score_weighted": use_scores,
        "proteins_with_evidence": len(prot2pep),
        "selected_proteins": len(selected),
        "excluded_proteins": len(prot2pep) - len(selected),
        "unique_peptides": len(unique_peptides),
        "shared_peptides": len(pep2prot) - len(unique_peptides),
        "selected_with_unique_peptides": selected_with_unique,
        "info_score_percentiles": percentiles,
        "db_breakdown": {db: dict(s) for db, s in db_stats.items()},
    }

    return InferenceResult(selected_proteins=selected, stats=stats)


# ---------------------------------------------------------------------------
# Bayesian Protein Probability (evidence-driven parsimony)
# ---------------------------------------------------------------------------


def infer_bayesian(
    inp: InferenceInput,
    n_iterations: int = 10,
    min_probability: float = 0.01,
) -> InferenceResult:
    """
    Bayesian protein probability: iterative evidence redistribution.

    Models the probability that each protein is truly present in the sample,
    using de novo peptide evidence. Unlike razor (which forces a binary
    winner/loser decision for shared peptides), this distributes evidence
    proportionally and iteratively refines probabilities.

    Algorithm (Expectation-Maximization)
    ------------------------------------

    **Initialization:**
    Each protein starts with a prior probability proportional to its
    total peptide evidence (number of matching peptides, weighted by
    de novo confidence scores). This gives proteins with more evidence
    a head start, but doesn't commit to winners yet.

    **E-step (expectation):**
    For each shared peptide, distribute its evidence across matching
    proteins proportional to their current probabilities:

        weight(protein, peptide) = P(protein) / sum(P(p) for p matching peptide)

    A protein with high probability "claims" more of each shared peptide.
    Unique peptides always contribute fully to their single protein.

    **M-step (maximization):**
    Update each protein's probability based on its total claimed evidence,
    normalized by protein length to avoid bias toward long proteins:

        P(protein) ∝ sum(denovo_score(pep) × weight(protein, pep)) / sqrt(length)

    The sqrt(length) normalization is gentler than dividing by length
    (which would over-penalize long proteins) but still corrects for
    the fact that longer proteins match more peptides by chance.

    **Convergence:**
    Repeat E/M steps until probabilities stabilize (typically 5-10
    iterations). Proteins below ``min_probability`` are excluded.

    Why this is better than info_score
    -----------------------------------
    1. **Iterative refinement:** info_score computes 1/n_matching once.
       Bayesian iterates — as weak proteins lose probability, strong
       proteins claim more of the shared peptides, amplifying signal.

    2. **Length normalization:** avoids bias toward long proteins that
       match more peptides by chance.

    3. **Convergent evidence:** if 5 peptides each match 50 proteins
       but all 5 share the same top protein, that protein accumulates
       evidence across iterations. Info_score misses this synergy.

    4. **Soft decisions:** doesn't force winner-take-all like razor.
       Two closely related strains with similar evidence both survive
       at moderate probability rather than one being eliminated.

    Parameters
    ----------
    inp : InferenceInput
        Pre-computed peptide-protein mappings and sequences.
    n_iterations : int
        Number of EM iterations (default: 10).
    min_probability : float
        Minimum probability threshold for inclusion (default: 0.01).

    Returns
    -------
    InferenceResult
        Selected protein IDs and statistics.
    """
    import math

    pep2prot, prot2pep = inp.pep2prot, inp.prot2pep
    proteins = inp.proteins
    scores = inp.peptide_scores
    use_scores = bool(scores)

    all_prots = list(prot2pep.keys())
    if not all_prots:
        return InferenceResult(selected_proteins=set(), stats={"strategy": "bayesian"})

    # Protein lengths (for normalization)
    prot_len = {}
    for pid in all_prots:
        seq = proteins.get(pid, "")
        prot_len[pid] = max(len(seq), 50)  # floor at 50 to avoid division issues

    # Initialize: prior proportional to total score-weighted peptide count
    prob: dict[str, float] = {}
    for pid in all_prots:
        # DETERMINISM GUARD: prot2pep[pid] is a set of peptide strings, so the
        # summation order depended on the per-process hash seed and float
        # addition is not associative. Sum in sorted peptide order.
        raw = (
            sum(scores.get(pep, 1.0) for pep in sorted(prot2pep[pid]))
            if use_scores
            else len(prot2pep[pid])
        )
        prob[pid] = raw / math.sqrt(prot_len[pid])

    # Normalize to sum to 1
    total = sum(prob.values())
    if total > 0:
        for pid in prob:
            prob[pid] /= total

    convergence_history = []

    # EM iterations
    for iteration in range(n_iterations):
        # E-step: compute weights for shared peptides
        new_evidence: dict[str, float] = defaultdict(float)

        for pep, matching_prots in pep2prot.items():
            denovo_score = scores.get(pep, 1.0) if use_scores else 1.0

            if len(matching_prots) == 1:
                # Unique peptide: full evidence to its protein
                pid = next(iter(matching_prots))
                new_evidence[pid] += denovo_score
            else:
                # Shared peptide: distribute proportional to current probabilities
                # DETERMINISM GUARD: matching_prots is a set of protein-ID
                # strings; the float sum below was therefore hash-seed ordered.
                # Sort once and reuse for the sum and both accumulation loops.
                ordered_prots = sorted(matching_prots)
                prob_sum = sum(prob.get(p, 0) for p in ordered_prots)
                if prob_sum == 0:
                    # Equal split if all probabilities are zero
                    share = denovo_score / len(ordered_prots)
                    for p in ordered_prots:
                        new_evidence[p] += share
                else:
                    for p in ordered_prots:
                        weight = prob.get(p, 0) / prob_sum
                        new_evidence[p] += denovo_score * weight

        # M-step: update probabilities with length normalization
        old_prob = dict(prob)
        for pid in all_prots:
            prob[pid] = new_evidence.get(pid, 0) / math.sqrt(prot_len[pid])

        # Normalize
        total = sum(prob.values())
        if total > 0:
            for pid in prob:
                prob[pid] /= total

        # Track convergence (max change)
        max_delta = max(abs(prob.get(pid, 0) - old_prob.get(pid, 0)) for pid in all_prots)
        convergence_history.append(max_delta)

        if max_delta < 1e-8:
            break

    # Select proteins above threshold
    # Scale probabilities so the maximum is 1.0 for interpretability
    max_prob = max(prob.values()) if prob else 1.0
    if max_prob > 0:
        scaled = {pid: prob[pid] / max_prob for pid in prob}
    else:
        scaled = prob

    selected = {pid for pid, p in scaled.items() if p >= min_probability}

    # Stats
    unique_peptides = sum(1 for pep, prots in pep2prot.items() if len(prots) == 1)
    selected_with_unique = sum(
        1 for p in selected if any(len(pep2prot[pep]) == 1 for pep in prot2pep[p])
    )

    # Database breakdown
    db_stats: dict[str, dict[str, int]] = defaultdict(lambda: {"selected": 0, "excluded": 0})
    for pid in prot2pep:
        db_type = get_db_type(pid)
        if pid in selected:
            db_stats[db_type]["selected"] += 1
        else:
            db_stats[db_type]["excluded"] += 1

    # Probability distribution of selected
    selected_probs = sorted([scaled[p] for p in selected], reverse=True)
    percentiles = {}
    if selected_probs:
        import statistics

        percentiles = {
            "median": round(statistics.median(selected_probs), 6),
            "p25": round(selected_probs[int(len(selected_probs) * 0.75)], 6),
            "p75": round(selected_probs[int(len(selected_probs) * 0.25)], 6),
            "max": round(selected_probs[0], 6),
            "min": round(selected_probs[-1], 6),
        }

    stats = {
        "strategy": "bayesian",
        "n_iterations": len(convergence_history),
        "min_probability": min_probability,
        "converged": convergence_history[-1] < 1e-8 if convergence_history else False,
        "final_max_delta": convergence_history[-1] if convergence_history else 0,
        "score_weighted": use_scores,
        "proteins_with_evidence": len(prot2pep),
        "selected_proteins": len(selected),
        "excluded_proteins": len(prot2pep) - len(selected),
        "unique_peptides": unique_peptides,
        "shared_peptides": len(pep2prot) - unique_peptides,
        "selected_with_unique_peptides": selected_with_unique,
        "probability_percentiles": percentiles,
        "db_breakdown": {db: dict(s) for db, s in db_stats.items()},
    }

    return InferenceResult(selected_proteins=selected, stats=stats)


# ---------------------------------------------------------------------------
# Cross-Sample Bayesian (lake-level evidence scoring)
# ---------------------------------------------------------------------------


def score_lake_bayesian(
    lake_proteins: dict[str, str],
    sample_peptides: dict[str, dict[str, float]],
    n_iterations: int = 10,
) -> dict[str, float]:
    """
    Score every protein in the evidence lake using Bayesian EM across ALL samples.

    This runs at the LAKE level, not per-sample. It uses de novo predictions
    from every sample to compute a single probability for each protein in the
    lake. Proteins consistently supported across multiple samples get boosted;
    proteins supported only by highly-shared peptides in one sample get
    downweighted.

    Parameters
    ----------
    lake_proteins : dict[str, str]
        Protein ID → sequence for the entire evidence lake.
    sample_peptides : dict[str, dict[str, float]]
        sample_id → {peptide → denovo_score}.
        Each sample's top de novo predictions.
    n_iterations : int
        Number of EM iterations.

    Returns
    -------
    dict[str, float]
        Protein ID → probability (0 to 1, where 1 = most confident).
        Can be used to threshold the lake before per-sample extraction.

    Algorithm
    ---------
    1. Build a unified peptide→protein map across the entire lake.
    2. For each sample, each peptide contributes evidence weighted by:
       - its de novo confidence score
       - 1 / n_proteins matching (uniqueness)
       - current protein probability (from previous iteration)
    3. Cross-sample evidence is SUMMED — a protein supported by unique
       peptides in 10 samples gets 10x the evidence of one in 1 sample.
    4. Normalize by sqrt(protein length).
    5. After convergence, the probability reflects:
       - How many samples provide evidence
       - How strong and specific that evidence is
       - How consistently the protein appears across the experiment
    """
    import math

    all_prot_ids = list(lake_proteins.keys())
    if not all_prot_ids:
        return {}

    # Build peptide → protein map for entire lake
    logger.info("Building peptide-protein map for %d proteins...", len(all_prot_ids))
    from fasta_lake.helpers import map_peptides_to_proteins

    # Collect all unique peptides across all samples
    all_peptides: dict[str, float] = {}
    for sample_id, peps in sample_peptides.items():
        for pep, score in peps.items():
            if pep not in all_peptides or score > all_peptides[pep]:
                all_peptides[pep] = score

    pep2prot, prot2pep = map_peptides_to_proteins(lake_proteins, all_peptides)

    # Protein lengths
    prot_len = {pid: max(len(seq), 50) for pid, seq in lake_proteins.items()}

    # Count per-sample evidence for cross-sample boosting
    # sample_pep2prot: for each sample, which peptides map to which proteins
    sample_evidence: dict[str, dict[str, float]] = defaultdict(lambda: defaultdict(float))
    for sample_id, peps in sample_peptides.items():
        for pep, score in peps.items():
            if pep in pep2prot:
                for pid in pep2prot[pep]:
                    sample_evidence[sample_id][pid] += score / len(pep2prot[pep])

    # Count in how many samples each protein has evidence
    prot_sample_count = Counter()
    for sample_id, prot_scores in sample_evidence.items():
        for pid in prot_scores:
            prot_sample_count[pid] += 1

    # Initialize probabilities
    prob: dict[str, float] = {}
    for pid in all_prot_ids:
        if pid in prot2pep:
            # DETERMINISM GUARD: prot2pep[pid] is a set of peptide strings, so
            # the summation order was hash-seed dependent and float addition is
            # not associative. Sum in sorted peptide order.
            raw = sum(
                all_peptides.get(pep, 1.0) / len(pep2prot.get(pep, {pid}))
                for pep in sorted(prot2pep[pid])
            )
            # Boost by sample frequency
            sample_boost = math.sqrt(prot_sample_count.get(pid, 1))
            prob[pid] = raw * sample_boost / math.sqrt(prot_len[pid])
        else:
            prob[pid] = 0.0

    total = sum(prob.values())
    if total > 0:
        for pid in prob:
            prob[pid] /= total

    # EM iterations
    for iteration in range(n_iterations):
        new_evidence: dict[str, float] = defaultdict(float)

        for pep, matching_prots in pep2prot.items():
            denovo_score = all_peptides.get(pep, 1.0)

            if len(matching_prots) == 1:
                pid = next(iter(matching_prots))
                # Boost unique peptides by number of samples that see them
                n_samples_with_pep = sum(1 for sid, sp in sample_peptides.items() if pep in sp)
                new_evidence[pid] += denovo_score * math.sqrt(n_samples_with_pep)
            else:
                # DETERMINISM GUARD: matching_prots is a set of protein-ID
                # strings; sort once and reuse for the float sum and both
                # accumulation loops.
                ordered_prots = sorted(matching_prots)
                prob_sum = sum(prob.get(p, 0) for p in ordered_prots)
                if prob_sum == 0:
                    share = denovo_score / len(ordered_prots)
                    for p in ordered_prots:
                        new_evidence[p] += share
                else:
                    for p in ordered_prots:
                        weight = prob.get(p, 0) / prob_sum
                        new_evidence[p] += denovo_score * weight

        old_prob = dict(prob)
        for pid in all_prot_ids:
            sample_boost = math.sqrt(prot_sample_count.get(pid, 1))
            prob[pid] = new_evidence.get(pid, 0) * sample_boost / math.sqrt(prot_len[pid])

        total = sum(prob.values())
        if total > 0:
            for pid in prob:
                prob[pid] /= total

        max_delta = max(abs(prob.get(pid, 0) - old_prob.get(pid, 0)) for pid in all_prot_ids)
        if max_delta < 1e-8:
            logger.info("Converged at iteration %d (delta=%.2e)", iteration, max_delta)
            break

    # Scale to 0-1
    max_prob = max(prob.values()) if prob else 1.0
    if max_prob > 0:
        return {pid: prob[pid] / max_prob for pid in prob}
    return prob


# ---------------------------------------------------------------------------
# Strategy registry
# ---------------------------------------------------------------------------

STRATEGIES = {
    "species_budget": infer_species_budget,
    "razor": infer_razor,
    "uniform": infer_uniform,
    "uniform_2pep": infer_uniform_2pep,
    "info_score": infer_info_score,
    "bayesian": infer_bayesian,
}
