"""Entrapment generation and conditional FDP estimates.

For N_E entrapment and N_T original-target discoveries, the combined estimate
is N_E * (1 + 1/r) / (N_T + N_E), where r is the effective entrapment-to-target
search-space ratio at the counting level. Under the equal-chance-type
assumptions in Wen et al. (2025), this is an upper estimate in expectation.
It is not a per-experiment confidence bound. A zero observed count does not
prove that the underlying error rate is zero.

The ratio N_E / (N_T + N_E) is a lower bound and cannot validate successful
FDR control. The paired estimator requires an explicit one-to-one pairing
of unique target and entrapment entities, scored on a comparable scale.
Holding K/R positions while shuffling whole proteins does not establish
such a peptide pairing, nor does partitioning an ordinary reference into
two random halves establish that one half is absent from the sample.

Empty target or entrapment partitions and zero total discoveries are
undefined. Callers must inspect ``evaluable`` before quoting an estimate.
The injection stage, counting unit, conserved matches, ratio definition,
and null comparability must accompany any biological interpretation.

Reference: Wen B, Freestone J, Riffle M, MacCoss MJ, Noble WS, Keich U.
Assessment of false discovery rate control in tandem mass spectrometry
analysis using entrapment. Nature Methods (2025).
https://doi.org/10.1038/s41592-025-02719-x
"""

from __future__ import annotations

import logging
import random
from pathlib import Path

from fasta_lake.config import EntrapmentConfig
from fasta_lake.helpers import load_fasta_proteins

logger = logging.getLogger(__name__)


def _validate_counts(**counts: int) -> None:
    """Reject negative counts, non-integers and booleans before FDP accounting."""
    for name, value in counts.items():
        if isinstance(value, bool) or not isinstance(value, int) or value < 0:
            raise ValueError(f"{name} must be a nonnegative integer; got {value!r}")


def generate_shuffled_entrapment(
    source_proteins: dict[str, str],
    n_proteins: int = 1000,
    seed: int = 42,
    prefix: str = "ENTRAP_SHUFFLED",
    preserve_cleavage_sites: bool = True,
) -> dict[str, str]:
    """
    Generate shuffled entrapment proteins from source sequences.

    Each shuffled sequence preserves the amino acid composition of its
    source protein but randomizes the order. This ensures:

    - Same mass distribution as real proteins
    - A reproducible sequence perturbation, conditional on seed and input order

    This does not guarantee distinct or known-absent peptides. Target overlap,
    tryptic constraints, homology, and null comparability require validation.

    .. note::
       ``preserve_cleavage_sites`` defaults to ``True``, which holds K and R
       in place so each entrapment protein's tryptic peptide lengths match
       its source protein's K/R position pattern. This matches the reference
       implementation used for the manuscript
       (``projects/masters_circularity_F1/scripts/make_shuffled_entrapment.py``).

       Earlier versions of this function did a naive full shuffle. That
       preserves the K/R *count*, so the aggregate peptide-length
       distribution is about right, but not the per-protein pairing:
       measured on the 512-protein P. furiosus panel, the naive shuffle
       reproduced a protein's exact tryptic length multiset for 1 of 512
       proteins (0.2%) against 512 of 512 (100%) when K/R are held.
       Pass ``preserve_cleavage_sites=False`` for the old behaviour.

    For lake-scale work use
    :func:`fasta_lake.entrapment_workflow.generate_entrapment_fasta`, which
    streams instead of building a dict in memory.

    Parameters
    ----------
    source_proteins : dict[str, str]
        Source proteins to shuffle (protein_id -> sequence).
    n_proteins : int
        Number of entrapment proteins to generate. If larger than
        source pool, all source proteins are used. Default: 1000.
    seed : int
        Random seed for reproducibility. Default: 42.
    prefix : str
        Prefix for entrapment protein IDs. Default: ``'ENTRAP_SHUFFLED'``.
    preserve_cleavage_sites : bool
        Hold K/R residues in place while shuffling. Default: ``True``.

    Returns
    -------
    dict[str, str]
        Entrapment proteins (``'{prefix}_{i:05d}'`` -> shuffled sequence).
    """
    rng = random.Random(seed)

    # Sample source proteins
    source_ids = sorted(source_proteins.keys())
    if n_proteins < len(source_ids):
        source_ids = rng.sample(source_ids, n_proteins)

    entrapment: dict[str, str] = {}
    for i, source_id in enumerate(source_ids):
        source_seq = source_proteins[source_id]
        if preserve_cleavage_sites:
            movable = [j for j, c in enumerate(source_seq) if c not in ("K", "R")]
            chars = [source_seq[j] for j in movable]
            rng.shuffle(chars)
            out = list(source_seq)
            for j, c in zip(movable, chars):
                out[j] = c
            seq = out
        else:
            seq = list(source_seq)
            rng.shuffle(seq)
        entrap_id = f"{prefix}_{i:05d}"
        entrapment[entrap_id] = "".join(seq)

    logger.info(
        "Generated %d shuffled entrapment proteins from %d sources",
        len(entrapment),
        len(source_proteins),
    )

    return entrapment


def load_entrapment_fasta(
    fasta_path: str | Path,
    prefix: str = "ENTRAP",
) -> dict[str, str]:
    """
    Load user-supplied entrapment proteins from a FASTA file.

    Optionally prefixes protein IDs with ``prefix`` to mark them as
    entrapment in downstream analysis.

    .. warning::
       This builds the whole file into a dict and keys it on the **first
       header token**, so descriptions are discarded and memory scales with
       the input. It is fine for an entrapment panel of a few thousand
       proteins; it is not usable at lake scale. For streaming, header-
       preserving work use
       :func:`fasta_lake.entrapment_workflow.generate_entrapment_fasta` and
       :func:`fasta_lake.entrapment_workflow.iter_fasta`.

    Parameters
    ----------
    fasta_path : str or Path
        Path to entrapment FASTA file.
    prefix : str
        Prefix to add to protein IDs. Set to empty string to keep
        original IDs. Default: ``'ENTRAP'``.

    Returns
    -------
    dict[str, str]
        Entrapment proteins (optionally prefixed ID -> sequence).
    """
    raw = load_fasta_proteins(fasta_path)

    if prefix:
        entrapment = {f"{prefix}_{pid}": seq for pid, seq in raw.items()}
    else:
        entrapment = raw

    logger.info("Loaded %d entrapment proteins from %s", len(entrapment), fasta_path)
    return entrapment


def generate_noble_paired(
    proteins: dict[str, str],
    seed: int = 42,
) -> tuple[dict[str, str], dict[str, str]]:
    """
    Noble paired entrapment: split proteins into target and entrapment sets.

    The Noble method randomly partitions proteins into two equal-sized
    groups. Search both against the sample. The entrapment group hit rate
    estimates FDR for the target group.

    Parameters
    ----------
    proteins : dict[str, str]
        Full protein database to partition.
    seed : int
        Random seed for reproducibility.

    Returns
    -------
    target : dict[str, str]
        Target proteins (half of input).
    entrapment : dict[str, str]
        Entrapment proteins (other half, prefixed with ``'ENTRAP_NOBLE_'``).
    """
    rng = random.Random(seed)
    ids = sorted(proteins.keys())
    rng.shuffle(ids)

    mid = len(ids) // 2
    target_ids = ids[:mid]
    entrap_ids = ids[mid:]

    target = {pid: proteins[pid] for pid in target_ids}
    entrapment = {f"ENTRAP_NOBLE_{pid}": proteins[pid] for pid in entrap_ids}

    logger.info(
        "Noble paired: %d target, %d entrapment",
        len(target),
        len(entrapment),
    )

    return target, entrapment


def estimate_fdp(
    entrapment_hits: int,
    target_hits: int,
    target_db_size: int,
    entrapment_db_size: int,
    method: str = "combined",
) -> dict:
    """
    Estimate False Discovery Proportion using entrapment hits.

    Implements the methods from Wen, Keich & Noble (2025) Nature Methods.

    Parameters
    ----------
    entrapment_hits : int
        Number of PSMs matching entrapment proteins (N_E).
    target_hits : int
        Number of PSMs matching target proteins (N_T).
    target_db_size : int
        Number of proteins in the target database.
    entrapment_db_size : int
        Number of proteins in the entrapment database.
    method : str
        FDP estimation method:
        - ``'combined'``: Upper estimate in expectation under the null assumptions (Eq 1).
        - ``'lower_bound'``: Lower bound only (Eq 2). Can only prove
          failure, NOT success. Included for completeness but should
          never be the sole reported metric.
        - ``'paired'``: Requires r=1 and score-order pairing information,
          which this function does not take. Call
          :func:`estimate_fdp_paired` instead for the full Wen Eq 4 paired
          analysis; this function cannot compute it from scalar counts.
        Default: ``'combined'``.

    Returns
    -------
    dict
        ``n_entrapment``, ``n_target``, ``n_total``, ``r`` (db ratio),
        ``fdp_lower_bound``, ``fdp_combined``, ``method``.
    """
    _validate_counts(
        entrapment_hits=entrapment_hits,
        target_hits=target_hits,
        target_db_size=target_db_size,
        entrapment_db_size=entrapment_db_size,
    )
    if method not in {"combined", "lower_bound"}:
        raise ValueError("method must be combined or lower_bound; paired needs score pairing")
    n_total = target_hits + entrapment_hits
    # r = entrapment_db_size / target_db_size (Noble 2025 convention)
    # When r < 1, entrapment is smaller than target, so each entrapment
    # hit represents more false discoveries → (1 + 1/r) inflates upward.
    r = entrapment_db_size / target_db_size if target_db_size > 0 else None

    # An estimate is only defined when something was actually discovered AND
    # the searched database actually contained an entrapment partition to hit.
    # Zero discoveries, or an empty entrapment partition, means the test could
    # not run — NOT that it passed. See `evaluable` in the returned dict.
    evaluable = n_total > 0 and entrapment_db_size > 0 and target_db_size > 0

    # Eq 2: Lower bound — can only demonstrate failure, not success
    fdp_lower = entrapment_hits / n_total if evaluable else None

    # Eq 1: Combined — valid upper bound on FDP
    # FDP = N_E × (1 + 1/r) / (N_T + N_E)
    # When r=1 (equal DBs): reduces to 2 × N_E / (N_T + N_E), same as TDC
    fdp_combined = (entrapment_hits * (1 + 1 / r)) / n_total if (evaluable and r > 0) else None

    result = {
        "n_entrapment": entrapment_hits,
        "n_target": target_hits,
        "n_total": n_total,
        "r": round(r, 4) if r is not None else None,
        "evaluable": evaluable,
        "fdp_lower_bound": round(fdp_lower, 6) if fdp_lower is not None else None,
        "fdp_combined": round(fdp_combined, 6) if fdp_combined is not None else None,
        "method": method,
    }

    if evaluable:
        logger.info(
            "FDP estimate: %d/%d entrapment hits (r=%.4f), lower=%.4f, combined=%.4f",
            entrapment_hits,
            n_total,
            r,
            fdp_lower,
            fdp_combined,
        )
    else:
        logger.warning(
            "FDP UNDEFINED (not zero): %d discoveries, entrapment_db_size=%d. "
            "Nothing could be detected, so this is not a pass.",
            n_total,
            entrapment_db_size,
        )

    return result


def estimate_fdp_paired(
    n_target: int,
    n_entrapment: int,
    n_entrap_target_undiscovered: int,
    n_entrap_gt_target: int,
) -> dict:
    """
    Wen 2025 Eq 4 paired FDP estimator — a *tighter* valid upper bound.

    The paper's claim is that this is tighter than the COMBINED method, not that
    it is the tightest bound available: it also introduces a "k-matched"
    generalization (r = k) and proves both valid. Verbatim, §2.2 heading: "The
    paired estimation method yields a tighter upper bound on the FDP".

    Implements the paired estimation method of Wen, Freestone, Riffle, MacCoss,
    Noble & Keich (2025), Nature Methods 22, 1454-1463::

        FDP* = (N_E + N_{E>=s>T} + 2 * N_{E>T>=s}) / (N_T + N_E)

    Unlike the combined method (:func:`estimate_fdp`), the paired method uses
    the score-order relationship between each entrapment peptide and its
    *unique paired* target peptide to reduce conservative bias while remaining
    an upper estimate in expectation under the paper's assumptions. It requires
    ``r = 1`` and a verified one-to-one pairing of distinct peptide entities.
    A whole-protein K/R-preserving shuffle alone does not establish this pairing.
    Target and entrapment scores must be comparable and obtained in the
    declared common search experiment.

    Partitioning of the ``N_E`` *discovered* entrapment peptides:

    - ``n_entrap_target_undiscovered`` = ``N_{E>=s>T}``: entrapment peptide
      discovered, its paired target peptide NOT discovered (target score < s).
    - ``n_entrap_gt_target`` = ``N_{E>T>=s}``: entrapment peptide discovered
      AND its paired target also discovered, but the entrapment scored higher.
    - the remaining discovered entrapment peptides (paired target discovered
      and scored at least as high) contribute only their base count.

    Parameters
    ----------
    n_target : int
        ``N_T``, number of discovered pure-target peptides at the cutoff.
    n_entrapment : int
        ``N_E``, number of discovered entrapment peptides at the cutoff.
    n_entrap_target_undiscovered : int
        ``N_{E>=s>T}`` (see above). Must be <= ``n_entrapment``.
    n_entrap_gt_target : int
        ``N_{E>T>=s}`` (see above). Must be <= ``n_entrapment``.

    Returns
    -------
    dict
        ``n_target``, ``n_entrapment``, ``n_entrap_target_undiscovered``,
        ``n_entrap_gt_target``, ``numerator``, ``fdp_paired`` (the upper
        bound, or ``None`` when undefined), ``evaluable``, ``method``.

    Notes
    -----
    Undefined is not zero: when nothing was discovered (``N_T + N_E == 0``)
    the estimate is ``None`` and ``evaluable`` is ``False`` — the test could
    not run, which is not a pass. See the module docstring.
    """
    _validate_counts(
        n_target=n_target,
        n_entrapment=n_entrapment,
        n_entrap_target_undiscovered=n_entrap_target_undiscovered,
        n_entrap_gt_target=n_entrap_gt_target,
    )
    if n_entrap_target_undiscovered > n_entrapment or n_entrap_gt_target > n_entrapment:
        raise ValueError(
            "paired sub-counts cannot exceed n_entrapment "
            f"(N_E={n_entrapment}, N_{{E>=s>T}}={n_entrap_target_undiscovered}, "
            f"N_{{E>T>=s}}={n_entrap_gt_target})"
        )
    if (n_entrap_target_undiscovered + n_entrap_gt_target) > n_entrapment:
        raise ValueError(
            "N_{E>=s>T} and N_{E>T>=s} are disjoint subsets of the discovered "
            "entrapment peptides; their sum cannot exceed n_entrapment "
            f"({n_entrap_target_undiscovered} + {n_entrap_gt_target} > {n_entrapment})"
        )

    n_total = n_target + n_entrapment
    evaluable = n_total > 0
    numerator = n_entrapment + n_entrap_target_undiscovered + 2 * n_entrap_gt_target
    fdp_paired = (numerator / n_total) if evaluable else None

    result = {
        "n_target": n_target,
        "n_entrapment": n_entrapment,
        "n_entrap_target_undiscovered": n_entrap_target_undiscovered,
        "n_entrap_gt_target": n_entrap_gt_target,
        "numerator": numerator,
        "evaluable": evaluable,
        "fdp_paired": round(fdp_paired, 6) if fdp_paired is not None else None,
        "method": "paired",
    }

    if evaluable:
        logger.info(
            "Paired FDP (Wen Eq 4): num=%d / (N_T=%d + N_E=%d) = %.6f",
            numerator,
            n_target,
            n_entrapment,
            fdp_paired,
        )
    else:
        logger.warning("Paired FDP UNDEFINED (not zero): 0 discoveries. Not a pass.")
    return result


def validate_entrapment_fdr(
    entrapment_hits: int,
    target_hits: int,
    target_db_size: int,
    entrapment_db_size: int,
    expected_fdr: float = 0.01,
) -> dict:
    """
    Validate FDR calibration using entrapment analysis.

    Uses the combined method (Eq 1, Keich & Noble 2025) which provides
    a valid upper bound on FDP. If this upper bound is at or below the
    nominal FDR threshold, the FDR control is validated.

    Parameters
    ----------
    entrapment_hits : int
        Number of PSMs matching entrapment proteins.
    target_hits : int
        Number of PSMs matching target proteins.
    target_db_size : int
        Number of proteins in the target database.
    entrapment_db_size : int
        Number of proteins in the entrapment database.
    expected_fdr : float
        Nominal FDR threshold to validate against. Default: 0.01.

    Returns
    -------
    dict
        Includes ``fdp_combined`` (the valid upper bound), ``passed``
        (whether FDR is well-calibrated), and all fields from
        ``estimate_fdp()``.
    """
    if not 0 < expected_fdr <= 1:
        raise ValueError("expected_fdr must lie in (0, 1]")
    result = estimate_fdp(
        entrapment_hits,
        target_hits,
        target_db_size,
        entrapment_db_size,
        method="combined",
    )
    result["expected_fdr"] = expected_fdr
    result["interpretation"] = (
        "Comparison of an entrapment FDP point estimate with a nominal threshold. "
        "PASS does not establish a confidence bound or prove FDR control. Validity "
        "depends on the null construction, searched hypothesis ratio and counting unit."
    )

    if not result["evaluable"]:
        # No discoveries, or no entrapment in the searched database. There was
        # no test, so there can be no pass. Returning True here would be the
        # most favourable verdict standing in for a measurement never made.
        result["passed"] = None
        result["verdict"] = "UNDEFINED"
        result["warning"] = (
            "FDR calibration is UNDEFINED, not validated: "
            f"{result['n_total']} discoveries against an entrapment database of "
            f"{entrapment_db_size}. Nothing could be detected, so this is not a "
            "pass. Report it as unevaluable."
        )
        logger.warning(result["warning"])
        return result

    result["passed"] = result["fdp_combined"] <= expected_fdr
    result["verdict"] = "PASS" if result["passed"] else "FAIL"

    if not result["passed"]:
        result["warning"] = (
            f"FDP upper bound ({result['fdp_combined']:.2%}) exceeds "
            f"nominal FDR ({expected_fdr:.2%}). FDR control may be "
            f"poorly calibrated at this database size."
        )
        logger.warning(result["warning"])

    return result


def prepare_entrapment(
    config: EntrapmentConfig,
    source_proteins: dict[str, str] | None = None,
) -> dict[str, str]:
    """
    Prepare entrapment proteins based on configuration.

    Dispatcher that calls the appropriate generation method based on
    ``config.method``.

    Parameters
    ----------
    config : EntrapmentConfig
        Entrapment configuration.
    source_proteins : dict[str, str] or None
        Source proteins for shuffled method. Required for ``'shuffled'``
        method. For ``'phylogenetic'`` and ``'user'``, the FASTA is loaded
        from ``config.source``.

    Returns
    -------
    dict[str, str]
        Entrapment proteins ready to merge into the search database.

    Raises
    ------
    ValueError
        If method is not recognized or required inputs are missing.
    """
    if not config.enabled:
        return {}

    if config.method == "shuffled":
        if source_proteins is None:
            raise ValueError("source_proteins required for 'shuffled' entrapment method")
        # The seed is validated at config construction (nonzero); it must reach
        # the generator, or every config-driven shuffle silently uses 42.
        return generate_shuffled_entrapment(
            source_proteins,
            n_proteins=config.n_proteins,
            seed=config.seed,
        )

    elif config.method in ("phylogenetic", "user"):
        if not config.source:
            raise ValueError(f"config.source (FASTA path) required for '{config.method}' method")
        proteins = load_entrapment_fasta(config.source)
        # Limit to n_proteins if specified
        if config.n_proteins and len(proteins) > config.n_proteins:
            ids = sorted(proteins.keys())[: config.n_proteins]
            proteins = {pid: proteins[pid] for pid in ids}
        return proteins

    else:
        raise ValueError(
            f"Unknown entrapment method '{config.method}'. "
            "Must be one of: 'phylogenetic', 'shuffled', 'user'"
        )
