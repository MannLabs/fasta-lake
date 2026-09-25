"""
Inference runner: dispatches strategies and handles file I/O.

This module ties together the loading, mapping, strategy execution,
and output writing into a single ``run_inference()`` function.

When ``kmer_rescue_db`` is provided, the runner identifies peptides that
failed exact matching against the Stage 2 FASTA and rescues them via
mass k-mer anchoring against the given database. Rescued proteins are
injected into the protein set before inference runs, closing the gap
where per-sample exact extraction misses near-match proteins.
"""

from __future__ import annotations

import csv
import json
import logging
import subprocess
import tempfile
import time
from pathlib import Path

from fasta_lake.helpers import (
    collect_peptides_with_scores,
    load_fasta_proteins,
    map_peptides_to_proteins,
)
from fasta_lake.inference.base import InferenceInput, InferenceResult
from fasta_lake.inference.strategies import STRATEGIES
from fasta_lake.provenance import write_provenance

logger = logging.getLogger(__name__)

# Default path to mass_kmer_anchor binary (override via KMER_BINARY env var)
_KMER_BINARY_DEFAULT = (
    Path(__file__).resolve().parents[2]
    / "rust"
    / "fasta_extractor_v2"
    / "target"
    / "release"
    / "mass_kmer_anchor"
)


def _rescue_failed_peptides(
    peptide_scores: dict[str, float],
    pep2prot: dict[str, set[str]],
    rescue_db: Path,
    min_votes: int = 5,
    threads: int = 4,
) -> dict[str, str]:
    """
    Run mass k-mer anchoring on peptides that failed exact matching.

    Parameters
    ----------
    peptide_scores : dict[str, float]
        All peptides with their scores (I/L-normalized).
    pep2prot : dict[str, set[str]]
        Peptides that DID match proteins (from exact mapping).
    rescue_db : Path
        FASTA database to rescue against (typically the evidence lake).
    min_votes : int
        Minimum k-mer votes for a rescue hit.
    threads : int
        Threads for mass_kmer_anchor.

    Returns
    -------
    dict[str, str]
        Mapping of rescued protein_id -> sequence, loaded from rescue_db.
    """
    import os

    kmer_binary = os.environ.get("KMER_BINARY", str(_KMER_BINARY_DEFAULT))
    if not Path(kmer_binary).exists():
        # Rescue was requested. Skipping it would change the database with only a
        # log line to show for it, so a missing binary is a hard failure.
        raise FileNotFoundError(
            f"mass_kmer_anchor not found at {kmer_binary}; k-mer rescue was requested. "
            "Build rust/fasta_extractor_v2 or set KMER_BINARY, or run without "
            "--kmer-rescue-db."
        )

    # Identify failed peptides
    failed = {pep: score for pep, score in peptide_scores.items() if pep not in pep2prot}
    if not failed:
        logger.info("K-mer rescue: no failed peptides — nothing to rescue")
        return {}

    logger.info("K-mer rescue: %d failed peptides (of %d total)", len(failed), len(peptide_scores))

    with tempfile.TemporaryDirectory() as tmpdir:
        # Write failed peptides CSV
        failed_csv = Path(tmpdir) / "failed_peptides.csv"
        with open(failed_csv, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["peptide_prediction_detokenized_unmodified", "score"])
            for pep, score in sorted(failed.items(), key=lambda x: -x[1]):
                writer.writerow([pep, f"{score:.6f}"])

        # Run mass_kmer_anchor
        rescued_tsv = Path(tmpdir) / "rescued.tsv"
        cmd = [
            kmer_binary,
            "--database",
            str(rescue_db),
            "--predictions",
            str(failed_csv),
            "--output",
            str(rescued_tsv),
            "--min-votes",
            str(min_votes),
            "--skip-exact",
            "-t",
            str(threads),
        ]
        logger.info("K-mer rescue: running %s", " ".join(cmd[-6:]))
        # A failed or timed-out rescue must not degrade into a smaller database
        # with exit 0: the output would differ from a run on an idle machine.
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        except subprocess.TimeoutExpired as e:
            raise RuntimeError(
                f"mass_kmer_anchor timed out after {e.timeout:.0f} s on {len(failed)} "
                "failed peptides; the rescue result is undefined, so the run stops "
                "rather than continuing with a partial database."
            ) from e
        if result.returncode != 0:
            raise RuntimeError(
                f"mass_kmer_anchor failed (exit {result.returncode}): "
                f"{result.stderr.strip()[:2000]}"
            )

        # Read rescued protein IDs. mass_kmer_anchor exits 0 without writing the
        # file when it has nothing to report; that is a legitimate zero.
        rescued_ids: set[str] = set()
        if rescued_tsv.exists():
            with open(rescued_tsv) as f:
                reader = csv.DictReader(f, delimiter="\t")
                for row in reader:
                    pid = row.get("protein_id", "").strip()
                    if pid:
                        rescued_ids.add(pid)

    if not rescued_ids:
        logger.info("K-mer rescue: no proteins rescued")
        return {}

    # Load rescued protein sequences from the rescue database
    rescued_proteins: dict[str, str] = {}
    current_id: str | None = None
    current_seq: list[str] = []
    writing = False

    with open(rescue_db) as f:
        for line in f:
            if line.startswith(">"):
                if writing and current_id:
                    rescued_proteins[current_id] = "".join(current_seq)
                current_id = line[1:].strip().split()[0]
                current_seq = []
                writing = current_id in rescued_ids
            elif writing:
                current_seq.append(line.strip())
        if writing and current_id:
            rescued_proteins[current_id] = "".join(current_seq)

    logger.info(
        "K-mer rescue: %d proteins rescued from %d failed peptides",
        len(rescued_proteins),
        len(failed),
    )
    return rescued_proteins


def run_inference(
    sample_id: str,
    stage2_fasta: str | Path,
    alphanovo_csv: str | Path,
    output_dir: str | Path,
    strategy: str = "species_budget",
    min_peptides: int = 2,
    min_length: int = 9,  # v5: raised from 8 → 9. Matches SAGE min_len=9 default.
    max_length: int = 50,
    min_score: float = 0.0,
    taxonomy_map_path: str | Path | None = None,
    taxonomy_pruning: bool = False,
    kmer_rescue_db: str | Path | None = None,
    kmer_min_votes: int = 5,
    kmer_threads: int = 4,
) -> dict:
    """
    Run protein inference on a single sample.

    Parameters
    ----------
    sample_id : str
        Sample identifier (used in output file names).
    stage2_fasta : str or Path
        Path to Stage 2 sample-specific FASTA file.
    alphanovo_csv : str or Path
        Path to AlphaNovo predictions CSV file.
    output_dir : str or Path
        Directory for output files (FASTA + stats JSON).
    strategy : str
        Inference strategy name. One of:
        ``'species_budget'``, ``'razor'``, ``'uniform'``, ``'uniform_2pep'``.
    min_peptides : int
        Minimum peptides for uniform_2pep strategy (default: 2).
    min_length : int
        Minimum peptide length filter.
    max_length : int
        Maximum peptide length filter.
    min_score : float
        Minimum AlphaNovo confidence score (default: 0.0).
    kmer_rescue_db : str or Path, optional
        Path to FASTA database for k-mer rescue of failed peptides.
        When set, peptides that don't exact-match any Stage 2 protein
        are rescued via mass k-mer anchoring against this database.
        Rescued proteins are injected before inference runs.
    kmer_min_votes : int
        Minimum k-mer votes for rescue (default: 5).
    kmer_threads : int
        Threads for mass_kmer_anchor (default: 4).

    Returns
    -------
    dict
        Statistics dictionary with inference results.

    Raises
    ------
    ValueError
        If strategy name is not recognized.
    FileNotFoundError
        If input files do not exist.
    """
    stage2_fasta = Path(stage2_fasta)
    alphanovo_csv = Path(alphanovo_csv)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if strategy not in STRATEGIES:
        raise ValueError(f"Unknown strategy '{strategy}'. Available: {', '.join(STRATEGIES)}")

    if not stage2_fasta.exists():
        raise FileNotFoundError(f"Stage 2 FASTA not found: {stage2_fasta}")
    if not alphanovo_csv.exists():
        raise FileNotFoundError(f"AlphaNovo CSV not found: {alphanovo_csv}")

    start_time = time.time()

    # Load data
    peptide_scores = collect_peptides_with_scores(
        alphanovo_csv, min_length=min_length, max_length=max_length
    )

    # Apply min_score filter
    if min_score > 0.0:
        before = len(peptide_scores)
        peptide_scores = {pep: score for pep, score in peptide_scores.items() if score >= min_score}
        logger.info(
            "Score filter (>= %.2f): %d -> %d peptides",
            min_score,
            before,
            len(peptide_scores),
        )

    stage2_proteins = load_fasta_proteins(stage2_fasta)

    # Build initial mappings (exact match)
    pep2prot, prot2pep = map_peptides_to_proteins(stage2_proteins, peptide_scores)

    # K-mer rescue: inject rescued proteins for failed peptides
    kmer_stats: dict = {}
    kmer_rescued_ids: set[str] = set()  # v5.1: track for provenance
    if kmer_rescue_db:
        kmer_rescue_db = Path(kmer_rescue_db)
        if not kmer_rescue_db.exists():
            raise FileNotFoundError(f"K-mer rescue DB not found: {kmer_rescue_db}")
        n_before = len(stage2_proteins)
        rescued = _rescue_failed_peptides(
            peptide_scores,
            pep2prot,
            kmer_rescue_db,
            min_votes=kmer_min_votes,
            threads=kmer_threads,
        )
        # Add rescued proteins (skip duplicates)
        new_count = 0
        for pid, seq in rescued.items():
            if pid not in stage2_proteins:
                stage2_proteins[pid] = seq
                kmer_rescued_ids.add(pid)
                new_count += 1
        if new_count:
            logger.info(
                "K-mer rescue: injected %d new proteins (%d total now, was %d)",
                new_count,
                len(stage2_proteins),
                n_before,
            )
            # Re-map peptides with expanded protein set
            pep2prot, prot2pep = map_peptides_to_proteins(stage2_proteins, peptide_scores)

        # Recorded whenever rescue was requested, so the stats file says that it
        # ran and what it found, including a zero.
        kmer_stats = {
            "kmer_rescue_db": str(kmer_rescue_db),
            "kmer_proteins_rescued": len(rescued),
            "kmer_proteins_new": new_count,
            "kmer_proteins_already_present": len(rescued) - new_count,
            "kmer_min_votes": kmer_min_votes,
        }

    # Load taxonomy map if provided
    tax_map: dict[str, str] = {}
    if taxonomy_map_path:
        taxonomy_map_path = Path(taxonomy_map_path)
        if taxonomy_map_path.exists():
            with open(taxonomy_map_path) as f:
                for line in f:
                    parts = line.strip().split("\t")
                    if len(parts) >= 2 and not parts[0].startswith("#"):
                        tax_map[parts[0]] = parts[1]
            logger.info("Loaded taxonomy map: %d entries", len(tax_map))

    # Create inference input
    inp = InferenceInput(
        pep2prot=pep2prot,
        prot2pep=prot2pep,
        proteins=stage2_proteins,
        peptide_scores=peptide_scores,
        taxonomy_map=tax_map,
        taxonomy_pruning=taxonomy_pruning,
    )

    # Run strategy
    strategy_fn = STRATEGIES[strategy]
    if strategy == "uniform_2pep":
        result: InferenceResult = strategy_fn(inp, min_peptides=min_peptides)
    elif strategy == "info_score":
        result = strategy_fn(inp, min_info=min_score if min_score > 0 else 0.1)
    elif strategy == "bayesian":
        result = strategy_fn(inp, min_probability=min_score if min_score > 0 else 0.01)
    else:
        result = strategy_fn(inp)

    elapsed = time.time() - start_time

    # Determine output file name suffix
    suffix_map = {
        "species_budget": "species_budget",
        "razor": "razor",
        "uniform": "uniform",
        "uniform_2pep": "uniform2pep",
        "info_score": "info_score",
        "bayesian": "bayesian",
    }
    suffix = suffix_map.get(strategy, strategy)

    # Write output FASTA
    output_fasta = output_dir / f"{sample_id}_{suffix}.fasta"
    with open(output_fasta, "w") as f:
        for protein_id in sorted(result.selected_proteins):
            if protein_id in stage2_proteins:
                f.write(f">{protein_id}\n")
                seq = stage2_proteins[protein_id]
                for i in range(0, len(seq), 60):
                    f.write(seq[i : i + 60] + "\n")

    file_size_mb = output_fasta.stat().st_size / (1024**2)

    # Build complete stats
    stats = {
        "sample_id": sample_id,
        **result.stats,
        **kmer_stats,
        "runtime_seconds": elapsed,
        "stage2_proteins": len(stage2_proteins),
        "reduction_factor": (
            len(stage2_proteins) / len(result.selected_proteins) if result.selected_proteins else 0
        ),
        "file_size_mb": file_size_mb,
    }

    # Write stats JSON
    stats_json = output_dir / f"{sample_id}_{suffix}_stats.json"
    with open(stats_json, "w") as f:
        json.dump(stats, f, indent=2)

    # v5.1: write per-protein provenance JSONL
    provenance_path = output_dir / f"{sample_id}_{suffix}_provenance.jsonl"
    n_prov = write_provenance(
        out_path=provenance_path,
        selected_proteins=result.selected_proteins,
        prot2pep=prot2pep,
        peptide_scores=peptide_scores,
        strategy=strategy,
        kmer_rescued_proteins=kmer_rescued_ids,
    )
    stats["provenance_records"] = n_prov

    logger.info(
        "Inference complete: %s proteins selected (%.1fx reduction) in %.1fs",
        f"{len(result.selected_proteins):,}",
        stats["reduction_factor"],
        elapsed,
    )

    return stats
