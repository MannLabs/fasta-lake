"""
Rust binary discovery and execution.

Locates compiled Rust binaries (``fasta_extractor``, ``lake_builder``)
and provides Python wrappers for invoking them with proper arguments.
"""

from __future__ import annotations

import json
import logging
import math
import os
import shutil
import subprocess
from pathlib import Path

logger = logging.getLogger(__name__)

# Search order for Rust binaries relative to the package
_SEARCH_PATHS = [
    "rust/fasta_extractor_v2/target/release",
    "rust/parsimony_engine/target/release",
    "rust/aggregation_engine/target/release",
    "rust/fasta_export/target/release",
]


def find_binary(
    name: str,
    search_dir: str | Path | None = None,
    explicit_path: str | Path | None = None,
) -> Path:
    """
    Locate a compiled Rust binary.

    Parameters
    ----------
    name : str
        Binary name (e.g., ``'fasta_extractor'`` or ``'lake_builder'``).
    search_dir : str or Path, optional
        Base directory to search from. Defaults to the FastaLake repo root.
    explicit_path : str or Path, optional
        If provided, use this path directly (from config).

    Returns
    -------
    Path
        Absolute path to the binary.

    Raises
    ------
    FileNotFoundError
        If the binary cannot be found.
    """
    if explicit_path:
        p = Path(explicit_path)
        if p.is_file() and os.access(p, os.X_OK):
            return p.resolve()
        raise FileNotFoundError(f"Explicit binary path is missing or not executable: {p}")

    configured_dir = os.environ.get("FASTALAKE_BIN_DIR")
    if configured_dir:
        return find_binary(name, explicit_path=Path(configured_dir) / name)

    if search_dir is None:
        # Default: repo root is 2 levels up from this file
        search_dir = Path(__file__).resolve().parent.parent.parent

    search_dir = Path(search_dir)

    for rel_path in _SEARCH_PATHS:
        candidate = search_dir / rel_path / name
        if candidate.is_file() and os.access(candidate, os.X_OK):
            logger.info("Found %s at %s", name, candidate)
            return candidate.resolve()

    on_path = shutil.which(name)
    if on_path:
        return Path(on_path).resolve()

    raise FileNotFoundError(
        f"Rust binary '{name}' not found. "
        f"Searched in {search_dir}. "
        f"Build with: cd rust/fasta_extractor_v2 && cargo build --release"
    )


def run_fasta_extractor(
    database: str | Path,
    predictions: str | Path,
    output: str | Path,
    source_tag: str = "REF",
    normalize_il: bool = True,
    min_length: int = 9,
    max_length: int = 50,
    top_percent: float | None = None,
    output_hits: str | Path | None = None,
    threads: int | None = None,
    binary_path: str | Path | None = None,
) -> dict:
    """
    Run the Rust ``fasta_extractor`` binary.

    Parameters
    ----------
    database : str or Path
        Reference FASTA database.
    predictions : str or Path
        Predictions CSV file or directory.
    output : str or Path
        Output FASTA path.
    source_tag : str
        Source tag for output headers.
    normalize_il : bool
        Enable I/L normalization.
    min_length : int
        Minimum peptide length.
    max_length : int
        Maximum peptide length.
    top_percent : float, optional
        Keep top N% peptides per sample by score.
    output_hits : str or Path, optional
        Path for peptide-protein hits TSV.
    threads : int, optional
        Number of threads.
    binary_path : str or Path, optional
        Explicit path to binary.

    Returns
    -------
    dict
        Stats from the ``.stats.json`` sidecar file.
    """
    binary = find_binary("fasta_extractor", explicit_path=binary_path)
    if min_length < 1 or max_length < min_length:
        raise ValueError("Peptide lengths must satisfy 1 <= min_length <= max_length")
    if top_percent is not None and not (math.isfinite(top_percent) and 0 < top_percent <= 1):
        raise ValueError("top_percent is a fraction in (0, 1]; use 0.30 for 30%")
    if threads is not None and threads < 1:
        raise ValueError("threads must be at least 1")
    for input_path in (database, predictions):
        if not Path(input_path).exists():
            raise FileNotFoundError(f"Input does not exist: {input_path}")
    for output_path in (output, output_hits, Path(output).with_suffix(".stats.json")):
        if output_path and Path(output_path).exists():
            raise FileExistsError(f"Refusing to overwrite existing output: {output_path}")
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    if output_hits:
        Path(output_hits).parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        str(binary),
        "-d",
        str(database),
        "-p",
        str(predictions),
        "-o",
        str(output),
        "-s",
        source_tag,
        "--min-length",
        str(min_length),
        "--max-length",
        str(max_length),
    ]

    # The fasta_extractor_v2 binary defines --normalize-il with
    # `action = clap::ArgAction::Set`, so it REQUIRES an explicit value
    # (`--normalize-il true`). A bare `--normalize-il` is rejected by clap
    # with exit status 2. Always pass the value form. (v1 has no such flag;
    # binaries.py only ever targets v2 — see _SEARCH_PATHS ordering.)
    cmd.extend(["--normalize-il", "true" if normalize_il else "false"])
    if top_percent is not None:
        cmd.extend(["--top-percent", str(top_percent)])
    if output_hits:
        cmd.extend(["--output-hits", str(output_hits)])
    if threads:
        cmd.extend(["-t", str(threads)])

    logger.info("Running: %s", " ".join(cmd))

    result = subprocess.run(cmd, capture_output=True, text=True, check=True)

    if result.stdout:
        for line in result.stdout.strip().split("\n"):
            logger.info("  [rust] %s", line)

    # Read stats sidecar
    stats_path = Path(output).with_suffix(".stats.json")
    if stats_path.exists():
        with open(stats_path) as f:
            return json.load(f)

    return {}


def run_lake_builder(
    sources: list[tuple[str, int, str]],
    output: str | Path,
    threads: int | None = None,
    binary_path: str | Path | None = None,
) -> None:
    """
    Run the Rust ``lake_builder`` binary.

    Parameters
    ----------
    sources : list[tuple[str, int, str]]
        List of (tag, priority, path) tuples for each source database.
    output : str or Path
        Output lake FASTA path.
    threads : int, optional
        Number of threads.
    binary_path : str or Path, optional
        Explicit path to binary.
    """
    binary = find_binary("lake_builder", explicit_path=binary_path)
    if not sources:
        raise ValueError("At least one source FASTA is required")
    if threads is not None and threads < 1:
        raise ValueError("threads must be at least 1")
    if Path(output).exists():
        raise FileExistsError(f"Refusing to overwrite existing output: {output}")
    for tag, _priority, path in sources:
        if not tag or ":" in tag:
            raise ValueError("Source tags must be nonempty and cannot contain ':'")
        if not Path(path).is_file():
            raise FileNotFoundError(f"Source FASTA does not exist: {path}")
    Path(output).parent.mkdir(parents=True, exist_ok=True)

    cmd = [str(binary), "-o", str(output)]

    for tag, priority, path in sources:
        cmd.extend(["-s", f"{tag}:{priority}:{path}"])

    if threads:
        cmd.extend(["-t", str(threads)])

    logger.info("Running: %s", " ".join(cmd))
    subprocess.run(cmd, capture_output=True, text=True, check=True)
