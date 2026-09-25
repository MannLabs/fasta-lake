"""
Cross-language hash equality: Python ``fasta_lake.hashing.sequence_hash``
vs Rust ``rust/fasta_extractor_v2/src/bin/hash_fasta.rs``.

This is the **load-bearing correctness test** for the v5.3 hashing contract.
If it fails, a stage of the pipeline disagrees about which sequences are
identical, and dedup is silently broken — a category of bug that otherwise
goes unnoticed for years.

The test is parameterised over a variety of real FASTA fixtures (canonical
human SwissProt, splice isoforms, the UP000005640 reference proteome, a
personal-variant hybrid, and a tiny synthetic fixture). For every sequence
in every file, the Rust binary's hex digest must exactly equal Python's.

Skipped when:
- The Rust binary is not built (build with
  ``cargo build --release --bin hash_fasta``).
- The fixture FASTAs are not on the filesystem (they live on the MPI GPFS
  pool and are not redistributed with the repo).
"""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

import pytest

from fasta_lake.hashing import sequence_hash

REPO_ROOT = Path(__file__).resolve().parent.parent
RUST_BIN = REPO_ROOT / "rust" / "fasta_extractor_v2" / "target" / "release" / "hash_fasta"
RUST_BIN = Path(os.environ.get("FASTALAKE_BIN_DIR", str(RUST_BIN.parent))) / RUST_BIN.name

ISOFORM_DIR = Path(
    os.environ.get(
        "FASTALAKE_TEST_HASHING_CROSSLANG_1", "tests/optional_data/test_hashing_crosslang_1"
    )
)

FIXTURES: list[tuple[str, Path]] = [
    ("tiny", REPO_ROOT / "tests" / "data" / "tiny.fasta"),
    ("swissprot_human", ISOFORM_DIR / "uniprot" / "homo_sapiens_swissprot.fasta"),
    ("alt_splicing_isoforms", ISOFORM_DIR / "uniprot" / "alternative_splicing_isoforms.fasta"),
    ("UP000005640_canonical", ISOFORM_DIR / "uniprot" / "UP000005640_canonical.fasta"),
    (
        "personal_hybrid",
        ISOFORM_DIR / "personal_proteome" / "NBR66A34_hybrid_swissprot_personal.fasta",
    ),
]


def _read_fasta(path: Path) -> list[tuple[str, str]]:
    """Return [(header, sequence), ...]. Matches the Rust binary's parsing."""
    out: list[tuple[str, str]] = []
    with open(path) as f:
        header = None
        parts: list[str] = []
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    out.append((header, "".join(parts)))
                header = line[1:]
                parts = []
            elif line:
                parts.append(line)
        if header is not None:
            out.append((header, "".join(parts)))
    return out


def _rust_hashes(fasta: Path) -> list[tuple[str, str, int]]:
    """Run the Rust hash_fasta binary and return [(hex, accession, length), ...]."""
    proc = subprocess.run(
        [str(RUST_BIN), "-i", str(fasta)],
        capture_output=True,
        check=True,
        text=True,
    )
    rows: list[tuple[str, str, int]] = []
    for line in proc.stdout.splitlines():
        hex_, accession, length = line.split("\t")
        rows.append((hex_, accession, int(length)))
    return rows


@pytest.fixture(scope="module")
def rust_binary_available() -> bool:
    return RUST_BIN.exists()


@pytest.mark.parametrize("name,path", FIXTURES, ids=[n for n, _ in FIXTURES])
def test_python_rust_hash_equality(
    name: str,
    path: Path,
    rust_binary_available: bool,
) -> None:
    if not rust_binary_available:
        pytest.skip(
            f"Rust binary not built at {RUST_BIN}. "
            "Build with: cargo build --release --bin hash_fasta "
            "--manifest-path rust/fasta_extractor_v2/Cargo.toml"
        )
    if not path.exists():
        pytest.skip(f"Fixture not on this filesystem: {path}")

    # Python side — hash every sequence in the file.
    py_records = _read_fasta(path)
    py_rows = [(sequence_hash(seq), header.split()[0]) for header, seq in py_records]

    # Rust side — one invocation, one line per sequence.
    rust_rows = _rust_hashes(path)

    # Same row count. If not, parsing diverged.
    assert len(py_rows) == len(rust_rows), (
        f"[{name}] record count differs: python={len(py_rows)} rust={len(rust_rows)}"
    )

    # Per-record equality. Report the first mismatch clearly; the tail of
    # mismatches is usually a copy of the root cause.
    mismatches: list[str] = []
    for i, ((py_hex, py_acc), (r_hex, r_acc, r_len)) in enumerate(zip(py_rows, rust_rows)):
        if py_acc != r_acc:
            mismatches.append(f"  row {i}: accession differs — python={py_acc} rust={r_acc}")
        if py_hex != r_hex:
            mismatches.append(f"  row {i} ({py_acc}): hash differs — python={py_hex} rust={r_hex}")
        if len(mismatches) >= 5:
            mismatches.append("  (further mismatches suppressed)")
            break

    if mismatches:
        pytest.fail(
            f"[{name}] Python/Rust hash divergence on {path.name}:\n" + "\n".join(mismatches)
        )


def test_rust_binary_exists(rust_binary_available: bool) -> None:
    """Fast fail if someone forgot to build. Not a correctness gate — every
    per-file test already skips on missing binary — but it surfaces the
    'you forgot to cargo build' case with a clear message."""
    if not rust_binary_available:
        pytest.skip(f"Rust binary not built at {RUST_BIN}")
    assert RUST_BIN.is_file() and RUST_BIN.stat().st_mode & 0o111, (
        f"{RUST_BIN} exists but is not executable"
    )
