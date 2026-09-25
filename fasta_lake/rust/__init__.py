"""Rust binary integration for performance-critical stages."""

from fasta_lake.rust.binaries import find_binary, run_fasta_extractor, run_lake_builder

__all__ = ["find_binary", "run_fasta_extractor", "run_lake_builder"]
