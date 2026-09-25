"""
Protein inference engine for FASTA Lake.

Provides 4 inference strategies consolidated from the e2e-tested
standalone scripts:

- **species_budget**: Cross-species aware; budgets shared peptides by species
- **razor**: MaxQuant-style; shared peptides → protein with most unique peptides
- **uniform**: Baseline; keep all proteins with any evidence
- **uniform_2pep**: Conservative; require >=2 peptides per protein
"""

from fasta_lake.inference.base import InferenceInput, InferenceResult
from fasta_lake.inference.runner import run_inference
from fasta_lake.inference.strategies import STRATEGIES

__all__ = [
    "InferenceInput",
    "InferenceResult",
    "STRATEGIES",
    "run_inference",
]
