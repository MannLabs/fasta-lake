"""De novo protein-database curation and study analysis for metaproteomics.

The supported manifest workflow uses Rust for reference matching and compact
database selection, Sage for searching, and Python for orchestration, study
grouping and quantitative evidence. Matched metaG/metaT can add exact proteins
after selection. Sum/directLFQ quantities feed automatic AlphaPeptTools QC/PCA.

See ``docs/reference/function-map.md`` for the public stage functions and their tests.
Optional annotation and model dependencies are loaded by their own entry points.
Older JSON configuration and Python-inference helpers remain separate interfaces.
"""

__version__ = "1.1.0rc2"

from fasta_lake.complexity import estimate_complexity
from fasta_lake.config import FastaLakeConfig, load_config
from fasta_lake.entrapment import (
    generate_noble_paired,
    generate_shuffled_entrapment,
    load_entrapment_fasta,
    prepare_entrapment,
    validate_entrapment_fdr,
)
from fasta_lake.helpers import (
    clean_sequence,
    collect_peptides_with_scores,
    get_db_type,
    get_species_id,
    load_fasta_proteins,
    map_peptides_to_proteins,
    normalize_il,
)
from fasta_lake.quality import (
    compute_sequence_entropy,
    filter_sequences_by_quality,
    has_poly_run,
    is_fragment,
)

__all__ = [
    "__version__",
    # Helpers (verbatim from e2e-tested scripts)
    "normalize_il",
    "clean_sequence",
    "get_db_type",
    "get_species_id",
    "load_fasta_proteins",
    "collect_peptides_with_scores",
    "map_peptides_to_proteins",
    # Config
    "FastaLakeConfig",
    "load_config",
    # Quality filtering (v3)
    "has_poly_run",
    "compute_sequence_entropy",
    "is_fragment",
    "filter_sequences_by_quality",
    # Complexity estimation (v3)
    "estimate_complexity",
    # Entrapment (v3)
    "generate_shuffled_entrapment",
    "load_entrapment_fasta",
    "generate_noble_paired",
    "prepare_entrapment",
    "validate_entrapment_fdr",
]
