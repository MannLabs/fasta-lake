"""
FastaLake parameter presets (v5).

Legacy entrapment observations refer to the historical 592-sample MicrobPredict sweep;
they do not calibrate the current whole workflow.
See DESIGN_V5.md for the empirical justification.

Use via CLI:
    fasta-lake infer --preset stringent ...
    fasta-lake export --preset default ...

Or programmatically:
    from fasta_lake.presets import apply_preset
    cfg = apply_preset(cfg, "stringent")
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Preset:
    """Immutable preset definition — a named bundle of parameter values."""

    name: str
    description: str

    # Stage 1-3 (de novo peptide filter)
    denovo_min_length: int
    denovo_max_length: int = 50
    min_score: float = 0.0
    top_percent: float = 0.30

    # Stage 4 (SAGE search)
    search_min_length: int = 9
    search_max_length: int = 50
    search_missed_cleavages: int = 2
    search_max_variable_mods: int = 2
    search_precursor_ppm: float = 10.0
    search_fdr: float = 0.01

    # Pipeline features
    kmer_rescue: bool = False
    cluster_enabled: bool = True
    cluster_identity: float = 0.97
    infer_strategy: str = "species_budget"
    infer_min_peptides: int = 2

    def as_dict(self) -> dict:
        """Return a dataclass dictionary of this preset's settings."""
        from dataclasses import asdict

        return asdict(self)


# ── Preset definitions (validated parameters) ──────────────────────────────

STRINGENT = Preset(
    name="stringent",
    description=(
        "Conservative historical parameter bundle. Historical shuffled-entrapment "
        "FDP <0.1% at q≤0.001; not whole-workflow calibration."
    ),
    denovo_min_length=9,
    search_min_length=10,
    search_missed_cleavages=1,
    search_max_variable_mods=1,
    search_precursor_ppm=5.0,
    search_fdr=0.001,
    kmer_rescue=False,
    cluster_identity=0.95,
    infer_strategy="razor",
    infer_min_peptides=2,
)

DEFAULT = Preset(
    name="default",
    description=(
        "Historical default parameter bundle. Historical shuffled-entrapment "
        "FDP ~1.4% at q≤0.01; not whole-workflow calibration."
    ),
    denovo_min_length=9,
    search_min_length=9,
    search_missed_cleavages=2,
    search_max_variable_mods=2,
    search_precursor_ppm=10.0,
    search_fdr=0.01,
    kmer_rescue=False,
    cluster_identity=0.97,
    infer_strategy="species_budget",
    infer_min_peptides=2,
)

LENIENT = Preset(
    name="lenient",
    description=(
        "Exploratory / high-sensitivity mode. More permissive — acceptable only for "
        "single-species or hypothesis-generation. FDP NOT guaranteed at 1%."
    ),
    denovo_min_length=7,
    search_min_length=7,
    search_missed_cleavages=2,
    search_max_variable_mods=2,
    search_precursor_ppm=20.0,
    search_fdr=0.05,
    kmer_rescue=True,  # OK in exploratory mode
    cluster_identity=0.99,  # retain strain-level variation
    infer_strategy="uniform",
    infer_min_peptides=1,
)


PROTEOMICS = Preset(
    name="proteomics",
    description=(
        "Single-organism (host/human) proteomics — NOT metaproteomics. "
        "Personal-genome-anchored extension over a compact canonical reference; "
        "no population pooling, no clustering (reference is already small). "
        "Value is proteoform specificity (variants/IG/HLA), not protein-count depth. "
        "See fasta_lake/proteomics/ARCHITECTURE.md."
    ),
    denovo_min_length=7,  # DIA de novo (DiaNovo/AlphaNovo) peptides run short
    search_min_length=7,
    search_missed_cleavages=2,
    search_max_variable_mods=2,
    search_precursor_ppm=10.0,  # DIA-NN default regime
    search_fdr=0.01,
    kmer_rescue=False,  # never in proteomics mode — it adds boundary IDs
    cluster_enabled=False,  # compact human reference — clustering is unnecessary/harmful
    cluster_identity=1.0,
    infer_strategy="razor",  # canonical Tier-1 backbone; variant class handled by cascade
    infer_min_peptides=1,  # DIA precursor-level
)


PRESETS: dict[str, Preset] = {
    "stringent": STRINGENT,
    "default": DEFAULT,
    "lenient": LENIENT,
    "proteomics": PROTEOMICS,
}


def get_preset(name: str) -> Preset:
    """Return the named preset or raise ValueError listing the available names."""
    if name not in PRESETS:
        choices = ", ".join(sorted(PRESETS.keys()))
        raise ValueError(f"Unknown preset '{name}'. Choose from: {choices}")
    return PRESETS[name]


def apply_preset(cfg, preset_name: str):
    """Apply preset fields to an existing RunConfig-like object."""
    preset = get_preset(preset_name)
    for field_name, value in preset.as_dict().items():
        if field_name in ("name", "description"):
            continue
        if hasattr(cfg, field_name):
            setattr(cfg, field_name, value)
    return cfg


def preset_sage_config(base_config: dict, preset_name: str) -> dict:
    """Apply preset to a SAGE JSON config dict. Returns modified copy."""
    import json

    preset = get_preset(preset_name)
    cfg = json.loads(json.dumps(base_config))  # deep copy
    cfg["database"]["enzyme"]["min_len"] = preset.search_min_length
    cfg["database"]["enzyme"]["max_len"] = preset.search_max_length
    cfg["database"]["enzyme"]["missed_cleavages"] = preset.search_missed_cleavages
    cfg["database"]["max_variable_mods"] = preset.search_max_variable_mods
    cfg["precursor_tol"] = {"ppm": [-preset.search_precursor_ppm, preset.search_precursor_ppm]}
    return cfg


def describe_preset(name: str) -> str:
    """Human-readable description of what a preset does."""
    p = get_preset(name)
    lines = [
        f"PRESET: {p.name.upper()}",
        p.description,
        "",
        "  Stage 1-3 (de novo evidence filter):",
        f"    min_length      : {p.denovo_min_length}",
        f"    max_length      : {p.denovo_max_length}",
        f"    top_percent     : {p.top_percent}",
        f"    min_score       : {p.min_score}",
        "",
        "  Stage 4 (SAGE search):",
        f"    min_length      : {p.search_min_length}",
        f"    max_length      : {p.search_max_length}",
        f"    missed_cleavages: {p.search_missed_cleavages}",
        f"    max_var_mods    : {p.search_max_variable_mods}",
        f"    precursor_ppm   : ±{p.search_precursor_ppm}",
        f"    q_threshold     : {p.search_fdr}",
        "",
        "  Pipeline:",
        f"    clustering      : {p.cluster_identity * 100:.0f}% identity",
        f"    inference       : {p.infer_strategy} (min {p.infer_min_peptides} peptides)",
        f"    kmer_rescue     : {p.kmer_rescue}",
    ]
    return "\n".join(lines)
