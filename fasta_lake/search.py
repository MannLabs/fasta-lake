"""Sage DDA settings: cleavage and precursor mass window are independent choices.

Upstream search engine: https://github.com/lazear/sage; Lazear (2023),
https://doi.org/10.1021/acs.jproteome.3c00486. FastaLake writes explicit search
settings; Sage performs spectrum matching, scoring and its native confidence
estimation. Search q-values are not whole adaptive-workflow FDR calibration.
"""

DIGESTIONS = ("tryptic", "semi-tryptic", "unspecific")
SEARCH_MODES = ("closed", "open")


def _dda_config(fasta, mzml, output):
    """The explicitly versioned DDA settings used in the release acceptance runs."""
    return {
        "database": {
            "bucket_size": 32768,
            "enzyme": {
                "missed_cleavages": 2,
                "min_len": 5,
                "max_len": 50,
                "cleave_at": "KR",
                "restrict": "P",
            },
            "fragment_min_mz": 200,
            "fragment_max_mz": 2000,
            "peptide_min_mass": 500,
            "peptide_max_mass": 5000,
            "ion_kinds": ["b", "y"],
            "min_ion_index": 2,
            "static_mods": {"C": 57.0215},
            "variable_mods": {"M": [15.9949]},
            "max_variable_mods": 2,
            "decoy_tag": "rev_",
            "generate_decoys": True,
            "fasta": str(fasta),
        },
        "quant": {
            "tmt": None,
            "lfq": True,
            "lfq_settings": {"peak_scoring": "Hybrid", "integration": "Sum", "spectral_angle": 0.7},
        },
        "precursor_tol": {"ppm": [-20, 20]},
        "fragment_tol": {"ppm": [-20, 20]},
        "isotope_errors": [-1, 3],
        "deisotope": False,
        "chimera": False,
        "predict_rt": True,
        "wide_window": False,
        "min_peaks": 15,
        "max_peaks": 150,
        "min_matched_peaks": 4,
        "max_fragment_charge": 2,
        "report_psms": 1,
        "output_directory": str(output),
        "mzml_paths": [str(mzml)],
    }


def sage_config(
    fasta, mzml, output, *, digestion="tryptic", search_mode="closed", min_length=5, max_length=50
):
    """Build the complete configuration for one Sage 0.14.6 acquisition.

    Parameters
    ----------
    fasta, mzml, output : path-like
        Selected target database, one measured acquisition, and its new search
        directory. Sage generates decoys; input FASTA must contain targets only.
    digestion : {"tryptic", "semi-tryptic", "unspecific"}
        Require two, one, or no tryptic termini. Protein termini are valid
        enzymatic boundaries. Permissiveness does not establish protease activity.
    search_mode : {"closed", "open"}
        Closed uses +/-20 ppm. Open uses [-500, +100] Da, no variable
        modifications, no isotope-error expansion, and no RT prediction.
    min_length, max_length : int
        Search peptide lengths, independent of the de novo evidence filter.

    Returns
    -------
    dict
        A fresh JSON-compatible config. LFQ and target/decoy generation are
        retained. Every actual search saves these settings in its JSON.
    """
    if digestion not in DIGESTIONS:
        raise ValueError(f"Unknown digestion: {digestion}; choose from {DIGESTIONS}")
    if search_mode not in SEARCH_MODES:
        raise ValueError(f"Unknown search mode: {search_mode}; choose from {SEARCH_MODES}")
    if (
        type(min_length) is not int
        or type(max_length) is not int
        or not 1 <= min_length <= max_length
    ):
        raise ValueError("Search peptide lengths require integers: 1 <= min <= max")
    cfg = _dda_config(fasta, mzml, output)
    enzyme = cfg["database"]["enzyme"]
    enzyme.update(min_len=min_length, max_len=max_length)
    if digestion == "semi-tryptic":
        enzyme["semi_enzymatic"] = True
    elif digestion == "unspecific":
        enzyme.update(cleave_at="", restrict=None, missed_cleavages=0, semi_enzymatic=False)
    if search_mode == "open":
        cfg["precursor_tol"] = {"da": [-500, 100]}
        cfg["isotope_errors"] = [0, 0]
        cfg["predict_rt"] = False
        cfg["database"]["variable_mods"] = {}
        cfg["database"]["max_variable_mods"] = 0
    return cfg
