"""Selectable quantification of an existing, fixed study dictionary.

The original study outputs remain the source of assignments and accepted features.
A new destination contains the selected matrix and method-specific diagnostics.


Optional backend: MannLabs directLFQ, https://github.com/MannLabs/directlfq.
Ammar et al. (2023), Accurate Label-Free Quantification by directLFQ to Compare
Unlimited Numbers of Proteomes. This adapter builds and checks the fixed-group
ion input; directLFQ performs the quantification.
"""

from __future__ import annotations

import csv
import hashlib
import importlib.metadata
import json
import math
import shutil
import tempfile
from pathlib import Path

from fasta_lake.gzip_io import open_gzip_text
from fasta_lake.study import canonical_peptide

DIRECTLFQ_VERSION = "0.3.3"


def stamp(path):
    """Return basename, byte size and a streamed SHA256 for a pathlib.Path input."""
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return {"file": path.name, "sha256": digest.hexdigest(), "bytes": path.stat().st_size}


def require_directlfq():
    """Fail before running a workflow if the optional backend is unavailable."""
    try:
        version = importlib.metadata.version("directlfq")
        if version != DIRECTLFQ_VERSION:
            raise RuntimeError(f"Expected directlfq {DIRECTLFQ_VERSION}; found {version}")
        from directlfq import lfq_manager  # noqa: F401
    except ImportError as error:
        raise RuntimeError(
            "Install the optional backend: pip install 'fasta-lake[directlfq]'"
        ) from error
    return version


def _baseline(path):
    """Validate the sum matrix without requiring numerical packages for sum mode."""
    with path.open() as stream:
        reader = csv.reader(stream, delimiter="\t")
        header = next(reader, [])
        if len(header) < 2 or header[0] != "study_group" or len(set(header)) != len(header):
            raise ValueError("Invalid study matrix header or duplicate acquisition names")
        counts, groups = [0] * (len(header) - 1), set()
        for row in reader:
            if len(row) != len(header) or not row[0] or row[0] in groups:
                raise ValueError("Invalid or duplicate study matrix group")
            groups.add(row[0])
            for i, cell in enumerate(row[1:]):
                value = float(cell) if cell else 0.0
                if not math.isfinite(value) or value < 0:
                    raise ValueError("Nonfinite or negative study intensity")
                counts[i] += value > 0
    return header[1:], counts


def quantify_study(groups, output, *, method="sum", threads=1):
    """Quantify an existing study dictionary without assigning peptides again.

    Parameters
    ----------
    groups : path-like
        Completed sequence-checked grouping bundle from group_searches.
    output : path-like
        New directory for the selected-method matrix, coverage and summary.
    method : {"sum", "directlfq"}
        Sum copies the original matrix unchanged. directLFQ estimates quantities
        from assigned ions using the required backend version.
    threads : int
        Positive worker count passed to directLFQ.

    Returns
    -------
    dict
        Method, acquisition roster, input/output hashes and coverage information
        also written to summary.json.

    Raises
    ------
    FileExistsError
        If the output already exists.
    ValueError
        For invalid settings, inconsistent group/feature data or changed inputs.
    RuntimeError
        If the requested directLFQ backend is unavailable or has the wrong version.

    Notes
    -----
    The fixed peptide-to-group dictionary and acquisition order are preserved.
    Unquantified output cells remain blank; no missing values are imputed.
    """
    if method not in {"sum", "directlfq"}:
        raise ValueError("method must be sum or directlfq")
    if not isinstance(threads, int) or threads < 1:
        raise ValueError("threads must be a positive integer")
    source, out = Path(groups).resolve(), Path(output).absolute()
    if out.exists() or out.is_symlink():
        raise FileExistsError(f"Output exists: {out}; choose a new directory")
    out = out.resolve()
    if source == out or out in source.parents:
        raise ValueError("Output must not contain the source study directory")
    version = require_directlfq() if method == "directlfq" else None
    names = [
        "summary.json",
        "study_group_matrix.tsv",
        "peptide_to_study_group.tsv",
        "feature_assignments.tsv.gz",
        "acquisition_summary.tsv",
        "input_manifest.json",
    ]
    names += [
        n for n in ["peptide_evidence.tsv.gz", "peptide_ambiguity.tsv"] if (source / n).is_file()
    ]
    inputs = [stamp(source / name) for name in names]
    summary = json.loads((source / "summary.json").read_text())
    if summary.get("schema_version") != 1 or not summary.get(
        "representative_sequence_compatibility_checked"
    ):
        raise ValueError("Expected a completed, sequence-checked study grouping bundle")
    samples, counts = _baseline(source / "study_group_matrix.tsv")
    if samples != summary.get("acquisition_names"):
        raise ValueError("Study matrix acquisition roster differs from summary")
    out.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=f".{out.name}-", dir=out.parent) as directory:
        temp = Path(directory)
        if method == "sum":
            shutil.copyfile(source / "study_group_matrix.tsv", temp / "study_group_matrix.tsv")
            details = {"sample_normalization": False, "quantified_groups": counts}
        else:
            details = _directlfq(source, temp, samples, threads)
        with (temp / "coverage.tsv").open("w", newline="") as stream:
            writer = csv.writer(stream, delimiter="\t")
            writer.writerow(["sample", "sum_groups", "selected_method_groups", "groups_lost"])
            writer.writerows(
                (s, a, b, a - b)
                for s, a, b in zip(samples, counts, details.pop("quantified_groups"))
            )
        # Detect changes to source artifacts during export/estimation.
        if inputs != [stamp(source / name) for name in names]:
            raise ValueError("Study inputs changed during quantification")
        record = {
            "schema_version": 1,
            "method": method,
            "directlfq_version": version,
            "acquisition_names": samples,
            "source_study_directory": str(source),
            "matrix": "study_group_matrix.tsv",
            "inputs": inputs,
            "implementation": stamp(Path(__file__)),
            **details,
            "missing_values": "blank in selected matrix; zero in directLFQ tables; no imputation",
            "scope": "Relative study-group quantities; no new identification or group FDR. "
            "Representative labels do not establish unique protein or organism origin.",
            "outputs": [stamp(p) for p in sorted(temp.iterdir())],
        }
        (temp / "summary.json").write_text(json.dumps(record, indent=2) + "\n")
        out.mkdir(exist_ok=False)
        for path in sorted(temp.iterdir(), key=lambda p: (p.name == "summary.json", p.name)):
            path.rename(out / path.name)
    return record


def _read_frame(path, **kwargs):
    """Read a TSV with pandas while preserving identifier-like NA strings."""
    import pandas as pd

    return pd.read_csv(path, sep="\t", keep_default_na=False, **kwargs)


def _directlfq(source, temp, samples, threads):
    """Run directLFQ on ions assigned to the fixed study dictionary.

    Write inputs and backend outputs in the caller-owned temporary directory.
    Preserve the acquisition roster and validate the resulting quantities;
    this step does not infer new protein groups.
    """
    import numpy as np
    import pandas as pd
    from directlfq import lfq_manager, utils

    if {"protein", "ion"} & set(samples):
        raise ValueError("Acquisition names protein and ion are reserved by directLFQ")
    baseline = _read_frame(source / "study_group_matrix.tsv", dtype={"study_group": str})
    # Empty cells are the text "" (keep_default_na=False); replace them with the
    # text "0" before the numeric cast, so pandas never has to downcast objects.
    baseline = baseline.set_index("study_group").replace("", "0").astype(float)
    mapping = _read_frame(source / "peptide_to_study_group.tsv", dtype=str)
    if mapping.canonical_peptide.duplicated().any():
        raise ValueError("A canonical peptide maps to multiple study groups")
    mapping = mapping.set_index("canonical_peptide").study_group
    columns = [
        "sample",
        "input_line",
        "peptide",
        "canonical_peptide",
        "charge",
        "study_group",
        "intensity",
    ]
    features = _read_frame(
        source / "feature_assignments.tsv.gz",
        usecols=columns,
        dtype={**{c: str for c in columns if c != "intensity"}, "intensity": float},
    )
    if features.duplicated(["sample", "input_line"]).any():
        raise ValueError("A physical source feature was counted more than once")
    if not set(features["sample"]) <= set(samples):
        raise ValueError("Feature acquisition outside the study roster")
    if not features.charge.str.fullmatch(r"(?:[1-9][0-9]*|-1)").all():
        raise ValueError("Expected a positive integer charge or Sage combined-charge sentinel -1")
    if (
        not features.input_line.str.fullmatch(r"[0-9]+").all()
        or (features.input_line.astype(int) < 2).any()
    ):
        raise ValueError("Invalid source feature line number")
    if not np.isfinite(features.intensity).all() or not (features.intensity > 0).all():
        raise ValueError("Expected finite positive accepted feature intensities")
    if not features.peptide.map(canonical_peptide).equals(features.canonical_peptide):
        raise ValueError("Feature canonical peptide is inconsistent")
    if not (features.canonical_peptide.map(mapping) == features.study_group).all():
        raise ValueError("Feature assignment differs from the fixed study dictionary")
    acquisitions = _read_frame(source / "acquisition_summary.tsv", dtype={"sample": str})
    if acquisitions["sample"].tolist() != samples:
        raise ValueError("Acquisition summary roster differs from the study matrix")
    expected_counts = acquisitions.set_index("sample").lfq_features
    if not np.array_equal(
        features.groupby("sample").size().reindex(samples, fill_value=0), expected_counts
    ):
        raise ValueError("Accepted feature counts differ from the acquisition summary")
    modified = features.peptide.str.replace("I", "L", regex=False)
    charge_modes = pd.DataFrame({"peptide": modified, "combined": features.charge == "-1"})
    if (charge_modes.groupby("peptide").combined.nunique() > 1).any():
        raise ValueError("A modified peptide mixes combined-charge and charge-resolved quantities")
    features["ion"] = (
        features.peptide.str.replace("I", "L", regex=False) + "|charge=" + features.charge
    )
    keys = features[["ion", "study_group", "canonical_peptide"]].drop_duplicates()
    if keys.ion.duplicated().any():
        raise ValueError("An ion belongs to more than one study group or canonical peptide")
    ion_codes, ions = pd.factorize(features.ion, sort=True)
    values = np.zeros((len(ions), len(samples)), dtype=np.float64)
    sample_codes = pd.Categorical(features["sample"], categories=samples).codes
    np.add.at(values, (ion_codes, sample_codes), features.intensity.to_numpy())
    keys = keys.set_index("ion").reindex(ions)
    matrix = pd.DataFrame(
        values,
        index=pd.MultiIndex.from_arrays([keys.study_group, ions], names=["protein", "ion"]),
        columns=samples,
    ).sort_index()
    sums = matrix.groupby(level="protein").sum()
    if set(sums.index) != set(baseline.index) or not np.allclose(
        sums.reindex(baseline.index), baseline, rtol=1e-10, atol=0
    ):
        raise ValueError("Exported ions do not reproduce every study-group sum")
    path = temp / "ions.aq_reformat.tsv"
    matrix.to_csv(path, sep="\t")
    settings = dict(
        num_cores=threads,
        min_nonan=1,
        number_of_quadratic_samples=50,
        maximum_number_of_quadratic_ions_to_use_per_protein=10,
        deactivate_normalization=False,
        log_processed_proteins=False,
        compile_normalized_ion_table=True,
    )
    if len(matrix):
        imported = utils.import_data(str(path))
        if list(imported.columns) != ["protein", "ion"] + samples or len(imported) != len(matrix):
            raise ValueError("directLFQ importer changed the ion table or acquisition roster")
        imported = imported.set_index(["protein", "ion"])
        if not imported.index.equals(matrix.index) or not np.allclose(
            imported, matrix, rtol=1e-13, atol=0
        ):
            raise ValueError("directLFQ importer changed ion identities or intensities")
        del imported, features, values
        lfq_manager.run_lfq(str(path), **settings)
        protein = _read_frame(Path(str(path) + ".protein_intensities.tsv"), dtype={"protein": str})
        protein = protein.set_index("protein")
        ion_output = _read_frame(
            Path(str(path) + ".ion_intensities.tsv"), dtype={"protein": str, "ion": str}
        ).set_index(["protein", "ion"])
    else:
        protein = pd.DataFrame(columns=samples, dtype=float)
        ion_output = matrix.copy()
    if (
        list(protein.columns) != samples
        or not protein.index.is_unique
        or not set(protein.index) <= set(baseline.index)
        or not np.isfinite(protein.to_numpy()).all()
        or (protein < 0).any().any()
    ):
        raise ValueError("Invalid directLFQ protein output")
    if (
        list(ion_output.columns) != samples
        or not ion_output.index.is_unique
        or not ion_output.index.isin(matrix.index).all()
        or not np.isfinite(ion_output.to_numpy()).all()
        or (ion_output < 0).any().any()
    ):
        raise ValueError("Invalid directLFQ ion output")
    retained = ion_output.reindex(matrix.index, fill_value=0) > 0
    if (retained & ~(matrix > 0)).any().any():
        raise ValueError("directLFQ created a new positive ion cell")
    protein = protein.reindex(baseline.index, fill_value=0)
    if ((protein > 0) & ~(baseline > 0)).any().any():
        raise ValueError("directLFQ created a new positive group cell")
    # Keep all source groups and acquisitions, including entirely lost/empty ones.
    ion_index = matrix.index.get_level_values("ion")
    diagnostics = keys.reindex(ion_index).copy()
    diagnostics.index.name = "ion"
    diagnostics["input_acquisitions"] = (matrix > 0).sum(axis=1).to_numpy()
    diagnostics["retained_acquisitions"] = retained.sum(axis=1).to_numpy()
    diagnostics["in_output_ion_roster"] = matrix.index.isin(ion_output.index).astype(int)
    diagnostics.reset_index().to_csv(temp / "ion_mapping.tsv", sep="\t", index=False)
    protein.index.name = "study_group"
    protein.where(protein > 0).to_csv(temp / "study_group_matrix.tsv", sep="\t")
    _peptide_support(temp, matrix, retained, keys, protein)
    return {
        "settings": settings,
        "backend_implementation": [
            stamp(Path(module.__file__))
            for module in [lfq_manager, lfq_manager.lfqnorm, lfq_manager.lfqprot_estimation, utils]
        ],
        "dependency_versions": {
            name: importlib.metadata.version(name)
            for name in ["numpy", "pandas", "numba", "multiprocess", "pyarrow"]
        },
        "sample_normalization": True,
        "quantified_groups": (protein > 0).sum().tolist(),
        "input_ions": len(matrix),
        "input_positive_ion_cells": int((matrix > 0).sum().sum()),
        "retained_positive_ion_cells": int(retained.sum().sum()),
        "backend_ran": bool(len(matrix)),
        "baseline_sums_reproduced": True,
        "ion_selection": "directLFQ 0.3.3 caps groups at 100 ions; sparse/disconnected "
        "profiles may lose quantities. min_nonan=1 is not a two-peptide rule.",
    }


def _peptide_support(temp, matrix, retained, keys, protein):
    """Write observed peptide support before/after estimation, one acquisition at a time."""

    import pandas as pd

    key = keys.reindex(matrix.index.get_level_values("ion"))
    index = pd.MultiIndex.from_arrays(
        [matrix.index.get_level_values("protein"), key.canonical_peptide],
        names=["study_group", "canonical_peptide"],
    )
    with open_gzip_text(temp / "peptide_support.tsv.gz") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(
            [
                "sample",
                "study_group",
                "canonical_peptide",
                "input_ions",
                "retained_ions",
                "input_intensity",
                "group_quantified",
            ]
        )
        for sample in matrix.columns:
            frame = pd.DataFrame(
                {
                    "input_ions": (matrix[sample] > 0).to_numpy(),
                    "retained_ions": retained[sample].to_numpy(),
                    "input_intensity": matrix[sample].to_numpy(),
                },
                index=index,
            )
            frame = frame.groupby(level=[0, 1]).sum()
            frame = frame[frame.input_ions > 0]
            for (group, peptide), row in frame.iterrows():
                writer.writerow(
                    [
                        sample,
                        group,
                        peptide,
                        int(row.input_ions),
                        int(row.retained_ions),
                        row.input_intensity,
                        int(protein.loc[group, sample] > 0),
                    ]
                )
