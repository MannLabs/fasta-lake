#!/usr/bin/env python3
"""Portable, explicit Rust workflow for a declared acquisition manifest.

Inputs and prior outputs are preserved. Chunked extraction reclaims its own
consumed scratch pieces. A run requires a new output directory; failed runs remain
available for diagnosis. Run this script from an installed FastaLake environment.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

from fasta_lake import __version__
from fasta_lake.resources import run_recorded
from fasta_lake.rust.binaries import find_binary
from fasta_lake.search import DIGESTIONS, SEARCH_MODES, sage_config

ID = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]*\Z")


def sha256(path):
    """Stream file bytes into a SHA-256 digest for the run's provenance."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_manifest(path):
    """Validate acquisitions and resolve inputs relative to the manifest.

    Parameters
    ----------
    path : path-like
        TSV with sample, predictions and mzml columns. Molecular rows additionally
        need paired specimen and molecular_source fields.

    Returns
    -------
    list[dict]
        Rows with absolute input paths and verified molecular-source identities.

    Raises
    ------
    ValueError
        For duplicate/invalid sample IDs, missing files, unsupported input
        schemas or an inconsistent specimen/source pairing.

    Notes
    -----
    One sample ID names one acquisition. The function reads and checks inputs;
    it does not create workflow outputs or infer biological pairing from names.
    """
    path = Path(path).resolve()
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        fields = reader.fieldnames or []
        if len(fields) != len(set(fields)):
            raise ValueError("Manifest column names must be unique")
        required = {"sample", "predictions", "mzml"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError("Manifest requires tab-separated sample, predictions, mzml columns")
        rows = list(reader)
    if not rows:
        raise ValueError("Manifest contains no acquisitions")
    seen = set()
    for row in rows:
        sample = row["sample"]
        if not sample or not ID.fullmatch(sample) or sample in seen:
            raise ValueError(f"Invalid or duplicate sample ID: {sample!r}")
        seen.add(sample)
        for key in ("predictions", "mzml"):
            if not row.get(key):
                raise ValueError(f"{sample}: missing {key}")
            source = Path(row[key]).expanduser()
            if not source.is_absolute():
                source = path.parent / source
            if not source.is_file() or not source.stat().st_size:
                raise ValueError(f"{sample}: missing or empty {key} file: {source}")
            row[key] = str(source.resolve())
        source = row.get("molecular_source", "")
        specimen = row.get("specimen", "")
        if bool(source) != bool(specimen):
            raise ValueError(f"{sample}: molecular_source and specimen must be supplied together")
        if source:
            from fasta_lake.multi_omics.sources import load_reference

            source = Path(source)
            if not source.is_absolute():
                source = path.parent / source
            source = source.resolve(strict=True)
            reference = load_reference(source)
            if reference.manifest["specimen"] != specimen:
                raise ValueError(f"{sample}: acquisition/source specimen mismatch")
            row["molecular_source"] = str(source)
            row["molecular_source_sha256"] = sha256(source)
        if Path(row["predictions"]).suffix.lower() != ".csv":
            raise ValueError(f"{sample}: this workflow requires a predictions CSV")
        if Path(row["mzml"]).suffix.lower() != ".mzml":
            raise ValueError(f"{sample}: this workflow has been validated with mzML inputs")
        with open(row["predictions"], newline="") as stream:
            header = next(csv.reader(stream), [])
        if not (
            {"score", "sequence"}.issubset(header)
            or {"score", "peptide_prediction_detokenized_unmodified"}.issubset(header)
        ):
            raise ValueError(f"{sample}: expected AlphaNovo or generic sequence,score CSV columns")
    return rows


def main(argv=None):
    """Validate and execute the acquisition-manifest workflow.

    Parameters
    ----------
    argv : list[str] or None
        Command-line arguments; None reads the process arguments.

    Returns
    -------
    int
        Zero for a completed requested scope or successful validation.

    Notes
    -----
    Preflight checks precede output creation. The normal route extracts pooled
    evidence, then candidates, selection and Sage for each acquisition before
    study grouping, quantification and QC. Optional molecular additions enter
    after selection. --validate-only does not execute analyses; --databases-only
    omits searching and downstream reporting. A failed command stops the run,
    preserving logs; COMPLETE.json marks only the scope that actually finished.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--lake", type=Path, help="Required reservoir in augment mode")
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--sage", help="Explicit path to SAGE 0.14.6 (otherwise PATH)")
    parser.add_argument(
        "--threads", type=int, help="Rust threads (default: 4, automatic in laptop mode)"
    )
    parser.add_argument(
        "--laptop", action="store_true", help="Detect resources and enable bounded reference pieces"
    )
    parser.add_argument(
        "--memory-budget-gib",
        type=float,
        help="Laptop planning budget; use Docker/Slurm for a hard RAM cap",
    )
    parser.add_argument(
        "--memory-fraction",
        type=float,
        default=0.35,
        help="Laptop share of visible RAM (default: 0.35)",
    )
    parser.add_argument(
        "--cpu-fraction",
        type=float,
        default=0.5,
        help="Laptop share of visible CPU slots (default: 0.5)",
    )
    parser.add_argument(
        "--reference-chunk-mib",
        type=int,
        help="Extract reference pieces of this target size, then merge before razor (e.g. 256)",
    )
    parser.add_argument("--top-fraction", type=float, default=0.30)
    parser.add_argument("--min-length", type=int, default=9)
    parser.add_argument("--max-length", type=int, default=50)
    parser.add_argument(
        "--digestion",
        choices=DIGESTIONS,
        default="tryptic",
        help="Search cleavage specificity; does not change de novo selection",
    )
    parser.add_argument(
        "--search-mode",
        choices=SEARCH_MODES,
        default="closed",
        help="Closed +/-20 ppm or open [-500,+100] Da precursor search",
    )
    parser.add_argument("--search-min-length", type=int, default=5)
    parser.add_argument("--search-max-length", type=int, default=50)
    parser.add_argument("--databases-only", action="store_true")
    parser.add_argument(
        "--skip-qc",
        action="store_true",
        help="Explicitly omit automatic QC/PCA for a core-only installation",
    )
    parser.add_argument("--validate-only", action="store_true")
    parser.add_argument(
        "--quantification",
        choices=["sum", "directlfq"],
        help="Study quantification after searching (default: sum); QC/PCA follows automatically",
    )
    parser.add_argument(
        "--molecular-reference",
        choices=["augment", "shared"],
        default="augment",
        help="Augment a reservoir-derived baseline, or select from each shared gene reference",
    )
    parser.add_argument("--molecular-dna-percent", type=float, default=1.0)
    parser.add_argument("--molecular-rna-percent", type=float, default=0.0)
    parser.add_argument("--molecular-dna-absolute", type=float)
    parser.add_argument("--molecular-rna-absolute", type=float)
    parser.add_argument("--molecular-metabolites", action="store_true")
    parser.add_argument(
        "--background",
        type=Path,
        help="Common target-only host/contaminant FASTA for molecular runs",
    )
    parser.add_argument(
        "--annotation-settings",
        type=Path,
        help="Optional JSON for local eggNOG/ESM-C annotation after grouping",
    )
    args = parser.parse_args(argv)
    # Validate scientific options before loading dependencies or creating output.
    sage_config(
        "targets.fasta",
        "sample.mzML",
        "search",
        digestion=args.digestion,
        search_mode=args.search_mode,
        min_length=args.search_min_length,
        max_length=args.search_max_length,
    )
    if args.databases_only and (
        args.digestion != "tryptic"
        or args.search_mode != "closed"
        or args.search_min_length != 5
        or args.search_max_length != 50
    ):
        parser.error("Search options require searches; incompatible with --databases-only")
    if args.annotation_settings:
        from fasta_lake.annotation import read_annotation_settings

        if args.databases_only:
            parser.error("Study annotation requires searches and grouping")
        read_annotation_settings(args.annotation_settings)
    resource_plan = None
    if args.laptop:
        from fasta_lake.capacity import plan_laptop

        resource_plan = plan_laptop(
            memory_fraction=args.memory_fraction,
            cpu_fraction=args.cpu_fraction,
            memory_gib=args.memory_budget_gib,
            threads=args.threads,
            chunk_mib=args.reference_chunk_mib,
        )
        args.threads = resource_plan["threads_per_process"]
        args.reference_chunk_mib = resource_plan["reference_chunk_mib"]
        # Set helper limits before preflight imports NumPy/AlphaPeptTools.
        for name in (
            "OMP_NUM_THREADS",
            "OPENBLAS_NUM_THREADS",
            "MKL_NUM_THREADS",
            "NUMBA_NUM_THREADS",
        ):
            os.environ[name] = "1"
        print(
            f"Laptop plan: {args.threads} Rust threads, "
            f"{resource_plan['memory_budget_bytes'] / 1024**3:.2f} GiB RAM planning target, "
            f"{args.reference_chunk_mib} MiB reference pieces; sequential acquisitions. "
            "Use Docker/Slurm for an enforced memory limit.",
            flush=True,
        )
    else:
        if (
            args.memory_budget_gib is not None
            or args.memory_fraction != 0.35
            or args.cpu_fraction != 0.5
        ):
            parser.error("Resource budget options require --laptop")
        args.threads = 4 if args.threads is None else args.threads
    if args.threads < 1 or not 1 <= args.min_length <= args.max_length:
        parser.error("Require positive threads and 1 <= min-length <= max-length")
    if not math.isfinite(args.top_fraction) or not 0 < args.top_fraction <= 1:
        parser.error("--top-fraction must be in (0, 1]; 0.30 means 30%")
    if args.out.exists() or args.out.is_symlink():
        parser.error("Output already exists; choose a new folder to preserve previous results")
    if args.molecular_reference == "augment":
        if args.lake is None or not args.lake.is_file() or not args.lake.stat().st_size:
            parser.error("Lake FASTA is missing or empty")
    elif args.lake is not None:
        parser.error("Shared-reference mode uses specimen sources; omit --lake")
    if args.reference_chunk_mib is not None and args.reference_chunk_mib < 1:
        parser.error("--reference-chunk-mib must be positive")
    rows = load_manifest(args.manifest)
    from fasta_lake.multi_omics.selection import validate_rules
    from fasta_lake.multi_omics.sources import load_reference, write_reference

    validate_rules(
        args.molecular_dna_percent,
        args.molecular_rna_percent,
        args.molecular_dna_absolute,
        args.molecular_rna_absolute,
    )
    molecular_rows = [r for r in rows if r.get("molecular_source")]
    if args.molecular_reference == "shared" and len(molecular_rows) != len(rows):
        parser.error("Shared-reference mode requires molecular_source and specimen for every row")
    if (
        args.background
        or args.molecular_metabolites
        or args.molecular_rna_percent
        or args.molecular_rna_absolute is not None
        or args.molecular_dna_absolute is not None
    ) and not molecular_rows:
        parser.error("Molecular options require specimen-qualified source columns")
    if args.background:
        from fasta_lake.multi_omics.selection import _targets

        _targets(args.background, reject_duplicates=False)
    if args.molecular_metabolites:
        for row in molecular_rows:
            reference = load_reference(row["molecular_source"])
            if not {"annotations", "metabolites", "compound_links"} <= set(reference.inputs):
                parser.error("Metabolite selection requires all three source annotation inputs")
    binaries = {name: find_binary(name) for name in ("fasta_extractor", "parsimony_engine")}
    if not args.databases_only:
        binaries["sage"] = find_binary(
            "sage", explicit_path=args.sage or shutil.which("sage") or "/missing/sage"
        )
        binaries["aggregation_engine"] = find_binary("aggregation_engine")
        version = subprocess.check_output([binaries["sage"], "--version"], text=True).strip()
        if version.split()[-1] != "0.14.6":
            parser.error(f"This configuration was validated with SAGE 0.14.6; found {version}")
    if args.quantification and args.databases_only:
        parser.error("--quantification requires searches; incompatible with --databases-only")
    if not args.databases_only:
        args.quantification = args.quantification or "sum"
        if not args.skip_qc:
            from fasta_lake.downstream import require_analysis

            require_analysis()
    if args.quantification == "directlfq":
        from fasta_lake.quantification import require_directlfq

        require_directlfq()
    if args.validate_only:
        print(f"Validated {len(rows)} acquisitions, required executables, and output destination")
        return 0
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    env = {**os.environ, "RAYON_NUM_THREADS": str(args.threads)}
    if resource_plan is not None:
        for name in (
            "OMP_NUM_THREADS",
            "OPENBLAS_NUM_THREADS",
            "MKL_NUM_THREADS",
            "NUMBA_NUM_THREADS",
        ):
            env[name] = "1"
        (out / "RESOURCE_PLAN.json").write_text(json.dumps(resource_plan, indent=2) + "\n")
    commands = []

    def execute(name, command):
        """Run one command with resource measurements and append its result to commands.json."""
        command = list(map(str, command))
        begin = time.monotonic()
        with (
            (out / (name + ".stdout")).open("x") as stdout,
            (out / (name + ".stderr")).open("x") as stderr,
        ):
            resources = run_recorded(command, env=env, stdout=stdout, stderr=stderr)
            status = resources["exit_code"]
        commands.append(
            {
                "stage": name,
                "argv": command,
                "exit_code": status,
                "seconds": time.monotonic() - begin,
                "resources": resources,
            }
        )
        (out / "commands.json").write_text(json.dumps(commands, indent=2) + "\n")
        if status:
            raise RuntimeError(
                f"{name} failed with exit {status}; inspect {out / (name + '.stderr')}"
            )

    provenance = {
        "version": __version__,
        "parameters": vars(args),
        "samples": rows,
        "binary_sha256": {name: sha256(path) for name, path in binaries.items()},
        "runner_sha256": sha256(__file__),
        "manifest_sha256": sha256(args.manifest),
        "lake_sha256": sha256(args.lake) if args.lake else None,
        "background_sha256": sha256(args.background) if args.background else None,
        "input_sha256": {
            r["sample"]: {key: sha256(r[key]) for key in ("predictions", "mzml")} for r in rows
        },
        "inference": {
            "strategy": "razor",
            "tiebreak": "hash-acc",
            "top_fraction": None,
            "normalize_il": True,
            "clustering": False,
        },
    }
    (out / "PROVENANCE.json").write_text(json.dumps(provenance, indent=2, default=str) + "\n")
    pool = out / "pooled_predictions"
    pool.mkdir()
    for row in rows:
        (pool / (row["sample"] + ".csv")).symlink_to(row["predictions"])
    evidence = out / "evidence.fasta"

    def extract(name, database, predictions, output):
        """Run whole-reference or chunked exact lookup with the same declared evidence settings."""
        if args.reference_chunk_mib is not None:
            pieces = output.parent / (output.stem + "_pieces")
            execute(
                name,
                [
                    sys.executable,
                    "-m",
                    "fasta_lake.cli",
                    "extract-chunked",
                    "--database",
                    database,
                    "--predictions",
                    predictions,
                    "--out",
                    pieces,
                    "--chunk-mib",
                    args.reference_chunk_mib,
                    "--top-fraction",
                    args.top_fraction,
                    "--min-length",
                    args.min_length,
                    "--max-length",
                    args.max_length,
                    "--threads",
                    args.threads,
                ],
            )
            # Linking the already merged product avoids a second large evidence copy.
            # Keep relative links so a completed output bundle remains relocatable.
            output.symlink_to(os.path.relpath(pieces / "evidence.fasta", output.parent))
            output.with_suffix(".stats.json").symlink_to(
                os.path.relpath(pieces / "evidence.stats.json", output.parent)
            )
            stats = json.loads(output.with_suffix(".stats.json").read_text())
            if not stats["proteins_matched"]:
                raise RuntimeError(f"{name}: no matched proteins across the reference pieces")
            return
        execute(
            name,
            [
                binaries["fasta_extractor"],
                "--database",
                database,
                "--peptides-dir",
                predictions,
                "--output",
                output,
                "--top-percent",
                args.top_fraction,
                "--min-length",
                args.min_length,
                "--max-length",
                args.max_length,
                "--normalize-il",
                "true",
                "-t",
                args.threads,
            ],
        )
        stats = json.loads(output.with_suffix(".stats.json").read_text())
        if not stats["proteins_matched"]:
            raise RuntimeError(
                f"{name}: no matched proteins; inspect predictions and reference coverage"
            )

    if args.molecular_reference == "augment":
        extract("evidence", args.lake.resolve(), pool, evidence)
    for row in rows:
        sample = row["sample"]
        sample_dir = out / "samples" / sample
        sample_dir.mkdir(parents=True)
        drawn = sample_dir / "sample_lake.fasta"
        if row.get("molecular_source"):
            if sha256(row["molecular_source"]) != row["molecular_source_sha256"]:
                raise RuntimeError("Molecular source manifest changed after preflight")
        database = evidence
        if args.molecular_reference == "shared":
            database = sample_dir / "shared_reference.fasta"
            write_reference(load_reference(row["molecular_source"]), database)
        extract(sample + "_extract", database, row["predictions"], drawn)
        inferred = sample_dir / "inference"
        execute(
            sample + "_infer",
            [
                binaries["parsimony_engine"],
                "--fasta",
                drawn,
                "--peptides",
                row["predictions"],
                "--output-dir",
                inferred,
                "--sample",
                sample,
                "--strategy",
                "razor",
                "--tiebreak",
                "hash-acc",
                "--min-length",
                args.min_length,
                "--max-length",
                args.max_length,
                "--normalize-il",
                "-t",
                args.threads,
            ],
        )
        fasta = inferred / (sample + "_razor.fasta")
        if not fasta.is_file() or not fasta.stat().st_size:
            raise RuntimeError(f"{sample}: inferred database is empty")
        if row.get("molecular_source"):
            destination = sample_dir / "molecular"
            command = [
                sys.executable,
                "-m",
                "fasta_lake.cli",
                "molecular",
                "add",
                "--baseline",
                fasta,
                "--source",
                row["molecular_source"],
                "--specimen",
                row["specimen"],
                "--out",
                destination,
                "--dna-percent",
                args.molecular_dna_percent,
                "--rna-percent",
                args.molecular_rna_percent,
            ]
            for flag, value in (
                ("--dna-absolute", args.molecular_dna_absolute),
                ("--rna-absolute", args.molecular_rna_absolute),
                ("--background", args.background),
            ):
                if value is not None:
                    command.extend([flag, value])
            if args.molecular_metabolites:
                command.append("--metabolites")
            execute(sample + "_molecular", command)
            fasta = destination / "search.fasta"
        if not args.databases_only:
            search = out / "search" / sample
            search.mkdir(parents=True)
            config = search / "config.json"
            config.write_text(
                json.dumps(
                    sage_config(
                        fasta,
                        row["mzml"],
                        search,
                        digestion=args.digestion,
                        search_mode=args.search_mode,
                        min_length=args.search_min_length,
                        max_length=args.search_max_length,
                    ),
                    indent=2,
                )
                + "\n"
            )
            # Sage backend: https://github.com/lazear/sage; doi:10.1021/acs.jproteome.3c00486.
            execute(sample + "_search", [binaries["sage"], config, "--write-pin"])
            for name in ("results.sage.tsv", "lfq.tsv"):
                with (search / name).open() as stream:
                    next(stream, None)
                    if next(stream, None) is None:
                        raise RuntimeError(f"{sample}: {name} contains no data rows")
    if not args.databases_only:
        execute(
            "aggregate",
            [
                binaries["aggregation_engine"],
                "--input",
                out / "search",
                "--output",
                out / "aggregate",
                "--dataset",
                "Study",
                "--method",
                "FastaLake",
                "--data-type",
                "lfq",
                "--group-key-mode",
                "member-set",
                "--q-thresholds",
                "0.01",
                "--no-functional",
                "--threads",
                args.threads,
            ],
        )
    if args.quantification:
        execute(
            "study_quantification",
            [
                sys.executable,
                Path(__file__).with_name("group_study.py"),
                "--search",
                out / "search",
                "--out",
                out / "study_groups",
                "--quantification",
                args.quantification,
                "--threads",
                args.threads,
            ]
            + (["--skip-qc"] if args.skip_qc else []),
        )
    if args.annotation_settings:
        execute(
            "study_annotation",
            [
                sys.executable,
                "-m",
                "fasta_lake.cli",
                "annotate",
                "study",
                "--groups",
                out / "study_groups",
                "--settings",
                args.annotation_settings.resolve(),
                "--out",
                out / "annotation",
            ],
        )
    (out / "COMPLETE.json").write_text(
        json.dumps(
            {
                "seconds": time.monotonic() - started,
                "acquisitions": len(rows),
                "databases_only": args.databases_only,
                "quantification": args.quantification,
                "annotation": "annotation/STUDY_COMPLETE.json"
                if args.annotation_settings
                else None,
                "qc": (
                    "NOT_APPLICABLE: databases only"
                    if args.databases_only
                    else "SKIPPED: requested with --skip-qc"
                    if args.skip_qc
                    else "study_groups/qc/summary.json"
                ),
            },
            indent=2,
        )
        + "\n"
    )
    print(f"Completed {len(rows)} acquisitions: {out}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, RuntimeError, subprocess.SubprocessError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(1)
