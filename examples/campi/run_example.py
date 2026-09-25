#!/usr/bin/env python3
"""Run and verify two real CAMPI acquisition windows from the bundled inputs.

First use needs Python 3.12, venv/pip, Cargo, a C compiler and internet for Python
and Rust dependencies. A local Python environment and binaries are created inside
the new output folder. SAGE 0.14.6 for Linux x86_64 is bundled with its MIT licence.
Use --runtime PREVIOUS_OUTPUT/_runtime to reuse a completed, verified setup.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import platform
import shutil
import subprocess
import sys
import tarfile
import time
from pathlib import Path

from check_results import sha, validate

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CRATES = ("fasta_extractor_v2", "parsimony_engine", "aggregation_engine")
BINARIES = ("lake_builder", "fasta_extractor", "parsimony_engine", "aggregation_engine")
SAGE_ARCHIVE = "sage-v0.14.6-x86_64-unknown-linux-gnu.tar.gz"
SAGE_SHA256 = "3492d04b9c922d47ee44fd7fd4248b946d8afb7e88893d1d97c403bb3760be30"


def source_identity():
    paths = [
        ROOT / "pyproject.toml",
        ROOT / "rust-toolchain.toml",
        ROOT / "requirements/release.txt",
        ROOT / "requirements/container.txt",
    ]
    paths += [
        p for p in (ROOT / "fasta_lake").rglob("*") if p.is_file() and "__pycache__" not in p.parts
    ]
    paths += [p for p in (ROOT / "rust/common").rglob("*.rs") if p.is_file()]
    for crate in CRATES:
        directory = ROOT / "rust" / crate
        paths += [
            p
            for p in directory.rglob("*")
            if p.is_file() and "target" not in p.relative_to(directory).parts
        ]
    manifest = "".join(str(p.relative_to(ROOT)) + "\t" + sha(p) + "\n" for p in sorted(paths))
    return hashlib.sha256(manifest.encode()).hexdigest()


def verify_inputs():
    directory = HERE / "inputs"
    seen = set()
    for line in (directory / "SHA256SUMS").read_text().splitlines():
        digest, name = line.split("  ", 1)
        if name in seen:
            raise ValueError("Duplicate input checksum entry: " + name)
        seen.add(name)
        if Path(name).name != name or sha(directory / name) != digest:
            raise ValueError("Input checksum mismatch: " + name)
    actual = {p.name for p in directory.iterdir() if p.is_file() and p.name != "SHA256SUMS"}
    if not seen or seen != actual:
        raise ValueError("Input checksum manifest must cover every bundled input exactly once")
    if sha(ROOT / "examples/campi/runtime" / SAGE_ARCHIVE) != SAGE_SHA256:
        raise ValueError("Bundled SAGE archive checksum mismatch")
    if not (HERE / "expected/results.json").is_file():
        raise ValueError("Expected results are missing from the example")


def main(example_dir=None, cohort="CAMPI"):
    global HERE
    if example_dir is not None:
        HERE = Path(example_dir).resolve()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--out", type=Path, required=True, help="New output directory; never reused or deleted"
    )
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--runtime", type=Path, help="Reuse a completed PREVIOUS_OUTPUT/_runtime")
    args = parser.parse_args()
    if sys.version_info[:2] != (3, 12):
        parser.error(
            "Python 3.12 is required: the example installs the hash-locked package set "
            "compiled for Python 3.12 on Linux x86_64 (requirements/container.txt)"
        )
    if platform.system() != "Linux" or platform.machine() not in ("x86_64", "AMD64"):
        parser.error(
            "This bundled runtime was tested on Linux x86_64. "
            "Use the portable manual workflow on other platforms."
        )
    if args.threads < 1:
        parser.error("--threads must be positive")
    if args.out.exists() or args.out.is_symlink():
        parser.error("Output already exists; choose a new directory")
    if not args.runtime and not shutil.which("cargo"):
        parser.error(
            "Cargo was not found. Install the Rust toolchain, "
            "or pass a completed --runtime directory."
        )
    verify_inputs()
    identity = source_identity()
    runtime = args.runtime.resolve() if args.runtime else args.out.resolve() / "_runtime"
    if args.runtime:
        saved = json.loads((runtime / "RUNTIME.json").read_text())
        if saved["source_sha256"] != identity:
            parser.error("The runtime was built from different source files; create a new runtime")
        for name, digest in saved["binaries_sha256"].items():
            if sha(runtime / "bin" / name) != digest:
                parser.error("Runtime binary changed: " + name)
        if not (runtime / "venv/bin/python").is_file():
            parser.error("Runtime Python environment is unavailable")
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    env = {
        **os.environ,
        "PYTHONPATH": "",
        "PYTHONDONTWRITEBYTECODE": "1",
        "RAYON_NUM_THREADS": str(args.threads),
        "CARGO_BUILD_JOBS": str(args.threads),
    }
    commands = []

    def execute(name, command, extra=None, cwd=out):
        print(name.replace("_", " ") + "…", flush=True)
        begin = time.monotonic()
        command = list(map(str, command))
        with (out / (name + ".log")).open("x") as log:
            status = subprocess.run(
                command, cwd=cwd, env={**env, **(extra or {})}, stdout=log, stderr=subprocess.STDOUT
            ).returncode
        commands.append(
            {"stage": name, "argv": command, "status": status, "seconds": time.monotonic() - begin}
        )
        (out / "commands.json").write_text(json.dumps(commands, indent=2) + "\n")
        if status:
            raise RuntimeError(name + " failed; inspect " + str(out / (name + ".log")))

    if not args.runtime:
        runtime.mkdir()
        (runtime / "bin").mkdir()
        execute("create_python_environment", [sys.executable, "-m", "venv", runtime / "venv"])
        # Build/install from a retained copy so the shared source tree stays clean.
        source = runtime / "python_source"
        source.mkdir()
        for name in ("pyproject.toml", "README.md", "LICENSE", "MANIFEST.in"):
            shutil.copy2(ROOT / name, source / name)
        shutil.copytree(
            ROOT / "fasta_lake",
            source / "fasta_lake",
            ignore=shutil.ignore_patterns("__pycache__", "*.pyc"),
        )
        # The locked dependency set, identical to the container and CI, then the
        # package itself without resolving anything further from PyPI.
        execute(
            "install_locked_dependencies",
            [
                runtime / "venv/bin/python",
                "-m",
                "pip",
                "install",
                "--no-clean",
                "--require-hashes",
                "-r",
                ROOT / "requirements/container.txt",
            ],
        )
        execute(
            "install_python_package",
            [
                runtime / "venv/bin/python",
                "-m",
                "pip",
                "install",
                "--no-clean",
                "--no-deps",
                str(source) + "[analysis]",
            ],
        )
        for crate in CRATES:
            execute(
                "build_" + crate,
                [
                    "cargo",
                    "build",
                    "--release",
                    "--locked",
                    "--manifest-path",
                    ROOT / "rust" / crate / "Cargo.toml",
                ],
                extra={
                    "CARGO_TARGET_DIR": str(runtime / "cargo"),
                    "FASTALAKE_BUILD_REV": "example-" + identity[:16],
                },
                cwd=ROOT,
            )
        for name in BINARIES:
            shutil.copy2(runtime / "cargo/release" / name, runtime / "bin" / name)
        with tarfile.open(ROOT / "examples/campi/runtime" / SAGE_ARCHIVE) as archive:
            member = archive.getmember("sage-v0.14.6-x86_64-unknown-linux-gnu/sage")
            if not member.isfile():
                raise ValueError("SAGE archive contains an invalid executable entry")
            with archive.extractfile(member) as source, (runtime / "bin/sage").open("xb") as target:
                shutil.copyfileobj(source, target)
        (runtime / "bin/sage").chmod(0o755)
        execute("check_sage_version", [runtime / "bin/sage", "--version"])
        if (out / "check_sage_version.log").read_text().strip() != "sage 0.14.6":
            raise ValueError("Unexpected SAGE version")
        (runtime / "RUNTIME.json").write_text(
            json.dumps(
                {
                    "source_sha256": identity,
                    "platform": platform.platform(),
                    "binaries_sha256": {
                        name: sha(runtime / "bin" / name) for name in (*BINARIES, "sage")
                    },
                },
                indent=2,
            )
            + "\n"
        )
    python = runtime / "venv/bin/python"
    env["FASTALAKE_BIN_DIR"] = str(runtime / "bin")
    execute(
        "check_python_runtime",
        [
            python,
            "-c",
            "from importlib.metadata import version; import fasta_lake; "
            "assert version('click') == '8.3.3'; assert version('pyahocorasick') == '2.2.0'; "
            "assert fasta_lake.__version__ == '1.1.0rc2'; print('Runtime versions verified')",
        ],
    )
    workflow_start = time.monotonic()
    reference = HERE / "inputs/reference.fasta"
    compressed = HERE / "inputs/reference.fasta.gz"
    parts_manifest = HERE / "inputs/reference_parts.json"
    if sum(p.is_file() for p in (reference, compressed, parts_manifest)) != 1:
        raise ValueError("Bundle must declare exactly one reference representation")
    parts = [compressed] if compressed.is_file() else []
    if parts_manifest.is_file():
        names = json.loads(parts_manifest.read_text())["parts"]
        if not names or len(set(names)) != len(names):
            raise ValueError("Reference parts must be nonempty and distinct")
        if any(Path(name).name != name for name in names):
            raise ValueError("Reference parts must be local input filenames")
        parts = [HERE / "inputs" / name for name in names]
    if parts:
        reference = out / "reference.fasta"
        with reference.open("xb") as target:
            for part in parts:
                with gzip.open(part, "rb") as source:
                    shutil.copyfileobj(source, target)
        expected_sha = json.loads((HERE / "inputs/input_provenance.json").read_text())["reference"][
            "sha256"
        ]
        if sha(reference) != expected_sha:
            raise ValueError("Reconstructed reference checksum mismatch")
    execute(
        "build_reference_lake",
        [
            runtime / "bin/lake_builder",
            "--source",
            cohort + ":1:" + str(reference),
            "--uniform-tag",
            cohort,
            "--track-all-sources",
            "--output",
            out / "lake.fasta",
            "--output-headers",
            out / "source_headers.tsv",
            "--output-stats",
            out / "lake_stats.txt",
            "--threads",
            args.threads,
        ],
    )
    execute(
        "run_acquisitions",
        [
            python,
            ROOT / "tools/run_manifest.py",
            "--manifest",
            HERE / "inputs/samples.tsv",
            "--lake",
            out / "lake.fasta",
            "--out",
            out / "workflow",
            "--sage",
            runtime / "bin/sage",
            "--threads",
            args.threads,
        ],
    )
    # The normal runner now groups and plots automatically. Keep the documented
    # example path as a relative alias without repeating grouping or quantification.
    (out / "study").symlink_to("workflow/study_groups", target_is_directory=True)
    qc = out / "study/qc"
    qc_record = json.loads((qc / "summary.json").read_text())
    if qc_record["status"] != "PASS" or qc_record["samples"] != 2:
        raise RuntimeError("Automatic two-acquisition QC did not complete")
    if not qc_record["pca"]["status"].startswith("SKIPPED:"):
        raise RuntimeError("Two-acquisition example must explicitly skip PCA")
    for ext in ("pdf", "png"):
        if (qc / ("QC_overview." + ext)).stat().st_size <= 1000:
            raise RuntimeError("Automatic QC plot is missing or empty")
    result = validate(out, HERE / "expected")
    result.update(
        workflow_seconds=time.monotonic() - workflow_start,
        total_seconds=time.monotonic() - started,
        runtime=str(runtime),
        source_sha256=identity,
        example_runner_sha256=sha(__file__),
        expected_results_sha256=sha(HERE / "expected/results.json"),
        automatic_qc={
            "status": qc_record["status"],
            "pca": qc_record["pca"],
            "summary": "study/qc/summary.json",
        },
    )
    (out / "VALIDATION.json").write_text(json.dumps(result, indent=2) + "\n")
    if result["status"] != "passed":
        raise RuntimeError("Expected-result checks failed: " + "; ".join(result["errors"]))
    (out / "COMPLETE.json").write_text(
        json.dumps(
            {
                "status": "passed",
                "validation": "VALIDATION.json",
                "scope": (
                    "Two real acquisition windows; database construction, "
                    "search, LFQ and study reporting"
                ),
            },
            indent=2,
        )
        + "\n"
    )
    print("\nPASS: both real-data acquisitions and expected-result checks completed.", flush=True)
    for sample, row in result["measurements"]["samples"].items():
        print(
            f"{sample}: {row['selected_sequences']} selected sequences; "
            f"{row['canonical_target_peptides_peptide_q_001']} accepted canonical peptides"
        )
    print("Inferred study groups:", result["measurements"]["study"]["inferred_study_groups"])
    print("Results:", out)
    print("Reuse this setup with --runtime", runtime)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, RuntimeError, subprocess.SubprocessError) as error:
        print("ERROR:", error, file=sys.stderr)
        raise SystemExit(1)
