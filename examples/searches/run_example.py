"""Exercise six Sage search choices on the same measured CAMPI windows.

This validates execution, preserved candidate sequences and reporting. A larger
search space is not assumed to improve identification or establish protease activity.
"""

import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path

from fasta_lake.multi_omics.exact import sha256
from fasta_lake.search import DIGESTIONS, SEARCH_MODES, sage_config

ROOT = Path(__file__).resolve().parents[2]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campi", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--search-max-length", type=int, default=30)
    args = parser.parse_args()
    source, out = args.campi.resolve(), args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    manifest = ROOT / "examples/campi/inputs/samples.tsv"
    with manifest.open() as stream:
        samples = list(csv.DictReader(stream, delimiter="\t"))
    checks = []
    for digestion in DIGESTIONS:
        for mode in SEARCH_MODES:
            label = digestion + "_" + mode
            destination = out / label
            argv = [
                sys.executable,
                str(ROOT / "tools/run_manifest.py"),
                "--manifest",
                str(manifest),
                "--lake",
                str(source / "lake.fasta"),
                "--out",
                str(destination),
                "--threads",
                str(args.threads),
                "--digestion",
                digestion,
                "--search-mode",
                mode,
                "--search-min-length",
                "7",
                "--search-max-length",
                str(args.search_max_length),
            ]
            with (out / (label + ".log")).open("x") as log:
                result = subprocess.run(argv, stdout=log, stderr=subprocess.STDOUT)
            if result.returncode:
                raise RuntimeError(f"{label} failed; inspect {out / (label + '.log')}")
            assert (destination / "evidence.fasta").read_bytes() == (
                source / "workflow/evidence.fasta"
            ).read_bytes()
            for row in samples:
                sample = row["sample"]
                relative = f"samples/{sample}/inference/{sample}_razor.fasta"
                assert (destination / relative).read_bytes() == (
                    source / "workflow" / relative
                ).read_bytes()
                search = destination / "search" / sample
                config = json.loads((search / "config.json").read_text())
                assert config == sage_config(
                    destination / relative,
                    manifest.parent / row["mzml"],
                    search,
                    digestion=digestion,
                    search_mode=mode,
                    min_length=7,
                    max_length=args.search_max_length,
                )
                with (search / "results.sage.tsv").open() as stream:
                    psms = list(csv.DictReader(stream, delimiter="\t"))
                assert psms, f"No PSM output for {label}/{sample}"
                assert all(Path(p["filename"]).name == row["mzml"] for p in psms)
                checks.append(
                    {
                        "digestion": digestion,
                        "mode": mode,
                        "sample": sample,
                        "psms_reported": len(psms),
                        "fasta_sha256": sha256(destination / relative),
                        "config_sha256": sha256(search / "config.json"),
                        "unchanged_de_novo_fasta": True,
                    }
                )
            qc = json.loads((destination / "study_groups/qc/summary.json").read_text())
            assert qc["status"] == "PASS"
            assert (destination / "study_groups/qc/QC_overview.png").stat().st_size > 1000
    (out / "VALIDATION.json").write_text(
        json.dumps(
            {
                "status": "PASS",
                "searches": checks,
                "scope": "Cropped CAMPI spectra; six choices, fixed candidate FASTAs, automatic QC",
                "biological_or_fdr_calibration": False,
            },
            indent=2,
        )
        + "\n"
    )
    print("PASS: six search choices, identical de novo databases, measured spectra and QC")


if __name__ == "__main__":
    main()
