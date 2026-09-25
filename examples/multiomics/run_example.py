"""Real CAMPI spectra with synthetic molecular inputs: an additive workflow check.

TPM, KO and compound values below are artificial test inputs, not cohort measurements.
"""

import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path

from fasta_lake.multi_omics.exact import read_fasta, sha256

ROOT = Path(__file__).resolve().parents[2]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campi", type=Path, required=True, help="Completed bundled CAMPI example")
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--reference-chunk-mib", type=int)
    parser.add_argument("--laptop", action="store_true")
    args = parser.parse_args()
    campi, out = args.campi.resolve(), args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    commands = []

    def command(name, argv):
        argv = list(map(str, argv))
        with (out / (name + ".log")).open("x") as log:
            completed = subprocess.run(argv, stdout=log, stderr=subprocess.STDOUT)
        commands.append(dict(stage=name, argv=argv, exit_code=completed.returncode))
        (out / "commands.json").write_text(json.dumps(commands, indent=2) + "\n")
        if completed.returncode:
            raise RuntimeError(f"{name} failed; inspect {out / (name + '.log')}")

    module = [sys.executable, "-m", "fasta_lake.cli", "molecular"]
    candidates = {r[2]: r[1] for r in read_fasta(campi / "workflow/evidence.fasta").values()}
    with (ROOT / "examples/campi/inputs/samples.tsv").open() as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    expected = {}
    manifest_rows = []
    for row in rows:
        sample = row["sample"]
        raw = out / "inputs" / sample
        raw.mkdir(parents=True)
        baseline = campi / "workflow/samples" / sample / "inference" / (sample + "_razor.fasta")
        original = read_fasta(baseline)
        base_hashes = {r[2] for r in original.values()}
        extras = sorted(set(candidates) - base_hashes)[:3]
        if len(extras) != 3:
            raise RuntimeError(
                "Expected three real candidate proteins outside the de novo baseline"
            )
        # g0, g1 and g2 provide DNA, RNA and metabolite additions, respectively.
        ordered = extras + sorted(base_hashes)
        with (raw / "proteins.faa").open("x") as fasta, (raw / "tpm.tsv").open("x") as tpm:
            tpm.write("prodigal_protein\tmetaG_tpm\tmetaT_tpm\n")
            for i, digest in enumerate(ordered):
                fasta.write(f">g{i}\n{candidates[digest]}\n")
                tpm.write(f"g{i}\t{100 if i == 0 else 0}\t{100 if i == 1 else 'NA'}\n")
            tpm.write("unavailable_test_gene\t10000\t10000\n")
        (raw / "annotations.tsv").write_text("prodigal_protein\tKEGG_ko\ng2\tK00003\n")
        (raw / "metabolites.tsv").write_text(f"specimen\tcompound\tvalue\n{sample}\tC00001\t1\n")
        (raw / "links.tsv").write_text("compound\treaction\tKO\nC00001\tR00001\tK00003\n")
        (raw / "provenance.json").write_text(
            json.dumps(
                dict(
                    kind="synthetic molecular values and annotations over real CAMPI proteins",
                    spectra="measured bundled CAMPI acquisition window",
                    biological_performance_claim=False,
                )
            )
            + "\n"
        )
        destination = out / "sources" / sample
        command(
            sample + "_prepare",
            [
                *module,
                "prepare-source",
                "--out",
                destination,
                "--specimen",
                sample,
                "--reference-id",
                "synthetic_molecular_campi",
                "--reference-scope",
                "available",
                "--tpm",
                raw / "tpm.tsv",
                "--proteins",
                raw / "proteins.faa",
                "--provenance",
                raw / "provenance.json",
                "--annotations",
                raw / "annotations.tsv",
                "--metabolites",
                raw / "metabolites.tsv",
                "--compound-links",
                raw / "links.tsv",
            ],
        )
        manifest_rows.append(
            dict(
                sample=sample,
                predictions=str(ROOT / "examples/campi/inputs" / row["predictions"]),
                mzml=str(ROOT / "examples/campi/inputs" / row["mzml"]),
                molecular_source=str(destination / "source.json"),
                specimen=sample,
            )
        )
        expected[sample] = dict(
            baseline=baseline, hashes=base_hashes | set(extras), baseline_records=original
        )
    manifest = out / "samples.tsv"
    with manifest.open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(manifest_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(manifest_rows)
    common = [
        sys.executable,
        ROOT / "tools/run_manifest.py",
        "--manifest",
        manifest,
        "--lake",
        campi / "lake.fasta",
        "--out",
        out / "workflow",
        "--molecular-dna-percent",
        "100",
        "--molecular-rna-percent",
        "100",
        "--molecular-metabolites",
        "--quantification",
        "directlfq",
        "--threads",
        args.threads,
    ]
    if args.reference_chunk_mib is not None:
        common += ["--reference-chunk-mib", args.reference_chunk_mib]
    if args.laptop:
        common += ["--laptop"]
    command("preflight", [*common, "--validate-only"])
    if (out / "workflow").exists():
        raise RuntimeError("Preflight created a workflow output")
    command("workflow", common)
    verified = []
    for sample, record in expected.items():
        sample_dir = out / "workflow/samples" / sample
        original = sample_dir / "inference" / (sample + "_razor.fasta")
        if original.read_bytes() != record["baseline"].read_bytes():
            raise RuntimeError("Molecular options altered the de novo inference baseline")
        final = read_fasta(sample_dir / "molecular/search.fasta")
        if {r[2] for r in final.values()} != record["hashes"]:
            raise RuntimeError("Final database does not contain baseline plus the three additions")
        if any(final[a] != r for a, r in record["baseline_records"].items()):
            raise RuntimeError("A de novo protein or header changed")
        config = json.loads((out / "workflow/search" / sample / "config.json").read_text())
        if Path(config["database"]["fasta"]) != sample_dir / "molecular/search.fasta":
            raise RuntimeError("Sage did not search the molecular union")
        command(
            sample + "_compare",
            [
                *module,
                "compare",
                "--left",
                campi / "workflow/search" / sample,
                "--right",
                out / "workflow/search" / sample,
                "--left-label",
                "denovo",
                "--right-label",
                "denovo_plus_molecular",
                "--out",
                out / "comparisons" / sample,
            ],
        )
        command(
            sample + "_increments",
            [
                *module,
                "increments",
                "--left",
                campi / "workflow/search" / sample,
                "--right",
                out / "workflow/search" / sample,
                "--reference",
                campi / "lake.fasta",
                "--out",
                out / "increments" / sample,
            ],
        )
        comparison = json.loads((out / "comparisons" / sample / "COMPLETE.json").read_text())
        increments = json.loads((out / "increments" / sample / "summary.json").read_text())
        for key, expected_count in (
            ("peptides_gained", comparison["peptide_right_only"]),
            ("peptides_lost", comparison["peptide_left_only"]),
            ("peptides_retained", comparison["peptide_shared"]),
            ("lfq_left", comparison["lfq_features"]["denovo"]),
            ("lfq_right", comparison["lfq_features"]["denovo_plus_molecular"]),
            ("candidate_sequences_added", 3),
            ("candidate_sequences_removed", 0),
        ):
            if increments[key] != expected_count:
                raise RuntimeError(f"Reference-wide increment accounting disagrees for {key}")
        verified.append(
            dict(
                sample=sample,
                baseline_sequences=len(record["baseline_records"]),
                additions=3,
                final_sequences=len(final),
                fasta_sha256=sha256(sample_dir / "molecular/search.fasta"),
            )
        )
    command(
        "shared_reference_databases",
        [
            sys.executable,
            ROOT / "tools/run_manifest.py",
            "--manifest",
            manifest,
            "--molecular-reference",
            "shared",
            "--out",
            out / "shared_reference",
            "--molecular-dna-percent",
            "100",
            "--molecular-rna-percent",
            "100",
            "--molecular-metabolites",
            "--databases-only",
            "--threads",
            args.threads,
        ]
        + (["--laptop"] if args.laptop else [])
        + (
            ["--reference-chunk-mib", args.reference_chunk_mib]
            if args.reference_chunk_mib is not None
            else []
        ),
    )
    for sample in expected:
        folder = out / "shared_reference/samples" / sample
        baseline = folder / "inference" / (sample + "_razor.fasta")
        final = folder / "molecular/search.fasta"
        records = read_fasta(baseline)
        final_records = read_fasta(final)
        if any(final_records[a] != r for a, r in records.items()):
            raise RuntimeError("Shared-reference mode dropped a de novo protein/header")
        if not final.read_bytes().startswith(baseline.read_bytes()):
            raise RuntimeError("Shared-reference baseline bytes changed")
        selected = json.loads((folder / "molecular/COMPLETE.json").read_text())
        if selected["reference_coverage"]["unavailable_genes"] != 1:
            raise RuntimeError("Shared-reference mode lost the missing-gene ledger")
    receipt = dict(
        status="PASS",
        scope=__doc__,
        samples=verified,
        real_sage_searches=len(verified),
        quantification="directlfq",
        shared_reference_database_checks=len(expected),
        reference_wide_increment_reports=len(verified),
        preserves_every_denovo_sequence_and_header=True,
        unavailable_genes_excluded_from_both_tpm_layers=True,
    )
    (out / "VALIDATION.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    main()
