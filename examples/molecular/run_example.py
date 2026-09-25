"""Synthetic format/selection example; these are not measured biological data."""

import argparse
import json
import subprocess
import sys
from pathlib import Path

from fasta_lake.multi_omics.exact import read_fasta, sha256


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    source = args.out / "inputs"
    source.mkdir()
    (source / "genes.faa").write_text(
        ">gene_A original gene A\nMACDEK\n>gene_B same sequence\nMACDEK\n"
        ">gene_C RNA evidence\nMLLLK\n>gene_D unselected\nMVVVVK\n"
    )
    (source / "tpm.tsv").write_text(
        "prodigal_protein\tmetaG_tpm\tmetaT_tpm\n"
        "gene_A\t1\tNA\ngene_B\t10\tNA\ngene_C\t0\t2\ngene_D\t0\t0\n"
    )
    (source / "provenance.json").write_text(
        json.dumps(
            {
                "kind": "synthetic example; no measured TPM or biological provenance",
                "reference": "the four original protein records bundled by this script",
            }
        )
        + "\n"
    )
    manifest = dict(
        schema="fastalake.molecular-source.v1",
        specimen="synthetic_A",
        reference_id="synthetic_reference_v1",
        reference_status="verified_original_quantification_reference",
    )
    for key, name in [
        ("tpm", "tpm.tsv"),
        ("proteins", "genes.faa"),
        ("provenance", "provenance.json"),
    ]:
        manifest[key] = dict(path=name, sha256=sha256(source / name))
    (source / "source.json").write_text(json.dumps(manifest, indent=2) + "\n")
    (source / "baseline.fasta").write_text(">baseline preserve this header\nMQQQQK\n")
    command = [
        sys.executable,
        "-m",
        "fasta_lake.cli",
        "molecular",
        "add",
        "--baseline",
        str(source / "baseline.fasta"),
        "--source",
        str(source / "source.json"),
        "--specimen",
        "synthetic_A",
        "--out",
        str(args.out / "union"),
        "--dna-percent",
        "50",
        "--rna-absolute",
        "1",
    ]
    completed = subprocess.run(command, check=True, capture_output=True, text=True)
    result = json.loads(completed.stdout)
    observed = {r[1] for r in read_fasta(args.out / "union/search.fasta").values()}
    if observed != {"MQQQQK", "MACDEK", "MLLLK"} or result["additions"] != 2:
        raise RuntimeError("Synthetic example did not reproduce its expected union")
    print(
        json.dumps(
            {
                "status": "PASS",
                "scope": "synthetic exact-union example",
                "baseline_sequences": 1,
                "added_sequences": 2,
                "final_sequences": 3,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
