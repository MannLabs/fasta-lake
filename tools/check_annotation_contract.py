#!/usr/bin/env python3
"""Test the post-search annotation wiring using an explicitly synthetic mapper double.

Rust, Sage, grouping and QC execute normally on the measured CAMPI example.
The tiny mapper double validates orchestration and annotation accounting only;
real eggNOG/ESM inference is covered by separate optional-resource acceptance.
"""

import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

MAPPER = """import sys
from pathlib import Path
if "--version" in sys.argv:
    print("emapper-2.1.12 / Installed eggNOG DB version: 5.0.2")
    raise SystemExit(0)
args=sys.argv
source=Path(args[args.index("-i")+1])
out=Path(args[args.index("--output_dir")+1])
queries=[line[1:].split()[0] for line in source.read_text().splitlines() if line.startswith(">")]
assert "--dbmem" not in args and args[args.index("--block_size")+1]=="0.5"
# One synthetic assignment; all other queries deliberately have no hit.
(out/"annotation.emapper.annotations").write_text(
    "#query\\tseed_ortholog\\tKEGG_ko\\tEC\\n"+queries[0]+"\\tTEST_ONLY\\tko:K00001\\t-\\n")
"""


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campi", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    tool = out / "synthetic_contract_mapper.py"
    tool.write_text(MAPPER)
    data = out / "synthetic_contract_databases"
    data.mkdir()
    for name in (
        "eggnog.db",
        "eggnog.taxa.db",
        "eggnog.taxa.db.traverse.pkl",
        "eggnog_proteins.dmnd",
    ):
        (data / name).write_text("Synthetic contract fixture; not a biological database.\n")
    settings = out / "settings.json"
    settings.write_text(
        json.dumps(
            {
                "emapper": str(tool),
                "emapper_python": sys.executable,
                "data_dir": str(data),
                "threads": 1,
            }
        )
    )
    workflow = out / "workflow"
    argv = [
        sys.executable,
        str(ROOT / "tools/run_manifest.py"),
        "--manifest",
        str(ROOT / "examples/campi/inputs/samples.tsv"),
        "--lake",
        str(args.campi.resolve() / "lake.fasta"),
        "--out",
        str(workflow),
        "--threads",
        "2",
        "--annotation-settings",
        str(settings),
    ]
    with (out / "workflow.log").open("x") as log:
        subprocess.run(argv, stdout=log, stderr=subprocess.STDOUT, check=True)
    completed = json.loads((workflow / "COMPLETE.json").read_text())
    assert completed["annotation"] == "annotation/STUDY_COMPLETE.json"
    annotation = json.loads((workflow / completed["annotation"]).read_text())
    assert annotation["status"] == "PASS"
    assert annotation["unknown_sequences"] == annotation["unique_sequences"] - 1
    with (workflow / "annotation/annotations.tsv").open() as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    assert sum(row["annotation_status"] == "assigned" for row in rows) == 1
    assert sum(row["annotation_status"] == "no_hit" for row in rows) == len(rows) - 1
    # The extra stage must not change the search/grouping result.
    assert (workflow / "study_groups/study_group_matrix.tsv").read_bytes() == (
        args.campi / "workflow/study_groups/study_group_matrix.tsv"
    ).read_bytes()
    (out / "VALIDATION.json").write_text(
        json.dumps(
            {
                "status": "PASS",
                "scope": "Actual CAMPI searches; synthetic mapper contract only",
                "real_eggnog_inference": False,
                "biological_result": False,
                "post_search_annotation": True,
                "unchanged_study_matrix": True,
            },
            indent=2,
        )
        + "\n"
    )
    print("PASS: automatic annotation wiring; synthetic mapper, unchanged real search matrix")


if __name__ == "__main__":
    main()
