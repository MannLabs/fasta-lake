#!/usr/bin/env python3
"""Create the CAMPI teaching inputs from original files, preserving measured arrays.

The manifest is tab-separated: sample, predictions, mzml. This preparation tool
uses only the Python standard library. The runnable example already includes its
outputs; recreating them requires the original acquisition files and predictions.
"""

from __future__ import annotations

import argparse
import base64
import csv
import hashlib
import json
import math
import re
import xml.etree.ElementTree as ET
from pathlib import Path

URI = "http://psi.hupo.org/ms/mzml"
NS = {"m": URI}
ET.register_namespace("", URI)
ET.register_namespace("xsi", "http://www.w3.org/2001/XMLSchema-instance")


def sha(path, algorithm="sha256"):
    h = hashlib.new(algorithm)
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def info(spectrum):
    native = spectrum.attrib["id"]
    match = re.search(r"\bscan=(\d+)", native)
    if not match:
        raise ValueError("Expected Thermo scan= native identifiers")
    level = int(spectrum.find("m:cvParam[@accession='MS:1000511']", NS).get("value"))
    node = spectrum.find("m:scanList/m:scan/m:cvParam[@accession='MS:1000016']", NS)
    rt = float(node.get("value"))
    if node.get("unitAccession") == "UO:0000010":
        rt /= 60
    elif node.get("unitAccession") != "UO:0000031":
        raise ValueError("Unknown retention-time unit")
    precursor = spectrum.find("m:precursorList/m:precursor", NS)
    ion = spectrum.find(".//m:selectedIon/m:cvParam[@accession='MS:1000744']", NS)
    return {
        "native_id": native,
        "scan": int(match.group(1)),
        "ms_level": level,
        "rt_min": rt,
        "precursor_ref": precursor.get("spectrumRef", "") if precursor is not None else "",
        "precursor_mz": float(ion.get("value")) if ion is not None else None,
    }


def spectra(path):
    """Stream spectra without retaining cleared spectrum shells in their parent."""
    parent = None
    for event, element in ET.iterparse(path, events=("start", "end")):
        if event == "start" and element.tag == f"{{{URI}}}spectrumList":
            parent = element
        if event == "end" and element.tag == f"{{{URI}}}spectrum":
            yield element
            if parent is not None:
                parent.remove(element)
            element.clear()


def signal_sha(spectrum):
    digest = hashlib.sha256()
    for node in spectrum.findall("m:binaryDataArrayList/m:binaryDataArray/m:binary", NS):
        data = base64.b64decode(node.text or "")
        digest.update(len(data).to_bytes(8, "big"))
        digest.update(data)
    return digest.hexdigest()


def subset(mzml, predictions, out, sample, start, end):
    selected = {}
    references = set()
    for spectrum in spectra(mzml):
        row = info(spectrum)
        if row["ms_level"] in (1, 2) and start <= row["rt_min"] < end:
            selected[row["native_id"]] = row
            if row["precursor_ref"]:
                references.add(row["precursor_ref"])
    wanted = set(selected) | references
    if not selected:
        raise ValueError("The requested retention-time interval contains no spectra")

    # Preserve the original metadata and write valid non-indexed mzML. Native scan
    # identifiers stay unchanged; the sequential XML index is regenerated.
    with mzml.open() as stream:
        head = stream.read(65536)
    first = re.search(r"<mzML\b", head)
    marker = re.search(r"<spectrumList\b[^>]*>", head)
    if first is None or marker is None:
        raise ValueError("Could not locate the mzML header and spectrumList")
    prefix = head[first.start() : marker.start()]
    if 'xmlns="' not in prefix.split(">", 1)[0]:
        prefix = prefix.replace("<mzML ", f'<mzML xmlns="{URI}" ', 1)
    opening = re.sub(r'\bcount="\d+"', f'count="{len(wanted)}"', marker.group())
    name = sample + f"_{start:g}-{end:g}min.mzML"
    destination = out / name
    written = []
    with destination.open("x") as stream:
        stream.write('<?xml version="1.0" encoding="utf-8"?>\n' + prefix + opening + "\n")
        for spectrum in spectra(mzml):
            if spectrum.attrib["id"] not in wanted:
                continue
            row = info(spectrum)
            row["signal_sha256"] = signal_sha(spectrum)
            row["inside_requested_window"] = start <= row["rt_min"] < end
            spectrum.set("index", str(len(written)))
            stream.write(ET.tostring(spectrum, encoding="unicode"))
            written.append(row)
        stream.write("\n</spectrumList>\n</run>\n</mzML>\n")
    if {r["native_id"] for r in written} != wanted:
        raise ValueError("A precursor reference did not resolve to a retained spectrum")
    observed = [info(s) for s in spectra(destination)]
    if len(observed) != len(written):
        raise ValueError("Reparsed spectrum count differs from the written count")
    if not all(
        signal_sha(s) == row["signal_sha256"] for s, row in zip(spectra(destination), written)
    ):
        raise ValueError("Measured binary-array payload changed during subsetting")

    ms2 = {r["scan"]: r for r in written if r["ms_level"] == 2 and r["inside_requested_window"]}
    fields = ["sequence", "score", "spec_idx", "prec_charge", "prec_mz", "beam_rank"]
    prediction_name = sample + "_predictions.csv"
    scans = set()
    errors = []
    with (
        predictions.open(newline="") as source,
        (out / prediction_name).open("x", newline="") as target,
    ):
        writer = csv.DictWriter(target, fieldnames=fields)
        writer.writeheader()
        for row in csv.DictReader(source):
            scan = int(row["spec_idx"])
            if scan not in ms2:
                continue
            if scan in scans:
                raise ValueError("This example expects one prediction row per MS2 scan")
            scans.add(scan)
            mass = ms2[scan]["precursor_mz"]
            error = abs(float(row["prec_mz"]) - mass) / mass * 1e6
            if error > 0.2:
                raise ValueError("Prediction and mzML precursor m/z disagree")
            errors.append(error)
            writer.writerow(
                {
                    "sequence": row["peptide_prediction_detokenized_unmodified"],
                    **{key: row[key] for key in fields if key != "sequence"},
                }
            )
    if scans != set(ms2):
        raise ValueError("Prediction rows do not cover every selected MS2 scan exactly once")
    with (out / (sample + "_spectra.tsv")).open("x", newline="") as target:
        writer = csv.DictWriter(target, fieldnames=list(written[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(written)
    raw_sha = re.search(r'accession="MS:1000569" value="([a-f0-9]+)"', head)
    return {
        "sample": sample,
        "mzml": name,
        "predictions": prediction_name,
        "window_minutes": [start, end],
        "window_interval": "start inclusive, end exclusive",
        "ms1_spectra": sum(r["ms_level"] == 1 for r in written),
        "ms2_spectra": len(ms2),
        "extra_precursor_spectra": sum(not r["inside_requested_window"] for r in written),
        "prediction_rows": len(scans),
        "max_precursor_pairing_error_ppm": max(errors),
        "source_mzml_name": mzml.name,
        "source_mzml_sha256": sha(mzml),
        "source_predictions_name": predictions.name,
        "source_predictions_sha256": sha(predictions),
        "source_raw_sha1_from_mzml": raw_sha.group(1) if raw_sha else None,
        "output_mzml_bytes": destination.stat().st_size,
        "output_mzml_sha256": sha(destination),
        "prediction_sha256": sha(out / prediction_name),
        "modifications": (
            "Retained measured binary arrays and native scan IDs; rebuilt sequential indices; "
            "omitted original chromatogram and index sections. Added any MS1 precursor "
            "referenced from the selected interval."
        ),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--start-minutes", type=float, default=45)
    parser.add_argument("--end-minutes", type=float, default=50)
    args = parser.parse_args()
    if not 0 <= args.start_minutes < args.end_minutes or not math.isfinite(args.end_minutes):
        parser.error("Require finite 0 <= start-minutes < end-minutes")
    with args.manifest.open() as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if not {"sample", "predictions", "mzml"}.issubset(reader.fieldnames or []):
            parser.error("Manifest requires sample, predictions and mzml columns")
        rows = list(reader)
    if len(rows) != 2 or {row["sample"] for row in rows} != {"S01", "S02"}:
        parser.error("This preparation script requires exactly the CAMPI S01 and S02 inputs")
    for row in rows:
        for key in ("predictions", "mzml"):
            path = Path(row[key])
            if not path.is_absolute():
                path = args.manifest.resolve().parent / path
            if not path.is_file():
                parser.error(f"Missing {key} for {row['sample']}: {path}")
            row[key] = path
    if not args.reference.is_file():
        parser.error("Reference FASTA is missing")
    args.out.mkdir(parents=True, exist_ok=False)
    results = [
        subset(
            Path(r["mzml"]),
            Path(r["predictions"]),
            args.out,
            r["sample"],
            args.start_minutes,
            args.end_minutes,
        )
        for r in rows
    ]
    targets = decoys = 0
    include = False
    with args.reference.open() as source, (args.out / "reference.fasta").open("x") as target:
        for line in source:
            if line.startswith(">"):
                include = not line[1:].startswith("DECOY_")
                targets += int(include)
                decoys += int(not include)
            if include:
                target.write(line)
    with (args.out / "samples.tsv").open("x", newline="") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=["sample", "predictions", "mzml"], delimiter="\t"
        )
        writer.writeheader()
        writer.writerows({key: r[key] for key in writer.fieldnames} for r in results)
    provenance = {
        "project": "PXD023217",
        "samples": results,
        "reference": {
            "source_name": args.reference.name,
            "source_sha256": sha(args.reference),
            "source_sha1": sha(args.reference, "sha1"),
            "target_records": targets,
            "excluded_DECOY_records": decoys,
            "sha256": sha(args.out / "reference.fasta"),
            "policy": (
                "All original target records retained; only deposited DECOY_ records excluded. "
                "The example runner normalizes and hashes sequences before selection; "
                "SAGE generates new decoys."
            ),
        },
    }
    (args.out / "input_provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    paths = sorted(p for p in args.out.iterdir() if p.is_file())
    (args.out / "SHA256SUMS").write_text("".join(sha(p) + "  " + p.name + "\n" for p in paths))
    print(json.dumps(provenance, indent=2))


if __name__ == "__main__":
    main()
