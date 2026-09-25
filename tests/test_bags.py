"""Native AlphaNovo matching parity and original-evidence selection accounting."""

import csv
import gzip
import json
import os
import subprocess
from pathlib import Path

import pytest

from fasta_lake.bags import (
    bag_plan,
    compile_queries,
    control_form,
    fnv1a,
    norm,
    select_databases,
    tokens,
    upstream,
)


@pytest.fixture
def native():
    source = os.environ.get("FASTALAKE_ALPHANOVO_SOURCE")
    if not source:
        pytest.skip("Pinned external AlphaNovo checkout required")
    return source, upstream(source)


def test_streaming_matches_native_all_protein_sets(tmp_path, native):
    source, (mass, bag) = native
    engine = Path(os.environ["FASTALAKE_BAG_EXTRACT"])
    sequences = ["PEPTIDERK", "PEPNTIDERK", "PEPQTIDERK", "PEPCMTIDERK"]
    patterns = ["PE[PT]IDERK", "PEP[NT]IDERK", "PEPQTIDERK", "PEPCMTIDERK"]
    modified = sequences[:3] + ["PEP[Carbamidomethyl@C][Oxidation@M]TIDERK"]
    fields = [
        "peptide_prediction_detokenized_unmodified",
        "peptide_prediction_detokenized",
        "score",
        "beam_rank",
        "spec_idx",
        "sequence_with_bags",
    ]
    path = tmp_path / "predictions.csv"
    with path.open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for i, (seq, pattern, mod) in enumerate(zip(sequences, patterns, modified)):
            writer.writerow(dict(zip(fields, [seq, mod, 0.99, 0, i + 1, pattern])))
    receipt = compile_queries(path, source, tmp_path)
    proteins = {
        "exact": "PEPTIDERK",
        "swapped": "PETPIDERK",
        "longer": "XXPEPGGTIDERKXX",
        "q_to_ag": "PEPAGTIDERK",
        "not_a_reference_explanation": "WWWWWWWWWW",
        "modified": "PEPCGMTIDERK",
        "wrong_order": "PEDTIPERK",
        "unknown": "PEPXTIDERK",
        "boundary_left": "PEP",
        "boundary_right": "TIDERK",
    }
    fasta = tmp_path / "reservoir.fasta"
    fasta.write_text("".join(">" + name + "\n" + seq + "\n" for name, seq in proteins.items()))
    run = subprocess.run(
        [str(engine), str(tmp_path / "queries.jsonl"), str(fasta), str(tmp_path / "scan.json")],
        capture_output=True,
        text=True,
        check=True,
    )
    observed = [json.loads(line) for line in run.stdout.splitlines()]
    lookup = {record["header"]: dict(record["evidence"]) for record in observed}
    table = {
        key: tuple(v for v in values if (key, v) not in {("D", "N"), ("E", "Q")})
        for key, values in mass.build_swap_table().items()
    }
    for index, (seq, pattern, mod) in enumerate(zip(sequences, patterns, modified)):
        uid = index * 3
        expected_exact = {name for name, protein in proteins.items() if norm(seq) in norm(protein)}
        expected_bag = bag.anchor_bag_search([pattern], proteins, find_any=False)[1].get(
            pattern, set()
        )
        alternatives = mass.generate_variants(seq.replace("I", "L").replace("L", "I"), table, 2)
        if "[Carbamidomethyl@C]" in mod:
            alternatives |= mass.generate_modified_mass_equivalent_variants(mod)
        alternatives = {norm(v) for v in alternatives} - {norm(seq)}
        expected_mass = {
            name
            for name, protein in proteins.items()
            if any(v in norm(protein) for v in alternatives)
        }
        for bit, expected in [(1, expected_exact), (2, expected_bag), (4, expected_mass)]:
            actual = {name for name, evidence in lookup.items() if evidence.get(uid, 0) & bit}
            assert actual == expected
    with gzip.open(tmp_path / "candidates.jsonl.gz", "wt") as handle:
        handle.write(run.stdout)
    selected = select_databases(tmp_path / "candidates.jsonl.gz", tmp_path, receipt["query_units"])
    assert selected["arms"]["combined"]["assigned_original_queries"] <= len(sequences)
    assert json.loads((tmp_path / "scan.json").read_text())["reservoir_records"] == len(proteins)


def test_variant_spelling_does_not_multiply_support(tmp_path):
    records = [
        {"header": "A", "sequence": "PEPTIDERK", "lookup_routes": 7, "evidence": [[0, 7], [3, 4]]},
        {"header": "B", "sequence": "PETPIDERK", "lookup_routes": 7, "evidence": [[0, 7]]},
    ]
    path = tmp_path / "candidates.jsonl.gz"
    with gzip.open(path, "wt") as handle:
        for record in records:
            handle.write(json.dumps(record) + "\n")
    result = select_databases(path, tmp_path, 6)
    assert result["arms"]["combined"]["selected_proteins"] == 1
    assert result["arms"]["combined"]["assigned_original_queries"] == 2
    assert (tmp_path / "combined.fasta").read_text().startswith(">A\n")


def test_bag_span_and_fnv_contract(native):
    _, (_, bag) = native
    with pytest.raises(ValueError, match="does not describe"):
        bag_plan("PE[PP]IDERK", "PEPTLDERK", bag)
    assert bag_plan("[PEPTLDERK]", "PEPTLDERK", bag) is None
    assert fnv1a("") == 14695981039346656037
    assert fnv1a("a") == 0xAF63DC4C8601EC8C


def test_terminal_modifications_do_not_become_residues_or_move(native):
    _, (_, bag) = native
    modified = "[Acetyl@Protein_N-term][Oxidation@M]AHVEADYEK"
    units, prefix, suffix = tokens(modified, "MAHVEADYEK")
    assert len(units) == 10
    assert prefix == "[Acetyl@Protein_N-term]" and suffix == ""
    sequence, pattern, shuffled = control_form("MAHVEADYEK", "MAHVEADYEK", modified, 3, bag)
    assert shuffled.startswith(prefix)
    shuffled_units, _, _ = tokens(shuffled, sequence)
    assert sorted(shuffled_units) == sorted(units)
    assert pattern == sequence


def test_native_cam_expansion_retains_terminal_input_semantics(tmp_path, native):
    source, (mass, _) = native
    fields = [
        "peptide_prediction_detokenized_unmodified",
        "peptide_prediction_detokenized",
        "score",
        "beam_rank",
        "spec_idx",
        "sequence_with_bags",
    ]
    modified = "[Acetyl@Protein_N-term]PEP[Carbamidomethyl@C]TIDERK"
    path = tmp_path / "predictions.csv"
    with path.open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerow(dict(zip(fields, ["PEPCTIDERK", modified, 0.99, 0, 1, "PEPCTIDERK"])))
    compile_queries(path, source, tmp_path)
    seeds = {
        r["seed"]
        for r in map(json.loads, (tmp_path / "queries.jsonl").read_text().splitlines())
        if r["query"] == 0 and r["route"] == 4
    }
    assert {norm(s) for s in mass.generate_modified_mass_equivalent_variants(modified)} <= seeds
