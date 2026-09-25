"""Unit tests for fasta_lake.provenance (v5.1)."""

import json
from pathlib import Path

from fasta_lake.provenance import summarise_provenance, write_provenance


def test_write_provenance_basic(tmp_path: Path):
    out = tmp_path / "prov.jsonl"
    selected = {"sp|P04637|P53_HUMAN", "MGYG000000001_00001"}
    prot2pep = {
        "sp|P04637|P53_HUMAN": {"MEEPQSDPSVEPPLSQ", "ETPPLGE"},
        "MGYG000000001_00001": {"AAAAAAAAAAK"},
    }
    peptide_scores = {
        "MEEPQSDPSVEPPLSQ": 0.97,
        "ETPPLGE": 0.88,
        "AAAAAAAAAAK": 0.55,
    }
    n = write_provenance(
        out,
        selected,
        prot2pep,
        peptide_scores,
        strategy="species_budget",
    )
    assert n == 2

    lines = out.read_text().strip().split("\n")
    assert len(lines) == 2
    records = [json.loads(ln) for ln in lines]
    # Sorted by protein_id
    assert records[0]["protein_id"] == "MGYG000000001_00001"
    assert records[1]["protein_id"] == "sp|P04637|P53_HUMAN"

    p53 = records[1]
    assert p53["n_peptides"] == 2
    assert set(p53["peptides"]) == {"MEEPQSDPSVEPPLSQ", "ETPPLGE"}
    assert p53["best_score"] == 0.97
    assert p53["admission"] == "exact"
    assert p53["strategy"] == "species_budget"


def test_write_provenance_kmer_rescue(tmp_path: Path):
    out = tmp_path / "prov.jsonl"
    selected = {"PROT_A", "PROT_B"}
    prot2pep = {"PROT_A": {"PEPTIDE"}, "PROT_B": {"OTHERPEP"}}
    scores = {"PEPTIDE": 0.8, "OTHERPEP": 0.7}
    n = write_provenance(
        out,
        selected,
        prot2pep,
        scores,
        strategy="razor",
        kmer_rescued_proteins={"PROT_B"},
    )
    assert n == 2

    records = [json.loads(ln) for ln in out.read_text().strip().split("\n")]
    by_id = {r["protein_id"]: r for r in records}
    assert by_id["PROT_A"]["admission"] == "exact"
    assert by_id["PROT_B"]["admission"] == "kmer_rescue"


def test_write_provenance_missing_peptides(tmp_path: Path):
    """Protein selected but not in prot2pep map — defensive path."""
    out = tmp_path / "prov.jsonl"
    n = write_provenance(
        out,
        {"ORPHAN"},
        {},
        {},
        strategy="uniform",
    )
    assert n == 1
    rec = json.loads(out.read_text().strip())
    assert rec["n_peptides"] == 0
    assert rec["peptides"] == []
    assert rec["best_score"] is None


def test_write_provenance_truncates_large_peptide_lists(tmp_path: Path):
    out = tmp_path / "prov.jsonl"
    peps = {f"PEPTIDE{i:03d}K" for i in range(50)}
    prot2pep = {"BIG": peps}
    scores = {p: 0.5 + i * 0.01 for i, p in enumerate(sorted(peps))}
    write_provenance(
        out,
        {"BIG"},
        prot2pep,
        scores,
        strategy="razor",
    )
    rec = json.loads(out.read_text().strip())
    assert rec["n_peptides"] == 50
    assert len(rec["peptides"]) == 20
    assert rec["peptides_truncated"] is True


def test_write_provenance_contaminant_flag(tmp_path: Path):
    """Trypsin (P00761) should be flagged as contaminant."""
    out = tmp_path / "prov.jsonl"
    write_provenance(
        out,
        {"sp|P00761|TRYP_PIG", "sp|P04637|P53_HUMAN"},
        {"sp|P00761|TRYP_PIG": {"PEP"}, "sp|P04637|P53_HUMAN": {"OTHER"}},
        {"PEP": 0.9, "OTHER": 0.9},
        strategy="razor",
    )
    records = [json.loads(ln) for ln in out.read_text().strip().split("\n")]
    by_id = {r["protein_id"]: r for r in records}
    assert by_id["sp|P00761|TRYP_PIG"]["contaminant"] is True
    assert by_id["sp|P04637|P53_HUMAN"]["contaminant"] is False


def test_summarise_provenance(tmp_path: Path):
    out = tmp_path / "prov.jsonl"
    write_provenance(
        out,
        {"A", "B", "C"},
        {"A": {"P1", "P2"}, "B": {"P3"}, "C": {"P4", "P5", "P6"}},
        {"P1": 0.9, "P2": 0.8, "P3": 0.5, "P4": 0.7, "P5": 0.6, "P6": 0.85},
        strategy="razor",
        kmer_rescued_proteins={"B"},
    )
    summary = summarise_provenance(out)
    assert summary["n_proteins"] == 3
    assert summary["n_kmer_rescued"] == 1
    assert summary["proteins_with_1_peptide"] == 1
    assert summary["proteins_with_2plus_peptides"] == 2
