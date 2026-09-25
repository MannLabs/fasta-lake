"""
End-to-end test for byte-identical-collapse tracing.

Flow:
  1. Write two tiny FASTA files — SwissProt-style and TrEMBL-style — where
     ONE sequence is byte-identical across both.
  2. Run build_lake() on them → lake.fasta + lake.provenance.jsonl.
  3. Read the sidecar back:
       - the shared sequence must show up once (dedup correct)
       - its record must have n_sources=2, alt_accessions populated
  4. Run `fasta-lake lake-inspect --duplicates` on the sidecar. The filter
     must return exactly the shared record.

This is the "user asks 'where did my paralog go?'" path end-to-end. If any
link breaks — dedup, sidecar derivation, or CLI filter — this test fails.
"""

from __future__ import annotations

import json
from pathlib import Path

from click.testing import CliRunner

from fasta_lake.cli import cli
from fasta_lake.lake_builder import build_lake
from fasta_lake.lake_merge import read_provenance_sidecar

SHARED_SEQ = "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAP"
SP_UNIQUE = "MASEFKKKLSVDLTAVEKACR"
TR_UNIQUE = "MKQSTIALALLPLLFTPVTKA"


def _write_fasta(path: Path, entries: list[tuple[str, str]]) -> None:
    with open(path, "w") as f:
        for header, seq in entries:
            f.write(f">{header}\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i : i + 60] + "\n")


def test_duplicate_trace_end_to_end(tmp_path: Path) -> None:
    # ---- 1. Write two sources — SP and TR — with ONE byte-identical sequence.
    sp_fasta = tmp_path / "sp.fasta"
    tr_fasta = tmp_path / "tr.fasta"

    _write_fasta(
        sp_fasta,
        [
            (
                "sp|P04637|P53_HUMAN Cellular tumor antigen p53 "
                "OS=Homo sapiens OX=9606 GN=TP53 PE=1 SV=4",
                SHARED_SEQ,
            ),
            (
                "sp|P00000|UNIQ_HUMAN Something unique to SwissProt "
                "OS=Homo sapiens OX=9606 GN=UNIQSP PE=1 SV=1",
                SP_UNIQUE,
            ),
        ],
    )
    _write_fasta(
        tr_fasta,
        [
            # Same sequence as SP|P04637, different TrEMBL accession.
            (
                "tr|H8LKR3|H8LKR3_HUMAN Cellular tumor antigen p53 "
                "OS=Homo sapiens OX=9606 GN=TP53 PE=2 SV=1",
                SHARED_SEQ,
            ),
            (
                "tr|Q9ZZZ0|Q9ZZZ0_HUMAN Something unique to TrEMBL "
                "OS=Homo sapiens OX=9606 GN=UNIQTR PE=2 SV=1",
                TR_UNIQUE,
            ),
        ],
    )

    # ---- 2. Build the lake.
    out_fasta = tmp_path / "lake.fasta"
    out_prov = tmp_path / "lake.provenance.jsonl"
    stats = build_lake(
        sources={"SP": sp_fasta, "TR": tr_fasta},
        out_fasta=out_fasta,
        out_provenance=out_prov,
    )

    assert stats["n_input_sequences"] == 4
    assert stats["n_unique_sequences"] == 3  # SP_UNIQUE, TR_UNIQUE, SHARED
    assert stats["n_shared_across_sources"] == 1

    # ---- 3. Sidecar records trace the collapse.
    recs = read_provenance_sidecar(out_prov)
    shared_rec = next(r for r in recs if len(r.sources) > 1)
    assert shared_rec.primary_tag == "SP"
    assert shared_rec.primary_accession == "P04637"

    # Re-open the raw JSONL to check the derived fields are emitted.
    lines = [json.loads(ln) for ln in out_prov.read_text().splitlines() if ln.strip()]
    shared_json = next(d for d in lines if d["n_sources"] > 1)
    assert shared_json["n_sources"] == 2
    assert shared_json["primary_accession"] == "P04637"
    assert shared_json["alt_accessions"] == ["H8LKR3"]

    # Singleton records must have empty alt_accessions.
    for d in lines:
        if d["n_sources"] == 1:
            assert d["alt_accessions"] == []

    # ---- 4. CLI filter `--duplicates` returns exactly the collapsed record.
    runner = CliRunner()
    result = runner.invoke(cli, ["lake-inspect", str(out_prov), "--duplicates"])
    assert result.exit_code == 0, result.output
    # The output is a JSON blob; the shared primary must appear, the unique
    # ones must not.
    assert "P04637" in result.output
    assert '"n_sources": 2' in result.output
    assert "H8LKR3" in result.output
    assert "UNIQ_HUMAN" not in result.output
    assert "Q9ZZZ0" not in result.output

    # `--summary` flag shows the aggregate count the user is asked to grep.
    result = runner.invoke(cli, ["lake-inspect", str(out_prov), "--summary"])
    assert result.exit_code == 0, result.output
    assert "Shared (≥2 src):" in result.output
