"""PROVENANCE.json records the settings the run actually used."""

from argparse import Namespace

from fasta_lake.workflow import provenance_record


def test_provenance_records_the_evidence_top_fraction(tmp_path):
    manifest = tmp_path / "samples.tsv"
    manifest.write_text("sample\tpredictions\tmzml\n")
    args = Namespace(
        manifest=manifest, lake=None, background=None, lake_headers=None, top_fraction=0.25
    )
    record = provenance_record(args, [], {})
    assert record["inference"]["top_fraction"] == 0.25
    assert record["parameters"]["top_fraction"] == 0.25
