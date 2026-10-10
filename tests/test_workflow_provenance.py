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


def test_manifest_accepts_gzipped_predictions_and_spectra(tmp_path):
    import gzip

    import pytest

    from fasta_lake.workflow import load_manifest

    with gzip.open(tmp_path / "S.csv.gz", "wt") as stream:
        stream.write("sequence,score\nMPEPTIDEAK,0.9\n")
    (tmp_path / "S.mzML.gz").write_bytes(gzip.compress(b"<mzML/>"))
    manifest = tmp_path / "samples.tsv"
    manifest.write_text("sample\tpredictions\tmzml\nS\tS.csv.gz\tS.mzML.gz\n")
    row = load_manifest(manifest)[0]
    assert row["predictions"].endswith("S.csv.gz") and row["mzml"].endswith("S.mzML.gz")
    # A gzipped file must still carry the prediction columns.
    with gzip.open(tmp_path / "S.csv.gz", "wt") as stream:
        stream.write("peptide,confidence\nMPEPTIDEAK,0.9\n")
    with pytest.raises(ValueError, match="expected AlphaNovo or generic"):
        load_manifest(manifest)
    (tmp_path / "S.txt").write_text("sequence,score\n")
    manifest.write_text("sample\tpredictions\tmzml\nS\tS.txt\tS.mzML.gz\n")
    with pytest.raises(ValueError, match="predictions CSV"):
        load_manifest(manifest)
