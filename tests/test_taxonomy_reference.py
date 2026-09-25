"""Keep custom Unipept references sequence-exact and taxonomy-explicit."""

import json

import pytest

from fasta_lake.taxonomy_reference import prepare_unipept_reference


def test_aliases_and_il_variants_stay_distinct(tmp_path):
    fasta = tmp_path / "proteins.fasta"
    fasta.write_text(">UHGG_A\nACDEILK\n>ALIAS_A\nACDEILK\n>UHGG_B\nACDELLK\n")
    mapping = tmp_path / "taxonomy.tsv"
    mapping.write_text("accession\tncbi_taxon_id\nUHGG_A\t2\n")
    out = tmp_path / "out"
    r = prepare_unipept_reference(fasta, out, taxonomy=mapping)
    assert r["proteins"] == 3 and r["mapped_proteins"] == 1 and r["unresolved_proteins"] == 2
    assert (out / "proteins.tsv").read_text().splitlines() == [
        "UHGG_A\t2\tACDEILK\t",
        "ALIAS_A\t0\tACDEILK\t",
        "UHGG_B\t0\tACDELLK\t",
    ]
    assert json.loads((out / "COMPLETE.json").read_text())["experimental"]


@pytest.mark.parametrize(
    "rows",
    [
        "UHGG_A\t2\nUHGG_A\t2\n",
        "foreign\t2\n",
        "UHGG_A\t0\n",
        "UHGG_A\tGTDB:s__Species\n",
        "UHGG_A\t4294967296\n",
    ],
)
def test_no_duplicate_foreign_or_fabricated_taxonomy(tmp_path, rows):
    fasta = tmp_path / "p.fasta"
    fasta.write_text(">UHGG_A\nACDEILK\n")
    mapping = tmp_path / "taxa.tsv"
    mapping.write_text("accession\tncbi_taxon_id\n" + rows)
    with pytest.raises(ValueError):
        prepare_unipept_reference(fasta, tmp_path / "out", taxonomy=mapping)
    assert not (tmp_path / "out").exists()


def test_unmapped_proteins_are_retained_without_an_invented_taxon(tmp_path):
    fasta = tmp_path / "p.fasta"
    fasta.write_text(">MGYG000000001_1\nACDEILK\n")
    out = tmp_path / "out"
    r = prepare_unipept_reference(fasta, out)
    assert r["unresolved_proteins"] == 1
    assert "\t0\t" in (out / "proteins.tsv").read_text()
    with pytest.raises(FileExistsError):
        prepare_unipept_reference(fasta, out)
