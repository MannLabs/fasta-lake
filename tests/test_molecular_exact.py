import csv
import json

import pytest

from fasta_lake.multi_omics.complementarity import compare
from fasta_lake.multi_omics.exact import build_union, read_fasta, read_tpm, sha256, top_cutoff


def _tsv_rows(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def source(tmp_path):
    (tmp_path / "genes.faa").write_text(
        ">g1 original\nMACDEK\n>g2 duplicate\nMACDEK\n>g3 distinct I\nMILK\n>g4 distinct L\nMLLK\n"
    )
    (tmp_path / "tpm.tsv").write_text(
        "prodigal_protein\tmetaG_tpm\tmetaT_tpm\ng1\t1\tNA\ng2\t10\tNA\ng3\t0\t1\ng4\t0\t1.0001\n"
    )
    (tmp_path / "provenance.json").write_text('{"source":"synthetic test fixture"}\n')
    obj = dict(
        schema="fastalake.molecular-source.v1",
        specimen="s1",
        reference_id="r1",
        reference_status="verified_original_quantification_reference",
    )
    for key, name in [
        ("tpm", "tpm.tsv"),
        ("proteins", "genes.faa"),
        ("provenance", "provenance.json"),
    ]:
        obj[key] = dict(path=name, sha256=sha256(tmp_path / name))
    path = tmp_path / "source.json"
    path.write_text(json.dumps(obj))
    (tmp_path / "baseline.fasta").write_bytes(b">base preserve this header\r\nMQQQQK\r\n")
    return path


def test_gene_threshold_precedes_dedup_and_rna_is_strict(tmp_path):
    manifest = source(tmp_path)
    result = build_union(
        tmp_path / "baseline.fasta",
        manifest,
        "s1",
        tmp_path / "out",
        dna_percent=50,
        rna_absolute=1,
    )
    sequences = {r[1] for r in read_fasta(tmp_path / "out/search.fasta").values()}
    assert sequences == {"MQQQQK", "MACDEK", "MLLK"}
    assert result["additions"] == 2
    rows = _tsv_rows(tmp_path / "out/gene_evidence.tsv")
    assert len(rows) == 4 and rows[0]["dna_pass"] == "0" and rows[1]["dna_pass"] == "1"
    assert rows[2]["rna_pass"] == "0" and rows[3]["rna_pass"] == "1"
    assert rows[0]["included"] == "1"
    assert rows[0]["final_accession"] == rows[1]["final_accession"]
    assert (
        (tmp_path / "out/search.fasta")
        .read_bytes()
        .startswith((tmp_path / "baseline.fasta").read_bytes())
    )


def test_disabled_layers_preserve_baseline_bytes(tmp_path):
    manifest = source(tmp_path)
    result = build_union(
        tmp_path / "baseline.fasta", manifest, "s1", tmp_path / "out", dna_percent=0
    )
    assert result["additions"] == 0
    assert (tmp_path / "out/search.fasta").read_bytes() == (
        tmp_path / "baseline.fasta"
    ).read_bytes()


def test_all_zero_and_rank_boundary():
    assert top_cutoff([0, None, 0], 100) is None
    assert top_cutoff(range(1, 102), 1, "floor") == 101
    assert top_cutoff(range(1, 102), 1, "ceil") == 100
    assert top_cutoff([1, 2, 2, 3], 50) == 2
    assert top_cutoff([3, 2, 2, 1], 100) == 1
    assert top_cutoff([1], 0) is None


@pytest.mark.parametrize("value", [-1, 101, float("nan"), float("inf")])
def test_invalid_percent(value):
    with pytest.raises(ValueError):
        top_cutoff([1], value)


def test_foreign_specimen_fails_before_output(tmp_path):
    manifest = source(tmp_path)
    with pytest.raises(ValueError, match="specimen mismatch"):
        build_union(tmp_path / "baseline.fasta", manifest, "s2", tmp_path / "out")
    assert not (tmp_path / "out").exists()


def test_changed_source_is_rejected(tmp_path):
    manifest = source(tmp_path)
    (tmp_path / "genes.faa").write_text(">g1\nMCHANGEDK\n")
    with pytest.raises(ValueError, match="checksum mismatch"):
        build_union(tmp_path / "baseline.fasta", manifest, "s1", tmp_path / "out")


def test_rehashed_missing_gene_still_rejected(tmp_path):
    manifest = source(tmp_path)
    (tmp_path / "genes.faa").write_text(">g1\nMACDEK\n")
    obj = json.loads(manifest.read_text())
    obj["proteins"]["sha256"] = sha256(tmp_path / "genes.faa")
    manifest.write_text(json.dumps(obj))
    with pytest.raises(ValueError, match="All-gene reference mismatch"):
        build_union(tmp_path / "baseline.fasta", manifest, "s1", tmp_path / "out", dna_percent=0)


def test_acc_collision_is_rejected(tmp_path):
    manifest = source(tmp_path)
    digest = read_fasta(tmp_path / "genes.faa")["g2"][2]
    (tmp_path / "baseline.fasta").write_text(">FMO_SHA256_" + digest + "\nMQQQQK\n")
    with pytest.raises(ValueError, match="collides"):
        build_union(tmp_path / "baseline.fasta", manifest, "s1", tmp_path / "out")


def test_fasta_duplicates_and_tpm_nonfinite_are_rejected(tmp_path):
    p = tmp_path / "f.fasta"
    p.write_text(">x\nMILK\n>x\nMLLK\n")
    with pytest.raises(ValueError, match="duplicate"):
        read_fasta(p)
    p = tmp_path / "t.tsv"
    p.write_text("prodigal_protein\tmetaG_tpm\tmetaT_tpm\nx\tnan\tNA\n")
    with pytest.raises(ValueError):
        read_tpm(p)


def searches(tmp_path):
    mzml = tmp_path / "spectra.mzML"
    mzml.write_text("synthetic spectrum fixture")
    for arm in ["left", "right"]:
        folder = tmp_path / arm
        folder.mkdir()
        seqs = ">a\nMACDEK\n>b\nMQQQQK\n>c\nMACDEKMQQQQK\n"
        if arm == "right":
            seqs = seqs.replace(">a", ">alias_a")
        (folder / "db.fasta").write_text(seqs)
        config = dict(
            database=dict(fasta="db.fasta", generate_decoys=True),
            mzml_paths=[str(mzml)],
            output_directory=str(folder),
        )
        (folder / "config.json").write_text(json.dumps(config))
        a = "a" if arm == "left" else "alias_a"
        rows = [["ACDEK", a + ";c", "1", ".001", "spectra.mzML"]]
        if arm == "right":
            rows.append(["QQQQK", "b;c", "1", ".001", "spectra.mzML"])
        rows.append(["QQQQK", "rev_b", "-1", ".0001", "spectra.mzML"])
        with (folder / "results.sage.tsv").open("w") as f:
            w = csv.writer(f, delimiter="\t", lineterminator="\n")
            w.writerow(["peptide", "proteins", "label", "peptide_q", "filename"])
            w.writerows(rows)
        lfq = [["ACDEK", "2", a, "100" if arm == "left" else "200"]]
        if arm == "right":
            lfq.append(["QQQQK", "2", "b", "300"])
        with (folder / "lfq.tsv").open("w") as f:
            w = csv.writer(f, delimiter="\t", lineterminator="\n")
            w.writerow(["peptide", "charge", "proteins", "spectra.mzML"])
            w.writerows(lfq)
    return tmp_path / "left", tmp_path / "right"


def test_exact_aliases_shared_peptides_and_quantitative_gains(tmp_path):
    a, b = searches(tmp_path)
    summary = compare(a, b, tmp_path / "comparison")
    assert summary["protein_cells"] == {"both": 2, "metaG_only": 0, "de_novo_only": 1, "neither": 0}
    assert summary["peptide_left_only"] == 0 and summary["peptide_right_only"] == 1
    assert summary["lfq_shared_features"] == 1
    assert summary["lfq_median_log2_right_over_left"] == 1
    rows = _tsv_rows(tmp_path / "comparison/protein_overlap.tsv")
    assert all(r["metaG_single_sequence_candidate_peptides"] == "0" for r in rows)


def test_different_search_settings_fail(tmp_path):
    a, b = searches(tmp_path)
    p = b / "config.json"
    cfg = json.loads(p.read_text())
    cfg["database"]["generate_decoys"] = False
    p.write_text(json.dumps(cfg))
    with pytest.raises(ValueError, match="settings differ"):
        compare(a, b, tmp_path / "out")


def test_false_peptide_protein_membership_fails(tmp_path):
    a, b = searches(tmp_path)
    p = b / "results.sage.tsv"
    p.write_text(p.read_text().replace("alias_a;c", "b"))
    with pytest.raises(ValueError, match="not contained"):
        compare(a, b, tmp_path / "out")


def test_duplicate_quant_feature_fails(tmp_path):
    a, b = searches(tmp_path)
    with (b / "lfq.tsv").open("a") as f:
        f.write("ACDEK\t2\talias_a\t200\n")
    with pytest.raises(ValueError, match="Duplicate accepted LFQ"):
        compare(a, b, tmp_path / "out")


def test_il_variant_lfq_rows_are_retained_and_excluded_from_agreement(tmp_path):
    a, b = searches(tmp_path)
    for folder in [a, b]:
        with (folder / "db.fasta").open("a") as f:
            f.write(">ile\nMILK\n>leu\nMLLK\n")
        with (folder / "results.sage.tsv").open("a") as f:
            f.write("ILK\tile;leu\t1\t.001\tspectra.mzML\n")
        with (folder / "lfq.tsv").open("a") as f:
            f.write("ILK\t2\tile\t100\nLLK\t2\tleu\t101\n")
    summary = compare(a, b, tmp_path / "out")
    assert summary["lfq_shared_exact_rows"] == 3
    assert summary["lfq_shared_features"] == 1
    assert summary["lfq_shared_rows_excluded_for_il_ambiguity"] == 2
    rows = _tsv_rows(tmp_path / "out/lfq_overlap.tsv")
    ambiguous = [r for r in rows if r["il_family_multiple_rows"] == "1"]
    assert len(ambiguous) == 2 and all(r["log2_right_over_left"] == "" for r in ambiguous)


@pytest.mark.parametrize("value", ["-0.1", "1.1", "nan"])
def test_invalid_protein_q_is_rejected(tmp_path, value):
    a, b = searches(tmp_path)
    for folder in (a, b):
        p = folder / "results.sage.tsv"
        lines = p.read_text().splitlines()
        lines[0] += "\tprotein_q"
        lines[1:] = [line + "\t" + value for line in lines[1:]]
        p.write_text("\n".join(lines) + "\n")
    with pytest.raises(ValueError):
        compare(a, b, tmp_path / "out", protein_q=0.01)


def test_custom_decoy_tag_and_novel_peptide_export(tmp_path):
    a, b = searches(tmp_path)
    for folder in (a, b):
        p = folder / "config.json"
        cfg = json.loads(p.read_text())
        cfg["database"]["decoy_tag"] = "DECOY_"
        p.write_text(json.dumps(cfg))
        with (folder / "db.fasta").open("a") as f:
            f.write(">DECOY_x\nMDECOYK\n")
    summary = compare(a, b, tmp_path / "out")
    assert summary["protein_universe"] == 3
    rows = _tsv_rows(tmp_path / "out/protein_overlap.tsv")
    row = next(r for r in rows if r["category"] == "de_novo_only")
    assert row["protein_sequence"] == "MQQQQK"
    assert row["de_novo_support_peptides_not_accepted_by_other_arm"] == "QQQQK"


def test_compressed_tables_and_receipt(tmp_path):
    import gzip

    a, b = searches(tmp_path)
    result = compare(a, b, tmp_path / "out", compress_tables=True)
    for name in result["tables"].values():
        path = tmp_path / "out" / name
        assert name.endswith(".tsv.gz")
        assert sha256(path) == result["output_sha256"][name]
        with gzip.open(path, "rt") as stream:
            reader = csv.DictReader(stream, delimiter="\t")
            assert reader.fieldnames
            if "excluded" not in name:
                assert list(reader)


def test_sage_combined_charge_sentinel_is_retained(tmp_path):
    a, b = searches(tmp_path)
    for folder in (a, b):
        p = folder / "lfq.tsv"
        p.write_text(p.read_text().replace("\t2\t", "\t-1\t"))
    result = compare(a, b, tmp_path / "out")
    assert result["lfq_shared_features"] == 1
    assert result["lfq_charge_modes"] == {"metaG": ["combined"], "de_novo": ["combined"]}


@pytest.mark.parametrize("charge", [0, -2])
def test_invalid_lfq_charge_is_rejected(tmp_path, charge):
    a, b = searches(tmp_path)
    p = b / "lfq.tsv"
    p.write_text(p.read_text().replace("\t2\t", "\t" + str(charge) + "\t"))
    with pytest.raises(ValueError, match="LFQ charge"):
        compare(a, b, tmp_path / "out")


def test_nonpositive_lfq_is_excluded_with_raw_value_ledger(tmp_path):
    a, b = searches(tmp_path)
    p = a / "lfq.tsv"
    p.write_text(p.read_text().replace("\t100\n", "\t-0.7079\n"))
    result = compare(a, b, tmp_path / "out")
    assert result["lfq_features"]["metaG"] == 0
    assert result["lfq_excluded_counts"]["metaG"] == {"nonpositive": 1}
    rows = _tsv_rows(tmp_path / "out/lfq_excluded.tsv")
    assert rows[0]["raw_intensity"] == "-0.7079" and rows[0]["reason"] == "nonpositive"
