"""Annotation gaps, sequence identity and model pooling must be explicit."""

import json

import numpy as np
import pytest

from fasta_lake.annotation import (
    annotate_fasta,
    prepare_annotation_targets,
    read_eggnog_annotations,
)
from fasta_lake.embeddings import embed_unknown, residue_embedding_sum, sequence_windows
from fasta_lake.multi_omics.exact import read_fasta, sha256


def test_exact_aliases_deduplicate_but_il_proteins_remain_distinct(tmp_path):
    fasta = tmp_path / "selected.faa"
    fasta.write_text(">a full header\nPEPTIDEK\n>b alias\nPEPTIDEK\n>c\nPEPTLDEK\n")
    sequences, original = prepare_annotation_targets(fasta, tmp_path / "annotation")
    assert len(sequences) == 2 and original == sha256(fasta)
    assert len(read_fasta(tmp_path / "annotation/queries.fasta")) == 2
    assert len((tmp_path / "annotation/aliases.tsv").read_text().splitlines()) == 4


@pytest.mark.parametrize("text", [">rev_x\nPEPTIDEK\n", ">a\nPEPTIDEK\n>a\nPEPTLDEK\n", ""])
def test_invalid_targets_fail_before_output(tmp_path, text):
    fasta = tmp_path / "input.faa"
    fasta.write_text(text)
    with pytest.raises(ValueError):
        prepare_annotation_targets(fasta, tmp_path / "output")
    assert not (tmp_path / "output").exists()


def test_record_limit_stops_before_whole_catalogue_is_loaded(tmp_path):
    fasta = tmp_path / "input.faa"
    # The malformed later record must never be consumed after reaching the cap.
    fasta.write_text(">a\nPEPTIDEK\n>b\nPEPTLDEK\n>\nINVALID\n")
    with pytest.raises(ValueError, match="record limit"):
        prepare_annotation_targets(fasta, tmp_path / "output", max_proteins=1)
    assert not (tmp_path / "output").exists()


def test_no_hit_and_hit_without_requested_terms_are_kept(tmp_path):
    table = tmp_path / "annotations"
    table.write_text(
        "## tool version\n#query\tseed_ortholog\tKEGG_ko\tEC\n"
        "q1\tseed1\tko:K00001\t-\nq2\tseed2\t-\t-\n"
    )
    rows, unknown = read_eggnog_annotations(table, {"q1", "q2", "q3"})
    assert unknown == ["q2", "q3"]
    assert [r["annotation_status"] for r in rows] == [
        "assigned",
        "hit_without_requested_terms",
        "no_hit",
    ]
    assert rows[1]["seed_ortholog"] == "seed2"


@pytest.mark.parametrize(
    "text",
    [
        "#query\tEC\nq1\t1.1.1.1\n",
        "#query\tEC\tKEGG_ko\nforeign\t-\t-\n",
        "#query\tEC\tKEGG_ko\nq1\t-\t-\nq1\t-\t-\n",
        "#query\tEC\tKEGG_ko\nq1\t-\n",
        "q1\t-\t-\n",
    ],
)
def test_malformed_or_foreign_annotations_fail(tmp_path, text):
    table = tmp_path / "annotations"
    table.write_text(text)
    with pytest.raises(ValueError):
        read_eggnog_annotations(table, {"q1"})


def test_windows_cover_every_residue_and_pool_by_residue_count():
    assert list(sequence_windows("ABCDEFGHI", 4)) == [(0, "ABCD"), (4, "EFGH"), (8, "I")]
    # Distinct special-token values expose accidental BOS/EOS averaging.
    first = np.array([[100, 200], [1, 2], [3, 4], [100, 200]], dtype=float)
    last = np.array([[100, 200], [8, 9], [100, 200]], dtype=float)
    total = residue_embedding_sum(first, 2) + residue_embedding_sum(last, 1)
    np.testing.assert_array_equal(total / 3, [4, 5])


@pytest.mark.parametrize(
    "values,length", [(np.ones((2, 3)), 2), (np.ones(4), 2), (np.array([[0], [np.nan], [0]]), 1)]
)
def test_invalid_model_output_is_rejected(values, length):
    with pytest.raises(ValueError):
        residue_embedding_sum(values, length)


def test_empty_unknown_roster_needs_no_model_and_checks_source(tmp_path):
    source = tmp_path / "annotation"
    source.mkdir()
    fasta = source / "unknown.fasta"
    fasta.write_text("")
    (source / "COMPLETE.json").write_text(
        json.dumps(
            {
                "status": "PASS",
                "unknown_sequences": 0,
                "unknown_fields": ["KEGG_ko", "EC"],
                "files": {"unknown.fasta": sha256(fasta)},
            }
        )
    )
    result = embed_unknown(source, tmp_path / "vectors")
    assert result["embedded_sequences"] == 0 and result["backend_ran"] is False
    fasta.write_text(">changed\nPEPTIDEK\n")
    with pytest.raises(ValueError, match="changed"):
        embed_unknown(source, tmp_path / "bad_vectors")
    assert not (tmp_path / "bad_vectors").exists()


def test_failed_eggnog_keeps_logs_and_never_claims_completion(tmp_path, monkeypatch):
    data = tmp_path / "database"
    data.mkdir()
    for name in (
        "eggnog.db",
        "eggnog.taxa.db",
        "eggnog.taxa.db.traverse.pkl",
        "eggnog_proteins.dmnd",
    ):
        (data / name).write_text("test database")
    tool = tmp_path / "emapper.py"
    tool.write_text("# fake external tool used only for failure-path regression")
    fasta = tmp_path / "selected.faa"
    fasta.write_text(">a\nPEPTIDEK\n")
    monkeypatch.setattr(
        "fasta_lake.annotation.subprocess.check_output",
        lambda *a, **k: "emapper-2.1.12 / Installed eggNOG DB version: 5.0.2",
    )

    def fail(argv, **kwargs):
        assert "--dbmem" not in argv and "--override" not in argv
        assert argv[argv.index("--tax_scope") + 1] == "auto"
        kwargs["stderr"].write("controlled failure\n")
        return {"exit_code": 2}

    monkeypatch.setattr("fasta_lake.annotation.run_recorded", fail)
    out = tmp_path / "annotation"
    with pytest.raises(RuntimeError, match="eggNOG failed"):
        annotate_fasta(fasta, out, emapper=tool, data_dir=data)
    assert (out / "eggnog.stderr").read_text() == "controlled failure\n"
    assert not (out / "COMPLETE.json").exists()


def test_dangling_output_symlink_is_preserved(tmp_path):
    destination = tmp_path / "output"
    destination.symlink_to(tmp_path / "missing")
    with pytest.raises(FileExistsError):
        annotate_fasta("irrelevant", destination, emapper="unused", data_dir="unused")
    assert destination.is_symlink()


def test_study_annotation_uses_only_verified_representatives(tmp_path, monkeypatch):
    from fasta_lake.annotation import annotate_study

    groups = tmp_path / "groups"
    groups.mkdir()
    (groups / "representatives.fasta").write_text(">rep\nPEPTIDEK\n")
    (groups / "study_groups.tsv").write_text("study_group\tselection_rank\nrep\t0\n")
    (groups / "summary.json").write_text(
        json.dumps(
            {
                "annotation_targets": {
                    "fasta": "representatives.fasta",
                    "sha256": sha256(groups / "representatives.fasta"),
                }
            }
        )
    )
    tool = tmp_path / "emapper.py"
    tool.write_text("# locally installed mapper")
    (tmp_path / "data").mkdir()
    settings = tmp_path / "settings.json"
    settings.write_text(json.dumps({"emapper": "emapper.py", "data_dir": "data"}))
    called = []

    def annotate(fasta, out, **kwargs):
        from pathlib import Path

        called.append((read_fasta(fasta), kwargs))
        Path(out).mkdir()
        return {"status": "PASS", "unique_sequences": 1, "unknown_sequences": 0}

    monkeypatch.setattr("fasta_lake.annotation.annotate_fasta", annotate)
    result = annotate_study(groups, settings, tmp_path / "out")
    assert set(called[0][0]) == {"rep"}
    assert called[0][1]["emapper"] == tool
    assert result["study"]["groups_sha256"] == sha256(groups / "study_groups.tsv")
    assert (tmp_path / "out/STUDY_COMPLETE.json").is_file()
    (groups / "representatives.fasta").write_text(">rep\nPEPTLDEK\n")
    with pytest.raises(ValueError, match="missing or changed"):
        annotate_study(groups, settings, tmp_path / "changed")
    assert len(called) == 1


@pytest.mark.parametrize(
    "extra",
    [
        {"threads": 0},
        {"unknown_fields": ["invented"]},
        {"generate_embeddings": "yes"},
        {"unknown_option": True},
        {"model": "esmc_unknown"},
    ],
)
def test_annotation_settings_reject_bad_options_before_work(tmp_path, extra):
    from fasta_lake.annotation import read_annotation_settings

    (tmp_path / "emapper.py").write_text("# tool")
    (tmp_path / "data").mkdir()
    config = tmp_path / "settings.json"
    config.write_text(json.dumps({"emapper": "emapper.py", "data_dir": "data", **extra}))
    with pytest.raises(ValueError):
        read_annotation_settings(config)
