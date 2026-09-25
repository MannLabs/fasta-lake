"""Regression checks for preselection controls and conservative error accounting."""

import json

import pytest

from tools.run_bag_entrapment import (
    augment,
    membership,
    read_psms,
    summarize,
    validate_absence,
)


def test_absence_must_be_confirmed_for_the_acquisition():
    evidence = dict(
        source="fixture only",
        reviewed_by="test",
        species=["fixture species"],
        contaminant_policy="include known contaminants",
        samples=["s1"],
        absence_confirmed=True,
    )
    validate_absence(evidence, "s1")
    with pytest.raises(ValueError, match="No confirmed absence"):
        validate_absence(evidence, "s2")
    evidence["absence_confirmed"] = False
    with pytest.raises(ValueError, match="No confirmed absence"):
        validate_absence(evidence, "s1")


def test_controls_enter_the_reservoir_with_original_identifiers(tmp_path):
    target, control, out = (tmp_path / name for name in ("t.fa", "e.fa", "all.fa"))
    target.write_text(">t\nPEPTIDERK\n")
    control.write_text(">e\nPETPIDERK\n")
    controls, count = augment(target, control, out)
    assert count == 2 and controls == {"e": "PETPIDERK"}
    assert out.read_text() == target.read_text() + control.read_text()
    with pytest.raises(FileExistsError):
        augment(target, control, out)
    control.write_text(">t\nPETPIDERK\n")
    with pytest.raises(ValueError, match="collision"):
        augment(target, control, tmp_path / "collision.fa")


def test_full_sources_and_contaminants_exclude_shared_peptides(tmp_path):
    target, control, possible = (tmp_path / name for name in ("t.fa", "e.fa", "p.fa"))
    target.write_text(">not_selected\nPEPTIDERK\n")
    control.write_text(">e\nPEPTLDERKACDEFGHKMNPQRSTK\n")
    possible.write_text(">contaminant\nACDEFGHK\n")
    classes = membership({"PEPTLDERK", "ACDEFGHK", "MNPQRSTK"}, target, control, possible)
    assert classes == {"PEPTLDERK": "shared", "ACDEFGHK": "shared", "MNPQRSTK": "entrapment_only"}
    with pytest.raises(ValueError, match="missing from input"):
        membership({"WWWWWWWK"}, target, control, possible)
    assert membership(set(), target, control, possible) == {}


def test_threshold_levels_gains_losses_and_empty_denominators():
    classes = {"A": "target_compatible", "B": "entrapment_only", "C": "shared"}
    exact = [("A", "file", "1", 0.009, 0.02)]
    relaxed = [("B", "file", "1", 0.009, 0.02), ("C", "file", "2", 0.02, 0.009)]
    rows = summarize(relaxed, classes, exact)

    def row(level, population, q=0.01):
        return next(
            r for r in rows if (r["level"], r["population"], r["q"]) == (level, population, q)
        )

    assert row("canonical_peptide", "gained_vs_exact")["entrapment_only"] == 1
    assert row("canonical_peptide", "lost_vs_exact")["target_compatible"] == 1
    assert row("PSM", "all")["shared"] == 1
    assert row("PSM", "all")["entrapment_only"] == 0
    assert row("PSM", "all", 0.001)["observed_control_fraction"] is None
    assert "NaN" not in json.dumps(rows, allow_nan=False)


def test_psm_reader_excludes_decoys_and_rejects_duplicate_winners(tmp_path):
    path = tmp_path / "sage.tsv"
    header = "label\trank\tpeptide\tfilename\tscannr\tpeptide_q\tspectrum_q\n"
    target = "1\t1\tPEPTIDERK\tf\t1\t0.001\t0.002\n"
    path.write_text(header + target + "-1\t1\tPETPIDERK\tf\t2\t0.001\t0.002\n")
    assert len(read_psms(path)) == 1
    path.write_text(header + target + target)
    with pytest.raises(ValueError, match="Multiple rank-one"):
        read_psms(path)


def test_collection_does_not_promote_execution_to_fdr_validation(tmp_path):
    from tools.run_bags import collect

    collect({"tasks": [], "output": str(tmp_path / "outputs")})
    result = json.loads((tmp_path / "COLLECTION.json").read_text())
    assert result["workflow_fdr"] == "UNVERIFIED"
    assert result["same_fdr_sensitivity_claim_supported"] is False


def test_native_preselection_control_reaches_search_and_diagnostic(tmp_path, monkeypatch):
    """Real native matching/selection; simulated search rows test accounting only."""
    import csv
    import os
    from pathlib import Path

    from fasta_lake.bags import ARMS, metadata
    from tools import run_bag_entrapment as runner

    source = os.environ.get("FASTALAKE_ALPHANOVO_SOURCE")
    engine = os.environ.get("FASTALAKE_BAG_EXTRACT")
    if not source or not engine:
        pytest.skip("Pinned AlphaNovo and streaming engine required")
    target = tmp_path / "target.fa"
    control = tmp_path / "control.fa"
    possible = tmp_path / "possible.fa"
    # Only the control is compatible with the bag; no forced post-selection addition.
    target.write_text(">target\nWWWWWWWWK\n")
    control.write_text(">control\nPETPIDERK\n")
    possible.write_text(">known\nWWWWWWWWK\n")
    predictions = tmp_path / "predictions.csv"
    fields = [
        "peptide_prediction_detokenized_unmodified",
        "peptide_prediction_detokenized",
        "score",
        "beam_rank",
        "spec_idx",
        "sequence_with_bags",
    ]
    with predictions.open("w") as handle:
        writer = csv.writer(handle)
        writer.writerow(fields)
        writer.writerow(["PEPTIDERK", "PEPTIDERK", 0.99, 0, "1", "PE[PT]IDERK"])
    config = tmp_path / "config.json"
    config.write_text('{"database": {"generate_decoys": true}}')
    absence = tmp_path / "absence.json"
    absence.write_text(
        json.dumps(
            dict(
                source="SYNTHETIC SOFTWARE FIXTURE ONLY",
                reviewed_by="test",
                species=["invented fixture"],
                contaminant_policy="fixture",
                samples=["fixture"],
                absence_confirmed=True,
            )
        )
    )
    root = Path(runner.__file__).resolve().parents[1]
    manifest = dict(
        tasks=[
            dict(
                cohort="fixture",
                sample="fixture",
                candidates=metadata(target),
                predictions=metadata(predictions),
                search_config=metadata(config),
                mzml=[],
                expected_reservoir_records=1,
            )
        ],
        code=[
            metadata(root / name)
            for name in ("fasta_lake/bags.py", "tools/run_bags.py", "tools/run_bag_entrapment.py")
        ],
        bag_engine=metadata(engine),
        sage=metadata(engine),
        alphanovo_source=source,
    )
    controls = dict(
        entrapment_fasta=metadata(control),
        possible_present_fasta=metadata(possible),
        absence_evidence=metadata(absence),
    )

    def search_fixture(local, index):
        dest = runner.case(local, index)
        assert "control" in (dest / "bags.fasta").read_text()
        assert (dest / "exact.fasta").read_text() == ""
        for arm in ARMS:
            if (dest / (arm + ".fasta")).stat().st_size:
                out = dest / "searches" / arm
                out.mkdir(parents=True)
                (out / "results.sage.tsv").write_text(
                    "label\trank\tpeptide\tfilename\tscannr\tpeptide_q\tspectrum_q\n"
                    "1\t1\tPETPIDERK\tfixture\t1\t0.001\t0.001\n"
                )

    monkeypatch.setattr(runner, "search", search_fixture)
    runner.run(manifest, 0, controls, tmp_path / "output")
    report = json.loads((tmp_path / "output" / "ENTRAPMENT_DIAGNOSTIC.json").read_text())
    assert report["stage_counts"]["bags"]["retrieved_entrapment"] == 1
    assert report["stage_counts"]["bags"]["selected_entrapment"] == 1
    accepted = next(
        r
        for r in report["arms"]["bags"]
        if r["q"] == 0.01 and r["level"] == "canonical_peptide" and r["population"] == "all"
    )
    assert accepted["entrapment_only"] == 1
    assert report["workflow_fdr"] == "UNVERIFIED"
