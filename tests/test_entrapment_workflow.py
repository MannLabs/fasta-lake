"""Tests for the user-facing entrapment workflow (generate / inject / fdp).

These cover the correctness requirements that came out of real bugs:

- protein counts come from header lines only (a ``grep -cv`` idiom inflated them ~8x)
- an undefined FDP is ``None``, never ``0.0``
- decoys are read from the engine's ``label`` column, not a guessed header prefix
- generation is deterministic: same seed, byte-identical output
- results from different injection stages or classes are not pooled
- conserved (genuinely shared) peptides are separated from class-specific ones
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from click.testing import CliRunner

from fasta_lake.cli import cli
from fasta_lake.config import EntrapmentConfig, FastaLakeConfig
from fasta_lake.entrapment_classes import ENTRAPMENT_CLASSES, class_names, get_class
from fasta_lake.entrapment_stages import (
    CONSTRUCTION_STAGES,
    ENTRAPMENT_STAGES,
    VALID_ENTRAPMENT_STAGES,
    get_stage,
    normalize_stage,
    stages_comparable,
)
from fasta_lake.entrapment_workflow import (
    ENTRAP_PREFIX,
    STAGE_APPENDED,
    STAGE_PRESELECT,
    PrefixMismatchError,
    bootstrap_median_ci,
    classify_psms,
    compute_fdp,
    conserved_peptide_scan,
    count_fasta_headers,
    detect_prefixes,
    expected_results_subdir,
    fdp_estimates,
    find_sample_results,
    generate_entrapment_fasta,
    inject_entrapment,
    load_manifests,
    partition_counts,
    prepare_from_config,
    shuffle_sequence,
    strip_modifications,
    summarize,
    verify_prefix_matches,
    write_fdp_reports,
)

# ── fixtures ─────────────────────────────────────────────────────────────────

# Deliberately wrapped over multiple sequence lines: this is what broke the
# `grep -cv '^>ENTRAP_'` counting idiom.
WRAPPED_FASTA = """\
>prot1 first protein
MKVLATGSED
KRAAWQFTGE
DPLMNK
>prot2 second protein
MSTNPKVYFD
IAVDGEPLGR
VSFELFADKV
PK
>prot3 third protein
MEEKLKAQAR
EAGVSFDELM
K
"""


@pytest.fixture
def lake(tmp_path: Path) -> Path:
    p = tmp_path / "lake.fasta"
    p.write_text(WRAPPED_FASTA)
    return p


def _sage_tsv(path: Path, rows: list[dict]) -> Path:
    cols = ["peptide", "proteins", "peptide_q", "peptide_len", "label"]
    lines = ["\t".join(cols)]
    for r in rows:
        lines.append("\t".join(str(r.get(c, "")) for c in cols))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")
    return path


# ── header-line-only counting (the ~8x inflation bug) ─────────────────────────


class TestHeaderCounting:
    def test_counts_headers_not_lines(self, lake: Path):
        """The file has 3 proteins across 11 sequence lines."""
        assert count_fasta_headers(lake) == 3
        assert lake.read_text().count("\n") > 3  # would-be inflated count

    def test_partition_excludes_sequence_lines(self, tmp_path: Path):
        p = tmp_path / "mixed.fasta"
        p.write_text(WRAPPED_FASTA + ">ENTRAP_SHUF_x\nMKVL\nATGS\n")
        n_target, n_entrap = partition_counts(p, ENTRAP_PREFIX)
        assert (n_target, n_entrap) == (3, 1)

    def test_target_count_is_not_line_count(self, tmp_path: Path):
        """Reproduces the shape of the real bug: n_target must be 3, not 14."""
        p = tmp_path / "mixed.fasta"
        p.write_text(WRAPPED_FASTA + ">ENTRAP_SHUF_x\nMKVL\nATGS\n")
        n_target, _ = partition_counts(p, ENTRAP_PREFIX)
        buggy = sum(1 for line in p.read_text().splitlines() if not line.startswith(">ENTRAP_"))
        assert n_target == 3
        assert buggy > n_target


# ── determinism ──────────────────────────────────────────────────────────────


class TestDeterminism:
    def test_same_seed_byte_identical(self, lake: Path, tmp_path: Path):
        a = tmp_path / "a.fasta"
        b = tmp_path / "b.fasta"
        generate_entrapment_fasta(a, entrap_class="shuffled", source_fasta=lake, seed=42)
        generate_entrapment_fasta(b, entrap_class="shuffled", source_fasta=lake, seed=42)
        assert a.read_bytes() == b.read_bytes()

    def test_different_seed_differs(self, lake: Path, tmp_path: Path):
        a = tmp_path / "a.fasta"
        b = tmp_path / "b.fasta"
        generate_entrapment_fasta(a, entrap_class="shuffled", source_fasta=lake, seed=42)
        generate_entrapment_fasta(b, entrap_class="shuffled", source_fasta=lake, seed=7)
        assert a.read_bytes() != b.read_bytes()

    def test_seed_zero_rejected(self, lake: Path, tmp_path: Path):
        with pytest.raises(ValueError, match="seed 0 is rejected"):
            generate_entrapment_fasta(
                tmp_path / "x.fasta",
                entrap_class="shuffled",
                source_fasta=lake,
                seed=0,
            )

    def test_subsample_deterministic(self, lake: Path, tmp_path: Path):
        a = tmp_path / "a.fasta"
        b = tmp_path / "b.fasta"
        for out in (a, b):
            generate_entrapment_fasta(
                out,
                entrap_class="shuffled",
                source_fasta=lake,
                seed=42,
                n_proteins=2,
            )
        assert a.read_bytes() == b.read_bytes()
        assert count_fasta_headers(a) == 2

    def test_bootstrap_deterministic(self):
        vals = [0.01 * i for i in range(20)]
        assert bootstrap_median_ci(vals, seed=42) == bootstrap_median_ci(vals, seed=42)


# ── generation ───────────────────────────────────────────────────────────────


class TestGenerate:
    def test_tags_with_class_prefix(self, lake: Path, tmp_path: Path):
        out = tmp_path / "e.fasta"
        stats = generate_entrapment_fasta(
            out,
            entrap_class="shuffled",
            source_fasta=lake,
            seed=42,
        )
        assert stats.prefix == "ENTRAP_SHUF_"
        assert stats.n_written == 3
        headers = [ln for ln in out.read_text().splitlines() if ln.startswith(">")]
        assert all(h.startswith(">ENTRAP_SHUF_") for h in headers)

    def test_shuffle_preserves_composition(self, lake: Path, tmp_path: Path):
        out = tmp_path / "e.fasta"
        generate_entrapment_fasta(out, entrap_class="shuffled", source_fasta=lake, seed=42)
        src = "MKVLATGSEDKRAAWQFTGEDPLMNK"
        got = out.read_text().splitlines()[1]
        assert sorted(got) == sorted(src)

    def test_shuffle_preserves_cleavage_sites(self, lake: Path, tmp_path: Path):
        """K/R stay put, so tryptic peptide lengths match the source."""
        out = tmp_path / "e.fasta"
        generate_entrapment_fasta(out, entrap_class="shuffled", source_fasta=lake, seed=42)
        src = "MKVLATGSEDKRAAWQFTGEDPLMNK"
        got = out.read_text().splitlines()[1]
        src_kr = [i for i, c in enumerate(src) if c in "KR"]
        got_kr = [i for i, c in enumerate(got) if c in "KR"]
        assert src_kr == got_kr

    def test_naive_shuffle_moves_cleavage_sites(self):
        import random

        seq = "MKVLATGSEDKRAAWQFTGEDPLMNK"
        naive = shuffle_sequence(seq, random.Random(1), preserve_cleavage_sites=False)
        kept = shuffle_sequence(seq, random.Random(1), preserve_cleavage_sites=True)
        assert [i for i, c in enumerate(kept) if c in "KR"] == [
            i for i, c in enumerate(seq) if c in "KR"
        ]
        assert sorted(naive) == sorted(seq)

    def test_shuffled_requires_source(self, tmp_path: Path):
        with pytest.raises(ValueError, match="derives its null from the target"):
            generate_entrapment_fasta(tmp_path / "e.fasta", entrap_class="shuffled")

    def test_user_class_copies_sequences_verbatim(self, lake: Path, tmp_path: Path):
        out = tmp_path / "e.fasta"
        generate_entrapment_fasta(out, entrap_class="user", source_fasta=lake, seed=42)
        assert out.read_text().splitlines()[1] == "MKVLATGSEDKRAAWQFTGEDPLMNK"

    def test_unknown_class_rejected(self, lake: Path, tmp_path: Path):
        with pytest.raises(ValueError, match="Unknown entrapment class"):
            generate_entrapment_fasta(
                tmp_path / "e.fasta",
                entrap_class="nope",
                source_fasta=lake,
            )


# ── class registry ───────────────────────────────────────────────────────────


class TestClassRegistry:
    def test_ships_required_classes(self):
        for name in ("shuffled", "thermophilic", "plant", "user"):
            assert name in ENTRAPMENT_CLASSES

    def test_shuffled_is_default_and_headline(self):
        assert class_names()[0] == "shuffled"
        assert get_class("shuffled").recommended_headline is False

    def test_every_class_declares_difficulty_and_caveat(self):
        for name, c in ENTRAPMENT_CLASSES.items():
            assert c.difficulty, name
            assert c.estimates, name
            assert len(c.caveat) > 40, name

    def test_non_self_is_not_an_fdp_and_not_poolable(self):
        c = get_class("non_self")
        assert "NOT FDP" in c.estimates or "NOT FDP" in c.estimates.upper()
        assert c.poolable is False

    def test_plant_declares_homology_pressure(self):
        assert get_class("plant").difficulty == "cross-kingdom-homologous"

    def test_prefixes_are_unique_and_namespaced(self):
        prefixes = [c.prefix for c in ENTRAPMENT_CLASSES.values()]
        assert len(prefixes) == len(set(prefixes))
        assert all(p.startswith(ENTRAP_PREFIX) for p in prefixes)


# ── injection + manifest ─────────────────────────────────────────────────────


class TestInject:
    def _entrap(self, lake: Path, tmp_path: Path) -> Path:
        out = tmp_path / "e.fasta"
        generate_entrapment_fasta(out, entrap_class="shuffled", source_fasta=lake, seed=42)
        return out

    def test_manifest_records_stage_class_seed_and_ratio(self, lake: Path, tmp_path: Path):
        e = self._entrap(lake, tmp_path)
        m = inject_entrapment(
            lake,
            e,
            tmp_path / "aug.fasta",
            stage=STAGE_PRESELECT,
            entrap_class="shuffled",
            seed=42,
        )
        assert m["stage"] == STAGE_PRESELECT
        assert m["class"] == "shuffled"
        assert m["seed"] == 42
        assert m["n_target"] == 3
        assert m["n_entrapment"] == 3
        assert m["ratio_entrapment_to_target"] == 1.0
        assert "construction" in m["measures"]
        assert Path(m["manifest_path"]).exists()

    def test_appended_stage_says_search_error(self, lake: Path, tmp_path: Path):
        e = self._entrap(lake, tmp_path)
        m = inject_entrapment(
            lake,
            e,
            tmp_path / "aug.fasta",
            stage=STAGE_APPENDED,
            entrap_class="shuffled",
            seed=42,
        )
        assert "NOT construction error" in m["measures"]

    def test_output_contains_both_partitions(self, lake: Path, tmp_path: Path):
        e = self._entrap(lake, tmp_path)
        out = tmp_path / "aug.fasta"
        inject_entrapment(lake, e, out, entrap_class="shuffled")
        assert partition_counts(out, ENTRAP_PREFIX) == (3, 3)

    def test_double_injection_refused(self, lake: Path, tmp_path: Path):
        e = self._entrap(lake, tmp_path)
        out = tmp_path / "aug.fasta"
        inject_entrapment(lake, e, out, entrap_class="shuffled")
        with pytest.raises(ValueError, match="already contains"):
            inject_entrapment(out, e, tmp_path / "aug2.fasta", entrap_class="shuffled")

    def test_untagged_entrapment_refused(self, lake: Path, tmp_path: Path):
        with pytest.raises(ValueError, match="without the"):
            inject_entrapment(lake, lake, tmp_path / "aug.fasta", entrap_class="shuffled")

    def test_bad_stage_refused(self, lake: Path, tmp_path: Path):
        e = self._entrap(lake, tmp_path)
        with pytest.raises(ValueError, match="Unknown entrapment stage"):
            inject_entrapment(lake, e, tmp_path / "aug.fasta", stage="whenever")


class TestRefusesToPool:
    def _manifest(self, path: Path, stage: str, cls: str) -> Path:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps({"stage": stage, "class": cls}))
        return path

    def test_mixed_stages_refused(self, tmp_path: Path):
        a = self._manifest(tmp_path / "a" / "m.json", STAGE_PRESELECT, "shuffled")
        b = self._manifest(tmp_path / "b" / "m.json", STAGE_APPENDED, "shuffled")
        with pytest.raises(ValueError, match="different injection"):
            load_manifests([a, b])

    def test_mixed_classes_refused(self, tmp_path: Path):
        a = self._manifest(tmp_path / "a" / "m.json", STAGE_PRESELECT, "shuffled")
        b = self._manifest(tmp_path / "b" / "m.json", STAGE_PRESELECT, "plant")
        with pytest.raises(ValueError, match="different classes"):
            load_manifests([a, b])

    def test_allow_mixed_overrides(self, tmp_path: Path):
        a = self._manifest(tmp_path / "a" / "m.json", STAGE_PRESELECT, "shuffled")
        b = self._manifest(tmp_path / "b" / "m.json", STAGE_APPENDED, "plant")
        resolved, manifests = load_manifests([a, b], allow_mixed=True)
        assert resolved["stage"] is None
        assert len(manifests) == 2

    def test_consistent_manifests_resolve(self, tmp_path: Path):
        a = self._manifest(tmp_path / "a" / "m.json", STAGE_PRESELECT, "shuffled")
        b = self._manifest(tmp_path / "b" / "m.json", STAGE_PRESELECT, "shuffled")
        resolved, _ = load_manifests([a, b])
        assert resolved == {"stage": STAGE_PRESELECT, "class": "shuffled"}


# ── undefined FDP is None, never 0.0 ─────────────────────────────────────────


class TestUndefinedFdp:
    def test_empty_entrapment_partition_is_none_not_zero(self):
        r = fdp_estimates(entrap=0, target=100, db_target=1000, db_entrap=0)
        assert r["fdp_elias_gygi"] is None
        assert r["fdp_combined_noble"] is None
        assert r["entrapment_in_db"] is False

    def test_unknown_db_sizes_give_none(self):
        r = fdp_estimates(entrap=1, target=100, db_target=None, db_entrap=None)
        assert r["fdp_elias_gygi"] is None
        assert r["entrapment_in_db"] is None

    def test_no_discoveries_is_none_not_zero(self):
        r = fdp_estimates(entrap=0, target=0, db_target=100, db_entrap=100)
        assert r["fdp_uncorrected"] is None
        assert r["evaluable"] is False

    def test_a_real_zero_is_still_zero(self):
        """Zero entrapment hits WITH a populated partition is a genuine 0.0."""
        r = fdp_estimates(entrap=0, target=100, db_target=100, db_entrap=100)
        assert r["fdp_uncorrected"] == 0.0
        assert r["evaluable"] is True

    def test_estimator_definitions(self):
        """Each estimator computes the equation it is named for.

        RENAMED AND CORRECTED 2026-08-02. This test previously asserted
        ``fdp_lower_bound == 1/100``, i.e. E/T. That is Wen equation 3, the
        "sample method", which the paper states "cannot be used to provide
        empirical evidence that a tool controls the FDR nor that the tool fails" —
        it is not a bound at all. The test was encoding the mislabelling.

        It also asserted the estimators were mutually distinct. They are not:
        ``fdp_uncorrected`` and ``fdp_lower_bound`` were always the SAME quantity
        under two names, E/(T+E). Asserting distinctness enforced the duplication.
        """
        r = fdp_estimates(entrap=1, target=100, db_target=1000, db_entrap=100)
        # eq 2, the lower bound on the FDP itself
        assert r["fdp_lower_bound_eq2"] == pytest.approx(1 / 101)
        assert r["fdp_uncorrected"] == pytest.approx(1 / 101)
        assert r["fdp_lower_bound"] == pytest.approx(1 / 101)
        # eq 3, the sample method — NOT a bound, kept only so it is nameable
        assert r["fdp_sample_method_eq3"] == pytest.approx(1 / 100)
        # Elias & Gygi concatenated, 2E/(T+E); equals eq 1 at r = 1
        assert r["fdp_elias_gygi"] == pytest.approx(2 / 101)
        assert r["fdp_elias_gygi_concatenated"] == pytest.approx(2 / 101)
        assert r["fdp_wen_equal_odds"] == pytest.approx(2 / 101)
        # eq 1, the combined upper bound on the FDR, at r = db_entrap/db_target
        assert r["fdp_combined_upper_eq1"] == pytest.approx((1 + 1 / 0.1) / 101)

    def test_lower_bound_is_not_the_sample_method(self):
        """The two must never collapse to the same number again.

        E/(T+E) and E/T differ, and the whole defect was reporting the second
        under the first's name.
        """
        r = fdp_estimates(entrap=5, target=100, db_target=1000, db_entrap=100)
        assert r["fdp_lower_bound_eq2"] == pytest.approx(5 / 105)
        assert r["fdp_sample_method_eq3"] == pytest.approx(5 / 100)
        assert r["fdp_lower_bound_eq2"] != r["fdp_sample_method_eq3"]

    def test_shared_sensitivity_bound(self):
        r = fdp_estimates(entrap=1, target=100, shared=9)
        assert r["fdp_shared_sensitivity"] == pytest.approx(10 / 110)
        assert r["fdp_uncorrected"] == pytest.approx(1 / 101)


class TestSummarize:
    def test_reports_n_median_iqr(self):
        s = summarize([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7])
        assert s["n"] == 7
        assert s["median"] == pytest.approx(0.4)
        assert s["iqr"] == pytest.approx(0.3)

    def test_ci_withheld_for_tiny_n(self):
        s = summarize([0.1, 0.2, 0.3])
        assert s["n"] == 3
        assert s["ci_lo"] is None

    def test_ci_present_for_adequate_n(self):
        s = summarize([0.01 * i for i in range(30)])
        assert s["ci_lo"] is not None and s["ci_hi"] is not None
        assert s["ci_lo"] <= s["median"] <= s["ci_hi"]

    def test_empty_is_none_not_zero(self):
        s = summarize([])
        assert s["n"] == 0 and s["median"] is None


# ── PSM classification / decoy handling ──────────────────────────────────────


class TestClassifyPsms:
    def test_uses_label_column_not_prefix(self, tmp_path: Path):
        """A decoy is a decoy because label says so, whatever the prefix looks like."""
        tsv = _sage_tsv(
            tmp_path / "results.sage.tsv",
            [
                {
                    "peptide": "PEPTIDEKAA",
                    "proteins": "prot1",
                    "peptide_q": 0.001,
                    "peptide_len": 10,
                    "label": 1,
                },
                {
                    "peptide": "DECOYPEPKA",
                    "proteins": "prot9",
                    "peptide_q": 0.001,
                    "peptide_len": 10,
                    "label": -1,
                },
                {
                    "peptide": "ENTRAPPEPK",
                    "proteins": "ENTRAP_SHUF_prot1",
                    "peptide_q": 0.001,
                    "peptide_len": 10,
                    "label": 1,
                },
            ],
        )
        c = classify_psms(tsv, prefix="ENTRAP_SHUF_")
        assert (c.target, c.entrap, c.decoy) == (1, 1, 1)
        assert c.decoy_source == "label column"

    def test_decoys_not_folded_into_targets(self, tmp_path: Path):
        tsv = _sage_tsv(
            tmp_path / "results.sage.tsv",
            [
                {
                    "peptide": f"PEPTIDEK{i:02d}",
                    "proteins": "rev_prot1",
                    "peptide_q": 0.001,
                    "peptide_len": 10,
                    "label": -1,
                }
                for i in range(5)
            ],
        )
        c = classify_psms(tsv, prefix="ENTRAP_SHUF_")
        assert c.target == 0 and c.decoy == 5

    def test_prefix_fallback_only_without_label(self, tmp_path: Path):
        tsv = tmp_path / "results.sage.tsv"
        tsv.write_text(
            "peptide\tproteins\tpeptide_q\tpeptide_len\n"
            "PEPTIDEKAA\trev_prot1\t0.001\t10\n"
            "PEPTIDEKBB\tprot1\t0.001\t10\n"
        )
        c = classify_psms(tsv, prefix="ENTRAP_SHUF_", decoy_prefix="rev_")
        assert (c.target, c.decoy) == (1, 1)
        assert "rev_" in c.decoy_source

    def test_shared_peptides_counted_separately(self, tmp_path: Path):
        tsv = _sage_tsv(
            tmp_path / "results.sage.tsv",
            [
                {
                    "peptide": "CONSERVEDK",
                    "proteins": "prot1;ENTRAP_SHUF_prot2",
                    "peptide_q": 0.001,
                    "peptide_len": 10,
                    "label": 1,
                },
            ],
        )
        c = classify_psms(tsv, prefix="ENTRAP_SHUF_")
        assert (c.target, c.entrap, c.shared) == (0, 0, 1)

    def test_filters_applied(self, tmp_path: Path):
        tsv = _sage_tsv(
            tmp_path / "results.sage.tsv",
            [
                {
                    "peptide": "SHORTK",
                    "proteins": "prot1",
                    "peptide_q": 0.001,
                    "peptide_len": 6,
                    "label": 1,
                },
                {
                    "peptide": "HIGHQPEPTK",
                    "proteins": "prot1",
                    "peptide_q": 0.5,
                    "peptide_len": 10,
                    "label": 1,
                },
                {
                    "peptide": "GOODPEPTIK",
                    "proteins": "prot1",
                    "peptide_q": 0.001,
                    "peptide_len": 10,
                    "label": 1,
                },
            ],
        )
        c = classify_psms(tsv, prefix="ENTRAP_SHUF_", q_threshold=0.01, min_length=9)
        assert c.target == 1

    def test_missing_q_column_is_loud(self, tmp_path: Path):
        tsv = tmp_path / "results.sage.tsv"
        tsv.write_text("peptide\tproteins\tpeptide_len\nPEPK\tprot1\t10\n")
        with pytest.raises(ValueError, match="has no 'peptide_q' column"):
            classify_psms(tsv)

    def test_strip_modifications(self):
        assert strip_modifications("M[+15.9949]PEPTIDEK") == "MPEPTIDEK"
        assert strip_modifications("PEPC[+57.0215]K") == "PEPCK"


# ── conserved-peptide analysis ───────────────────────────────────────────────


class TestConservedPeptides:
    def test_finds_peptide_present_in_target_space(self, lake: Path):
        scan = conserved_peptide_scan({"ATGSEDKR"}, lake)
        assert scan.n_conserved == 1
        assert "prot1" in scan.conserved["ATGSEDKR"]

    def test_class_specific_peptide_not_flagged(self, lake: Path):
        scan = conserved_peptide_scan({"WWWWWWWW"}, lake)
        assert scan.n_conserved == 0

    def test_il_normalisation(self, lake: Path):
        """prot2 contains ...IAVDGEPLGR...; the L/I variant must still match."""
        scan = conserved_peptide_scan({"LAVDGEPLGR"}, lake, normalize_il=True)
        assert scan.n_conserved == 1
        strict = conserved_peptide_scan({"LAVDGEPLGR"}, lake, normalize_il=False)
        assert strict.n_conserved == 0

    def test_reversed_negative_control_runs(self, lake: Path):
        scan = conserved_peptide_scan({"ATGSEDKR"}, lake)
        assert scan.n_control_tested == 1
        assert scan.n_control_conserved == 0
        assert scan.control_rate == 0.0

    def test_entrapment_records_excluded_from_target_space(self, tmp_path: Path):
        """A peptide must not be called conserved by matching entrapment itself."""
        f = tmp_path / "aug.fasta"
        f.write_text(">ENTRAP_SHUF_x\nWWWWWWWWWW\n")
        scan = conserved_peptide_scan({"WWWWWWWWWW"}, f)
        assert scan.n_conserved == 0

    def test_reports_which_proteins_share(self, lake: Path):
        scan = conserved_peptide_scan({"ATGSEDKR"}, lake)
        assert scan.conserved["ATGSEDKR"] == ["prot1"]


# ── end-to-end FDP over a results directory ──────────────────────────────────


def _results_tree(
    tmp_path: Path, lake: Path, stage: str = STAGE_PRESELECT, cls: str = "shuffled"
) -> Path:
    """Build a small results dir: 3 samples, one with an entrapment hit."""
    root = tmp_path / "results"
    root.mkdir()
    (root / "entrapment_manifest.json").write_text(
        json.dumps(
            {
                "schema": "fastalake.entrapment.manifest/1",
                "stage": stage,
                "class": cls,
                "seed": 42,
                "n_target": 3,
                "n_entrapment": 3,
            }
        )
    )
    for i, n_entrap in enumerate([0, 1, 2]):
        rows = [
            {
                "peptide": f"TARGETPEP{j}",
                "proteins": "prot1",
                "peptide_q": 0.001,
                "peptide_len": 10,
                "label": 1,
            }
            for j in range(20)
        ]
        # One entrapment hit is a genuinely conserved peptide from prot1.
        for j in range(n_entrap):
            pep = "ATGSEDKRAA" if j == 0 else f"WWWWWWWWW{j}"
            rows.append(
                {
                    "peptide": pep,
                    "proteins": "ENTRAP_SHUF_prot2",
                    "peptide_q": 0.001,
                    "peptide_len": 10,
                    "label": 1,
                }
            )
        rows.append(
            {
                "peptide": "DECOYPEPTK",
                "proteins": "rev_prot1",
                "peptide_q": 0.001,
                "peptide_len": 10,
                "label": -1,
            }
        )
        _sage_tsv(root / f"sample{i}" / "results.sage.tsv", rows)
    return root


class TestFindSampleResults:
    """A sample directory with several entrapment runs must never be resolved by glob order."""

    @staticmethod
    def _tree(tmp_path: Path, subdirs: list[str]) -> Path:
        root = tmp_path / "results"
        for sub in subdirs:
            f = root / "S1" / sub / "results.sage.tsv"
            f.parent.mkdir(parents=True)
            f.write_text("peptide\tproteins\n")
        return root

    def test_direct_layout(self, tmp_path: Path):
        root = tmp_path / "results"
        (root / "S1").mkdir(parents=True)
        (root / "S1" / "results.sage.tsv").write_text("peptide\tproteins\n")
        assert find_sample_results(root) == [("S1", root / "S1" / "results.sage.tsv")]

    def test_single_nested_run_is_taken(self, tmp_path: Path):
        root = self._tree(tmp_path, ["sage"])
        assert find_sample_results(root) == [("S1", root / "S1" / "sage" / "results.sage.tsv")]

    def test_class_subdir_is_preferred_over_glob_order(self, tmp_path: Path):
        root = self._tree(tmp_path, ["entrapment", "entrapment_archaea", "entrapment_plant"])
        found = find_sample_results(root, prefer_subdir=expected_results_subdir("plant"))
        assert found == [("S1", root / "S1" / "entrapment_plant" / "results.sage.tsv")]
        found = find_sample_results(root, prefer_subdir=expected_results_subdir("shuffled"))
        assert found == [("S1", root / "S1" / "entrapment" / "results.sage.tsv")]

    def test_ambiguous_without_preference_is_refused(self, tmp_path: Path):
        root = self._tree(tmp_path, ["entrapment", "entrapment_plant"])
        with pytest.raises(ValueError, match="will not be guessed"):
            find_sample_results(root)

    def test_preference_absent_is_refused(self, tmp_path: Path):
        root = self._tree(tmp_path, ["entrapment", "entrapment_plant"])
        with pytest.raises(ValueError, match="expected subdirectory 'entrapment_archaea'"):
            find_sample_results(root, prefer_subdir=expected_results_subdir("archaea"))

    def test_results_name_with_subdir_selects_explicitly(self, tmp_path: Path):
        root = self._tree(tmp_path, ["entrapment", "entrapment_plant"])
        found = find_sample_results(root, results_name="entrapment_plant/results.sage.tsv")
        assert found == [("S1", root / "S1" / "entrapment_plant" / "results.sage.tsv")]


class TestComputeFdpSelectsRunByClass:
    @staticmethod
    def _two_run_tree(tmp_path: Path) -> Path:
        """Three samples, each with a shuffled run (0 hits) and a plant run (1 hit)."""
        root = tmp_path / "results"
        root.mkdir()
        (root / "entrapment_manifest.json").write_text(
            json.dumps(
                {
                    "schema": "fastalake.entrapment.manifest/1",
                    "stage": STAGE_PRESELECT,
                    "class": "plant",
                    "seed": 42,
                    "n_target": 3,
                    "n_entrapment": 3,
                }
            )
        )
        target = [
            {
                "peptide": f"TARGETPEP{j}",
                "proteins": "prot1",
                "peptide_q": 0.001,
                "peptide_len": 10,
                "label": 1,
            }
            for j in range(20)
        ]
        plant_hit = {
            "peptide": "WWWWWWWWWQ",
            "proteins": "ENTRAP_PLANT_x",
            "peptide_q": 0.001,
            "peptide_len": 10,
            "label": 1,
        }
        for i in range(3):
            _sage_tsv(root / f"sample{i}" / "entrapment" / "results.sage.tsv", target)
            _sage_tsv(
                root / f"sample{i}" / "entrapment_plant" / "results.sage.tsv", [*target, plant_hit]
            )
        return root

    def test_plant_class_reads_the_plant_run(self, tmp_path: Path, lake: Path):
        root = self._two_run_tree(tmp_path)
        r = compute_fdp(root, target_fasta=lake)
        assert r["class"] == "plant"
        for row in r["per_sample"]:
            assert row["results_file"].endswith("entrapment_plant/results.sage.tsv")
            assert row["n_entrapment"] == 1
        assert r["coverage"]["n_evaluable"] == 3


class TestComputeFdp:
    def test_reads_stage_and_class_from_manifest(self, tmp_path: Path, lake: Path):
        root = _results_tree(tmp_path, lake)
        r = compute_fdp(root, target_fasta=lake)
        assert r["stage"] == STAGE_PRESELECT
        assert r["class"] == "shuffled"
        assert r["measures"] == "database-construction error"

    def test_refuses_when_stage_unknown(self, tmp_path: Path, lake: Path):
        root = _results_tree(tmp_path, lake)
        (root / "entrapment_manifest.json").unlink()
        with pytest.raises(ValueError, match="injection stage is unknown"):
            compute_fdp(root, entrap_class="shuffled", target_fasta=lake)

    def test_refuses_when_class_unknown(self, tmp_path: Path, lake: Path):
        root = _results_tree(tmp_path, lake)
        (root / "entrapment_manifest.json").unlink()
        with pytest.raises(ValueError, match="entrapment class is unknown"):
            compute_fdp(root, stage=STAGE_PRESELECT, target_fasta=lake)

    def test_explicit_stage_contradicting_manifest_refused(self, tmp_path: Path, lake: Path):
        root = _results_tree(tmp_path, lake)
        with pytest.raises(ValueError, match="contradicts"):
            compute_fdp(root, stage=STAGE_APPENDED, target_fasta=lake)

    def test_coverage_reported(self, tmp_path: Path, lake: Path):
        root = _results_tree(tmp_path, lake)
        r = compute_fdp(root, target_fasta=lake)
        assert r["coverage"]["n_attempted"] == 3
        assert r["coverage"]["n_evaluable"] == 3

    def test_raw_and_conserved_corrected_both_present(self, tmp_path: Path, lake: Path):
        root = _results_tree(tmp_path, lake)
        r = compute_fdp(root, target_fasta=lake)
        block = r["summary"]["fdp_uncorrected"]
        assert block["raw"]["n"] == 3
        assert block["conserved_corrected"]["n"] == 3
        # The conserved peptide (ATGSEDKRAA occurs in prot1) is subtracted, so the
        # corrected FDP must be strictly lower than the raw one.
        assert block["conserved_corrected"]["median"] < block["raw"]["median"]

    def test_conserved_peptide_identified(self, tmp_path: Path, lake: Path):
        root = _results_tree(tmp_path, lake)
        r = compute_fdp(root, target_fasta=lake)
        cons = r["conserved_peptides"]
        assert cons["performed"] is True
        assert cons["n_conserved"] == 1
        assert cons["negative_control"]["n_found"] == 0

    def test_shuffled_overlap_is_not_automatically_false_exculpation(
        self, tmp_path: Path, lake: Path
    ):
        """Whole-protein shuffling alone does not establish an absent peptide null."""
        root = _results_tree(tmp_path, lake)
        r = compute_fdp(root, target_fasta=lake)
        assert "false_exculpation_calibration" not in r["conserved_peptides"]
        assert r["conserved_peptides"]["n_conserved"] == 1

    def test_uncorrected_status_stated_when_skipped(self, tmp_path: Path, lake: Path):
        root = _results_tree(tmp_path, lake)
        r = compute_fdp(root, conserved_check=False)
        assert "NOT PERFORMED" in r["conserved_peptides"]["status"]
        assert r["summary"]["fdp_uncorrected"]["conserved_corrected"]["n"] == 0

    def test_homologous_class_requires_target_fasta(self, tmp_path: Path, lake: Path):
        root = _results_tree(tmp_path, lake, cls="plant")
        with pytest.raises(ValueError, match="homology pressure"):
            compute_fdp(root)

    def test_lake_level_db_sizes_are_flagged_as_such(self, tmp_path: Path, lake: Path):
        """Without --db-pattern the sizes come from the lake, and must say so."""
        root = _results_tree(tmp_path, lake)
        r = compute_fdp(root, target_fasta=lake)
        source = r["per_sample"][0]["db_size_source"]
        assert "injection manifest in" in source and "not the per-sample searched" in source

    def test_per_sample_db_used_when_available(self, tmp_path: Path, lake: Path):
        root = _results_tree(tmp_path, lake)
        for sample in ("sample0", "sample1", "sample2"):
            db = root / sample / "razor.fasta"
            db.write_text(WRAPPED_FASTA + ">ENTRAP_SHUF_prot2\nMKVL\nATGS\n")
        r = compute_fdp(root, target_fasta=lake, db_pattern="razor.fasta")
        row = r["per_sample"][0]
        assert row["db_target"] == 3 and row["db_entrapment"] == 1
        assert "per-sample" in row["db_size_source"]

    def test_class_metadata_carried_into_record(self, tmp_path: Path, lake: Path):
        root = _results_tree(tmp_path, lake)
        r = compute_fdp(root, target_fasta=lake)
        meta = r["class_metadata"]
        assert meta["difficulty"] and meta["estimates"] and meta["caveat"]

    def test_peptide_level_counting(self, tmp_path: Path, lake: Path):
        root = _results_tree(tmp_path, lake)
        r = compute_fdp(root, target_fasta=lake, peptide_level=True)
        assert r["counting"] == "peptide"

    def test_writes_record_and_per_sample(self, tmp_path: Path, lake: Path):
        root = _results_tree(tmp_path, lake)
        r = compute_fdp(root, dataset="testset", target_fasta=lake)
        written = write_fdp_reports(r, tmp_path / "out")
        assert written["per_sample"].exists()
        assert written["record"].exists()
        assert written["conserved"].exists()
        rec = json.loads(written["record"].read_text())
        for key in (
            "dataset",
            "class",
            "stage",
            "coverage",
            "summary",
            "conserved_peptides",
            "filters",
            "class_metadata",
        ):
            assert key in rec
        header = written["per_sample"].read_text().splitlines()[0]
        assert "fdp_uncorrected" in header
        assert "fdp_uncorrected_conserved_corrected" in header


# ── stage vocabulary (four distinct injection points) ────────────────────────


class TestStageVocabulary:
    def test_four_stages_in_injection_order(self):
        assert VALID_ENTRAPMENT_STAGES == (
            "prelake",
            "prebuild",
            "preselect",
            "appended",
        )

    def test_prebuild_is_legal(self):
        """prebuild is the arm the figures load; it must be expressible."""
        assert normalize_stage("prebuild") == "prebuild"
        assert EntrapmentConfig(stage="prebuild").stage == "prebuild"

    def test_prebuild_and_preselect_are_different_stages(self):
        """They differ by whether entrapment had to survive clustering."""
        pre, pb = get_stage("preselect"), get_stage("prebuild")
        assert pre.name != pb.name
        assert "cluster97" in pb.faces
        assert "cluster97" not in pre.faces
        assert not stages_comparable("prebuild", "preselect")

    def test_prelake_is_the_earliest_and_faces_everything(self):
        pl = get_stage("prelake")
        assert pl.order == 0
        assert "cohort evidence filter" in pl.faces
        assert len(pl.faces) > len(get_stage("prebuild").faces)

    def test_construction_vs_search_stages(self):
        assert CONSTRUCTION_STAGES == ("prelake", "prebuild", "preselect")
        assert get_stage("appended").is_construction_test is False
        assert "NOT construction" in get_stage("appended").measures

    def test_legacy_alias_normalises(self):
        assert normalize_stage("pre-selection") == "preselect"
        assert EntrapmentConfig(stage="pre-selection").stage == "preselect"

    def test_unknown_stage_rejected(self):
        with pytest.raises(ValueError, match="Unknown entrapment stage"):
            normalize_stage("whenever")

    def test_every_stage_names_its_reference_script(self):
        for name, s in ENTRAPMENT_STAGES.items():
            assert s.reference_script.endswith(".sh"), name
            assert s.injection_point and s.faces, name

    def test_stage_metadata_recorded_in_manifest(self, lake: Path, tmp_path: Path):
        e = tmp_path / "e.fasta"
        generate_entrapment_fasta(e, entrap_class="shuffled", source_fasta=lake, seed=42)
        m = inject_entrapment(
            lake, e, tmp_path / "aug.fasta", stage="prebuild", entrap_class="shuffled", seed=42
        )
        assert m["stage"] == "prebuild"
        assert m["is_construction_test"] is True
        assert "cluster97" in m["stage_metadata"]["faces"]

    def test_appended_manifest_says_search_error(self, lake: Path, tmp_path: Path):
        e = tmp_path / "e.fasta"
        generate_entrapment_fasta(e, entrap_class="shuffled", source_fasta=lake, seed=42)
        m = inject_entrapment(
            lake, e, tmp_path / "aug.fasta", stage="appended", entrap_class="shuffled"
        )
        assert m["is_construction_test"] is False
        assert "NOT construction error" in m["measures"]

    def test_fdp_refuses_to_pool_prebuild_with_preselect(self, tmp_path: Path, lake: Path):
        a = tmp_path / "a" / "m.json"
        b = tmp_path / "b" / "m.json"
        for p, st in ((a, "prebuild"), (b, "preselect")):
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_text(json.dumps({"stage": st, "class": "shuffled"}))
        with pytest.raises(ValueError, match="different injection"):
            load_manifests([a, b])


# ── prefix robustness ────────────────────────────────────────────────────────

NONSTANDARD_FASTA = """\
>ENTRAPHUMAN_P60709 actin
MKVLATGSEDKR
>ENTRAPNN_A0A123 near neighbour
MSTNPKVYFDIA
>prot1 real target
MEEKLKAQAREA
"""


class TestPrefixRobustness:
    def test_default_pattern_matches_all_conventions(self, tmp_path: Path):
        """ENTRAPHUMAN_ and ENTRAPNN_ do NOT match 'ENTRAP_' but must be caught."""
        f = tmp_path / "mixed.fasta"
        f.write_text(NONSTANDARD_FASTA)
        assert partition_counts(f, ENTRAP_PREFIX) == (1, 2)
        # The old convention would have scored both entrapment records as target.
        assert partition_counts(f, "ENTRAP_") == (3, 0)

    def test_detect_prefixes_reports_what_is_there(self, tmp_path: Path):
        f = tmp_path / "mixed.fasta"
        f.write_text(NONSTANDARD_FASTA)
        d = detect_prefixes(f)
        assert d["n_headers"] == 3
        assert set(d["observed"]) == {"ENTRAPHUMAN_", "ENTRAPNN_"}
        assert d["unknown"] == []

    def test_detect_prefixes_flags_unknown_convention(self, tmp_path: Path):
        f = tmp_path / "odd.fasta"
        f.write_text(">ENTRAPWEIRD_x1 y\nMKVL\n>prot1 t\nMEEK\n")
        d = detect_prefixes(f)
        assert "ENTRAPWEIRD_" in d["unknown"]

    def test_hyphen_family_landmine_is_caught_and_named(self, tmp_path: Path):
        """sihumi_gt_ceiling uses ENTRAPHUMAN- (hyphen); no ^ENTRAP_ matcher catches it."""
        f = tmp_path / "hyphen.fasta"
        f.write_text(">ENTRAPHUMAN-P60709 actin\nMKVL\n>prot1 t\nMEEK\n")
        # The old ENTRAP_ convention scores it as target -> silent ~0% FDP.
        assert partition_counts(f, "ENTRAP_") == (2, 0)
        # The registry catches it, and detect_prefixes names it verbatim as known.
        assert partition_counts(f, ENTRAP_PREFIX) == (1, 1)
        d = detect_prefixes(f)
        assert "ENTRAPHUMAN-" in d["observed"]
        assert "ENTRAPHUMAN-" in d["known"]

    def test_known_prefixes_registry_covers_all_four_families(self):
        from fasta_lake.entrapment_workflow import KNOWN_ENTRAP_PREFIXES

        for fam in ("ENTRAP_SHUF_", "ENTRAPHUMAN_", "ENTRAPNN_", "ENTRAPHUMAN-"):
            assert fam in KNOWN_ENTRAP_PREFIXES

    def test_verify_raises_on_mismatch(self, tmp_path: Path):
        f = tmp_path / "mixed.fasta"
        f.write_text(NONSTANDARD_FASTA)
        with pytest.raises(PrefixMismatchError, match="prefix mismatch"):
            verify_prefix_matches(f, "ENTRAP_SHUF_", declared_n_entrapment=2)

    def test_verify_error_names_the_prefixes_present(self, tmp_path: Path):
        f = tmp_path / "mixed.fasta"
        f.write_text(NONSTANDARD_FASTA)
        with pytest.raises(PrefixMismatchError, match="ENTRAPHUMAN_"):
            verify_prefix_matches(f, "ENTRAP_SHUF_", declared_n_entrapment=2)

    def test_verify_passes_when_prefix_correct(self, tmp_path: Path):
        f = tmp_path / "mixed.fasta"
        f.write_text(NONSTANDARD_FASTA)
        assert verify_prefix_matches(f, ENTRAP_PREFIX, declared_n_entrapment=2) == 2

    def test_zero_declared_never_raises(self, tmp_path: Path):
        """No entrapment declared means nothing to mismatch."""
        f = tmp_path / "mixed.fasta"
        f.write_text(NONSTANDARD_FASTA)
        assert verify_prefix_matches(f, "ENTRAP_SHUF_", declared_n_entrapment=0) == 0

    def test_genuine_zero_hits_is_not_a_prefix_error(self, tmp_path: Path, lake: Path):
        """Entrapment depleted to zero in the searched DB is the DESIRABLE result.

        Pre-selection curation is supposed to remove entrapment, so an empty
        per-sample partition must not be mistaken for a prefix mismatch.
        """
        root = _results_tree(tmp_path, lake)
        for sample in ("sample0", "sample1", "sample2"):
            db = root / sample / "razor.fasta"
            db.write_text(WRAPPED_FASTA)  # zero entrapment survived
        r = compute_fdp(root, target_fasta=lake, db_pattern="razor.fasta")
        assert r["coverage"]["n_empty_entrapment_partition"] == 3

    def test_fdp_raises_on_injected_db_prefix_mismatch(self, tmp_path: Path, lake: Path):
        """A declared partition invisible to the pattern must stop the run."""
        root = _results_tree(tmp_path, lake)
        aug = tmp_path / "aug.fasta"
        aug.write_text(NONSTANDARD_FASTA)
        m = json.loads((root / "entrapment_manifest.json").read_text())
        m.update(
            {
                "output_fasta": str(aug),
                "n_entrapment": 2,
                "entrapment_match_pattern": "ENTRAP_SHUF_",
            }
        )
        (root / "entrapment_manifest.json").write_text(json.dumps(m))
        with pytest.raises(PrefixMismatchError, match="prefix mismatch"):
            compute_fdp(root, target_fasta=lake)

    def test_fdp_records_pattern_and_its_source(self, tmp_path: Path, lake: Path):
        root = _results_tree(tmp_path, lake)
        r = compute_fdp(root, target_fasta=lake, prefix="ENTRAP_SHUF_")
        em = r["entrapment_matching"]
        assert em["pattern"] == "ENTRAP_SHUF_"
        assert "explicit" in em["pattern_source"]

    def test_fdp_uses_manifest_pattern_when_not_overridden(self, tmp_path: Path, lake: Path):
        root = _results_tree(tmp_path, lake)
        m = json.loads((root / "entrapment_manifest.json").read_text())
        m["entrapment_match_pattern"] = "ENTRAPNN_"
        (root / "entrapment_manifest.json").write_text(json.dumps(m))
        r = compute_fdp(root, target_fasta=lake)
        assert r["entrapment_matching"]["pattern"] == "ENTRAPNN_"
        assert "manifest" in r["entrapment_matching"]["pattern_source"]

    def test_unverified_pattern_is_flagged(self, tmp_path: Path, lake: Path):
        root = _results_tree(tmp_path, lake)
        r = compute_fdp(root, target_fasta=lake)
        assert r["entrapment_matching"]["verified_against_injected_database"] is False
        assert "NOT VERIFIED" in r["entrapment_matching"]["note"]

    def test_prefixes_cli(self, tmp_path: Path):
        f = tmp_path / "mixed.fasta"
        f.write_text(NONSTANDARD_FASTA)
        res = CliRunner().invoke(cli, ["entrapment", "prefixes", str(f)])
        assert res.exit_code == 0, res.output
        assert "ENTRAPHUMAN_" in res.output and "ENTRAPNN_" in res.output


# ── config threading ─────────────────────────────────────────────────────────


class TestConfigThreading:
    def test_disabled_is_a_no_op(self, lake: Path, tmp_path: Path):
        cfg = FastaLakeConfig(entrapment={"enabled": False})
        assert prepare_from_config(cfg, lake, tmp_path / "e") is None

    def test_enabled_generates_and_injects(self, lake: Path, tmp_path: Path):
        cfg = FastaLakeConfig(
            entrapment={
                "enabled": True,
                "entrap_class": "shuffled",
                "seed": 42,
            }
        )
        m = prepare_from_config(cfg, lake, tmp_path / "e")
        assert m["class"] == "shuffled"
        assert m["stage"] == STAGE_PRESELECT
        assert Path(m["output_fasta"]).exists()
        assert partition_counts(Path(m["output_fasta"]), ENTRAP_PREFIX) == (3, 3)

    def test_default_stage_is_preselect(self):
        assert EntrapmentConfig().stage == "preselect"

    def test_default_seed_is_42(self):
        assert EntrapmentConfig().seed == 42

    def test_invalid_stage_rejected(self):
        with pytest.raises(ValueError, match="entrapment.stage"):
            EntrapmentConfig(stage="whenever")

    def test_zero_seed_rejected(self):
        with pytest.raises(ValueError, match="nonzero"):
            EntrapmentConfig(seed=0)

    def test_bad_expected_fdr_rejected(self):
        with pytest.raises(ValueError, match="expected_fdr"):
            EntrapmentConfig(expected_fdr=1.5)

    def test_legacy_phylogenetic_maps_to_user(self, lake: Path, tmp_path: Path):
        cfg = FastaLakeConfig(
            entrapment={
                "enabled": True,
                "method": "phylogenetic",
                "source": str(lake),
            }
        )
        m = prepare_from_config(cfg, lake, tmp_path / "e")
        assert m["class"] == "user"

    def test_config_json_roundtrip(self, tmp_path: Path):
        from fasta_lake.config import load_config

        p = tmp_path / "c.json"
        p.write_text(
            json.dumps(
                {
                    "database": {},
                    "entrapment": {
                        "enabled": True,
                        "entrap_class": "plant",
                        "stage": "appended",
                        "seed": 7,
                    },
                    "output": {"directory": str(tmp_path)},
                }
            )
        )
        cfg = load_config(p)
        assert cfg.entrapment.entrap_class == "plant"
        assert cfg.entrapment.stage == "appended"
        assert cfg.entrapment.seed == 7
        assert cfg.has_entrapment() is True


# ── CLI ──────────────────────────────────────────────────────────────────────


class TestEntrapmentCli:
    def test_group_help(self):
        res = CliRunner().invoke(cli, ["entrapment", "--help"])
        assert res.exit_code == 0
        for cmd in ("generate", "inject", "fdp", "classes", "prepare"):
            assert cmd in res.output

    def test_classes_lists_all(self):
        res = CliRunner().invoke(cli, ["entrapment", "classes"])
        assert res.exit_code == 0
        for name in ("shuffled", "thermophilic", "plant", "user"):
            assert name in res.output
        assert "CAVEAT" in res.output

    def test_classes_json(self):
        res = CliRunner().invoke(cli, ["entrapment", "classes", "--json"])
        assert res.exit_code == 0
        data = json.loads(res.output)
        assert data["shuffled"]["recommended_headline"] is False

    def test_generate(self, lake: Path, tmp_path: Path):
        out = tmp_path / "e.fasta"
        res = CliRunner().invoke(
            cli,
            [
                "entrapment",
                "generate",
                "--class",
                "shuffled",
                "-i",
                str(lake),
                "-o",
                str(out),
                "--entrapment-seed",
                "42",
            ],
        )
        assert res.exit_code == 0, res.output
        assert count_fasta_headers(out) == 3
        assert "CAVEAT" in res.output

    def test_inject(self, lake: Path, tmp_path: Path):
        e = tmp_path / "e.fasta"
        generate_entrapment_fasta(e, entrap_class="shuffled", source_fasta=lake, seed=42)
        out = tmp_path / "aug.fasta"
        res = CliRunner().invoke(
            cli,
            [
                "entrapment",
                "inject",
                "-d",
                str(lake),
                "-e",
                str(e),
                "-o",
                str(out),
                "--class",
                "shuffled",
                "--stage",
                "preselect",
            ],
        )
        assert res.exit_code == 0, res.output
        assert "preselect" in res.output
        assert "database-construction error" in res.output
        assert (tmp_path / "entrapment_manifest.json").exists()

    def test_fdp(self, tmp_path: Path, lake: Path):
        root = _results_tree(tmp_path, lake)
        res = CliRunner().invoke(
            cli,
            [
                "entrapment",
                "fdp",
                "--results-dir",
                str(root),
                "-o",
                str(tmp_path / "out"),
                "--target-fasta",
                str(lake),
                "--dataset",
                "testset",
            ],
        )
        assert res.exit_code == 0, res.output
        assert "Coverage" in res.output
        assert "Conserved-peptide correction" in res.output
        assert "cons.corr" in res.output
        assert (tmp_path / "out" / "entrapment_fdp_record.json").exists()

    def test_fdp_refuses_mixed_stages(self, tmp_path: Path, lake: Path):
        root = _results_tree(tmp_path, lake)
        (root / "sample0" / "entrapment_manifest.json").write_text(
            json.dumps({"stage": "appended", "class": "shuffled"})
        )
        res = CliRunner().invoke(
            cli,
            [
                "entrapment",
                "fdp",
                "--results-dir",
                str(root),
                "-o",
                str(tmp_path / "out"),
                "--target-fasta",
                str(lake),
            ],
        )
        assert res.exit_code != 0
        assert "different injection" in str(res.exception)

    def test_prepare_from_config(self, tmp_path: Path, lake: Path):
        cfgp = tmp_path / "c.json"
        cfgp.write_text(
            json.dumps(
                {
                    "database": {},
                    "entrapment": {"enabled": True, "entrap_class": "shuffled", "seed": 42},
                    "output": {"directory": str(tmp_path)},
                }
            )
        )
        res = CliRunner().invoke(
            cli,
            [
                "entrapment",
                "prepare",
                "-c",
                str(cfgp),
                "-d",
                str(lake),
                "-o",
                str(tmp_path / "ent"),
            ],
        )
        assert res.exit_code == 0, res.output
        assert (tmp_path / "ent" / "augmented_database.fasta").exists()
        assert (tmp_path / "ent" / "entrapment_manifest.json").exists()

    def test_extract_sample_threading_disabled_passes_db_through(self, tmp_path: Path, lake: Path):
        from fasta_lake.cli import _maybe_inject_entrapment

        cfgp = tmp_path / "c.json"
        cfgp.write_text(
            json.dumps(
                {
                    "database": {},
                    "entrapment": {"enabled": False},
                    "output": {"directory": str(tmp_path)},
                }
            )
        )
        got = _maybe_inject_entrapment(str(cfgp), str(lake), str(tmp_path / "o.fasta"), None)
        assert got == str(lake)

    def test_extract_sample_threading_injects_when_enabled(self, tmp_path: Path, lake: Path):
        """`enabled: true` must change what stage 2 actually searches."""
        from fasta_lake.cli import _maybe_inject_entrapment

        cfgp = tmp_path / "c.json"
        cfgp.write_text(
            json.dumps(
                {
                    "database": {},
                    "entrapment": {"enabled": True, "entrap_class": "shuffled", "seed": 42},
                    "output": {"directory": str(tmp_path)},
                }
            )
        )
        got = _maybe_inject_entrapment(
            str(cfgp),
            str(lake),
            str(tmp_path / "o.fasta"),
            str(tmp_path / "ent"),
        )
        assert got != str(lake)
        assert partition_counts(Path(got), ENTRAP_PREFIX) == (3, 3)
        assert (tmp_path / "ent" / "entrapment_manifest.json").exists()

    def test_extract_sample_threading_refuses_appended_stage(self, tmp_path: Path, lake: Path):
        """`appended` happens after parsimony, not before extraction."""
        import click as _click

        from fasta_lake.cli import _maybe_inject_entrapment

        cfgp = tmp_path / "c.json"
        cfgp.write_text(
            json.dumps(
                {
                    "database": {},
                    "entrapment": {
                        "enabled": True,
                        "entrap_class": "shuffled",
                        "stage": "appended",
                    },
                    "output": {"directory": str(tmp_path)},
                }
            )
        )
        with pytest.raises(_click.ClickException, match="AFTER parsimony"):
            _maybe_inject_entrapment(
                str(cfgp),
                str(lake),
                str(tmp_path / "o.fasta"),
                None,
            )

    def test_prepare_disabled_is_noop(self, tmp_path: Path, lake: Path):
        cfgp = tmp_path / "c.json"
        cfgp.write_text(
            json.dumps(
                {
                    "database": {},
                    "entrapment": {"enabled": False},
                    "output": {"directory": str(tmp_path)},
                }
            )
        )
        res = CliRunner().invoke(
            cli,
            [
                "entrapment",
                "prepare",
                "-c",
                str(cfgp),
                "-d",
                str(lake),
                "-o",
                str(tmp_path / "ent"),
            ],
        )
        assert res.exit_code == 0
        assert "nothing to do" in res.output
