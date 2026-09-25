"""Integration tests: engine × mode × strategy matrix.

Tests the full pipeline flow from inference FASTA → engine-specific export,
verifying that each search engine gets correctly formatted input.

These tests use the tiny test data in tests/data/ and do NOT require
external tools (SAGE, DIA-NN, MSFragger) — they validate the FASTA output
format that each engine expects.
"""

import os
import subprocess
import sys
import tempfile

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

from fasta_lake.config import VALID_STRATEGIES, FastaLakeConfig
from fasta_lake.headers import (
    prepare_fasta_for_engine,
)
from fasta_lake.headers.validation import detect_experiment

# ═══════════════════════════════════════════════════════════════════════════════
# Test data: representative headers from each format
# ═══════════════════════════════════════════════════════════════════════════════

MULTI_FORMAT_FASTA = """\
>sp|P02768|ALBU_HUMAN Albumin OS=Homo sapiens OX=9606 GN=ALB PE=1 SV=2
MKWVTFISLLFLFSSAYSRGVFRRDAHK
>tr|A0A218KGR2|A0A218KGR2_HUMAN Amyloid-beta A4 protein OS=Homo sapiens GN=APP
MLPGLALLLLAAWTARALEVPTDGNAGLLAEPQIAMFCGRLNMHMNVQNG
>MGYG000001373_00501 Uncharacterized protein
MVKRTILGDTFKLPNMRVHPEDGQVSYVVHKNDEFGVIFDHGELYGEVLR
>GMGC10.050_359_667.DNAK DnaK protein
MGKIIGIDLGTTNSCVAVMEGKNPKVIENAEG
>var-nfe|P02768|ALB_Glu177Asp|4-73409403-A-T sp|P02768|ALBU_HUMAN Albumin GN=ALB
MKWVTFISLLFLFSSAYSRGVFRRDA
>k141_130475_1 # 1 # 1404 # 1 # ID=1_1;partial=00
MVKRTILGDTFKLPNMRVHPEDGQVS
"""


def _write_multi_format(path):
    with open(path, "w") as f:
        f.write(MULTI_FORMAT_FASTA)


def _count_proteins(path):
    with open(path) as f:
        return sum(1 for line in f if line.startswith(">"))


def _read_headers(path):
    with open(path) as f:
        return [line[1:].rstrip() for line in f if line.startswith(">")]


# ═══════════════════════════════════════════════════════════════════════════════
# Export tests: one per engine
# ═══════════════════════════════════════════════════════════════════════════════


class TestExportSAGE:
    """SAGE: GN= tags, NO decoys (SAGE generates internal decoys)."""

    def test_sage_adds_gn(self):
        with tempfile.TemporaryDirectory() as d:
            inp = os.path.join(d, "input.fasta")
            out = os.path.join(d, "sage.fasta")
            _write_multi_format(inp)
            stats = prepare_fasta_for_engine(inp, out, engine="sage")
            assert stats.target_proteins == 6
            assert stats.decoy_proteins == 0
            for hdr in _read_headers(out):
                assert "GN=" in hdr, f"Missing GN= in: {hdr}"

    def test_sage_no_decoys(self):
        with tempfile.TemporaryDirectory() as d:
            inp = os.path.join(d, "input.fasta")
            out = os.path.join(d, "sage.fasta")
            _write_multi_format(inp)
            prepare_fasta_for_engine(inp, out, engine="sage")
            for hdr in _read_headers(out):
                assert not hdr.startswith("rev_"), f"Unexpected decoy in SAGE output: {hdr}"

    def test_sage_swissprot_unchanged(self):
        with tempfile.TemporaryDirectory() as d:
            inp = os.path.join(d, "input.fasta")
            out = os.path.join(d, "sage.fasta")
            _write_multi_format(inp)
            prepare_fasta_for_engine(inp, out, engine="sage")
            headers = _read_headers(out)
            sp_header = [h for h in headers if h.startswith("sp|P02768")][0]
            # sp| headers should be unchanged (already have GN=)
            assert "Albumin" in sp_header
            assert "GN=ALB" in sp_header


class TestExportDIANN:
    """DIA-NN: gene as FIRST WORD of description for non-sp/tr."""

    def test_diann_gene_first_word_mgyg(self):
        with tempfile.TemporaryDirectory() as d:
            inp = os.path.join(d, "input.fasta")
            out = os.path.join(d, "diann.fasta")
            _write_multi_format(inp)
            prepare_fasta_for_engine(inp, out, engine="diann")
            headers = _read_headers(out)
            mgyg = [h for h in headers if h.startswith("MGYG")][0]
            parts = mgyg.split()
            assert parts[1] == "MGYG000001373_00501", (
                f"Gene not first word for MGYG: got {parts[1]}"
            )

    def test_diann_gene_first_word_gmgc(self):
        with tempfile.TemporaryDirectory() as d:
            inp = os.path.join(d, "input.fasta")
            out = os.path.join(d, "diann.fasta")
            _write_multi_format(inp)
            prepare_fasta_for_engine(inp, out, engine="diann")
            headers = _read_headers(out)
            gmgc = [h for h in headers if h.startswith("GMGC")][0]
            parts = gmgc.split()
            # Gene is placed as first word of description (may be full accession for GMGC)
            assert parts[1] != "DnaK", "Description should NOT start with raw description"
            assert "GN=" in gmgc, "GN= tag must be present"

    def test_diann_swissprot_not_modified(self):
        with tempfile.TemporaryDirectory() as d:
            inp = os.path.join(d, "input.fasta")
            out = os.path.join(d, "diann.fasta")
            _write_multi_format(inp)
            prepare_fasta_for_engine(inp, out, engine="diann")
            headers = _read_headers(out)
            sp = [h for h in headers if h.startswith("sp|P02768")][0]
            # sp| headers: DIA-NN reads GN= tag natively
            assert sp.startswith("sp|P02768|ALBU_HUMAN Albumin")

    def test_diann_no_decoys(self):
        with tempfile.TemporaryDirectory() as d:
            inp = os.path.join(d, "input.fasta")
            out = os.path.join(d, "diann.fasta")
            _write_multi_format(inp)
            stats = prepare_fasta_for_engine(inp, out, engine="diann")
            assert stats.decoy_proteins == 0


class TestExportMSFragger:
    """MSFragger: GN= tags + reversed decoy sequences with rev_ prefix."""

    def test_msfragger_has_decoys(self):
        with tempfile.TemporaryDirectory() as d:
            inp = os.path.join(d, "input.fasta")
            out = os.path.join(d, "msfragger.fasta")
            _write_multi_format(inp)
            stats = prepare_fasta_for_engine(inp, out, engine="msfragger")
            assert stats.target_proteins == 6
            assert stats.decoy_proteins == 6
            assert stats.total_proteins == 12

    def test_msfragger_decoy_count_equals_target(self):
        with tempfile.TemporaryDirectory() as d:
            inp = os.path.join(d, "input.fasta")
            out = os.path.join(d, "msfragger.fasta")
            _write_multi_format(inp)
            prepare_fasta_for_engine(inp, out, engine="msfragger")
            headers = _read_headers(out)
            targets = [h for h in headers if not h.startswith("rev_")]
            decoys = [h for h in headers if h.startswith("rev_")]
            assert len(targets) == len(decoys)

    def test_msfragger_decoy_sequences_reversed(self):
        with tempfile.TemporaryDirectory() as d:
            inp = os.path.join(d, "input.fasta")
            out = os.path.join(d, "msfragger.fasta")
            _write_multi_format(inp)
            prepare_fasta_for_engine(inp, out, engine="msfragger")

            seqs = {}
            with open(out) as f:
                current = None
                for line in f:
                    if line.startswith(">"):
                        current = line[1:].rstrip()
                        seqs[current] = ""
                    elif current:
                        seqs[current] += line.rstrip()

            # Check that ALBU decoy is reversed
            target_key = [
                k for k in seqs if k.startswith("sp|P02768") and not k.startswith("rev_")
            ][0]
            decoy_key = [k for k in seqs if k.startswith("rev_sp|P02768")][0]
            assert seqs[decoy_key] == seqs[target_key][::-1]


# ═══════════════════════════════════════════════════════════════════════════════
# Config validation tests
# ═══════════════════════════════════════════════════════════════════════════════


class TestConfigEngineValidation:
    """Test config validates acquisition + search_engine combinations."""

    def test_dda_sage_valid(self):
        cfg = FastaLakeConfig(acquisition="dda", search_engine="sage")
        assert cfg.acquisition == "dda"
        assert cfg.search_engine == "sage"

    def test_dia_diann_valid(self):
        cfg = FastaLakeConfig(acquisition="dia", search_engine="diann")
        assert cfg.acquisition == "dia"
        assert cfg.search_engine == "diann"

    def test_dda_msfragger_valid(self):
        cfg = FastaLakeConfig(acquisition="dda", search_engine="msfragger")
        assert cfg.search_engine == "msfragger"

    def test_invalid_acquisition(self):
        try:
            FastaLakeConfig(acquisition="maldi")
            assert False, "Should have raised ValueError"
        except ValueError as e:
            assert "maldi" in str(e)

    def test_invalid_engine(self):
        try:
            FastaLakeConfig(search_engine="maxquant")
            assert False, "Should have raised ValueError"
        except ValueError as e:
            assert "maxquant" in str(e)

    def test_all_strategies_valid_config(self):
        """All 6 inference strategies should be valid in config."""
        for strategy in VALID_STRATEGIES:
            cfg = FastaLakeConfig()
            cfg.inference.strategy = strategy


# ═══════════════════════════════════════════════════════════════════════════════
# Experiment auto-detection tests
# ═══════════════════════════════════════════════════════════════════════════════


class TestExperimentDetection:
    """Test auto-detection of experiment type from inputs."""

    def test_detect_alphanovo_dda(self):
        with tempfile.TemporaryDirectory() as d:
            csv_path = os.path.join(d, "sample.csv")
            with open(csv_path, "w") as f:
                f.write("peptide_prediction_detokenized_unmodified,score\n")
                f.write("PEPTIDEK,0.95\n")
            det = detect_experiment(predictions_dir=d)
            assert det.denovo_engine == "alphanovo"
            assert det.acquisition == "dda"
            assert det.search_engine == "sage"

    def test_detect_dianovo_dia(self):
        with tempfile.TemporaryDirectory() as d:
            csv_path = os.path.join(d, "sample.csv")
            with open(csv_path, "w") as f:
                f.write("graph_idx,pred_path,pred_prob,pred_seq\n")
                f.write("0,path,0.9 0.8,0.0 264.11 V V A\n")
            det = detect_experiment(predictions_dir=d)
            assert det.denovo_engine == "dianovo"
            assert det.acquisition == "dia"
            assert det.search_engine == "diann"
            assert det.needs_header_healing is True


# ═══════════════════════════════════════════════════════════════════════════════
# Rust binary integration tests
# ═══════════════════════════════════════════════════════════════════════════════


class TestRustFastaExport:
    """Test the Rust fasta_export binary."""

    BINARY = os.path.join(
        os.path.dirname(os.path.dirname(__file__)),
        "rust",
        "fasta_export",
        "target",
        "release",
        "fasta_export",
    )

    def test_binary_exists(self):
        """The binary is a BUILD ARTEFACT, so its absence is not a defect.

        This asserted unconditionally and was the single failure on a fresh clone
        of main (2026-08-03: 618 passed, 19 skipped, 1 failed). A suite that is red
        on a clean checkout for a reason that is not a bug trains everyone to
        ignore red, so this skips with the build command instead.
        """
        if not os.path.isfile(self.BINARY):
            pytest.skip(
                "fasta_export not built; run: cargo build --release "
                "--manifest-path rust/fasta_export/Cargo.toml"
            )
        assert os.path.isfile(self.BINARY)

    def test_export_sage(self):
        if not os.path.isfile(self.BINARY):
            return
        with tempfile.TemporaryDirectory() as d:
            inp = os.path.join(d, "input.fasta")
            out = os.path.join(d, "sage.fasta")
            _write_multi_format(inp)
            result = subprocess.run(
                [self.BINARY, "export", "-i", inp, "-o", out, "-e", "sage"],
                capture_output=True,
                text=True,
            )
            assert result.returncode == 0, f"STDERR: {result.stderr}"
            assert _count_proteins(out) == 6  # no decoys
            for hdr in _read_headers(out):
                assert "GN=" in hdr

    def test_export_msfragger(self):
        if not os.path.isfile(self.BINARY):
            return
        with tempfile.TemporaryDirectory() as d:
            inp = os.path.join(d, "input.fasta")
            out = os.path.join(d, "msfragger.fasta")
            _write_multi_format(inp)
            result = subprocess.run(
                [self.BINARY, "export", "-i", inp, "-o", out, "-e", "msfragger"],
                capture_output=True,
                text=True,
            )
            assert result.returncode == 0, f"STDERR: {result.stderr}"
            assert _count_proteins(out) == 12  # 6 targets + 6 decoys
            headers = _read_headers(out)
            assert sum(1 for h in headers if h.startswith("rev_")) == 6

    def test_export_diann(self):
        if not os.path.isfile(self.BINARY):
            return
        with tempfile.TemporaryDirectory() as d:
            inp = os.path.join(d, "input.fasta")
            out = os.path.join(d, "diann.fasta")
            _write_multi_format(inp)
            result = subprocess.run(
                [self.BINARY, "export", "-i", inp, "-o", out, "-e", "diann"],
                capture_output=True,
                text=True,
            )
            assert result.returncode == 0, f"STDERR: {result.stderr}"
            assert _count_proteins(out) == 6
            headers = _read_headers(out)
            mgyg = [h for h in headers if h.startswith("MGYG")][0]
            parts = mgyg.split()
            # Gene should be first word after accession
            assert parts[1] == "MGYG000001373_00501"
