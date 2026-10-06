"""The real-data checksum gate must reject incomplete integrity coverage."""

import hashlib
import importlib.util
from pathlib import Path

import pytest


def load_example(monkeypatch):
    directory = Path(__file__).resolve().parents[1] / "examples/campi"
    monkeypatch.syspath_prepend(str(directory))
    spec = importlib.util.spec_from_file_location("review_example", directory / "run_example.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize("defect", ["empty", "omitted", "duplicate"])
def test_incomplete_checksum_coverage_fails_before_execution(monkeypatch, tmp_path, defect):
    module = load_example(monkeypatch)
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    (inputs / "one.txt").write_bytes(b"one")
    line = hashlib.sha256(b"one").hexdigest() + "  one.txt\n"
    manifest = "" if defect == "empty" else line
    if defect == "duplicate":
        manifest += line
    if defect == "omitted":
        (inputs / "unverified.txt").write_bytes(b"unverified")
    (inputs / "SHA256SUMS").write_text(manifest)
    monkeypatch.setattr(module, "HERE", tmp_path)
    with pytest.raises(ValueError, match="checksum"):
        module.verify_inputs()


def test_real_campi_input_manifest_is_complete(monkeypatch):
    load_example(monkeypatch).verify_inputs()


def test_example_rejects_unsupported_qc_python_before_setup(monkeypatch, tmp_path, capsys):
    module = load_example(monkeypatch)
    destination = tmp_path / "result"
    monkeypatch.setattr(module.sys, "argv", ["run_example.py", "--out", str(destination)])
    # The locked package set is compiled for Python 3.12: older AND newer refuse.
    for version_info in ((3, 10, 14), (3, 11, 9), (3, 13, 2)):
        with monkeypatch.context() as version:
            version.setattr(module.sys, "version_info", version_info)
            with pytest.raises(SystemExit) as error:
                module.main()
        assert error.value.code == 2
        assert "Python 3.12 is required" in capsys.readouterr().err
    assert not destination.exists()


def test_bundled_sage_is_the_official_v0_14_7_release(monkeypatch):
    """v0.14.6 and v0.14.7 both report 'sage 0.14.6'; only checksums tell them apart."""
    import tarfile

    from fasta_lake.search import SAGE_EXECUTABLE_SHA256, SAGE_RELEASE, SAGE_REPORTED_VERSION

    module = load_example(monkeypatch)
    archive = module.ROOT / "examples/campi/runtime" / module.SAGE_ARCHIVE
    assert module.SAGE_RELEASE == SAGE_RELEASE == "v0.14.7"
    assert module.SAGE_ARCHIVE == "sage-v0.14.7-x86_64-unknown-linux-gnu.tar.gz"
    assert module.SAGE_SHA256 == (
        "e3dc6b41015cb167574f6c82525b75e946c094f30bd700271b05c051c30cbe8a"
    )
    assert (
        module.SAGE_EXECUTABLE_SHA256
        == SAGE_EXECUTABLE_SHA256
        == "d3820543a31bfa2a556e04f91719204ce0ac6f34dcca6670314a13df87064625"
    )
    assert module.SAGE_REPORTED_VERSION == "sage " + SAGE_REPORTED_VERSION
    assert module.sha(archive) == module.SAGE_SHA256
    with tarfile.open(archive) as bundle:
        executable = bundle.extractfile(bundle.getmember(module.SAGE_MEMBER)).read()
    assert hashlib.sha256(executable).hexdigest() == SAGE_EXECUTABLE_SHA256
    # The official v0.14.6 executable must not satisfy the pin.
    assert SAGE_EXECUTABLE_SHA256 != (
        "7367d53f5f25d21b8ef182f69a1dd14682a93d551acdf0184a28c50f84dc9f9f"
    )
