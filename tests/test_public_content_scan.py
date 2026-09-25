"""Prevent path leaks through metadata and archives missed by extension filters."""

import importlib.util
import io
import zipfile
from pathlib import Path


def scanner():
    path = Path(__file__).resolve().parents[1] / "tools/check_no_site_paths.py"
    spec = importlib.util.spec_from_file_location("content_scan", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_mzml_checks_metadata_without_matching_numeric_binary_payload():
    module = scanner()
    private = b"/" + b"home/example/private/run.raw"
    assert module.inspect_payload(b'<sourceFile location="' + private + b'"/>', "run.mzML")
    assert not module.inspect_payload(b"<binary>" + private + b"</binary>", "run.mzML")


def test_nested_archive_and_notebook_output_are_checked_without_echoing_secret():
    module = scanner()
    token = b"gh" + b"p_" + b"a" * 36
    inner = io.BytesIO()
    with zipfile.ZipFile(inner, "w") as archive:
        archive.writestr("result.ipynb", b'{"text": "' + token + b'"}')
    outer = io.BytesIO()
    with zipfile.ZipFile(outer, "w") as archive:
        archive.writestr("nested.zip", inner.getvalue())
    findings = module.inspect_payload(outer.getvalue(), "example.zip")
    assert any("nested.zip!result.ipynb" in hit for hit in findings)
    assert all(token.decode() not in hit for hit in findings)


def test_relative_and_container_paths_remain_portable():
    assert not scanner().inspect_payload(
        b"data/sample.mzML /opt/fastalake /usr/local/bin /tmp/example", "guide.md"
    )
