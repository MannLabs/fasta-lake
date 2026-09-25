"""Interruptions, wrong resources and hostile archives must never look complete."""

import gzip
import hashlib
import io
import json
import lzma
import tarfile
import threading
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

import pytest
from click.testing import CliRunner

from fasta_lake.download_cli import resources
from fasta_lake.downloads import (
    catalogue,
    download_file,
    download_resource,
    inspect_remote,
    resolve_resource,
    unpack_resource,
)
from fasta_lake.embeddings import resolve_checkpoint


@pytest.fixture
def server():
    state = {
        "body": b">a\nPEPTIDE\n" * 100,
        "etag": '"release1"',
        "ranges": [],
        "truncate": False,
        "ignore_range": False,
        "html": False,
        "head_405": False,
    }

    class Handler(BaseHTTPRequestHandler):
        def do_HEAD(self):
            if state["head_405"]:
                self.send_error(405)
                return
            self.send_response(200)
            self.send_header("Content-Length", str(len(state["body"])))
            self.send_header("ETag", state["etag"])
            self.send_header(
                "Content-Type", "text/html" if state["html"] else "application/octet-stream"
            )
            self.end_headers()

        def do_GET(self):
            body = state["body"]
            selected = self.headers.get("Range")
            state["ranges"].append(selected)
            start, stop = 0, len(body)
            if selected and not state["ignore_range"]:
                start_text, end_text = selected.removeprefix("bytes=").split("-")
                start = int(start_text)
                stop = int(end_text) + 1 if end_text else len(body)
                self.send_response(206)
                self.send_header("Content-Range", f"bytes {start}-{stop - 1}/{len(body)}")
            else:
                self.send_response(200)
            self.send_header("Content-Length", str(stop - start))
            self.send_header("ETag", state["etag"])
            self.send_header(
                "Content-Type", "text/html" if state["html"] else "application/octet-stream"
            )
            self.end_headers()
            self.wfile.write(body[start : stop // 2 if state["truncate"] else stop])

        def log_message(self, *args):
            pass

    httpd = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    thread = threading.Thread(target=httpd.serve_forever, daemon=True)
    thread.start()
    spec = {
        "url": f"http://127.0.0.1:{httpd.server_port}/proteins.gz",
        "name": "proteins.gz",
        "bytes": len(state["body"]),
        "checksum": "sha256:" + hashlib.sha256(state["body"]).hexdigest(),
    }
    yield state, spec
    httpd.shutdown()
    httpd.server_close()
    thread.join()


def test_stream_hash_and_no_overwrite(server, tmp_path):
    state, spec = server
    path = tmp_path / spec["name"]
    result = download_file(spec, path)
    assert path.read_bytes() == state["body"]
    assert result["upstream_checksum_verified"]
    assert result["hashes"]["sha256"] == hashlib.sha256(state["body"]).hexdigest()
    with pytest.raises(FileExistsError):
        download_file(spec, path)
    assert not Path(str(path) + ".lock").exists()


def test_truncation_resume_exact_bytes(server, tmp_path):
    state, spec = server
    state["truncate"] = True
    path = tmp_path / spec["name"]
    with pytest.raises(ValueError, match="Truncated"):
        download_file(spec, path)
    assert not path.exists()
    partial = Path(str(path) + ".partial")
    assert 0 < partial.stat().st_size < len(state["body"])
    offset = partial.stat().st_size
    state["truncate"] = False
    download_file(spec, path, resume=True)
    assert state["ranges"][-1] == f"bytes={offset}-"
    assert path.read_bytes() == state["body"]
    assert not partial.exists()


@pytest.mark.parametrize("change", ["etag", "weak", "ignore"])
def test_resume_refuses_changed_or_unsupported_resource(server, tmp_path, change):
    state, spec = server
    state["truncate"] = True
    path = tmp_path / spec["name"]
    if change == "weak":
        state["etag"] = 'W/"release1"'
    with pytest.raises(ValueError):
        download_file(spec, path)
    state["truncate"] = False
    if change == "etag":
        state["etag"] = '"release2"'
    if change == "ignore":
        state["ignore_range"] = True
    with pytest.raises(ValueError):
        download_file(spec, path, resume=True)
    assert not path.exists()


def test_wrong_checksum_never_publishes(server, tmp_path):
    _, spec = server
    spec["checksum"] = "sha256:" + "0" * 64
    with pytest.raises(ValueError, match="Checksum mismatch"):
        download_file(spec, tmp_path / spec["name"])
    assert not (tmp_path / spec["name"]).exists()
    assert (tmp_path / (spec["name"] + ".partial")).is_file()


def test_html_and_changed_size_rejected(server):
    state, spec = server
    state["html"] = True
    with pytest.raises(ValueError, match="HTML"):
        inspect_remote(spec)
    state["html"] = False
    spec["bytes"] += 1
    with pytest.raises(ValueError, match="size changed"):
        inspect_remote(spec)


def test_head_fallback_uses_range_metadata(server):
    state, spec = server
    state["head_405"] = True
    assert inspect_remote(spec)["bytes"] == len(state["body"])
    assert state["ranges"] == ["bytes=0-0"]


def test_resource_receipts_resume_and_tamper_detection(server, tmp_path):
    state, spec = server
    spec["remote"] = inspect_remote(spec)
    plan = {"id": "test", "version": "v1", "destination": str(tmp_path), "files": [spec]}
    result = download_resource(plan)
    assert result["status"] == "PASS"
    assert (tmp_path / "test/DOWNLOAD_COMPLETE.json").is_file()
    assert download_resource(plan, resume=True) == result
    (tmp_path / "test/proteins.gz").write_bytes(state["body"] + b"changed")
    with pytest.raises(ValueError, match="Checksum"):
        download_resource(plan, resume=True)


def test_existing_unowned_file_not_accepted(server, tmp_path):
    _, spec = server
    out = tmp_path / "test"
    out.mkdir()
    plan = {"id": "test", "version": "v1", "destination": str(tmp_path), "files": [spec]}
    (out / "DOWNLOAD_PLAN.json").write_text(json.dumps(plan))
    (out / spec["name"]).write_text("unverified")
    with pytest.raises(ValueError, match="completion receipt"):
        download_resource(plan, resume=True)


def test_uniprot_resolves_current_checksum(monkeypatch):
    xml = b"""<metalink xmlns="http://www.metalinker.org/"><version>release-test</version>
      <files><file name="uniprot_sprot.fasta.gz"><size>123</size><verification>
      <hash type="md5">12345678901234567890123456789012</hash></verification></file></files>
      </metalink>"""
    monkeypatch.setattr("fasta_lake.downloads._read_metadata", lambda url: xml)
    spec = resolve_resource("uniprot-swissprot")
    assert spec["version"] == "release-test"
    assert spec["files"][0]["bytes"] == 123
    assert spec["files"][0]["checksum"].startswith("md5:")


def test_gzip_and_xz_unpack_preserve_archives(tmp_path):
    a, b = tmp_path / "a.faa.gz", tmp_path / "b.tsv.xz"
    a.write_bytes(gzip.compress(b">p\nAAAA\n"))
    b.write_bytes(lzma.compress(b"gene\tKO\n"))
    before = a.read_bytes(), b.read_bytes()
    result = unpack_resource([a, b], tmp_path / "out", max_bytes=1000)
    assert result["expanded_bytes"] == 16
    assert (a.read_bytes(), b.read_bytes()) == before
    assert (tmp_path / "out/a.faa").read_bytes() == b">p\nAAAA\n"


def make_tar(path, name="data/proteins.faa", kind=tarfile.REGTYPE, duplicate=False):
    with tarfile.open(path, "w:gz") as archive:
        member = tarfile.TarInfo(name)
        member.type = kind
        if kind == tarfile.REGTYPE:
            member.size = 10
            archive.addfile(member, io.BytesIO(b">p\nAAAAAA\n"))
            if duplicate:
                archive.addfile(member, io.BytesIO(b">p\nAAAAAA\n"))
        else:
            member.linkname = "../../outside"
            archive.addfile(member)


@pytest.mark.parametrize("name", ["../../escape", "/absolute", "data/../../escape", "data\\escape"])
def test_tar_traversal_rejected(tmp_path, name):
    path = tmp_path / "bad.tar.gz"
    make_tar(path, name)
    with pytest.raises(ValueError, match="Unsafe archive"):
        unpack_resource([path], tmp_path / "out", max_bytes=1000)
    assert not (tmp_path / "out/UNPACK_COMPLETE.json").exists()


@pytest.mark.parametrize("kind", [tarfile.SYMTYPE, tarfile.LNKTYPE, tarfile.FIFOTYPE])
def test_tar_links_and_special_files_rejected(tmp_path, kind):
    path = tmp_path / "bad.tar.gz"
    make_tar(path, kind=kind)
    with pytest.raises(ValueError, match="link or special"):
        unpack_resource([path], tmp_path / "out", max_bytes=1000)


def test_tar_duplicate_and_budget_fail(tmp_path):
    path = tmp_path / "duplicate.tar.gz"
    make_tar(path, duplicate=True)
    with pytest.raises(FileExistsError):
        unpack_resource([path], tmp_path / "out", max_bytes=1000)
    with pytest.raises(OSError, match="budget"):
        unpack_resource([path], tmp_path / "small", max_bytes=5)
    assert not (tmp_path / "small/UNPACK_COMPLETE.json").exists()


def test_expansion_budget_stops_gzip_bomb(tmp_path):
    path = tmp_path / "big.gz"
    path.write_bytes(gzip.compress(b"A" * 100000))
    with pytest.raises(OSError, match="budget"):
        unpack_resource([path], tmp_path / "out", max_bytes=100)
    assert not (tmp_path / "out/UNPACK_COMPLETE.json").exists()


def test_resource_menu_is_offline_and_includes_requested_sources(monkeypatch):
    def unexpected(*args, **kwargs):
        raise AssertionError("offline menu tried the network")

    monkeypatch.setattr("fasta_lake.downloads.urlopen", unexpected)
    result = CliRunner().invoke(resources, ["list"])
    assert result.exit_code == 0, result.output
    for name in ("gmsc-90", "gmgc-95nr", "uhgp-95", "uniprot-human", "eggnog-5.0.2", "esmc-300m"):
        assert name in result.output
    assert len(catalogue()) >= 17


def test_local_checkpoint_required_before_any_model_download(tmp_path):
    with pytest.raises(FileNotFoundError, match="resources download esmc-300m"):
        resolve_checkpoint(None, "esmc_300m")
    path = tmp_path / "wrong.pth"
    path.write_bytes(b"wrong model")
    with pytest.raises(ValueError, match="published version"):
        resolve_checkpoint(path, "esmc_300m")


def test_multifile_space_check_prevents_any_transfer(tmp_path, monkeypatch):
    plan = {
        "id": "two-files",
        "version": "1",
        "destination": str(tmp_path),
        "files": [{"name": name, "bytes": 12} for name in ("first.gz", "second.gz")],
    }
    monkeypatch.setattr("fasta_lake.downloads.free_bytes", lambda path: 20)

    def unexpected(*args, **kwargs):
        raise AssertionError("Started a download that cannot fit")

    monkeypatch.setattr("fasta_lake.downloads.download_file", unexpected)
    with pytest.raises(OSError, match="complete remaining resource"):
        download_resource(plan)
    assert not (tmp_path / "two-files").exists()
