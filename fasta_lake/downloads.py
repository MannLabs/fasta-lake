"""Explicit, resumable resource downloads with local integrity records.

Downloads never run as a side effect of database construction or annotation.
Archive expansion is a separate bounded operation. Public functions accept plain
file specifications so callers can test them without fetching large catalogues.
"""

from __future__ import annotations

import copy
import gzip
import hashlib
import json
import lzma
import os
import re
import shutil
import tarfile
import time
import xml.etree.ElementTree as ET
from contextlib import contextmanager
from datetime import datetime, timezone
from importlib.resources import files
from pathlib import Path, PurePosixPath
from urllib.error import HTTPError
from urllib.parse import urlsplit
from urllib.request import Request, urlopen

BLOCK = 1024 * 1024
HEADERS = {"User-Agent": "FastaLake-resource-download/1", "Accept-Encoding": "identity"}


def catalogue():
    """Return the installed, dated resource catalogue without making network calls."""
    return json.loads(files("fasta_lake").joinpath("configs/download_resources.json").read_text())[
        "resources"
    ]


def _read_metadata(url):
    """Fetch at most 2 MiB of resource metadata, raising for a larger response."""
    with urlopen(Request(url, headers=HEADERS), timeout=30) as response:
        data = response.read(2 * BLOCK + 1)
    if len(data) > 2 * BLOCK:
        raise ValueError("Resource metadata exceeds 2 MiB")
    return data


def resolve_resource(name):
    """Resolve a catalogue choice, including live release identity for UniProt.

    No sequence, archive or model data is fetched here. Fixed-version resources
    retain their registered identities. UniProt's moving release is frozen in
    the returned download plan using its upstream metadata and checksum.
    """
    known = catalogue()
    if name not in known:
        raise ValueError(f"Unknown resource {name!r}; run 'fasta-lake resources list'")
    spec = copy.deepcopy(known[name])
    spec["id"] = name
    if "metalink" in spec:
        root = ET.fromstring(_read_metadata(spec["metalink"]))
        ns = {"m": "http://www.metalinker.org/"}
        spec["version"] = root.findtext("m:version", namespaces=ns)
        for item in spec["files"]:
            matches = [
                x for x in root.findall("m:files/m:file", ns) if x.get("name") == item["name"]
            ]
            if len(matches) != 1 or not spec["version"]:
                raise ValueError("UniProt release metadata does not identify the requested file")
            item["bytes"] = int(matches[0].findtext("m:size", namespaces=ns))
            checksum = matches[0].find("m:verification/m:hash[@type='md5']", ns)
            if checksum is None or not re.fullmatch(r"[a-fA-F0-9]{32}", checksum.text or ""):
                raise ValueError("UniProt release metadata lacks its MD5 digest")
            item["checksum"] = "md5:" + checksum.text.lower()
    elif "release_url" in spec:
        spec["version"] = _read_metadata(spec["release_url"]).decode().strip()
        # Reference proteome sizes can change independently of the saved listing.
        for item in spec["files"]:
            item["bytes"] = None
    return spec


def inspect_remote(item):
    """Read transfer metadata without downloading the resource body.

    Unknown size stays unknown. A refused HEAD request may fall back to a one-byte
    range request; servers ignoring Range are closed without reading the body.
    """
    request = Request(item["url"], headers=HEADERS, method="HEAD")
    try:
        response = urlopen(request, timeout=10)
    except (HTTPError, TimeoutError) as error:
        if isinstance(error, HTTPError) and error.code not in {403, 405, 501}:
            raise
        response = urlopen(
            Request(item["url"], headers={**HEADERS, "Range": "bytes=0-0"}), timeout=30
        )
    with response:
        headers = response.headers
        if headers.get_content_type() == "text/html":
            raise ValueError(f"Server returned an HTML page instead of {item['name']}")
        length = headers.get("Content-Length")
        if response.status == 206:
            match = re.fullmatch(r"bytes 0-0/(\d+)", headers.get("Content-Range", ""))
            if not match:
                raise ValueError("Invalid byte-range metadata")
            length = match[1]
        size = int(length) if length is not None else None
        if item.get("bytes") is not None and size is not None and size != item["bytes"]:
            raise ValueError(
                f"Remote size changed for {item['name']}; update the resource catalogue"
            )
        return {
            "bytes": size if size is not None else item.get("bytes"),
            "etag": headers.get("ETag"),
            "last_modified": headers.get("Last-Modified"),
        }


def free_bytes(destination):
    """Read free disk space from the closest existing destination parent."""
    path = Path(destination).absolute()
    while not path.exists():
        if path.is_symlink():
            raise ValueError(f"Dangling destination link: {path}")
        path = path.parent
    return shutil.disk_usage(path).free


def plan_download(name, destination):
    """Resolve files, transfer bytes, disk availability and expanded-size uncertainty."""
    spec = resolve_resource(name)
    for item in spec["files"]:
        item["remote"] = inspect_remote(item)
    sizes = [item["remote"]["bytes"] for item in spec["files"]]
    transfer = sum(sizes) if all(x is not None for x in sizes) else None
    spec.update(
        destination=str(Path(destination).expanduser().absolute()),
        transfer_bytes=transfer,
        free_disk_bytes=free_bytes(Path(destination).expanduser()),
        planned_at=datetime.now(timezone.utc).isoformat(),
    )
    return spec


def _safe_name(name):
    """Reject empty names, directory traversal and embedded path separators."""
    if not name or Path(name).name != name or name in {".", ".."} or "\\" in name:
        raise ValueError(f"Unsafe resource filename: {name!r}")


def _hashes(path, expected=None):
    """Stream SHA-256 and an optional expected MD5/SHA-256; reject mismatches."""
    algorithms = {"sha256": hashlib.sha256()}
    if expected:
        algorithm, digest = expected.split(":", 1)
        lengths = {"sha256": 64, "md5": 32}
        if algorithm not in lengths or not re.fullmatch(
            r"[0-9a-fA-F]{%d}" % lengths[algorithm], digest
        ):
            raise ValueError("Invalid expected resource checksum")
        algorithms.setdefault(algorithm, hashlib.new(algorithm))
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(BLOCK), b""):
            for hasher in algorithms.values():
                hasher.update(block)
    values = {key: hasher.hexdigest() for key, hasher in algorithms.items()}
    if expected and values[algorithm] != digest.lower():
        raise ValueError(f"Checksum mismatch for {Path(path).name}; partial file retained")
    return values


@contextmanager
def _lock(destination):
    """Hold a directory lock for a download and release it on context exit."""
    lock = Path(str(destination) + ".lock")
    try:
        lock.mkdir()
    except FileExistsError as error:
        raise FileExistsError(
            f"Download lock exists: {lock}. If the previous process was killed, "
            "verify it has stopped before removing this empty lock directory."
        ) from error
    try:
        yield
    finally:
        lock.rmdir()


def download_file(item, destination, *, resume=False, progress=None):
    """Stream one planned file; publish only after size/checksum validation.

    Resume requires the same URL, expected digest and a strong server ETag. Range
    responses must start at the saved offset. Interrupted bytes remain .partial;
    an existing final file is never overwritten. SHA256 is always recorded, but
    it authenticates against upstream only when an upstream digest was supplied.
    """
    destination = Path(destination)
    _safe_name(destination.name)
    if urlsplit(item["url"]).scheme not in {"http", "https"}:
        raise ValueError("Resource downloads require an HTTP(S) URL")
    destination.parent.mkdir(parents=True, exist_ok=True)
    with _lock(destination):
        if destination.exists() or destination.is_symlink():
            raise FileExistsError(f"Resource already exists: {destination}")
        partial = Path(str(destination) + ".partial")
        state_path = Path(str(partial) + ".json")
        remote = item.get("remote") or inspect_remote(item)
        identity = {"url": item["url"], "checksum": item.get("checksum"), **remote}
        offset = 0
        if partial.exists() or partial.is_symlink() or state_path.exists():
            if not resume or partial.is_symlink() or state_path.is_symlink():
                raise FileExistsError("Partial download exists; use --resume or a new destination")
            saved = json.loads(state_path.read_text())
            etag = saved.get("etag")
            if saved != identity or not etag or etag.startswith("W/"):
                raise ValueError(
                    "Cannot safely resume: resource identity or strong ETag unavailable"
                )
            offset = partial.stat().st_size if partial.exists() else 0
        else:
            with state_path.open("x") as stream:
                json.dump(identity, stream, indent=2)
        total = remote.get("bytes")
        if total is not None and offset > total:
            raise ValueError("Partial download is larger than the declared resource")
        if total is not None and total - offset > free_bytes(destination):
            raise OSError("Insufficient free disk space for the remaining download")
        headers = dict(HEADERS)
        if offset:
            headers.update(Range=f"bytes={offset}-", **{"If-Range": remote["etag"]})
        if total is None or offset != total:
            with urlopen(Request(item["url"], headers=headers), timeout=60) as response:
                if response.headers.get_content_type() == "text/html":
                    raise ValueError("Server returned HTML instead of a resource")
                if response.headers.get("Content-Encoding", "identity") != "identity":
                    raise ValueError("Unexpected HTTP content encoding")
                if remote.get("etag") and response.headers.get("ETag") != remote["etag"]:
                    raise ValueError("Remote resource changed since planning")
                if offset:
                    match = re.fullmatch(
                        r"bytes (\d+)-(\d+)/(\d+)", response.headers.get("Content-Range", "")
                    )
                    if response.status != 206 or not match or int(match[1]) != offset:
                        raise ValueError("Server did not honour the requested resume range")
                    if total is not None and int(match[3]) != total:
                        raise ValueError("Resume response changed the total size")
                    total = int(match[3])
                elif response.status != 200:
                    raise ValueError(f"Unexpected download response: {response.status}")
                length = response.headers.get("Content-Length")
                expected_response = int(length) if length is not None else None
                if total is None and expected_response is not None:
                    total = offset + expected_response
                received, last_update = 0, time.monotonic()
                with partial.open("ab" if partial.exists() and resume else "xb") as stream:
                    while block := response.read(BLOCK):
                        if total is not None and offset + received + len(block) > total:
                            raise ValueError("Download exceeds its declared size")
                        if free_bytes(destination) < len(block) + BLOCK:
                            raise OSError("Disk filled during download; partial file retained")
                        stream.write(block)
                        received += len(block)
                        if progress and time.monotonic() - last_update >= 5:
                            progress(offset + received, total)
                            last_update = time.monotonic()
                if expected_response is not None and received != expected_response:
                    raise ValueError("Truncated HTTP response; partial file retained")
        if total is not None and partial.stat().st_size != total:
            raise ValueError("Downloaded size does not match the declared size")
        hashes = _hashes(partial, item.get("checksum"))
        result = {
            "name": destination.name,
            "bytes": partial.stat().st_size,
            "hashes": hashes,
            "upstream_checksum_verified": bool(item.get("checksum")),
            "source": identity,
        }
        # A hard link publishes without replacing a concurrently created destination.
        os.link(partial, destination)
        partial.unlink()
        state_path.unlink()
        if progress:
            progress(result["bytes"], total)
        return result


def download_resource(plan, *, resume=False, progress=None):
    """Download a resolved resource into its named local folder, with per-file receipts."""
    out = Path(plan["destination"]) / plan["id"]
    _safe_name(plan["id"])
    # Check the entire known remaining transfer before spending bandwidth on the
    # first file of a multi-file resource. Unknown sizes remain guarded per block.
    remaining = 0
    for item in plan["files"]:
        _safe_name(item["name"])
        size = item.get("remote", {}).get("bytes", item.get("bytes"))
        final = out / item["name"]
        partial = Path(str(final) + ".partial")
        if resume and final.is_file() and not final.is_symlink():
            continue  # Verified against its receipt below before being accepted.
        offset = partial.stat().st_size if resume and partial.is_file() else 0
        if size is not None:
            remaining += max(0, size - offset)
    if remaining > free_bytes(out):
        raise OSError("Insufficient free disk for the complete remaining resource download")
    manifest = out / "DOWNLOAD_PLAN.json"
    if out.exists() or out.is_symlink():
        if out.is_symlink() or not resume:
            raise FileExistsError(f"Destination exists: {out}; use --resume or a new location")
        saved = json.loads(manifest.read_text())
        # Resume uses the original plan, not a changed release behind a moving URL.
        if saved["id"] != plan["id"] or saved["files"] != plan["files"]:
            raise ValueError("The resource changed since this download began")
    else:
        out.mkdir(parents=True)
        with manifest.open("x") as stream:
            json.dump(plan, stream, indent=2)
    results = []
    for item in plan["files"]:
        _safe_name(item["name"])
        path = out / item["name"]
        receipt = out / (item["name"] + ".download.json")
        if path.exists() and resume:
            if path.is_symlink() or receipt.is_symlink() or not receipt.is_file():
                raise ValueError(f"Existing resource lacks its completion receipt: {path}")
            result = json.loads(receipt.read_text())
            if _hashes(path, item.get("checksum")) != result["hashes"]:
                raise ValueError(f"Existing resource changed: {path}")
        else:
            callback = (
                (lambda count, size: progress(item["name"], count, size)) if progress else None
            )
            result = download_file(item, path, resume=resume, progress=callback)
            with receipt.open("x") as stream:
                json.dump(result, stream, indent=2)
        results.append(result)
    result = {
        "status": "PASS",
        "resource": plan["id"],
        "version": plan["version"],
        "files": results,
        "directory": str(out),
    }
    complete = out / "DOWNLOAD_COMPLETE.json"
    if not complete.exists():
        with complete.open("x") as stream:
            json.dump(result, stream, indent=2)
    return result


def unpack_resource(archives, destination, *, max_bytes):
    """Expand gzip/xz/tar archives into a fresh folder with a hard output-byte budget.

    Only regular files and directories are accepted. Links, traversal paths,
    duplicates and special files fail. Original archives remain untouched. The
    budget counts expanded file bytes, not RAM; extraction uses 1 MiB buffers.
    """
    if type(max_bytes) is not int or max_bytes < 1:
        raise ValueError("Choose a positive expanded byte budget")
    out = Path(destination)
    if out.exists() or out.is_symlink():
        raise FileExistsError(f"Unpack destination exists: {out}")
    sources = [Path(p).resolve(strict=True) for p in archives]
    if not sources:
        raise ValueError("Choose at least one archive")
    if max_bytes > free_bytes(out):
        raise OSError("Requested expansion budget exceeds available disk space")
    before = {str(p): _hashes(p)["sha256"] for p in sources}
    out.mkdir(parents=True)
    count, outputs = 0, []

    def copy_stream(stream, relative):
        """Extract one safe archive member within the byte budget and record its hash."""
        nonlocal count
        pure = PurePosixPath(relative)
        if pure.is_absolute() or ".." in pure.parts or "\\" in relative or ":" in relative:
            raise ValueError(f"Unsafe archive path: {relative}")
        target = out.joinpath(*pure.parts)
        target.parent.mkdir(parents=True, exist_ok=True)
        hasher, size = hashlib.sha256(), 0
        with target.open("xb") as writer:
            while block := stream.read(BLOCK):
                count += len(block)
                if count > max_bytes or free_bytes(out) < len(block) + BLOCK:
                    raise OSError("Expanded data exceeded the selected disk budget")
                writer.write(block)
                hasher.update(block)
                size += len(block)
        outputs.append(
            {"file": str(target.relative_to(out)), "bytes": size, "sha256": hasher.hexdigest()}
        )

    for path in sources:
        if path.name.endswith((".tar.gz", ".tgz", ".tar.xz", ".tar")):
            with tarfile.open(path, mode="r|*") as archive:
                for member in archive:
                    pure = PurePosixPath(member.name)
                    if pure.is_absolute() or ".." in pure.parts or "\\" in member.name:
                        raise ValueError(f"Unsafe archive path: {member.name}")
                    if member.isdir():
                        continue
                    if not member.isfile():
                        raise ValueError(f"Archive contains a link or special file: {member.name}")
                    if count + member.size > max_bytes:
                        raise OSError("Archive exceeds the selected expanded disk budget")
                    with archive.extractfile(member) as stream:
                        copy_stream(stream, member.name)
        elif path.suffix in {".gz", ".xz"}:
            opener = gzip.open if path.suffix == ".gz" else lzma.open
            with opener(path, "rb") as stream:
                copy_stream(stream, path.stem)
        else:
            raise ValueError(f"Unsupported archive: {path.name}")
    if any(_hashes(p)["sha256"] != digest for p, digest in before.items()):
        raise ValueError("Source archive changed during extraction")
    result = {
        "status": "PASS",
        "expanded_bytes": count,
        "max_bytes": max_bytes,
        "sources": before,
        "files": outputs,
    }
    with (out / "UNPACK_COMPLETE.json").open("x") as stream:
        json.dump(result, stream, indent=2)
    return result
