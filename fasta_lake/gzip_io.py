"""Deterministic gzip writing.

``gzip.open`` stamps the current wall-clock time into the 4-byte MTIME field
of the gzip header and the file's basename into FNAME, so two runs that
produce identical content produce different bytes. Every receipt in this
package hashes its outputs, so that timestamp made identical reruns look
different. Writers here set MTIME to 0 and omit FNAME.
"""

from __future__ import annotations

import gzip
import io
from contextlib import contextmanager
from pathlib import Path

#: pandas ``compression=`` argument with the same guarantee.
PANDAS_GZIP = {"method": "gzip", "mtime": 0}


@contextmanager
def open_gzip_text(path: str | Path, newline: str = "", encoding: str = "utf-8"):
    """Open ``path`` for writing gzip-compressed text with a zero MTIME header."""
    with (
        open(path, "wb") as raw,
        gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed,
        io.TextIOWrapper(compressed, encoding=encoding, newline=newline) as text,
    ):
        yield text
