"""
Regression test for the cross-sample header-identity defect (v4, 2026-06).

Background
----------
`lake_builder` deduplicates by SHA-256 of the sequence, but historically wrote
each sequence's *original* header as its FASTA identifier. For reference
catalogues (UniProt / UHGP / GMGC / GMSC) those accessions are globally unique,
so that is fine. For metagenomic / assembly sources (MEGAHIT + Prodigal) the
native `k141_<contig>_<orf>` identifiers are unique only WITHIN one sample's
assembly and recur as *different* sequences across samples. Pooling protein
groups across samples by such an identifier silently merges distinct proteins,
inflating every cross-sample quantity.

The fix is `--reassign-header-tags <TAG[,TAG...]>`: every sequence from a tagged
source is relabelled `<TAG>_<sha256hex>`, so identical sequences share one id
across all samples/studies and distinct sequences never collide.

This test pins both halves of the contract:
  * WITHOUT the flag, colliding assembly headers are reproduced verbatim
    (documents the failure mode, so a future "always reassign" change is a
    deliberate decision, not an accident).
  * WITH the flag, every FASTA identifier is unique, maps 1:1 to a distinct
    sequence, is of the form `<TAG>_<64 hex>`, and an identical sequence shared
    across two "samples" collapses to exactly one identifier.
"""

import os
import subprocess
from pathlib import Path

import pytest

from fasta_lake.hashing import sequence_hash

REPO = Path(__file__).resolve().parents[1]
LAKE_BUILDER = REPO / "rust" / "fasta_extractor_v2" / "target" / "release" / "lake_builder"
LAKE_BUILDER = (
    Path(os.environ.get("FASTALAKE_BIN_DIR", str(LAKE_BUILDER.parent))) / LAKE_BUILDER.name
)

# Two "sample assemblies" with COLLIDING k141 names but different sequences,
# plus one sequence shared across both samples under different names.
SHARED = "MSHAREDSEQUENCEPEPTIDEKAAAACDEFGHIK"
SAMPLE_A = (
    ">k141_1_1 contigA\nMKTAYIAKQRQISFVKSHFSRQLEERAAAACDK\n>k141_2_1 contigB\n" + SHARED + "\n"
)
SAMPLE_B = (
    ">k141_1_1 contigA\nMDIFFERENTSEQUENCEQRQISWYTCDEFGHIK\n>k141_9_9 contigZ\n" + SHARED + "\n"
)


def _read_headers(fasta: Path):
    return [ln[1:].split()[0] for ln in fasta.read_text().splitlines() if ln.startswith(">")]


def _read_records(fasta: Path):
    recs, hid, seq = {}, None, []
    for ln in fasta.read_text().splitlines():
        if ln.startswith(">"):
            if hid is not None:
                recs[hid] = "".join(seq)
            hid, seq = ln[1:].split()[0], []
        else:
            seq.append(ln.strip())
    if hid is not None:
        recs[hid] = "".join(seq)
    return recs


def _build(tmp_path, reassign):
    src = tmp_path / "asm"
    src.mkdir()
    (src / "sampleA.faa").write_text(SAMPLE_A)
    (src / "sampleB.faa").write_text(SAMPLE_B)
    out = tmp_path / ("lake_reassign.fasta" if reassign else "lake_plain.fasta")
    cmd = [str(LAKE_BUILDER), "-s", f"CAPSCAN:1:{src}:*.faa", "-o", str(out)]
    if reassign:
        cmd += ["--reassign-header-tags", "CAPSCAN"]
    subprocess.run(cmd, check=True, capture_output=True, text=True)
    return out


@pytest.mark.skipif(not LAKE_BUILDER.exists(), reason="lake_builder release binary not built")
def test_reassign_gives_unique_sequence_keyed_ids(tmp_path):
    out = _build(tmp_path, reassign=True)
    recs = _read_records(out)
    headers = list(recs)
    # three distinct sequences across the two samples (A-unique, B-unique, shared)
    distinct_seqs = {sequence_hash(s) for s in recs.values()}
    assert len(distinct_seqs) == 3
    # every identifier is unique...
    assert len(headers) == len(set(headers)) == 3
    # ...of the form CAPSCAN_<64 hex> and equal to the sequence hash
    for hid, seq in recs.items():
        assert hid == f"CAPSCAN_{sequence_hash(seq)}"
    # the shared sequence collapsed to exactly one identifier
    shared_ids = [h for h, s in recs.items() if sequence_hash(s) == sequence_hash(SHARED)]
    assert len(shared_ids) == 1


@pytest.mark.skipif(not LAKE_BUILDER.exists(), reason="lake_builder release binary not built")
def test_without_reassign_headers_collide(tmp_path):
    # Documents the defect: the two different k141_1_1 sequences are both kept
    # (dedup is by sequence) but share the identifier `k141_1_1`, so pooling by
    # identifier downstream conflates them. This asserts the unsafe behaviour so
    # any change to it is intentional.
    out = _build(tmp_path, reassign=False)
    headers = _read_headers(out)
    assert headers.count("k141_1_1") == 2  # same id, two distinct sequences
    assert len(headers) != len(set(headers))  # collision present
