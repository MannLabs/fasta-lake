"""
End-to-end entrapment workflow: generate -> inject -> measure FDP.

This module is the file-level, streaming counterpart to :mod:`fasta_lake.entrapment`
(which holds the in-memory library primitives). It exists so that a user can run a
complete entrapment experiment from the command line::

    fasta-lake entrapment classes
    fasta-lake entrapment generate --class shuffled -i lake.fasta -o entrap.fasta
    fasta-lake entrapment inject -d lake.fasta -e entrap.fasta -o augmented.fasta
    # ... run the pipeline against augmented.fasta, then search ...
    fasta-lake entrapment fdp --results-dir results/ --target-fasta lake.fasta

Injection stage is the whole ballgame
-------------------------------------
Entrapment can enter at **four** distinct points, and they answer different
questions. The stage is explicit everywhere: a CLI flag, a field in the manifest,
and a column in the per-sample output. See :mod:`fasta_lake.entrapment_stages` for
the full definitions and the run script implementing each.

``prelake``    input lake, before the cohort evidence filter -- construction error
``prebuild``   clustered lake + entrapment, re-clustered -- construction error
``preselect``  clustered lake + entrapment, no re-clustering (DEFAULT) -- construction
``appended``   finished per-sample database, after parsimony -- SEARCH error

``prebuild`` and ``preselect`` are **not** the same stage under two names: they
differ by whether entrapment had to survive clustering. All are legitimate, none
are the same quantity; this module refuses to pool across stages (or across
entrapment classes) without an explicit override.

Prefix conventions are not uniform
----------------------------------
This repository tags entrapment at least three incompatible ways --
``ENTRAP_PFUR_``/``ENTRAP_SHUF_``/… , ``ENTRAPHUMAN_`` (SIHUMIx) and ``ENTRAPNN_``
(near-neighbour). The latter two do not match ``ENTRAP_``, so a matcher hard-coded
to that convention scores them as target and drives the FDP toward zero *silently*.
The matcher here defaults to the common denominator ``ENTRAP``, is configurable,
records which pattern it used, and refuses to compute an FDP when a declared
entrapment partition is invisible to the configured pattern
(:class:`PrefixMismatchError`).

The conserved-peptide problem
-----------------------------
A peptide attributed to an entrapment protein may also genuinely occur in the
target sequence space -- universally conserved machinery such as ribosomal
proteins, elongation factors, chaperones and glycolytic enzymes is shared across
all of life. Matching one of those is not a false discovery, and counting it as one
inflates the FDP of any class with real homology pressure. Every FDP this module
reports therefore comes in a pair: **raw** and **conserved-corrected**, always
together, never one alone. See :func:`conserved_peptide_scan`.

Estimators
----------
Several FDP estimators are in circulation and their names collide across the
literature and across this codebase. Each is spelled out here, and each returns
``None`` -- never ``0.0`` -- when its denominator is undefined.

``fdp_uncorrected``        ``E / (E + T)``
``fdp_wen_equal_odds``     ``2E / (E + T)``              -- valid only when r ~ 1
``fdp_elias_gygi``         ``E(N_T/N_E) / (E(N_T/N_E) + T)``
``fdp_lower_bound``        ``E / T``                     -- proves failure, never success
``fdp_shared_sensitivity`` ``(E + S) / (E + S + T)``     -- S maps to both partitions
``fdp_combined_noble``     ``E(1 + 1/r) / (T + E)``, r = N_E/N_T

References
----------
- Wen B, Keich U, Noble WS et al. (2025). Nat Methods 22, 1454-1463.
- Elias JE & Gygi SP (2007). Nat Methods 4, 207-214.
"""

from __future__ import annotations

import csv
import hashlib
import json
import logging
import math
import random
import re
import statistics
import time
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterator

from fasta_lake.entrapment_classes import ENTRAPMENT_CLASSES, get_class
from fasta_lake.entrapment_stages import (
    ENTRAPMENT_STAGES,
    VALID_ENTRAPMENT_STAGES,
    get_stage,
    normalize_stage,
)
from fasta_lake.helpers import require_columns

logger = logging.getLogger(__name__)

#: Default pattern identifying an entrapment record by its accession token.
#:
#: NOT ``ENTRAP_``. The repository uses at least three incompatible conventions --
#: ``ENTRAP_PFUR_``/``ENTRAP_SHUF_``/``ENTRAP_PLANT_``/``ENTRAP_ARCH_`` (the
#: manuscript runs), ``ENTRAPHUMAN_`` (four SIHUMIx projects) and ``ENTRAPNN_``
#: (the near-neighbour run). The latter two do **not** match ``ENTRAP_``, so a
#: matcher assuming that convention scores them as TARGET and drives the FDP
#: toward zero silently, with the run completing and the output looking like a
#: pass. ``ENTRAP`` is the common denominator of all three. Override with
#: ``--entrap-prefix`` when a lake uses something else, and check with
#: :func:`detect_prefixes` before trusting a number.
ENTRAP_PREFIX = "ENTRAP"

#: Conventions observed in this repository, longest first so prefix matching is
#: unambiguous. This is the matcher ENTRAPMENT_PREFIX_CENSUS.md §8 recommends over
#: both a narrow literal (``ENTRAP_``, which silently misses ``ENTRAPHUMAN_`` and
#: ``ENTRAPNN_``) and a loose substring (which cannot tell a tag from a target
#: accession that merely contains the letters). ``ENTRAP_SHUFFLED_`` is the legacy
#: in-memory library prefix; ``ENTRAPHUMAN-`` (hyphen) is the sihumi_gt_ceiling
#: landmine the census found -- a real fourth family that no earlier matcher caught.
KNOWN_ENTRAP_PREFIXES = tuple(
    sorted(
        (
            "ENTRAP_SHUF_",
            "ENTRAP_SHUFFLED_",
            "ENTRAP_PFUR_",
            "ENTRAP_PLANT_",
            "ENTRAP_ARCH_",
            "ENTRAP_NOBLE_",
            "ENTRAP_NONSELF_",
            "ENTRAP_USER_",
            "ENTRAPHUMAN_",
            "ENTRAPHUMAN-",
            "ENTRAPNN_",
        ),
        key=len,
        reverse=True,
    )
)

# Injection stages live in entrapment_stages.py -- see that module for what each
# one measures and which run script implements it.
VALID_STAGES = VALID_ENTRAPMENT_STAGES

# Retained for callers that referenced these names; both are canonical stages.
STAGE_PRESELECT = "preselect"
STAGE_PREBUILD = "prebuild"
STAGE_PRELAKE = "prelake"
STAGE_APPENDED = "appended"

STAGE_MEANING = {
    name: (
        f"Entrapment injected at {s.injection_point}. It must survive: "
        f"{', '.join(s.faces)}. {s.description} Measures {s.measures.upper()}."
    )
    for name, s in ENTRAPMENT_STAGES.items()
}

MANIFEST_SCHEMA = "fastalake.entrapment.manifest/1"
RECORD_SCHEMA = "fastalake.entrapment.fdp/1"
MANIFEST_NAME = "entrapment_manifest.json"

# Amino acids held in place when ``preserve_cleavage_sites`` is on, so that the
# tryptic peptide-length distribution of the entrapment matches the source.
_CLEAVAGE_RESIDUES = frozenset("KR")

# Strips SAGE/DIA-NN modification annotations to a bare residue string.
_MOD_STRIP = re.compile(r"[^A-Z]")


def stage_measures(stage: str) -> str:
    """One-line statement of what an FDP at this injection stage estimates."""
    return get_stage(stage).measures


# ---------------------------------------------------------------------------
# Streaming FASTA primitives
# ---------------------------------------------------------------------------


def iter_fasta(path: str | Path) -> Iterator[tuple[str, str]]:
    """Stream ``(header, sequence)`` pairs from a FASTA file.

    The header is returned verbatim without the leading ``>``. Nothing is held in
    memory beyond one record, so this is safe on multi-hundred-million-sequence
    lakes.
    """
    header: str | None = None
    chunks: list[str] = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:].rstrip("\n\r")
                chunks = []
            else:
                chunks.append(line.strip())
    if header is not None:
        yield header, "".join(chunks)


def count_fasta_headers(path: str | Path, prefix: str | None = None) -> int:
    """Count proteins in a FASTA by counting **header lines only**.

    This is deliberately not ``grep -cv '^>ENTRAP_'``. That idiom counts every line
    that is not an entrapment header -- including all sequence lines -- and inflates
    protein counts by roughly the mean number of sequence lines per record
    (measured at 8.45x on a real per-sample database in this project, 122,189 lines
    against 14,465 proteins). Every count in this module goes through this function.
    """
    n = 0
    with open(path) as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            if prefix is None or line[1:].startswith(prefix):
                n += 1
    return n


def partition_counts(path: str | Path, prefix: str = ENTRAP_PREFIX) -> tuple[int, int]:
    """Return ``(n_target, n_entrapment)`` protein counts for a FASTA.

    Both counts come from header lines only; ``n_target`` never counts sequence
    lines.
    """
    n_total = 0
    n_entrap = 0
    with open(path) as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            n_total += 1
            if line[1:].startswith(prefix):
                n_entrap += 1
    return n_total - n_entrap, n_entrap


class PrefixMismatchError(ValueError):
    """The configured entrapment prefix matched nothing that was declared present.

    This is a distinct failure cause from "no entrapment was detected". A prefix
    mismatch silently reclassifies every entrapment record as target and drives
    the FDP toward zero while the run completes and the output looks like a pass.
    It must stop the pipeline, not produce a number.
    """


def detect_prefixes(
    path: str | Path,
    candidate: str = "ENTRAP",
    max_report: int = 20,
) -> dict:
    """Report the entrapment-like prefixes actually present in a FASTA.

    Scans header accession tokens for any that contain ``candidate`` (case
    insensitive) and reports the distinct prefixes found with their counts, so a
    convention mismatch is visible before it becomes a published number rather
    than after.

    Returns ``n_headers``, ``n_matching_default``, ``observed`` (prefix -> count)
    and ``known`` (which observed prefixes are conventions this repo already uses).
    """
    n_headers = 0
    n_default = 0
    observed: Counter = Counter()
    up = candidate.upper()

    with open(path) as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            n_headers += 1
            token = line[1:].split()[0] if len(line) > 1 and line[1:].split() else ""
            if not token:
                continue
            if token.startswith(candidate):
                n_default += 1
            if up in token.upper():
                # Prefer an exact known family (so ENTRAPHUMAN- with its hyphen
                # delimiter is reported verbatim, not collapsed). Fall back to the
                # underscore-cut heuristic only for genuinely novel forms.
                known = next((k for k in KNOWN_ENTRAP_PREFIXES if token.startswith(k)), None)
                if known is not None:
                    observed[known] += 1
                else:
                    idx = token.upper().index(up)
                    rest = token[idx:]
                    cut = rest.find("_", len(candidate))
                    observed[rest[: cut + 1] if cut != -1 else rest[: len(candidate)]] += 1

    return {
        "path": str(path),
        "n_headers": n_headers,
        "candidate": candidate,
        "n_matching_candidate": n_default,
        "observed": dict(observed.most_common(max_report)),
        "known": sorted(p for p in observed if p in KNOWN_ENTRAP_PREFIXES),
        "unknown": sorted(p for p in observed if p not in KNOWN_ENTRAP_PREFIXES),
    }


def verify_prefix_matches(
    path: str | Path,
    prefix: str,
    declared_n_entrapment: int,
    context: str = "",
) -> int:
    """Assert the configured prefix finds the entrapment partition that is declared.

    Returns the observed entrapment header count. Raises
    :class:`PrefixMismatchError` when ``declared_n_entrapment > 0`` but scanning
    with ``prefix`` finds none -- that combination can only mean the prefix is
    wrong for this file, and continuing would report a spuriously low FDP.

    Note what this does **not** flag: zero entrapment *hits* against a populated
    entrapment partition is a genuine, desirable 0% result and is left alone. The
    check is on the database partition, not on the discoveries.
    """
    _, observed = partition_counts(path, prefix)
    if declared_n_entrapment > 0 and observed == 0:
        diag = detect_prefixes(path)
        found = ", ".join(f"{k} ({v})" for k, v in diag["observed"].items()) or "none"
        raise PrefixMismatchError(
            f"entrapment prefix mismatch{' in ' + context if context else ''}: "
            f"{declared_n_entrapment} entrapment records were declared for "
            f"{path}, but scanning with prefix {prefix!r} found 0. "
            f"Entrapment-like prefixes actually present: {found}. "
            "Every entrapment record would have been counted as TARGET, driving "
            "the FDP toward zero. Set --entrap-prefix to the correct convention. "
            "Refusing to compute an FDP."
        )
    return observed


def _sha256_file(path: str | Path, chunk: int = 1 << 20) -> str:
    """SHA-256 of a file, streamed."""
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while True:
            block = fh.read(chunk)
            if not block:
                break
            h.update(block)
    return h.hexdigest()


# ---------------------------------------------------------------------------
# Generation
# ---------------------------------------------------------------------------


def shuffle_sequence(
    seq: str,
    rng: random.Random,
    preserve_cleavage_sites: bool = True,
) -> str:
    """Shuffle a sequence, optionally holding K/R residues in place.

    Holding lysine and arginine fixed keeps every tryptic cleavage site where it
    was, so each entrapment protein yields peptides whose lengths match its source
    protein's **one for one**, while sharing none of their order.

    A naive full shuffle preserves the K/R *count* per protein, so the aggregate
    peptide-length distribution comes out about the same (measured on the 512-protein
    P. furiosus panel: mean tryptic length 6.80 source vs 6.79 naive). What it does
    not preserve is the per-protein match: the same measurement found the naive
    shuffle reproduced a protein's exact tryptic length multiset for 1 of 512
    proteins (0.2%), against 512 of 512 (100%) when K/R are held. Per-protein
    matching is what makes each entrapment sequence a paired null for its source.
    """
    if not preserve_cleavage_sites:
        chars = list(seq)
        rng.shuffle(chars)
        return "".join(chars)

    movable = [i for i, c in enumerate(seq) if c not in _CLEAVAGE_RESIDUES]
    chars = [seq[i] for i in movable]
    rng.shuffle(chars)
    out = list(seq)
    for i, c in zip(movable, chars):
        out[i] = c
    return "".join(out)


@dataclass
class GenerateStats:
    """Summary of an entrapment-generation run."""

    entrap_class: str
    method: str
    seed: int
    prefix: str
    preserve_cleavage_sites: bool
    source: str
    n_source: int
    n_written: int
    n_duplicate_ids: int
    output: str


def generate_entrapment_fasta(
    output_fasta: str | Path,
    entrap_class: str = "shuffled",
    source_fasta: str | Path | None = None,
    seed: int = 42,
    n_proteins: int | None = None,
    prefix: str | None = None,
    preserve_cleavage_sites: bool = True,
) -> GenerateStats:
    """Write an entrapment FASTA for one entrapment class.

    Deterministic: the same class, seed and input file produce a byte-identical
    output file.

    Parameters
    ----------
    output_fasta
        Destination. Headers are written as ``>{prefix}{original_first_token}``.
    entrap_class
        A key of :data:`~fasta_lake.entrapment_classes.ENTRAPMENT_CLASSES`.
        ``'shuffled'`` (default) is the recommended headline class: it is the only
        class with no degenerate-coverage problem. Classes with
        ``generation='copy'`` tag a real organism's FASTA without altering the
        sequences.
    source_fasta
        For ``'shuffled'``, the target sequences to shuffle (required -- normally
        the lake being injected into). For copy classes, overrides the class's
        ``default_source``.
    seed
        Random seed. Used for shuffling and for subsampling when ``n_proteins`` is
        set. Seed 0 is rejected: it is a fixed point of the xorshift64 generator
        used by the Rust extractor, and rejecting it uniformly keeps one convention.
    n_proteins
        If set, write only this many proteins, chosen by seeded random subsample.
        If ``None`` (default) every source sequence is used.
    prefix
        Override the class's header prefix. Defaults to the class prefix.
    preserve_cleavage_sites
        Hold K/R in place while shuffling. Default ``True``. Ignored by copy classes.
    """
    cls = get_class(entrap_class)
    if seed == 0:
        raise ValueError(
            "entrapment seed 0 is rejected (fixed point of the xorshift64 generator "
            "used by the Rust extractor). Pick a nonzero seed, e.g. 42."
        )

    if cls.generation == "shuffle" and source_fasta is None:
        raise ValueError(
            f"class {cls.name!r} derives its null from the target sequences; pass "
            "--source pointing at the database you are going to inject into."
        )
    src = source_fasta or cls.default_source
    if src is None:
        raise ValueError(
            f"entrapment class {cls.name!r} needs a source FASTA (--source), and has no default."
        )
    src = Path(src)
    if not src.exists():
        raise ValueError(f"source FASTA for class {cls.name!r} does not exist: {src}")

    prefix = prefix or cls.prefix
    output_fasta = Path(output_fasta)
    output_fasta.parent.mkdir(parents=True, exist_ok=True)

    n_source = count_fasta_headers(src)
    if n_source == 0:
        raise ValueError(f"source FASTA has no sequences: {src}")

    # Seeded subsample of record indices: deterministic, and independent of how
    # the file happens to be chunked on disk.
    keep: set[int] | None = None
    if n_proteins is not None and n_proteins < n_source:
        keep = set(random.Random(seed).sample(range(n_source), n_proteins))

    rng = random.Random(seed)
    seen: set[str] = set()
    n_written = 0
    n_dup = 0

    with open(output_fasta, "w") as out:
        for idx, (header, seq) in enumerate(iter_fasta(src)):
            if keep is not None and idx not in keep:
                continue
            if not seq:
                continue
            tokens = header.split()
            pid = tokens[0] if tokens else f"seq{idx}"
            # A source that is already tagged (several shipped panels are) must
            # not end up double-prefixed.
            if pid.startswith(ENTRAP_PREFIX):
                pid = pid[len(ENTRAP_PREFIX) :]
                if pid.startswith(prefix[len(ENTRAP_PREFIX) :]):
                    pid = pid[len(prefix) - len(ENTRAP_PREFIX) :]
            if pid in seen:
                n_dup += 1
            seen.add(pid)
            payload = (
                shuffle_sequence(seq, rng, preserve_cleavage_sites)
                if cls.generation == "shuffle"
                else seq
            )
            out.write(f">{prefix}{pid}\n{payload}\n")
            n_written += 1

    if n_dup:
        logger.warning(
            "%d duplicate protein IDs in %s -- entrapment headers are not unique",
            n_dup,
            src,
        )
    logger.info(
        "Wrote %d entrapment sequences (class=%s, seed=%d) to %s",
        n_written,
        cls.name,
        seed,
        output_fasta,
    )
    return GenerateStats(
        entrap_class=cls.name,
        method=cls.generation,
        seed=seed,
        prefix=prefix,
        preserve_cleavage_sites=preserve_cleavage_sites and cls.generation == "shuffle",
        source=str(src),
        n_source=n_source,
        n_written=n_written,
        n_duplicate_ids=n_dup,
        output=str(output_fasta),
    )


# ---------------------------------------------------------------------------
# Injection
# ---------------------------------------------------------------------------


def inject_entrapment(
    database_fasta: str | Path,
    entrapment_fasta: str | Path,
    output_fasta: str | Path,
    stage: str = STAGE_PRESELECT,
    entrap_class: str | None = None,
    seed: int | None = None,
    prefix: str | None = None,
    manifest_path: str | Path | None = None,
    checksums: bool = True,
) -> dict:
    """Concatenate entrapment into a database and record what was done.

    Writes ``output_fasta`` (target records first, then entrapment) and a JSON
    manifest recording the class, the stage, the seed, the entrapment:target ratio
    and the counts. The manifest is what makes a downstream FDP number
    interpretable -- without it, a result file cannot be told apart from one
    produced at the other injection stage or from another class.

    Returns the manifest as a dict.
    """
    stage = normalize_stage(stage)
    stage_meta = get_stage(stage)
    cls = get_class(entrap_class) if entrap_class else None
    if prefix is None:
        prefix = cls.prefix if cls is not None else ENTRAP_PREFIX

    database_fasta = Path(database_fasta)
    entrapment_fasta = Path(entrapment_fasta)
    output_fasta = Path(output_fasta)
    output_fasta.parent.mkdir(parents=True, exist_ok=True)

    # Scan the target database with the generic pattern, so entrapment left over
    # under ANY convention is caught, not just the one we are about to write.
    db_target, db_entrap_preexisting = partition_counts(database_fasta, ENTRAP_PREFIX)
    e_nonentrap, e_entrap = partition_counts(entrapment_fasta, prefix)
    if e_nonentrap:
        diag = detect_prefixes(entrapment_fasta)
        found = ", ".join(f"{k} ({v})" for k, v in diag["observed"].items()) or "none"
        raise ValueError(
            f"{entrapment_fasta} contains {e_nonentrap} headers without the {prefix!r} "
            "prefix. Every entrapment sequence must be tagged, or hits cannot be "
            f"attributed. Entrapment-like prefixes present: {found}. "
            "Regenerate with `fasta-lake entrapment generate`, or pass the correct "
            "--entrap-prefix."
        )
    if e_entrap == 0:
        raise ValueError(f"entrapment FASTA has no {prefix!r}-tagged sequences: {entrapment_fasta}")
    if db_entrap_preexisting:
        raise ValueError(
            f"{database_fasta} already contains {db_entrap_preexisting} "
            f"{ENTRAP_PREFIX!r} headers. Injecting again would double-count and mix "
            "classes. Use the un-augmented database."
        )

    with open(output_fasta, "w") as out:
        for src in (database_fasta, entrapment_fasta):
            last = "\n"
            with open(src) as fh:
                for line in fh:
                    out.write(line)
                    last = line
            if not last.endswith("\n"):
                out.write("\n")

    manifest = {
        "schema": MANIFEST_SCHEMA,
        "created": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "stage": stage,
        "stage_meaning": STAGE_MEANING[stage],
        "stage_metadata": stage_meta.as_metadata(),
        "measures": stage_measures(stage),
        "is_construction_test": stage_meta.is_construction_test,
        "class": cls.name if cls else None,
        "class_metadata": cls.as_metadata() if cls else None,
        "seed": seed,
        "prefix": prefix,
        "entrapment_match_pattern": prefix,
        "n_target": db_target,
        "n_entrapment": e_entrap,
        "ratio_entrapment_to_target": (round(e_entrap / db_target, 6) if db_target else None),
        "database_fasta": str(database_fasta),
        "entrapment_fasta": str(entrapment_fasta),
        "output_fasta": str(output_fasta),
    }
    if checksums:
        manifest["database_sha256"] = _sha256_file(database_fasta)
        manifest["entrapment_sha256"] = _sha256_file(entrapment_fasta)

    mpath = Path(manifest_path) if manifest_path else output_fasta.parent / MANIFEST_NAME
    mpath.parent.mkdir(parents=True, exist_ok=True)
    mpath.write_text(json.dumps(manifest, indent=2) + "\n")
    manifest["manifest_path"] = str(mpath)

    logger.info(
        "Injected %d entrapment (class=%s) into %d target at stage=%s -> %s",
        e_entrap,
        cls.name if cls else "?",
        db_target,
        stage,
        output_fasta,
    )
    return manifest


def load_manifests(paths: list[Path], allow_mixed: bool = False) -> tuple[dict, list[dict]]:
    """Load entrapment manifests, refusing to mix injection stages or classes.

    Returns ``({'stage': ..., 'class': ...}, manifests)``. Raises ``ValueError``
    when the manifests disagree, unless ``allow_mixed`` is set: pre-selection and
    appended FDPs measure different quantities, and so do different entrapment
    classes, so pooling them produces a number that means nothing.
    """
    manifests = []
    for p in paths:
        try:
            manifests.append(json.loads(Path(p).read_text()))
        except (OSError, json.JSONDecodeError) as exc:
            raise ValueError(f"cannot read entrapment manifest {p}: {exc}") from exc

    stages = {m.get("stage") for m in manifests if m.get("stage")}
    classes = {m.get("class") for m in manifests if m.get("class")}

    if not allow_mixed:
        if len(stages) > 1:
            raise ValueError(
                "Refusing to combine entrapment results from different injection "
                f"stages: {sorted(stages)}. Pre-selection FDP measures "
                "database-construction error; appended FDP measures search error. "
                "They are different quantities. Pass --allow-mixed to override, and "
                "say so wherever the number is reported."
            )
        if len(classes) > 1:
            raise ValueError(
                "Refusing to combine entrapment results from different classes: "
                f"{sorted(classes)}. Classes differ in difficulty and in what their "
                "FDP estimates (see `fasta-lake entrapment classes`). Pass "
                "--allow-mixed to override."
            )
    elif len(stages) > 1 or len(classes) > 1:
        logger.warning(
            "POOLING INCOMPARABLE RESULTS: stages=%s classes=%s. The resulting FDP "
            "mixes different quantities and must not be reported as a single number.",
            sorted(stages),
            sorted(classes),
        )

    return (
        {
            "stage": next(iter(stages)) if len(stages) == 1 else None,
            "class": next(iter(classes)) if len(classes) == 1 else None,
        },
        manifests,
    )


# ---------------------------------------------------------------------------
# PSM classification
# ---------------------------------------------------------------------------


def strip_modifications(peptide: str) -> str:
    """Reduce a modified peptide string to bare uppercase residues.

    SAGE writes ``M[+15.9949]``-style annotations in the ``peptide`` column; the
    residue count and any sequence matching must be done on the stripped form.
    """
    return _MOD_STRIP.sub("", peptide.upper())


@dataclass
class PsmCounts:
    """Entrapment / target / shared PSM counts for one search result file."""

    target: int = 0
    entrap: int = 0
    shared: int = 0
    decoy: int = 0
    n_rows: int = 0
    decoy_source: str = ""
    peptides_target: set[str] = field(default_factory=set)
    peptides_entrap: set[str] = field(default_factory=set)
    peptides_shared: set[str] = field(default_factory=set)
    #: entrapment peptide -> number of spectra, so a conserved peptide can be
    #: subtracted at the same granularity it was counted.
    entrap_psm_counts: Counter = field(default_factory=Counter)


def classify_psms(
    results_tsv: str | Path,
    prefix: str = ENTRAP_PREFIX,
    q_threshold: float = 0.01,
    min_length: int = 9,
    q_column: str = "peptide_q",
    decoy_prefix: str = "rev_",
) -> PsmCounts:
    """Classify accepted PSMs into target / entrapment / shared.

    Decoys are identified from the search engine's own ``label`` column (SAGE
    writes ``1`` for target, ``-1`` for decoy) whenever that column exists. Only if
    it is absent do we fall back to a header-prefix test, and then a warning is
    emitted -- a hard-coded prefix assumption (``DECOY_`` where the engine actually
    writes ``rev_``) has silently corrupted results in this project before.

    Decoys are counted separately and are **never** folded into the target count;
    doing so inflates the FDP denominator and biases the estimate down.

    A PSM whose peptide maps to both an entrapment and a target protein is counted
    as ``shared`` rather than being silently dropped, so a sensitivity bound that
    treats shared as entrapment can be reported alongside the headline.
    """
    counts = PsmCounts()
    path = Path(results_tsv)
    if not path.exists():
        return counts

    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fields = reader.fieldnames or []
        has_label = "label" in fields
        counts.decoy_source = "label column" if has_label else f"{decoy_prefix!r} header prefix"
        if not has_label:
            logger.warning(
                "%s has no 'label' column; falling back to the %r header prefix to "
                "identify decoys. Verify this matches your search engine's decoy tag.",
                path,
                decoy_prefix,
            )
        if q_column not in fields:
            raise ValueError(
                f"{path} has no {q_column!r} column (found: {', '.join(fields)}). "
                "Pass --q-column to select the right FDR column."
            )
        require_columns(fields, ("proteins", "peptide"), path)

        for row in reader:
            counts.n_rows += 1
            prots = [p for p in (row.get("proteins", "") or "").replace(",", ";").split(";") if p]

            if has_label:
                try:
                    is_decoy = int(row["label"]) != 1
                except (KeyError, TypeError, ValueError):
                    continue
            else:
                is_decoy = any(p.startswith(decoy_prefix) for p in prots)

            if is_decoy:
                counts.decoy += 1
                continue

            peptide = strip_modifications(row.get("peptide", "") or "")
            try:
                if float(row[q_column]) > q_threshold:
                    continue
                plen = int(row.get("peptide_len", 0) or 0) or len(peptide)
                if plen < min_length:
                    continue
            except (KeyError, TypeError, ValueError):
                continue

            has_e = any(p.startswith(prefix) for p in prots)
            has_t = any(not p.startswith(prefix) for p in prots)
            if has_e and has_t:
                counts.shared += 1
                counts.peptides_shared.add(peptide)
            elif has_e:
                counts.entrap += 1
                counts.peptides_entrap.add(peptide)
                counts.entrap_psm_counts[peptide] += 1
            else:
                counts.target += 1
                counts.peptides_target.add(peptide)

    return counts


# ---------------------------------------------------------------------------
# Conserved-peptide analysis
# ---------------------------------------------------------------------------


@dataclass
class ConservedScan:
    """Which entrapment-attributed peptides also occur in the target space.

    ``conserved`` maps a peptide to the target protein IDs it was found in, so a
    user can see *what* is being shared rather than being handed a bare fraction.

    ``control_conserved`` is the negative control: the same test run on the
    **reversed** entrapment peptides. A reversed peptide has the same composition
    and length but is not a real conserved sequence, so it should almost never be
    found. If the control rate approaches the real rate, the matcher is picking up
    chance matches and the correction is not trustworthy.
    """

    n_tested: int = 0
    n_conserved: int = 0
    conserved: dict[str, list[str]] = field(default_factory=dict)
    n_control_tested: int = 0
    n_control_conserved: int = 0
    control_conserved: dict[str, list[str]] = field(default_factory=dict)
    warning: str = ""
    performed: bool = False
    target_fasta: str = ""

    @property
    def rate(self) -> float | None:
        """Fraction of entrapment peptides also present in the target space."""
        return self.n_conserved / self.n_tested if self.n_tested else None

    @property
    def control_rate(self) -> float | None:
        """Fraction of reversed entrapment peptides found -- should be ~0."""
        return self.n_control_conserved / self.n_control_tested if self.n_control_tested else None


def conserved_peptide_scan(
    peptides: set[str] | list[str],
    target_fasta: str | Path,
    normalize_il: bool = True,
    max_proteins_per_peptide: int = 5,
    run_control: bool = True,
    control_warn_ratio: float = 0.25,
) -> ConservedScan:
    """Test entrapment-attributed peptides for occurrence in the target space.

    Builds one Aho-Corasick automaton over the (small) peptide set and streams the
    (large) target FASTA past it, so cost is one pass over the database regardless
    of how many peptides are tested. Matching is exact substring on the stripped
    residue string, optionally with I/L normalised -- leucine and isoleucine are
    isobaric and indistinguishable by mass spectrometry, so treating them as
    distinct would call a genuinely shared peptide class-specific.

    A peptide found in the target space is a **conserved** hit: it is a real
    identification that happens also to occur in the entrapment organism, not a
    false discovery. Subtracting these gives the conserved-corrected FDP.

    The reversed-peptide negative control runs in the same pass.
    """
    import ahocorasick

    from fasta_lake.helpers import clean_sequence
    from fasta_lake.helpers import normalize_il as _norm_il

    def norm(s: str) -> str:
        """Apply I/L equivalence only when requested for this diagnostic scan."""
        return _norm_il(s) if normalize_il else s

    scan = ConservedScan(target_fasta=str(target_fasta))
    peps = sorted({p for p in peptides if p})
    scan.n_tested = len(peps)
    if not peps:
        scan.performed = True
        return scan

    automaton = ahocorasick.Automaton()
    # key -> list of (kind, original_peptide); a reversed peptide can coincide
    # with a real one (palindromes, or two peptides that reverse into each other).
    payload: dict[str, list[tuple[str, str]]] = {}
    for p in peps:
        payload.setdefault(norm(p), []).append(("real", p))
    if run_control:
        for p in peps:
            payload.setdefault(norm(p[::-1]), []).append(("control", p))
        scan.n_control_tested = len(peps)

    for key, entries in payload.items():
        automaton.add_word(key, (key, entries))
    automaton.make_automaton()

    hits_real: dict[str, list[str]] = {}
    hits_control: dict[str, list[str]] = {}

    for header, seq in iter_fasta(target_fasta):
        tokens = header.split()
        pid = tokens[0] if tokens else ""
        if pid.startswith(ENTRAP_PREFIX):
            continue  # target space only -- never match entrapment against itself
        seq_norm = norm(clean_sequence(seq))
        for _end, (_key, entries) in automaton.iter(seq_norm):
            for kind, pep in entries:
                bucket = hits_real if kind == "real" else hits_control
                prots = bucket.setdefault(pep, [])
                if len(prots) < max_proteins_per_peptide and pid not in prots:
                    prots.append(pid)

    scan.conserved = hits_real
    scan.n_conserved = len(hits_real)
    scan.control_conserved = hits_control
    scan.n_control_conserved = len(hits_control)
    scan.performed = True

    if run_control and scan.rate:
        if scan.control_rate and scan.control_rate > control_warn_ratio * scan.rate:
            scan.warning = (
                f"Negative control is not near-zero: {scan.n_control_conserved} of "
                f"{scan.n_control_tested} REVERSED entrapment peptides "
                f"({scan.control_rate:.2%}) were also found in the target space, "
                f"against {scan.rate:.2%} for the real peptides. Reversed peptides "
                "are not conserved sequence, so this indicates the matcher is "
                "finding chance matches -- most likely the target space is very "
                "large relative to the peptide length floor. Treat the "
                "conserved-corrected FDP as unreliable and raise --min-length."
            )
            logger.warning(scan.warning)

    return scan


def write_conserved_report(scan: ConservedScan, path: str | Path) -> Path:
    """Write the conserved peptides and their target proteins to TSV."""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["peptide", "kind", "n_target_proteins_shown", "target_proteins"])
        for pep, prots in sorted(scan.conserved.items()):
            w.writerow([pep, "conserved", len(prots), ";".join(prots)])
        for pep, prots in sorted(scan.control_conserved.items()):
            w.writerow([pep, "reversed_control", len(prots), ";".join(prots)])
    return p


# ---------------------------------------------------------------------------
# FDP estimators -- all return None (never 0.0) when undefined
# ---------------------------------------------------------------------------


def fdp_estimates(
    entrap: int,
    target: int,
    shared: int = 0,
    db_target: int | None = None,
    db_entrap: int | None = None,
) -> dict:
    """Compute every FDP estimator, returning ``None`` where undefined.

    An FDP of ``0.0`` reported from an undefined denominator is not an "honest
    zero" -- it is the most favourable value the statistic can take, standing in
    for a measurement that was never possible. In this project 52.6% of samples in
    one arm had **zero** entrapment sequences in the searched database and were
    reported as 0.000%. Every estimator here returns ``None`` in that situation,
    and callers report coverage separately.

    Parameters
    ----------
    entrap, target, shared
        Accepted PSM (or peptide) counts by class.
    db_target, db_entrap
        Protein counts of the two partitions of the **searched** database, from
        header lines only. Needed for the size-corrected estimators; pass ``None``
        if unknown, and those estimators return ``None``.
    """
    denom = entrap + target
    ratio = (db_entrap / db_target) if (db_target and db_entrap is not None) else None

    # NAMES CORRECTED 2026-08-02 against the version of record
    # (Nat Methods 22, 1454-1463, PMC12240826) + its Supplementary. Three of the
    # four estimators below were mislabelled, and a mislabelled bound is worse
    # than a missing one because it still prints a plausible number.
    #
    #   eq 2, LOWER BOUND        E / (T + E)
    #     This is what `uncorrected` already computed, unnamed. It is a rigorous
    #     bound on the FDP ITSELF (Claim 1), unlike eq 1 and eq 4 which bound only
    #     the FDR, i.e. the expectation. So it can DEMONSTRATE FAILURE of control
    #     but can never demonstrate control. Wen frames using it as evidence of
    #     control as the incorrect application of the combined method.
    #
    #   eq 3, SAMPLE method      E(1/r) / T   -- at r=1 this is E/T
    #     NOT a bound. Wen: it "cannot be used to provide empirical evidence that
    #     a tool controls the FDR nor that the tool fails". `lower` used to be
    #     E/target, i.e. this quantity, returned under the name lower_bound.
    #     That is the defect being fixed.
    #
    #   Elias & Gygi CONCATENATED   2E / (T + E)
    #     The factor of 2 is explicit in Elias & Gygi, Methods Mol Biol 604,
    #     55-71 (2010), eq 3 (FP = f*d, f = 2). The 2007 Nat Methods paper is
    #     closed access everywhere and must not be quoted. Wen confirms eq 1 at
    #     r = 1 "reduces to Elias and Gygi's original estimation of the FDR in the
    #     concatenated target-decoy database", so the identity holds from both
    #     ends. The old `fdp_elias_gygi` was a size-scaled E/(E+T), which equals
    #     E/(E+T) at r=1 -- the lower bound again, not the concatenated form.
    lower_bound = entrap / denom if denom else None  # eq 2
    sample_method = entrap / target if target else None  # eq 3 at r=1; NOT a bound

    # DEFINEDNESS GATE -- restored 2026-08-02 after the rename dropped it.
    # Elias-Gygi is only meaningful against a POPULATED entrapment partition. With
    # db_entrap = 0 or None the quantity is UNDEFINED, not zero. The first version
    # of this rename computed 2E/(T+E) unconditionally and so returned 0.0 for an
    # empty partition -- caught by TestUndefinedFdp, which exists because 52.6% of
    # samples in one arm had zero entrapment sequences in the searched database.
    # Reporting those as 0.000% is, in this module's own words, "the single most
    # misleading number this pipeline could emit".
    _partition_known = db_entrap is not None and db_entrap > 0
    elias_gygi_concat = (2 * entrap / denom) if (denom and _partition_known) else None
    # Retained under its old name so nothing downstream breaks, but it is the
    # same quantity as elias_gygi_concat: equal-odds == combined at r == 1.
    wen = elias_gygi_concat
    shared_sens = (
        (entrap + shared) / (entrap + shared + target) if (entrap + shared + target) else None
    )

    # eq 1, COMBINED UPPER BOUND    E(1 + 1/r) / (T + E)
    # Bounds the FDR (the expectation), not the FDP. Valid at any r.
    combined_noble = None
    size_scaled = None
    if db_entrap:  # zero or None entrapment partition -> undefined, not zero
        # Kept as a diagnostic under an honest name. It is NOT Elias & Gygi:
        # at r = 1 it collapses to E/(E+T), which is the eq-2 lower bound.
        scaled = entrap * (db_target / db_entrap)
        size_scaled = scaled / (scaled + target) if (scaled + target) else None
        if ratio and denom:
            combined_noble = entrap * (1 + 1 / ratio) / denom

    return {
        "n_entrapment": entrap,
        "n_target": target,
        "n_shared": shared,
        "db_target": db_target,
        "db_entrapment": db_entrap,
        "r": round(ratio, 6) if ratio is not None else None,
        # --- the suite, each named for the equation it implements -------------
        "fdp_lower_bound_eq2": lower_bound,  # E/(T+E)   bounds the FDP
        "fdp_combined_upper_eq1": combined_noble,  # E(1+1/r)/(T+E)  bounds the FDR
        "fdp_elias_gygi_concatenated": elias_gygi_concat,  # 2E/(T+E)
        "fdp_sample_method_eq3": sample_method,  # E/T  NOT A BOUND
        "fdp_size_scaled_diagnostic": size_scaled,  # diagnostic, not an estimator
        # --- retained keys so existing readers do not break -------------------
        # fdp_uncorrected and fdp_lower_bound were the SAME quantity under two
        # names, and fdp_lower_bound carried E/T, which is not a lower bound.
        "fdp_uncorrected": lower_bound,
        "fdp_wen_equal_odds": wen,
        "fdp_lower_bound": lower_bound,
        # CORRECTED VALUE under the old name. It used to hold the size-scaled
        # quantity, which at r=1 collapses to E/(E+T) -- the lower bound, not
        # Elias & Gygi. It now holds the concatenated form 2E/(T+E). Any report
        # comparing an old run to a new one will see this key roughly double,
        # and that change is the fix, not a regression.
        "fdp_elias_gygi": elias_gygi_concat,
        "fdp_shared_sensitivity": shared_sens,
        "fdp_combined_noble": combined_noble,
        "evaluable": denom > 0,
        "entrapment_in_db": None if db_entrap is None else db_entrap > 0,
    }


# The suite, one entry per equation. Extended 2026-08-02 when the estimators were
# reconciled against the version of record; the four legacy names are kept because
# existing reports and figure scripts read them by name.
ESTIMATOR_KEYS = (
    # --- correctly named, use these -----------------------------------------
    "fdp_lower_bound_eq2",  # E/(T+E)          bounds the FDP itself
    "fdp_combined_upper_eq1",  # E(1+1/r)/(T+E)   bounds the FDR
    "fdp_elias_gygi_concatenated",  # 2E/(T+E)
    "fdp_sample_method_eq3",  # E/T              NOT a bound
    "fdp_size_scaled_diagnostic",  # diagnostic only
    # --- legacy names, retained so existing readers keep working -------------
    "fdp_uncorrected",
    "fdp_wen_equal_odds",
    "fdp_elias_gygi",
    "fdp_lower_bound",
    "fdp_shared_sensitivity",
    "fdp_combined_noble",
)


def bootstrap_median_ci(
    values: list[float],
    seed: int = 42,
    n_draws: int = 2000,
    alpha: float = 0.05,
) -> tuple[float, float]:
    """Seeded percentile bootstrap CI on the median.

    Returns ``(nan, nan)`` when ``n <= 5`` -- a bootstrap CI on five or fewer
    observations is not informative and reporting one implies precision that is not
    there. Deterministic for a given ``seed``.
    """
    n = len(values)
    if n <= 5:
        return math.nan, math.nan
    rng = random.Random(seed)
    medians = []
    for _ in range(n_draws):
        draw = [values[rng.randrange(n)] for _ in range(n)]
        medians.append(statistics.median(draw))
    medians.sort()
    lo = medians[int((alpha / 2) * n_draws)]
    hi = medians[min(int((1 - alpha / 2) * n_draws), n_draws - 1)]
    return lo, hi


def summarize(values: list[float], seed: int = 42, n_draws: int = 2000) -> dict:
    """n, median, IQR and a seeded bootstrap CI for a list of per-sample values."""
    vals = sorted(v for v in values if v is not None and not math.isnan(v))
    if not vals:
        return {
            "n": 0,
            "median": None,
            "q1": None,
            "q3": None,
            "iqr": None,
            "ci_lo": None,
            "ci_hi": None,
        }
    median = statistics.median(vals)
    if len(vals) >= 4:
        q1, _, q3 = statistics.quantiles(vals, n=4, method="inclusive")
    else:
        q1 = q3 = median
    lo, hi = bootstrap_median_ci(vals, seed=seed, n_draws=n_draws)
    return {
        "n": len(vals),
        "median": median,
        "q1": q1,
        "q3": q3,
        "iqr": q3 - q1,
        "ci_lo": None if math.isnan(lo) else lo,
        "ci_hi": None if math.isnan(hi) else hi,
    }


# ---------------------------------------------------------------------------
# Per-sample FDP over a results directory
# ---------------------------------------------------------------------------


def expected_results_subdir(class_name: str) -> str:
    """The per-sample subdirectory the run templates write for an entrapment class.

    ``tune`` and ``diagnose`` document the layout ``SAMPLE/entrapment/`` for the
    shuffled class and ``SAMPLE/entrapment_<class>/`` for every other class.
    """
    return "entrapment" if class_name == "shuffled" else f"entrapment_{class_name}"


def find_sample_results(
    results_dir: str | Path,
    results_name: str = "results.sage.tsv",
    prefer_subdir: str | None = None,
) -> list[tuple[str, Path]]:
    """Find ``(sample_id, results_tsv)`` under ``results_dir``.

    Accepts either ``results_dir/SAMPLE/results.sage.tsv`` or the same file one
    level deeper (``results_dir/SAMPLE/sage/results.sage.tsv``), which is the layout
    the SLURM templates in this project write.

    A sample directory may hold several entrapment runs side by side
    (``entrapment/``, ``entrapment_plant/``, ``entrapment_archaea/``). When more
    than one nested candidate exists, the one in ``prefer_subdir`` is used; if
    that does not single out one file the sample is refused with the candidates
    listed, because silently taking the first would label one class's search
    results with another class's name.
    """
    root = Path(results_dir)
    found: list[tuple[str, Path]] = []
    for sample_dir in sorted(p for p in root.iterdir() if p.is_dir()):
        direct = sample_dir / results_name
        if direct.exists():
            found.append((sample_dir.name, direct))
            continue
        nested = sorted(sample_dir.glob(f"*/{results_name}"))
        if not nested:
            continue
        if len(nested) == 1:
            run_dir = nested[0].parent.name
            if prefer_subdir and run_dir.startswith("entrapment") and run_dir != prefer_subdir:
                raise ValueError(
                    f"sample {sample_dir.name!r}: the only {results_name} is in {run_dir!r}, "
                    f"which is a different entrapment run from the class being reported "
                    f"({prefer_subdir!r}). Pass --results-name '{run_dir}/{results_name}' "
                    "to report that run under its own class."
                )
            found.append((sample_dir.name, nested[0]))
            continue
        preferred = [p for p in nested if prefer_subdir and p.parent.name == prefer_subdir]
        if len(preferred) == 1:
            found.append((sample_dir.name, preferred[0]))
            continue
        candidates = ", ".join(p.parent.name for p in nested)
        raise ValueError(
            f"sample {sample_dir.name!r} holds {len(nested)} candidate {results_name} "
            f"files ({candidates})"
            + (
                f" and none is in the expected subdirectory {prefer_subdir!r}"
                if prefer_subdir
                else ""
            )
            + ". Which run the FDP is computed on will not be guessed: pass "
            f"--results-name '<subdir>/{results_name}' to choose one."
        )
    return found


def _find_search_db(sample_dir: Path, pattern: str | None) -> Path | None:
    """Return the first sorted matching search database, or None if absent."""
    if not pattern:
        return None
    hits = sorted(sample_dir.glob(pattern))
    return hits[0] if hits else None


def compute_fdp(
    results_dir: str | Path,
    dataset: str = "",
    entrap_class: str | None = None,
    stage: str | None = None,
    prefix: str | None = None,
    target_fasta: str | Path | None = None,
    q_threshold: float = 0.01,
    min_length: int = 9,
    q_column: str = "peptide_q",
    decoy_prefix: str = "rev_",
    db_pattern: str | None = None,
    results_name: str = "results.sage.tsv",
    seed: int = 42,
    n_draws: int = 2000,
    peptide_level: bool = False,
    allow_mixed: bool = False,
    conserved_check: bool = True,
    normalize_il: bool = True,
) -> dict:
    """Compute one FDP record for a (dataset, class, stage) triple.

    Reads any ``entrapment_manifest.json`` found beside the results to recover the
    class and injection stage, and refuses to pool samples that differ in either.
    Reports how many samples were **evaluable** and how many had an empty
    entrapment partition, so an undefined estimate can never be read as a zero.

    When ``target_fasta`` is given, every entrapment-attributed peptide is tested
    for occurrence in the target sequence space and both the raw and the
    conserved-corrected FDP are returned. For a class with real homology pressure
    the correction is mandatory: reporting the raw number alone would present
    genuinely conserved peptides as false discoveries.
    """
    root = Path(results_dir)
    if not root.is_dir():
        raise ValueError(f"results directory does not exist: {root}")

    manifest_paths = [m for m in sorted(root.glob(f"**/{MANIFEST_NAME}")) if m.is_file()]
    if manifest_paths:
        resolved, manifests = load_manifests(manifest_paths, allow_mixed=allow_mixed)
    else:
        resolved, manifests = {"stage": None, "class": None}, []

    stage = _reconcile("stage", stage, resolved["stage"])
    entrap_class = _reconcile("class", entrap_class, resolved["class"])

    if stage is None:
        raise ValueError(
            "injection stage is unknown: no single stage could be read from the "
            "manifests and --stage was not given. The stage determines what the FDP "
            "means (construction error vs search error) and will not be guessed."
            + (
                " You passed --allow-mixed, so you must also state the stage the "
                "pooled number is to be labelled with."
                if allow_mixed
                else ""
            )
        )
    if entrap_class is None:
        raise ValueError(
            "entrapment class is unknown: no single class could be read from the "
            "manifests and --class was not given. The class determines the "
            "difficulty of the test and what its FDP estimates, and will not be "
            "guessed."
            + (
                " You passed --allow-mixed, so you must also state the class the "
                "pooled number is to be labelled with."
                if allow_mixed
                else ""
            )
        )

    stage = normalize_stage(stage)
    stage_meta = get_stage(stage)
    cls = get_class(entrap_class)

    # Locate the per-sample results only now, so that a sample holding several
    # entrapment runs is resolved by the class being reported, never by glob order.
    samples = find_sample_results(
        root, results_name=results_name, prefer_subdir=expected_results_subdir(cls.name)
    )
    if not samples:
        raise ValueError(f"no {results_name} found under {root}")

    # An explicit --entrap-prefix wins; otherwise prefer what the manifest recorded
    # (the run may not use this class's default convention), then the class default.
    manifest_prefix = next(
        (
            m.get("entrapment_match_pattern") or m.get("prefix")
            for m in manifests
            if m.get("entrapment_match_pattern") or m.get("prefix")
        ),
        None,
    )
    prefix_source = "explicit" if prefix else "manifest" if manifest_prefix else "class default"
    prefix = prefix or manifest_prefix or cls.prefix

    if conserved_check and target_fasta is None:
        if cls.difficulty == "cross-kingdom-homologous":
            raise ValueError(
                f"class {cls.name!r} carries real homology pressure: its raw FDP is "
                "dominated by universally conserved peptides that are TRUE "
                "identifications. Pass --target-fasta so the conserved-corrected "
                "FDP can be computed, or --no-conserved-check to state explicitly "
                "that you are reporting an uncorrected number."
            )
        logger.warning(
            "No --target-fasta given: the conserved-peptide correction will not be "
            "computed and only an UNCORRECTED FDP will be reported."
        )

    # Partition sizes come from the manifest that governs a sample: the one in
    # its own directory or the nearest ancestor. Under per-sample (appended)
    # injection every sample has its own; applying the first manifest found
    # to all of them scaled the r-dependent estimators by one sample's r.
    manifest_db_target = manifests[0].get("n_target") if len(manifests) == 1 else None
    manifest_db_entrap = manifests[0].get("n_entrapment") if len(manifests) == 1 else None
    manifest_by_dir = {Path(p).resolve().parent: m for p, m in zip(manifest_paths, manifests)}

    def _governing_manifest(tsv: Path) -> tuple[dict | None, str]:
        """Find the nearest governing injection manifest, with a single-manifest fallback."""
        directory = tsv.resolve().parent
        stop = root.resolve()
        while True:
            if directory in manifest_by_dir:
                return manifest_by_dir[directory], f"injection manifest in {directory}"
            if directory == stop or directory.parent == directory:
                break
            directory = directory.parent
        if len(manifests) == 1:
            return manifests[0], "the single injection manifest"
        return None, "no manifest governs this sample"

    # Prefix-mismatch gate. The check is against the INJECTED database, never the
    # per-sample searched one: under pre-selection injection an empty per-sample
    # entrapment partition is the expected, desirable result (entrapment depleted
    # by curation), so it carries no information about the prefix. The injected
    # database, by contrast, must contain what the manifest says it contains.
    prefix_check = {"performed": False, "checked_fasta": "", "n_found": None}
    for m in manifests:
        injected = m.get("output_fasta")
        declared = m.get("n_entrapment") or 0
        if injected and declared and Path(injected).exists():
            n_found = verify_prefix_matches(
                injected,
                prefix,
                declared,
                context="injected database",
            )
            prefix_check = {
                "performed": True,
                "checked_fasta": injected,
                "n_declared": declared,
                "n_found": n_found,
            }
            break

    per_sample = []
    all_entrap_peptides: set[str] = set()
    entrap_psm_counts: Counter = Counter()

    for sample_id, tsv in samples:
        counts = classify_psms(
            tsv,
            prefix=prefix,
            q_threshold=q_threshold,
            min_length=min_length,
            q_column=q_column,
            decoy_prefix=decoy_prefix,
        )
        db = _find_search_db(tsv.parent, db_pattern) or _find_search_db(
            tsv.parent.parent, db_pattern
        )
        if db is not None:
            db_target, db_entrap = partition_counts(db, prefix)
            db_size_source = "per-sample searched database (header lines)"
        else:
            # Falling back to the injected lake's partition sizes. Under a
            # selection pipeline the per-sample database is far smaller than the
            # lake and its entrapment fraction is different, so the size-corrected
            # estimators computed from these are NOT per-sample quantities. Say so
            # rather than letting the number pass as if it were.
            governing, where = _governing_manifest(tsv)
            db_target = governing.get("n_target") if governing else None
            db_entrap = governing.get("n_entrapment") if governing else None
            db_size_source = (
                f"{where} (injected database, not the per-sample searched one -- "
                "size-corrected estimators are approximate; pass --db-pattern for "
                "exact values)"
                if db_target is not None
                else f"unknown ({where})"
            )

        if peptide_level:
            e, t, s = (
                len(counts.peptides_entrap),
                len(counts.peptides_target),
                len(counts.peptides_shared),
            )
        else:
            e, t, s = counts.entrap, counts.target, counts.shared

        row = fdp_estimates(e, t, s, db_target=db_target, db_entrap=db_entrap)
        row.update(
            {
                "sample": sample_id,
                "dataset": dataset,
                "class": cls.name,
                "stage": stage,
                "n_decoy": counts.decoy,
                "decoy_source": counts.decoy_source,
                "db_size_source": db_size_source,
                "search_db": str(db) if db else "",
                "results_file": str(tsv),
                "counting": "peptide" if peptide_level else "spectrum",
            }
        )
        row["_counts"] = counts
        per_sample.append(row)
        all_entrap_peptides |= counts.peptides_entrap
        entrap_psm_counts.update(counts.entrap_psm_counts)

    # ---- conserved-peptide correction -------------------------------------
    scan = ConservedScan()
    if target_fasta is not None and conserved_check:
        scan = conserved_peptide_scan(
            all_entrap_peptides,
            target_fasta,
            normalize_il=normalize_il,
        )
        conserved = set(scan.conserved)
        for row in per_sample:
            counts = row["_counts"]
            if peptide_level:
                e_cons = len(counts.peptides_entrap & conserved)
            else:
                e_cons = sum(n for p, n in counts.entrap_psm_counts.items() if p in conserved)
            e_spec = row["n_entrapment"] - e_cons
            row["n_entrapment_conserved"] = e_cons
            row["n_entrapment_class_specific"] = e_spec
            corrected = fdp_estimates(
                e_spec,
                row["n_target"],
                row["n_shared"],
                db_target=row["db_target"],
                db_entrap=row["db_entrapment"],
            )
            for key in ESTIMATOR_KEYS:
                row[f"{key}_conserved_corrected"] = corrected[key]
    else:
        for row in per_sample:
            row["n_entrapment_conserved"] = None
            row["n_entrapment_class_specific"] = None
            for key in ESTIMATOR_KEYS:
                row[f"{key}_conserved_corrected"] = None

    for row in per_sample:
        row.pop("_counts", None)

    n_total = len(per_sample)
    summaries = {}
    for key in ESTIMATOR_KEYS:
        summaries[key] = {
            "raw": summarize(
                [r[key] for r in per_sample if r[key] is not None], seed=seed, n_draws=n_draws
            ),
            "conserved_corrected": summarize(
                [
                    r[f"{key}_conserved_corrected"]
                    for r in per_sample
                    if r.get(f"{key}_conserved_corrected") is not None
                ],
                seed=seed,
                n_draws=n_draws,
            ),
        }

    pooled_e = sum(r["n_entrapment"] for r in per_sample)
    pooled_t = sum(r["n_target"] for r in per_sample)
    pooled_s = sum(r["n_shared"] for r in per_sample)
    pooled_e_spec = sum(r["n_entrapment_class_specific"] or 0 for r in per_sample)

    conserved_block = {
        "performed": scan.performed,
        "status": (
            "computed"
            if scan.performed
            else "NOT PERFORMED -- reported FDP is UNCORRECTED for conserved peptides"
        ),
        "target_fasta": scan.target_fasta,
        "normalize_il": normalize_il,
        "n_entrapment_peptides_tested": scan.n_tested,
        "n_conserved": scan.n_conserved,
        "conserved_rate": scan.rate,
        "negative_control": {
            "description": (
                "same shared-peptide test on the REVERSED entrapment peptides; should be near zero"
            ),
            "n_tested": scan.n_control_tested,
            "n_found": scan.n_control_conserved,
            "rate": scan.control_rate,
            "warning": scan.warning,
        },
    }
    if cls.is_definitional_null and scan.performed and scan.rate:
        conserved_block["false_exculpation_calibration"] = {
            "rate": scan.rate,
            "note": (
                f"{cls.name} is a null BY CONSTRUCTION -- no peptide from it is a "
                f"true identification. The {scan.rate:.1%} of its peptides that the "
                "conserved test exculpates is therefore a direct measurement of the "
                "correction's FALSE-EXCULPATION RATE, and bounds how much of any "
                "other class's exculpation is the same artifact."
            ),
        }

    return {
        "schema": RECORD_SCHEMA,
        "created": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "dataset": dataset,
        "class": cls.name,
        "class_metadata": cls.as_metadata(),
        "stage": stage,
        "stage_meaning": STAGE_MEANING[stage],
        "stage_metadata": stage_meta.as_metadata(),
        "measures": stage_measures(stage),
        "is_construction_test": stage_meta.is_construction_test,
        "counting": "peptide" if peptide_level else "spectrum",
        "counting_note": (
            "Counts are spectrum-level (one per PSM): a confidently repeated "
            "peptide contributes proportionally to its spectral count, so the "
            "effective n is below the nominal count."
            if not peptide_level
            else "Counts are peptide-level (deduplicated)."
        ),
        "filters": {
            "q_column": q_column,
            "q_threshold": q_threshold,
            "min_length": min_length,
            "prefix": prefix,
            "decoy_prefix": decoy_prefix,
        },
        # Which pattern classified a record as entrapment, where it came from, and
        # whether it was verified against the injected database. Recorded because
        # a wrong pattern silently reclassifies entrapment as target and drives the
        # FDP toward zero.
        "entrapment_matching": {
            "pattern": prefix,
            "pattern_source": (
                "explicit --entrap-prefix"
                if prefix_source == "explicit"
                else "injection manifest"
                if prefix_source == "manifest"
                else f"default for class {cls.name!r}"
            ),
            "verified_against_injected_database": prefix_check["performed"],
            "verification": prefix_check,
            "note": (
                ""
                if prefix_check["performed"]
                else "NOT VERIFIED: the injected database named in the manifest was not "
                "available, so a prefix mismatch could not be ruled out. If every "
                "sample reports zero entrapment, confirm the convention with "
                "`fasta-lake entrapment prefixes <fasta>` before trusting the FDP."
            ),
        },
        "coverage": {
            "n_attempted": n_total,
            "n_evaluable": sum(1 for r in per_sample if r["evaluable"]),
            "n_not_evaluable": sum(1 for r in per_sample if not r["evaluable"]),
            "n_empty_entrapment_partition": sum(
                1 for r in per_sample if r["entrapment_in_db"] is False
            ),
            "n_unknown_entrapment_partition": sum(
                1 for r in per_sample if r["entrapment_in_db"] is None
            ),
        },
        "pooled": {
            "raw": fdp_estimates(
                pooled_e,
                pooled_t,
                pooled_s,
                db_target=manifest_db_target,
                db_entrap=manifest_db_entrap,
            ),
            "conserved_corrected": (
                fdp_estimates(
                    pooled_e_spec,
                    pooled_t,
                    pooled_s,
                    db_target=manifest_db_target,
                    db_entrap=manifest_db_entrap,
                )
                if scan.performed
                else None
            ),
        },
        "conserved_peptides": conserved_block,
        "summary": summaries,
        "per_sample": per_sample,
        "seed": seed,
        "bootstrap_draws": n_draws,
        "_scan": scan,
    }


def _reconcile(what: str, given: str | None, from_manifest: str | None) -> str | None:
    """Prefer an explicit value, but never let it silently contradict a manifest."""
    if given and from_manifest and given != from_manifest:
        raise ValueError(
            f"--{what} {given!r} contradicts the {what} recorded in the entrapment "
            f"manifest ({from_manifest!r}). Refusing to mislabel the measurement."
        )
    return given or from_manifest


def write_fdp_reports(report: dict, output_dir: str | Path) -> dict[str, Path]:
    """Write the per-sample TSV, the (dataset, class, stage) record, and conserved peptides."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    written: dict[str, Path] = {}

    cols = [
        "dataset",
        "sample",
        "class",
        "stage",
        "counting",
        "n_entrapment",
        "n_entrapment_conserved",
        "n_entrapment_class_specific",
        "n_target",
        "n_shared",
        "n_decoy",
        "db_target",
        "db_entrapment",
        "r",
        "evaluable",
        "entrapment_in_db",
        *ESTIMATOR_KEYS,
        *[f"{k}_conserved_corrected" for k in ESTIMATOR_KEYS],
        "decoy_source",
        "db_size_source",
        "search_db",
        "results_file",
    ]
    tsv_path = out / "entrapment_fdp_per_sample.tsv"
    with open(tsv_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for row in report["per_sample"]:
            w.writerow({k: ("" if row.get(k) is None else row.get(k)) for k in cols})
    written["per_sample"] = tsv_path

    record = {k: v for k, v in report.items() if k not in ("per_sample", "_scan")}

    # 🔴 THE NEGATIVE CONTROL MUST REACH THE RECORD, OR THE CORRECTION HAS NO EVIDENCE.
    #
    # `_scan` was dropped wholesale -- correctly for its per-peptide dictionaries, which are
    # large and already written to entrapment_conserved_peptides.tsv, but it took the
    # reversed-peptide control with it. Measured 2026-08-08: the plant correction subtracted
    # 1,485 of 13,556 peptides (10.95%) and NOTHING in the record showed whether the control
    # had passed, or whether it had run at all. It had -- 170/13,556 (1.25%), a ratio of
    # 0.11 against the 0.25 warn threshold -- but that evidence existed only in a job log.
    #
    # A correction a reader cannot audit is a correction they must take on trust. The
    # SUMMARY is small; it goes in the record. The peptide lists stay in the TSV.
    _s = report.get("_scan")
    if _s is not None and getattr(_s, "performed", False):
        record["conserved_scan"] = {
            "target_fasta": _s.target_fasta,
            "n_tested": _s.n_tested,
            "n_conserved": _s.n_conserved,
            "rate": _s.rate,
            "n_control_tested": _s.n_control_tested,
            "n_control_conserved": _s.n_control_conserved,
            "control_rate": _s.control_rate,
            "control_warn_ratio": 0.25,
            "control_passed": (
                None
                if not _s.rate or _s.control_rate is None
                else _s.control_rate <= 0.25 * _s.rate
            ),
            "warning": _s.warning or "",
            "note": (
                "control_rate is the fraction of REVERSED entrapment peptides also found "
                "in the target space. Reversed peptides are not conserved sequence, so a "
                "control rate approaching the real rate means the matcher is finding "
                "chance matches and the conserved-corrected FDP is unreliable."
            ),
        }
    record_path = out / "entrapment_fdp_record.json"
    record_path.write_text(json.dumps(record, indent=2, default=str) + "\n")
    written["record"] = record_path

    scan = report.get("_scan")
    if scan is not None and scan.performed:
        written["conserved"] = write_conserved_report(
            scan, out / "entrapment_conserved_peptides.tsv"
        )
    return written


# ---------------------------------------------------------------------------
# Config threading
# ---------------------------------------------------------------------------


def prepare_from_config(
    cfg,
    database_fasta: str | Path,
    output_dir: str | Path,
) -> dict | None:
    """Honour ``entrapment.enabled`` in a FastaLake JSON config.

    Given a loaded :class:`~fasta_lake.config.FastaLakeConfig` (or anything with an
    ``.entrapment`` attribute) and the database about to be used, generate the
    configured entrapment and inject it, returning the manifest. Returns ``None``
    when entrapment is disabled, in which case the caller uses ``database_fasta``
    unchanged.

    This is the function that makes ``"entrapment": {"enabled": true}`` in a config
    JSON actually do something.
    """
    ec = getattr(cfg, "entrapment", None)
    if ec is None or not getattr(ec, "enabled", False):
        return None

    stage = normalize_stage(getattr(ec, "stage", STAGE_PRESELECT) or STAGE_PRESELECT)

    name = _config_class_name(ec)
    cls = get_class(name)
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    seed = getattr(ec, "seed", 42)

    entrap_fasta = out / "entrapment.fasta"
    generate_entrapment_fasta(
        entrap_fasta,
        entrap_class=cls.name,
        source_fasta=(database_fasta if cls.generation == "shuffle" else (ec.source or None)),
        seed=seed,
        n_proteins=(ec.n_proteins or None),
    )

    return inject_entrapment(
        database_fasta,
        entrap_fasta,
        out / "augmented_database.fasta",
        stage=stage,
        entrap_class=cls.name,
        seed=seed,
        manifest_path=out / MANIFEST_NAME,
    )


def _config_class_name(ec) -> str:
    """Map an :class:`EntrapmentConfig` onto a registered class name.

    ``method`` predates the class registry; ``'phylogenetic'`` is the historical
    name for "a real distant organism supplied by the user", which is the ``user``
    class once a source FASTA is given.
    """
    name = getattr(ec, "entrap_class", "") or ec.method
    if name in ENTRAPMENT_CLASSES:
        return name
    if name == "phylogenetic":
        return "user"
    raise ValueError(
        f"cannot map entrapment.method {name!r} onto an entrapment class. "
        f"Set entrapment.class to one of: {', '.join(sorted(ENTRAPMENT_CLASSES))}"
    )
