"""Group accepted study peptides and assign LFQ features directly by peptide.

The greedy cover follows the reviewed publication procedure. A representative
labels a group and explains its assigned peptides; it is not proof of their
unique biological origin. Pre-search selection and search confidence are unchanged.
"""

from __future__ import annotations

import csv
import hashlib
import heapq
import json
import math
import re
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

from fasta_lake.gzip_io import open_gzip_text

_MOD = re.compile(r"\[[^\]]*\]")
_LFQ_FIXED = {"peptide", "charge", "proteins", "q_value", "score", "spectral_angle"}


def canonical_peptide(peptide: str) -> str:
    """Remove Sage modifications and use the declared I/L equivalence.

    Sage writes residue mods in brackets (``PEPT[+15.9949]LDEK``) and N-terminal
    mods as a bracket followed by ``-`` (``[+42.0106]-PEPTLDEK``); both go.
    """
    return _MOD.sub("", peptide).replace("-", "").replace("I", "L")


def member_set_key(proteins: str) -> str:
    """Preserve the full distinct accession membership, without source-tag stripping."""
    return ";".join(sorted({p.strip() for p in proteins.split(";") if p.strip()}))


def greedy_set_cover(protein_peptides: dict[str, set[str]]) -> tuple[list, set[str]]:
    """Choose study representatives by maximum uncovered peptide gain.

    Parameters
    ----------
    protein_peptides : dict[str, set[str]]
        Observed canonical peptides for each compatible protein accession.

    Returns
    -------
    tuple
        Selection-ordered (accession, newly covered peptide count) pairs and the
        set of covered peptides. Equal gains resolve by lexical accession.

    Notes
    -----
    Gains are recalculated as evidence is covered. This is post-search reporting,
    distinct from pre-search razor's fixed evidence counts. Greedy cover does
    not guarantee the smallest possible protein set; inputs are not modified.
    """
    heap = [(-len(peptides), protein) for protein, peptides in protein_peptides.items()]
    heapq.heapify(heap)
    covered, selected = set(), []
    while heap:
        negative_gain, protein = heapq.heappop(heap)
        gain = len(protein_peptides[protein] - covered)
        if not gain:
            continue
        if gain != -negative_gain:
            heapq.heappush(heap, (-gain, protein))
            continue
        selected.append((protein, gain))
        covered.update(protein_peptides[protein])
    return selected, covered


def study_dictionary(protein_peptides: dict[str, set[str]], searched: set[str]) -> dict:
    """Assign accepted peptides to a fixed set of study representatives.

    Parameters
    ----------
    protein_peptides : dict[str, set[str]]
        Accepted canonical peptides associated with observed target accessions.
    searched : set[str]
        Accessions in the union of FASTAs actually searched for the study.

    Returns
    -------
    dict
        Selected representatives, covered peptides, selection ranks,
        peptide_to_group assignments and unobserved searched accessions.

    Raises
    ------
    ValueError
        If an observed accession is outside the searched universe.

    Notes
    -----
    Each peptide belongs to the first selected representative explaining it.
    A protein may be compatible with peptides assigned to several groups.
    This function checks graph membership; group_searches additionally verifies
    peptide containment in the actual representative sequences.
    """
    outside = set(protein_peptides) - searched
    if outside:
        raise ValueError(f"{len(outside)} observed proteins are outside the searched universe")
    selected, covered = greedy_set_cover(protein_peptides)
    peptide_to_group = {}
    for protein, _ in selected:
        for peptide in sorted(protein_peptides[protein]):
            peptide_to_group.setdefault(peptide, protein)
    return dict(
        selected=selected,
        covered=covered,
        rank={protein: i for i, (protein, _) in enumerate(selected)},
        peptide_to_group=peptide_to_group,
        unobserved=searched - set(protein_peptides),
    )


def _passing(value: str, threshold: float) -> bool:
    """Accept only a finite q value between zero and the inclusive threshold."""
    try:
        q = float(value)
        return math.isfinite(q) and 0 <= q <= threshold
    except (ValueError, TypeError):
        return False


def _sage_filename(value: str) -> str:
    """Sage emits the basename including mzML extension; accept recorded full paths."""
    return value.replace("\\", "/").rsplit("/", 1)[-1]


def _rows(path: Path, required: set[str], *, registry: list, acquisition=None, lfq=False):
    """Hash the actual bytes read, then parse strict, tab-delimited records."""
    digest = hashlib.sha256()
    with path.open("rb") as stream:

        def lines():
            """Decode TSV lines while hashing exactly the bytes consumed by the parser."""
            for raw in stream:
                digest.update(raw)
                yield raw.decode("utf-8")

        reader = csv.DictReader(lines(), delimiter="\t")
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f"{path}: missing columns {sorted(required)}")
        if lfq:
            columns = [c for c in reader.fieldnames if c not in _LFQ_FIXED]
            if len(columns) != 1 or _sage_filename(columns[0]) != acquisition:
                raise ValueError(f"{path}: LFQ acquisition identity does not match {acquisition}")
        for row in reader:
            if None in row or any(v is None for v in row.values()):
                raise ValueError(f"{path}: malformed row ending at line {reader.line_num}")
            if acquisition is not None and not lfq:
                if _sage_filename(row.get("filename", "")) != acquisition:
                    raise ValueError(
                        f"{path}: search acquisition identity does not match {acquisition}"
                    )
            yield reader.line_num, row, reader.fieldnames
    registry.append({"path": str(path), "sha256": digest.hexdigest()})


def _read_fasta(path: Path, sequences: dict[str, str], registry: list) -> set[str]:
    """Verify accession identity before merging acquisition evidence."""
    members, digest = set(), hashlib.sha256()
    accession, parts = None, []

    def retain():
        """Merge a nonempty FASTA record after checking local IDs and sequence identity."""
        if accession is None:
            return
        if accession in members:
            raise ValueError(f"{path}: duplicate FASTA accession {accession}")
        sequence = "".join(c for c in "".join(parts).upper() if "A" <= c <= "Z")
        if not sequence:
            raise ValueError(f"{path}: empty sequence for {accession}")
        if accession in sequences and sequences[accession] != sequence:
            raise ValueError(f"{path}: accession-to-sequence conflict for {accession}")
        sequences[accession] = sequence
        members.add(accession)

    with path.open("rb") as stream:
        for raw in stream:
            digest.update(raw)
            line = raw.decode("utf-8").strip()
            if line.startswith(">"):
                retain()
                fields = line[1:].split()
                if not fields:
                    raise ValueError(f"{path}: empty FASTA header")
                accession, parts = fields[0], []
            elif line:
                if accession is None:
                    raise ValueError(f"{path}: sequence precedes first header")
                parts.append(line)
    retain()
    if not members:
        raise ValueError(f"{path}: empty searched FASTA")
    registry.append({"path": str(path), "sha256": digest.hexdigest()})
    return members


def _members(field: str) -> set[str]:
    """Collect distinct nonempty accessions, excluding the rev_ decoy prefix."""
    return {p for p in field.split(";") if p and not p.startswith("rev_")}


def _write_table(path: Path, columns: list[str], rows):
    """Create a TSV exclusively, preserving the caller's column and row order."""
    with path.open("x", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(columns)
        writer.writerows(rows)


def group_searches(
    search_root: str | Path,
    output: str | Path,
    q_threshold: float = 0.01,
    *,
    config_name: str = "config.json",
    acquisition_names: list[str] | None = None,
) -> dict:
    """Build study groups and sum quantities from completed acquisition searches.

    Parameters
    ----------
    search_root : str or Path
        Directory containing one completed Sage search subdirectory per acquisition.
        Each contains a configuration, results.sage.tsv and single-acquisition lfq.tsv.
    output : str or Path
        New directory outside search_root. Existing paths are refused.
    q_threshold : float
        Inclusive target peptide-q threshold in (0, 1], applied within each search.
    config_name : {"config.json", "results.json"}
        Filename recording the exact searched FASTA and one mzML path.
    acquisition_names : list[str] or None
        Explicit acquisition directory names. None uses all direct subdirectories;
        the final roster is sorted and must contain unique, nonempty names.

    Returns
    -------
    dict
        The exported summary: roster, grouping rule, counts, intensity accounting
        and checks. Tables, representatives, peptide evidence and input hashes are
        written beside summary.json.

    Raises
    ------
    ValueError
        For inconsistent acquisitions, malformed evidence, accession/sequence
        conflicts, invalid confidence settings or incompatible representatives.
    FileExistsError
        If the requested output already exists.

    Notes
    -----
    Greedy cover uses only accepted evidence and actually searched proteins.
    Positive finite LFQ features require an accepted peptide in that acquisition.
    Each feature is assigned once; missing cells stay blank. Local members and
    representative presence are retained separately. This adds no study-wide
    or protein-group FDR estimate. Invalid inputs leave no completed bundle.
    """
    if not math.isfinite(q_threshold) or not 0 < q_threshold <= 1:
        raise ValueError("q_threshold must be in (0, 1]")
    if config_name not in {"config.json", "results.json"}:
        raise ValueError("config_name must be config.json or results.json")
    requested = Path(output).absolute()
    if requested.exists() or requested.is_symlink():
        raise FileExistsError(f"Output exists: {requested}; choose a new directory")
    root, out = Path(search_root).resolve(), requested.resolve()
    if root == out or root in out.parents:
        raise ValueError("Output must be outside the search input directory")
    if acquisition_names is None:
        samples = sorted(p for p in root.iterdir() if p.is_dir())
    else:
        if not acquisition_names or len(acquisition_names) != len(set(acquisition_names)):
            raise ValueError("Acquisition names must be nonempty and unique")
        if any(Path(name).name != name or name in {".", ".."} for name in acquisition_names):
            raise ValueError("Acquisition names must be directory names, not paths")
        samples = [root / name for name in sorted(acquisition_names)]
        if any(not sample.is_dir() for sample in samples):
            raise ValueError("A declared acquisition directory is missing")
    if not samples:
        raise ValueError("No acquisition directories found")
    graph, sequences, gates, searched_by_sample = defaultdict(set), {}, {}, {}
    registry = []
    acquisitions = {}
    for sample in samples:
        config = sample / config_name
        raw = config.read_bytes()
        registry.append({"path": str(config), "sha256": hashlib.sha256(raw).hexdigest()})
        cfg = json.loads(raw)
        paths = cfg.get("mzml_paths")
        if not isinstance(paths, list) or len(paths) != 1 or not isinstance(paths[0], str):
            raise ValueError(f"{sample}: config must record exactly one mzml_paths entry")
        acquisition = _sage_filename(paths[0])
        if not acquisition:
            raise ValueError(f"{sample}: empty acquisition identity")
        acquisitions[sample.name] = acquisition
        fasta = Path(cfg["database"]["fasta"])
        if not fasta.is_absolute():
            raise ValueError(f"{sample}: config must record an absolute searched FASTA path")
        searched = _read_fasta(fasta, sequences, registry)
        searched_by_sample[sample.name] = searched
        gate = set()
        for _, row, _ in _rows(
            sample / "results.sage.tsv",
            {"peptide", "proteins", "label", "peptide_q", "filename"},
            registry=registry,
            acquisition=acquisition,
        ):
            if row["label"] != "1" or not _passing(row["peptide_q"], q_threshold):
                continue
            peptide = canonical_peptide(row["peptide"])
            if not peptide:
                continue
            members = _members(row["proteins"])
            if members - searched:
                raise ValueError(f"{sample.name}: accepted proteins are outside its searched FASTA")
            gate.add(peptide)
            for protein in members:
                graph[protein].add(peptide)
        gates[sample.name] = gate
    dictionary = study_dictionary(graph, set(sequences))
    protein_counts = Counter(p for peptides in graph.values() for p in peptides)
    if not dictionary["covered"]:
        raise ValueError("No target peptide evidence passed the declared threshold")
    for peptide, group in dictionary["peptide_to_group"].items():
        if peptide not in sequences[group].replace("I", "L"):
            raise ValueError(f"Representative {group} does not contain assigned peptide {peptide}")

    # Build in a private temporary directory. Invalid input leaves no finished output.
    out.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=f".{out.name}-", dir=out.parent) as temporary:
        temp = Path(temporary)
        summary = _quantify(
            samples, dictionary, gates, searched_by_sample, temp, registry, acquisitions
        )
        _peptide_evidence(samples, dictionary, gates, protein_counts, temp, registry, q_threshold)
        summary.update(
            schema_version=1,
            acquisition_names=[sample.name for sample in samples],
            assignment_rule="direct canonical peptide to first selected explaining representative",
            grouping_rule="greedy uncovered peptide cover; lexical accession gain ties",
            inferred_study_groups=len(dictionary["selected"]),
            searched_universe=len(sequences),
            proteins_with_observed_evidence=len(graph),
            unobserved_searched_proteins=len(dictionary["unobserved"]),
            pooled_canonical_peptides=len(dictionary["covered"]),
            representative_sequence_compatibility_checked=True,
            acquisition_identity_checked=True,
            q_threshold=q_threshold,
            config_name=config_name,
            confidence="post-search reporting; no new protein-group or study-wide FDR estimate",
        )
        assigned = Counter(dictionary["peptide_to_group"].values())
        _write_table(
            temp / "study_groups.tsv",
            [
                "study_group",
                "selection_rank",
                "observed_support_peptides",
                "assigned_peptides",
            ],
            ((p, i, len(graph[p]), assigned[p]) for i, (p, _) in enumerate(dictionary["selected"])),
        )
        _write_table(
            temp / "peptide_to_study_group.tsv",
            ["canonical_peptide", "study_group"],
            sorted(dictionary["peptide_to_group"].items()),
        )
        # Export exactly the reporting representatives for optional annotation.
        # This does not collapse the retained peptide/member ambiguity ledgers.
        with (temp / "representatives.fasta").open("x") as stream:
            for protein, _ in dictionary["selected"]:
                stream.write(f">{protein}\n{sequences[protein]}\n")
        summary["annotation_targets"] = {
            "fasta": "representatives.fasta",
            "sha256": hashlib.sha256((temp / "representatives.fasta").read_bytes()).hexdigest(),
            "scope": "Selected reporting representatives; annotation does not prove unique origin",
        }
        module = Path(__file__)
        registry.append(
            {"path": str(module), "sha256": hashlib.sha256(module.read_bytes()).hexdigest()}
        )
        (temp / "input_manifest.json").write_text(json.dumps(registry, indent=2) + "\n")
        (temp / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
        # Exclusive creation also protects an output created concurrently after preflight.
        out.mkdir(exist_ok=False)
        for path in sorted(temp.iterdir(), key=lambda p: (p.name == "summary.json", p.name)):
            path.rename(out / path.name)
    return summary


def _quantify(samples, dictionary, gates, searched_by_sample, out, registry, acquisitions):
    """Assign positive finite LFQ features to the fixed peptide-to-group dictionary.

    Require accepted peptide evidence in the same acquisition. Write the feature
    ledger, group matrix and acquisition summaries, retaining missing cells as
    blanks. Return counts and intensity accounting; raise if mass balance fails.
    """
    matrix, long_rows, acquisition_rows = defaultdict(dict), [], []
    total_features, total_in = 0, 0.0
    off_local, off_members = 0, 0
    ledger_columns = [
        "sample",
        "input_line",
        "peptide",
        "canonical_peptide",
        "charge",
        "reported_proteins",
        "member_set",
        "study_group",
        "representative_in_local_fasta",
        "representative_in_local_members",
        "intensity",
    ]
    with open_gzip_text(out / "feature_assignments.tsv.gz") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(ledger_columns)
        for sample in samples:
            aggregates = defaultdict(lambda: [set(), 0, 0.0])
            sample_peptides, sample_features, sample_total = set(), 0, 0.0
            for line, row, fields in _rows(
                sample / "lfq.tsv",
                {"peptide", "proteins"},
                registry=registry,
                acquisition=acquisitions[sample.name],
                lfq=True,
            ):
                columns = [c for c in fields if c not in _LFQ_FIXED]
                if len(columns) != 1:
                    raise ValueError(f"{sample}: expected exactly one LFQ intensity column")
                peptide = canonical_peptide(row["peptide"])
                if peptide not in gates[sample.name]:
                    continue
                try:
                    intensity = float(row[columns[0]])
                except ValueError:
                    continue
                if not math.isfinite(intensity) or intensity <= 0:
                    continue
                members = _members(row["proteins"])
                if not members:
                    continue
                if members - searched_by_sample[sample.name]:
                    raise ValueError(
                        f"{sample.name}: quantified proteins are outside its searched FASTA"
                    )
                group = dictionary["peptide_to_group"].get(peptide)
                if group is None:
                    raise ValueError(f"{sample.name}: quantified peptide has no study group")
                in_fasta = group in searched_by_sample[sample.name]
                in_members = group in members
                off_local += not in_fasta
                off_members += not in_members
                writer.writerow(
                    [
                        sample.name,
                        line,
                        row["peptide"],
                        peptide,
                        row.get("charge", ""),
                        row["proteins"],
                        ";".join(sorted(members)),
                        group,
                        int(in_fasta),
                        int(in_members),
                        intensity,
                    ]
                )
                record = aggregates[group]
                record[0].add(peptide)
                record[1] += 1
                record[2] += intensity
                sample_peptides.add(peptide)
                sample_features += 1
                sample_total += intensity
            for group, (peptides, features, intensity) in sorted(aggregates.items()):
                long_rows.append((sample.name, group, len(peptides), features, intensity))
                matrix[group][sample.name] = intensity
            acquisition_rows.append(
                (sample.name, len(aggregates), len(sample_peptides), sample_features, sample_total)
            )
            total_features += sample_features
            total_in += sample_total
    total_out = sum(row[-1] for row in long_rows)
    if not math.isfinite(total_in) or not math.isfinite(total_out):
        raise ValueError("Non-finite aggregate intensity")
    if not math.isclose(total_in, total_out, rel_tol=1e-10, abs_tol=1e-6):
        raise RuntimeError("Intensity mass balance failed")
    _write_table(
        out / "study_group_long.tsv",
        [
            "sample",
            "study_group",
            "n_peptides",
            "n_features",
            "intensity",
        ],
        long_rows,
    )
    _write_table(
        out / "study_group_matrix.tsv",
        ["study_group"] + [p.name for p in samples],
        ((g, *(matrix[g].get(p.name, "") for p in samples)) for g in sorted(matrix)),
    )
    _write_table(
        out / "acquisition_summary.tsv",
        [
            "sample",
            "quantified_study_groups",
            "canonical_peptides",
            "lfq_features",
            "intensity",
        ],
        acquisition_rows,
    )
    return dict(
        acquisitions=len(samples),
        quantified_study_groups=len(matrix),
        lfq_features=total_features,
        intensity_in=total_in,
        intensity_out=total_out,
        features_representative_outside_local_fasta=off_local,
        features_representative_outside_local_members=off_members,
    )


def _peptide_evidence(samples, dictionary, gates, protein_counts, out, registry, q_threshold):
    """Retain accepted identifications, including peptides without positive LFQ features."""
    initial_hashes = {item["path"]: item["sha256"] for item in registry}
    with open_gzip_text(out / "peptide_evidence.tsv.gz") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(
            [
                "sample",
                "study_group",
                "canonical_peptide",
                "reported_proteins",
                "n_reported_proteins_sample",
                "n_reported_proteins_study",
            ]
        )
        for sample in samples:
            local = defaultdict(set)
            path = sample / "results.sage.tsv"
            for _, row, _ in _rows(
                path, {"peptide", "proteins", "label", "peptide_q"}, registry=registry
            ):
                if row["label"] == "1" and _passing(row["peptide_q"], q_threshold):
                    peptide = canonical_peptide(row["peptide"])
                    if peptide:
                        local[peptide].update(_members(row["proteins"]))
            if registry[-1]["sha256"] != initial_hashes[str(path)]:
                raise ValueError(f"Identification input changed during grouping: {path}")
            for peptide in sorted(gates[sample.name]):
                group = dictionary["peptide_to_group"].get(peptide)
                if group is not None:
                    writer.writerow(
                        [
                            sample.name,
                            group,
                            peptide,
                            ";".join(sorted(local[peptide])),
                            len(local[peptide]),
                            protein_counts[peptide],
                        ]
                    )
    _write_table(
        out / "peptide_ambiguity.tsv",
        ["n_reported_proteins_study", "n_canonical_peptides"],
        sorted(Counter(protein_counts.values()).items()),
    )
