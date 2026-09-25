"""Export accepted peptide memberships without conflating them with inference.

Nodes represent canonical peptides and the protein accessions reported by the
study's searches. They are not reference-wide sequence-specific identifications.
The assigned reporting representative is a separate attribute of each edge.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path


def _sha(path):
    """Return a streamed SHA-256 digest of the evidence file."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _rows(path, fields, limit):
    """Yield exact-schema TSV rows, enforcing the caller's row limit."""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if reader.fieldnames != fields:
            raise ValueError(f"{path}: expected columns {fields}")
        for number, row in enumerate(reader, 1):
            if number > limit:
                raise ValueError(f"{path}: exceeds {limit} rows; raise the explicit row limit")
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f"{path}: malformed row {number}")
            yield row


def _integer(value, label, minimum=0):
    """Parse an integer at least minimum, rejecting noncanonical numeric text."""
    try:
        parsed = int(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"Invalid integer for {label}: {value!r}") from error
    if str(parsed) != str(value) or parsed < minimum:
        raise ValueError(f"Invalid integer for {label}: {value!r}")
    return parsed


def _identity(value, label):
    """Validate a nonblank identifier without padding or control characters."""
    if (
        not isinstance(value, str)
        or not value
        or value != value.strip()
        or any(ord(c) < 32 for c in value)
    ):
        raise ValueError(f"Invalid {label}: {value!r}")
    return value


def read_evidence_network(groups, *, max_rows=1_000_000, max_edges=1_000_000):
    """Read and reconcile a schema-1 study's complete reported-membership graph.

    Each sample/peptide row must occur once. Edge counts refer to distinct
    acquisitions reporting that member, not PSM counts, quantities or unique
    proteins. Row/edge limits fail explicitly rather than truncate evidence.
    """
    max_rows = _integer(max_rows, "max_rows", 1)
    max_edges = _integer(max_edges, "max_edges", 1)
    root = Path(groups).resolve(strict=True)
    names = [
        "summary.json",
        "study_groups.tsv",
        "peptide_to_study_group.tsv",
        "peptide_evidence.tsv.gz",
        "input_manifest.json",
    ]
    identities = {name: _sha(root / name) for name in names}
    summary = json.loads((root / "summary.json").read_text())
    source_manifest = json.loads((root / "input_manifest.json").read_text())
    if summary.get("schema_version") != 1:
        raise ValueError("Network export requires a schema-1 study summary")
    samples = summary.get("acquisition_names")
    if not isinstance(samples, list) or any(not isinstance(s, str) for s in samples):
        raise ValueError("Study acquisition roster must be a list of distinct names")
    if len(samples) != len(set(samples)):
        raise ValueError("Study acquisition roster must be a list of distinct names")
    for sample in samples:
        _identity(sample, "acquisition name")
    if summary.get("acquisitions") != len(samples):
        raise ValueError("Study acquisition count disagrees with its roster")
    sample_set = set(samples)
    representatives = {}
    fields = ["study_group", "selection_rank", "observed_support_peptides", "assigned_peptides"]
    for row in _rows(root / "study_groups.tsv", fields, max_rows):
        group = _identity(row["study_group"], "study group")
        if group in representatives:
            raise ValueError(f"Duplicate study group: {group}")
        representatives[group] = {key: _integer(row[key], key) for key in fields[1:]}
    if summary.get("inferred_study_groups") != len(representatives):
        raise ValueError("Study group count disagrees with the dictionary")
    if sorted(r["selection_rank"] for r in representatives.values()) != list(
        range(len(representatives))
    ):
        raise ValueError("Study selection ranks must cover the dictionary exactly")
    assignment = {}
    for row in _rows(
        root / "peptide_to_study_group.tsv", ["canonical_peptide", "study_group"], max_rows
    ):
        peptide = _identity(row["canonical_peptide"], "canonical peptide")
        if (
            not peptide.isascii()
            or not peptide.isalpha()
            or not peptide.isupper()
            or "I" in peptide
        ):
            raise ValueError(
                f"Peptide is not canonical uppercase I/L-normalized sequence: {peptide}"
            )
        if peptide in assignment or row["study_group"] not in representatives:
            raise ValueError("Duplicate peptide or foreign group in the assignment dictionary")
        assignment[peptide] = row["study_group"]
    assigned_counts = Counter(assignment.values())
    if summary.get("pooled_canonical_peptides") != len(assignment):
        raise ValueError("Peptide count disagrees with the assignment dictionary")
    if any(assigned_counts[g] != r["assigned_peptides"] for g, r in representatives.items()):
        raise ValueError("Assigned peptide counts disagree with the study dictionary")

    memberships = defaultdict(set)
    peptide_samples = defaultdict(set)
    protein_samples = defaultdict(set)
    sample_peptides = defaultdict(set)
    sample_proteins = defaultdict(set)
    declared_degree = {}
    fields = [
        "sample",
        "study_group",
        "canonical_peptide",
        "reported_proteins",
        "n_reported_proteins_sample",
        "n_reported_proteins_study",
    ]
    for row in _rows(root / "peptide_evidence.tsv.gz", fields, max_rows):
        sample, peptide = row["sample"], row["canonical_peptide"]
        if sample not in sample_set or peptide not in assignment:
            raise ValueError("Peptide evidence contains a foreign acquisition or peptide")
        if row["study_group"] != assignment[peptide]:
            raise ValueError("Peptide evidence disagrees with the reporting assignment")
        if sample in peptide_samples[peptide]:
            raise ValueError(f"Duplicate acquisition/peptide evidence: {sample}, {peptide}")
        proteins = row["reported_proteins"].split(";")
        for protein in proteins:
            _identity(protein, "reported protein accession")
        if len(proteins) != len(set(proteins)) or len(proteins) != _integer(
            row["n_reported_proteins_sample"], "sample membership", 1
        ):
            raise ValueError("Reported sample membership count is inconsistent")
        degree = _integer(row["n_reported_proteins_study"], "study membership", 1)
        if declared_degree.setdefault(peptide, degree) != degree:
            raise ValueError("Study membership count changes between acquisitions")
        peptide_samples[peptide].add(sample)
        sample_peptides[sample].add(peptide)
        for protein in proteins:
            key = (peptide, protein)
            if key not in memberships and len(memberships) >= max_edges:
                raise ValueError(
                    f"Network exceeds {max_edges} edges; raise the explicit edge limit"
                )
            memberships[key].add(sample)
            protein_samples[protein].add(sample)
            sample_proteins[sample].add(protein)
    observed_degree = Counter(peptide for peptide, _ in memberships)
    support_degree = Counter(protein for _, protein in memberships)
    if set(peptide_samples) != set(assignment):
        raise ValueError("Peptide evidence does not cover the complete assignment dictionary")
    for peptide, group in assignment.items():
        if observed_degree[peptide] != declared_degree[peptide]:
            raise ValueError("Study membership total disagrees with the evidence rows")
        if (peptide, group) not in memberships:
            raise ValueError("Reporting representative lacks the declared study peptide membership")
    if any(support_degree[g] != r["observed_support_peptides"] for g, r in representatives.items()):
        raise ValueError("Observed support disagrees with the study dictionary")
    if summary.get("proteins_with_observed_evidence") != len(protein_samples):
        raise ValueError("Observed protein count disagrees with the evidence graph")
    if identities != {name: _sha(root / name) for name in names}:
        raise ValueError("Study inputs changed during network export")
    graph_id = hashlib.sha256(json.dumps(identities, sort_keys=True).encode()).hexdigest()

    def node_id(kind, label):
        """Derive a kind-prefixed node ID from the graph identity and node label."""
        return kind + ":" + hashlib.sha256((graph_id + "\0" + label).encode()).hexdigest()

    nodes = []
    for peptide in sorted(assignment):
        nodes.append(
            dict(
                id=node_id("peptide", peptide),
                node_type="peptide",
                label=peptide,
                reported_degree=observed_degree[peptide],
                observed_acquisitions=len(peptide_samples[peptide]),
                selected_representative=False,
                assigned_group=assignment[peptide],
            )
        )
    for protein in sorted(protein_samples):
        nodes.append(
            dict(
                id=node_id("protein", protein),
                node_type="reported_protein",
                label=protein,
                reported_degree=support_degree[protein],
                observed_acquisitions=len(protein_samples[protein]),
                selected_representative=protein in representatives,
                assigned_group="",
            )
        )
    edges = [
        dict(
            source=node_id("peptide", peptide),
            target=node_id("protein", protein),
            relation="reported_membership",
            observed_acquisitions=len(observed),
            assigned_representative=assignment[peptide] == protein,
            acquisitions=sorted(observed),
        )
        for (peptide, protein), observed in sorted(memberships.items())
    ]
    return dict(
        schema_version=1,
        graph_id=graph_id,
        source_files=identities,
        source_summary=summary,
        input_manifest=source_manifest,
        nodes=nodes,
        edges=edges,
        acquisitions=[
            dict(
                sample=sample,
                canonical_peptides=len(sample_peptides[sample]),
                reported_proteins=len(sample_proteins[sample]),
            )
            for sample in samples
        ],
        scope="Sage-reported peptide/protein memberships and separate study assignments; "
        "not reference-wide specificity, protein inference or physical interactions",
    )


def export_evidence_network(groups, output, *, max_rows=1_000_000, max_edges=1_000_000):
    """Write complete typed tables and a checked portable graph JSON in a new folder."""
    out = Path(output)
    if out.exists() or out.is_symlink():
        raise FileExistsError(f"Network output exists: {out}; choose a new directory")
    network = read_evidence_network(groups, max_rows=max_rows, max_edges=max_edges)
    out.mkdir(parents=True, exist_ok=False)
    for name, fields in [
        (
            "nodes",
            [
                "id",
                "node_type",
                "label",
                "reported_degree",
                "observed_acquisitions",
                "selected_representative",
                "assigned_group",
            ],
        ),
        (
            "edges",
            [
                "source",
                "target",
                "relation",
                "observed_acquisitions",
                "assigned_representative",
                "acquisitions",
            ],
        ),
        ("acquisitions", ["sample", "canonical_peptides", "reported_proteins"]),
    ]:
        with (out / (name + ".tsv")).open("w", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for original in network[name]:
                row = dict(original)
                if "acquisitions" in row:
                    row["acquisitions"] = json.dumps(row["acquisitions"], ensure_ascii=True)
                writer.writerow(row)
    (out / "network.json").write_text(json.dumps(network, indent=2, ensure_ascii=True) + "\n")
    files = ["network.json", "nodes.tsv", "edges.tsv", "acquisitions.tsv"]
    receipt = dict(
        status="PASS",
        schema_version=1,
        graph_id=network["graph_id"],
        nodes=len(network["nodes"]),
        edges=len(network["edges"]),
        scope=network["scope"],
        files={name: _sha(out / name) for name in files},
    )
    (out / "COMPLETE.json").write_text(json.dumps(receipt, indent=2) + "\n")
    return receipt


def load_evidence_network(bundle):
    """Read an unchanged, completed network export; reject partial or modified bundles."""
    root = Path(bundle).resolve(strict=True)
    receipt = json.loads((root / "COMPLETE.json").read_text())
    files = {"network.json", "nodes.tsv", "edges.tsv", "acquisitions.tsv"}
    if receipt.get("status") != "PASS" or receipt.get("schema_version") != 1:
        raise ValueError("Network bundle has no successful supported completion receipt")
    if set(receipt.get("files", {})) != files:
        raise ValueError("Network bundle manifest must name exactly its four data files")
    for name in files:
        if _sha(root / name) != receipt["files"][name]:
            raise ValueError(f"Network bundle checksum mismatch: {name}")
    data = json.loads((root / "network.json").read_text())
    if data.get("schema_version") != 1 or data.get("graph_id") != receipt.get("graph_id"):
        raise ValueError("Network identity disagrees with its completion receipt")
    if len(data["nodes"]) != receipt["nodes"] or len(data["edges"]) != receipt["edges"]:
        raise ValueError("Network counts disagree with its completion receipt")
    return data


def select_network_view(network, *, max_nodes=150, max_edges=300, preferred_proteins=()):
    """Select a deterministic bounded view, prioritizing shared peptide memberships.

    Every selected edge retains both endpoints. Limits apply only to the view;
    the complete export is never truncated. Reporting-assignment edges are
    visited first within each peptide's neighborhood.
    """
    max_nodes = _integer(max_nodes, "max_nodes", 2)
    max_edges = _integer(max_edges, "max_edges", 1)
    nodes = {node["id"]: node for node in network["nodes"]}
    preferred = set(preferred_proteins)
    if not preferred <= {
        key for key, node in nodes.items() if node["node_type"] == "reported_protein"
    }:
        raise ValueError("Preferred annotation proteins are not in the evidence graph")
    if len(nodes) != len(network["nodes"]):
        raise ValueError("Duplicate graph node identities")
    ordered = sorted(
        network["edges"],
        key=lambda edge: (
            edge["target"] not in preferred,
            -nodes[edge["source"]]["reported_degree"],
            nodes[edge["source"]]["label"],
            not edge["assigned_representative"],
            nodes[edge["target"]]["label"],
        ),
    )
    included = set()
    selected = []
    for edge in ordered:
        if len(selected) == max_edges:
            break
        endpoints = {edge["source"], edge["target"]}
        if len(included | endpoints) <= max_nodes:
            included.update(endpoints)
            selected.append(dict(edge))
    return dict(
        nodes=[dict(nodes[key]) for key in sorted(included)],
        edges=selected,
        total_nodes=len(nodes),
        total_edges=len(ordered),
        hidden_nodes=len(nodes) - len(included),
        hidden_edges=len(ordered) - len(selected),
        selection=("Annotated protein neighborhoods first; then " if preferred else "")
        + ("descending peptide membership degree, peptide label, assignment first, protein label"),
    )


def load_network_annotations(folder, graph_id):
    """Read an unchanged functional layer belonging to the requested evidence graph."""
    root = Path(folder).resolve(strict=True)
    receipt = json.loads((root / "COMPLETE.json").read_text())
    if receipt.get("status") != "PASS" or receipt.get("schema_version") != 1:
        raise ValueError("Functional layer has no completed supported receipt")
    if receipt.get("graph_id") != graph_id or "annotations.json" not in receipt.get("files", {}):
        raise ValueError("Functional annotations belong to another evidence graph")
    for name, expected in receipt["files"].items():
        if Path(name).name != name or _sha(root / name) != expected:
            raise ValueError(f"Functional annotation checksum mismatch: {name}")
    data = json.loads((root / "annotations.json").read_text())
    if data.get("graph_id") != graph_id:
        raise ValueError("Functional annotation graph identity mismatch")
    return data


def add_function_annotations(view, layer, *, max_nodes, max_edges, fields=("KEGG_ko", "EC")):
    """Extend a bounded view with typed orthology edges without changing assignments.

    Full exports retain all annotations. Display limits apply to the combined
    graph, and favour KO/EC over the larger GO/domain vocabularies.
    """
    view = {
        **view,
        "nodes": [dict(n) for n in view["nodes"]],
        "edges": [dict(e) for e in view["edges"]],
    }
    ids = {n["id"] for n in view["nodes"]}
    proteins = {r["node_id"]: r for r in layer["proteins"]}
    functions = {n["id"]: n for n in layer["nodes"]}
    for node in view["nodes"]:
        if node["id"] in proteins:
            row = proteins[node["id"]]
            node["annotation_status"] = row["annotation_status"]
            node["annotation_query"] = row["query"]
    priority = {"KEGG_ko": 0, "EC": 1, "COG_category": 2, "PFAMs": 3}
    for edge in sorted(
        layer["edges"],
        key=lambda e: (
            priority.get(e["namespace"], 4),
            functions[e["target"]]["label"],
            e["source"],
        ),
    ):
        if edge["namespace"] not in fields:
            continue
        if edge["source"] not in ids or len(view["edges"]) >= max_edges:
            continue
        if edge["target"] not in ids:
            if len(ids) >= max_nodes:
                continue
            view["nodes"].append(dict(functions[edge["target"]]))
            ids.add(edge["target"])
        view["edges"].append(dict(edge))
    view["annotation_scope"] = layer["scope"]
    view["visible_annotation_namespaces"] = list(fields)
    view["total_nodes"] += len(layer["nodes"])
    view["total_edges"] += len(layer["edges"])
    view["hidden_nodes"] = view["total_nodes"] - len(view["nodes"])
    view["hidden_edges"] = view["total_edges"] - len(view["edges"])
    return view


def render_ckg_network(
    bundle,
    output,
    *,
    ckg_python,
    ckg_source=None,
    max_nodes=150,
    max_edges=300,
    annotations=None,
):
    """Render a portable export with the user's separate compatible CKG environment.

    The complete graph stays in the export; only the interactive/static view is
    bounded. The worker reads files locally and does not connect to Neo4j.
    """
    import subprocess

    network = load_evidence_network(bundle)
    max_nodes = _integer(max_nodes, "max_nodes", 2)
    max_edges = _integer(max_edges, "max_edges", 1)
    destination = Path(output).resolve()
    if destination.exists() or Path(output).is_symlink():
        raise FileExistsError(f"CKG output exists: {output}; choose a new directory")
    from fasta_lake.resources import local_interpreter

    python = local_interpreter(ckg_python)
    worker = Path(__file__).with_name("_ckg_worker.py")
    command = [
        str(python),
        str(worker),
        "--bundle",
        str(Path(bundle).resolve(strict=True)),
        "--output",
        str(destination),
        "--max-nodes",
        str(max_nodes),
        "--max-edges",
        str(max_edges),
    ]
    if ckg_source is not None:
        command.extend(["--ckg-source", str(Path(ckg_source).resolve(strict=True))])
    if annotations is not None:
        load_network_annotations(annotations, network["graph_id"])
        command.extend(["--annotations", str(Path(annotations).resolve(strict=True))])
    process = subprocess.run(command, capture_output=True, text=True)
    destination.mkdir(parents=True, exist_ok=True)
    (destination / "worker.log").write_text(process.stdout + process.stderr)
    (destination / "PROCESS.json").write_text(
        json.dumps({"command": command, "exit_code": process.returncode}, indent=2) + "\n"
    )
    if process.returncode:
        raise RuntimeError(f"CKG rendering failed; inspect {destination / 'worker.log'}")
    receipt = json.loads((destination / "RENDER_COMPLETE.json").read_text())
    if receipt.get("status") != "PASS" or receipt.get("graph_id") != network["graph_id"]:
        raise ValueError("CKG renderer did not return a matching completed graph")
    for name, expected in receipt["files"].items():
        if Path(name).name != name or _sha(destination / name) != expected:
            raise ValueError(f"CKG output identity mismatch: {name}")
    return receipt
