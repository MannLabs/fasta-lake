"""Join eggNOG functions to a peptide/protein graph using verified exact sequences.

Functional edges describe orthology-based annotations. They do not turn peptide
membership into a unique protein identification, measured activity or flux.


Reads native eggNOG-mapper output: https://github.com/eggnogdb/eggnog-mapper.
Cantalapiedra et al. (2021), https://doi.org/10.1093/molbev/msab293.
FastaLake adds exact-sequence joins and explicit annotation coverage; the
orthology assignments are produced by eggNOG-mapper.
"""

from __future__ import annotations

import hashlib
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

from fasta_lake.annotation import MISSING, read_eggnog_annotations
from fasta_lake.multi_omics.exact import read_fasta, write_table
from fasta_lake.networks import (
    _identity,
    _rows,
    _sha,
    load_evidence_network,
)  # noqa: F401

FUNCTION_FIELDS = (
    "KEGG_ko",
    "EC",
    "GOs",
    "KEGG_Pathway",
    "KEGG_Module",
    "KEGG_Reaction",
    "CAZy",
    "PFAMs",
    "COG_category",
)


def load_eggnog_annotations(annotation, *, query_fasta=None, max_records=100_000):
    """Read a completed FastaLake annotation folder or native eggNOG v2 output.

    Native output requires the original query FASTA supplied by the user. A
    FastaLake folder supplies hashed queries, aliases and its completion record.
    Source/database metadata is retained; annotation databases are not reopened.
    """
    source = Path(annotation).resolve(strict=True)
    aliases = {}
    if source.is_dir():
        if query_fasta is not None:
            raise ValueError("A FastaLake annotation folder already supplies its query FASTA")
        receipt = json.loads((source / "COMPLETE.json").read_text())
        required = {
            "queries.fasta",
            "aliases.tsv",
            "annotation.emapper.annotations",
            "annotations.tsv",
            "PROVENANCE.json",
        }
        if receipt.get("status") != "PASS" or not required <= set(receipt.get("files", {})):
            raise ValueError("Annotation folder has no complete supported annotation receipt")
        for name, expected in receipt["files"].items():
            if Path(name).name != name or _sha(source / name) != expected:
                raise ValueError(f"Annotation checksum mismatch: {name}")
        provenance = json.loads((source / "PROVENANCE.json").read_text())
        query_path = source / "queries.fasta"
        table = source / "annotation.emapper.annotations"
        paths = [
            source / "COMPLETE.json",
            query_path,
            source / "aliases.tsv",
            table,
            source / "PROVENANCE.json",
        ]
        unknown_fields = tuple(provenance["unknown_fields"])
        linkage = "FastaLake sequence-hashed queries and checked annotation completion"
    else:
        if query_fasta is None:
            raise ValueError("Native eggNOG output requires --query-fasta from its annotation run")
        query_path = Path(query_fasta).resolve(strict=True)
        table = source
        paths = [table, query_path]
        unknown_fields = ("KEGG_ko", "EC")
        provenance = {"query_reference_linkage": "Original annotation query FASTA supplied by user"}
        linkage = "User-declared native query FASTA; identifiers and sequences checked locally"
    hashes = {str(path): _sha(path) for path in paths}
    queries = read_fasta(query_path, max_records=max_records)
    if source.is_dir():
        if provenance.get("query_sha256") != _sha(query_path):
            raise ValueError("Annotation provenance disagrees with its query FASTA")
        fields = ["accession", "query", "protein_sha256"]
        for row in _rows(source / "aliases.tsv", fields, max_records):
            accession = _identity(row["accession"], "annotation alias")
            query = row["query"]
            if accession in aliases or query not in queries:
                raise ValueError("Duplicate alias or foreign annotation query")
            if row["protein_sha256"] != queries[query][2] or query != "FLA_" + queries[query][2]:
                raise ValueError("Annotation alias disagrees with its exact query sequence")
            aliases[accession] = query
        if set(aliases.values()) != set(queries):
            raise ValueError("Annotation aliases do not cover the complete query roster")
    else:
        aliases = {query: query for query in queries}
    rows, _ = read_eggnog_annotations(table, queries, unknown_fields=unknown_fields)
    if hashes != {str(path): _sha(path) for path in paths}:
        raise ValueError("Annotation inputs changed while being loaded")
    return dict(
        queries=queries,
        aliases=aliases,
        rows={r["query"]: r for r in rows},
        source_files=hashes,
        provenance=provenance,
        linkage=linkage,
    )


def function_terms(value, field):
    """Split one annotation namespace without creating terms for missing values."""
    if field not in FUNCTION_FIELDS:
        raise ValueError(f"Unsupported functional namespace: {field}")
    if value.strip() in MISSING:
        return []
    terms = list(value.strip()) if field == "COG_category" else value.split(",")
    result = set()
    for term in terms:
        term = term.strip()
        if term in MISSING:
            continue
        if field == "KEGG_ko":
            term = term.removeprefix("ko:")
            if not re.fullmatch(r"K[0-9]{5}", term):
                raise ValueError(f"Invalid KEGG ortholog: {term}")
        elif field == "EC" and not re.fullmatch(r"[0-9]+(?:\.(?:[0-9]+|-)){3}", term):
            raise ValueError(f"Invalid EC number: {term}")
        elif field == "GOs" and not re.fullmatch(r"GO:[0-9]{7}", term):
            raise ValueError(f"Invalid GO identifier: {term}")
        _identity(term, "functional term")
        result.add(term)
    return sorted(result)


def annotate_evidence_network(
    bundle,
    annotation,
    fasta,
    output,
    *,
    query_fasta=None,
    fields=("KEGG_ko", "EC", "GOs", "PFAMs"),
    max_records=100_000,
):
    """Attach annotations to exact sequences from a recorded study FASTA.

    An accession match must also agree in sequence. Alternative names can match
    by full protein SHA-256; I/L variants remain distinct. Conflicting annotations
    for duplicate query sequences fail instead of selecting an arbitrary row.
    Only annotated reporting representatives enter the group-to-term exports.
    """
    out = Path(output)
    if out.exists() or out.is_symlink():
        raise FileExistsError(f"Annotation output exists: {out}; choose a new directory")
    fields = tuple(fields)
    if not fields or len(set(fields)) != len(fields) or not set(fields) <= set(FUNCTION_FIELDS):
        raise ValueError("Choose distinct supported annotation namespaces")
    network = load_evidence_network(bundle)
    fasta = Path(fasta).resolve(strict=True)
    fasta_hash = _sha(fasta)
    expected = {network["source_summary"].get("annotation_targets", {}).get("sha256")}
    manifest = network.get("input_manifest", [])
    if isinstance(manifest, list):
        expected.update(
            row.get("sha256")
            for row in manifest
            if isinstance(row, dict)
            and str(row.get("path", "")).endswith(
                (".fasta", ".fa", ".faa", ".fasta.gz", ".fa.gz", ".faa.gz")
            )
        )
    if fasta_hash not in expected:
        raise ValueError("Supply the study's representatives.fasta or a recorded search FASTA")
    proteins = read_fasta(fasta, max_records=max_records)
    loaded = load_eggnog_annotations(annotation, query_fasta=query_fasta, max_records=max_records)
    by_hash = defaultdict(list)
    for query, (_, _, digest) in loaded["queries"].items():
        by_hash[digest].append(query)
    nodes, edges, rows, used_queries = {}, [], [], set()
    maps = {field: [] for field in fields}
    for node in network["nodes"]:
        if node["node_type"] != "reported_protein":
            continue
        accession = node["label"]
        row = dict(
            protein_id=accession,
            node_id=node["id"],
            protein_sha256="",
            query="",
            match_rule="",
            annotation_status="sequence_not_supplied",
            selected_representative=node["selected_representative"],
            raw_annotation={},
        )
        if accession in proteins:
            digest = proteins[accession][2]
            row["protein_sha256"] = digest
            query = loaded["aliases"].get(accession)
            if query is not None:
                if loaded["queries"][query][2] != digest:
                    raise ValueError(
                        f"Annotation reuses accession for a different sequence: {accession}"
                    )
                row["match_rule"] = "accession_and_exact_sequence"
            else:
                candidates = sorted(by_hash.get(digest, []))
                signatures = {
                    json.dumps(
                        {k: v for k, v in loaded["rows"][q].items() if k != "query"}, sort_keys=True
                    )
                    for q in candidates
                }
                if len(signatures) > 1:
                    raise ValueError("Duplicate query sequences carry conflicting annotations")
                query = candidates[0] if candidates else None
                row["match_rule"] = "exact_sequence_alias" if query else ""
            row["annotation_status"] = "not_annotated"
            if query is not None:
                used_queries.add(query)
                values = loaded["rows"][query]
                if not set(fields) <= set(values):
                    raise ValueError("eggNOG output lacks a requested annotation namespace")
                row.update(
                    query=query,
                    annotation_status=values["annotation_status"],
                    raw_annotation=dict(values),
                )
                for field in fields:
                    for term in function_terms(values[field], field):
                        key = (
                            "function:" + hashlib.sha256((field + "\0" + term).encode()).hexdigest()
                        )
                        nodes[key] = dict(
                            id=key, node_type="annotated_function", namespace=field, label=term
                        )
                        edges.append(
                            dict(
                                source=node["id"],
                                target=key,
                                relation="orthology_annotation",
                                namespace=field,
                                query=query,
                                protein_sha256=digest,
                            )
                        )
                        if node["selected_representative"]:
                            maps[field].append({"study_group": accession, "term": term})
        rows.append(row)
    if not used_queries:
        raise ValueError(
            "No annotation query matches a graph protein with a supplied exact sequence"
        )
    if _sha(fasta) != fasta_hash or any(_sha(p) != h for p, h in loaded["source_files"].items()):
        raise ValueError("Functional annotation inputs changed during export")
    status_counts = dict(Counter(row["annotation_status"] for row in rows))
    scope = (
        "Orthology-based protein annotations; peptide memberships and reporting assignments "
        "remain unchanged. No peptide-to-function inference, activity or flux is asserted."
    )
    data = dict(
        schema_version=1,
        graph_id=network["graph_id"],
        proteins=rows,
        nodes=list(nodes.values()),
        edges=edges,
        status_counts=status_counts,
        source_fasta={"path": str(fasta), "sha256": fasta_hash},
        annotation_source_files=loaded["source_files"],
        annotation_provenance=loaded["provenance"],
        query_linkage=loaded["linkage"],
        unused_annotation_queries=sorted(set(loaded["queries"]) - used_queries),
        scope=scope,
    )
    out.mkdir(parents=True, exist_ok=False)
    (out / "annotations.json").write_text(json.dumps(data, indent=2) + "\n")
    columns = [
        "protein_id",
        "node_id",
        "protein_sha256",
        "query",
        "match_rule",
        "annotation_status",
        "selected_representative",
        *fields,
    ]
    write_table(
        out / "protein_annotations.tsv",
        columns,
        [
            {
                **{k: row[k] for k in columns if k not in fields},
                **{k: row["raw_annotation"].get(k, "") for k in fields},
            }
            for row in rows
        ],
    )
    write_table(
        out / "function_nodes.tsv", ["id", "node_type", "namespace", "label"], data["nodes"]
    )
    write_table(
        out / "function_edges.tsv",
        ["source", "target", "relation", "namespace", "query", "protein_sha256"],
        edges,
    )
    for field, pairs in maps.items():
        write_table(out / f"group_terms_{field}.tsv", ["study_group", "term"], pairs)
    files = {p.name: _sha(p) for p in out.iterdir() if p.is_file()}
    receipt = dict(
        status="PASS",
        schema_version=1,
        graph_id=network["graph_id"],
        protein_nodes=len(rows),
        function_nodes=len(nodes),
        annotation_edges=len(edges),
        status_counts=status_counts,
        files=files,
        scope=scope,
    )
    (out / "COMPLETE.json").write_text(json.dumps(receipt, indent=2) + "\n")
    return receipt
