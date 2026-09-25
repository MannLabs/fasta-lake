"""Optional exact-sequence annotation and explicitly defined unknown proteins.

eggNOG assigns annotations through orthology. An absent KO/EC assignment is an
annotation gap, not proof of a novel protein. ESM-C embeddings are separate
representations and never overwrite functional assignments.


Upstream: https://github.com/eggnogdb/eggnog-mapper.
Cantalapiedra et al. (2021), https://doi.org/10.1093/molbev/msab293.
FastaLake invokes the mapper and checks exact input/output identity; it does
not implement orthology assignment.
"""

from __future__ import annotations

import json
import re
import subprocess
from pathlib import Path

from fasta_lake.multi_omics.exact import read_fasta, sha256, write_table
from fasta_lake.resources import local_interpreter, run_recorded

MISSING = {"", "-", "NA", "N/A", "nan", "None"}


def prepare_annotation_targets(fasta, output, *, max_proteins=100000):
    """Deduplicate a selected FASTA by exact sequence; retain every input alias.

    Return query ID -> sequence, plus the original FASTA hash. Query IDs contain
    the full sequence SHA-256, so annotation cannot depend on accession spelling.
    Protein I/L variants remain distinct. This function expects a new directory;
    the count guard prevents accidentally annotating a worldwide catalogue.
    """
    before = sha256(fasta)
    records = read_fasta(fasta, max_records=max_proteins)
    if any(k.startswith("rev_") for k in records):
        raise ValueError("Annotation input must contain targets only, without rev_ decoys")
    if type(max_proteins) is not int or max_proteins < 1:
        raise ValueError("max_proteins must be a positive integer")
    sequences = {"FLA_" + row[2]: row[1] for row in records.values()}
    if not sequences or len(sequences) > max_proteins:
        raise ValueError("Select a nonempty protein subset within --max-proteins before annotation")
    if sha256(fasta) != before:
        raise ValueError("Annotation FASTA changed while being read")
    out = Path(output)
    out.mkdir(parents=True, exist_ok=False)
    with (out / "queries.fasta").open("x") as stream:
        for query, sequence in sorted(sequences.items()):
            stream.write(f">{query}\n{sequence}\n")
    write_table(
        out / "aliases.tsv",
        ["accession", "query", "protein_sha256"],
        [
            {"accession": accession, "query": "FLA_" + row[2], "protein_sha256": row[2]}
            for accession, row in sorted(records.items())
        ],
    )
    return sequences, before


def read_eggnog_annotations(path, queries, *, unknown_fields=("KEGG_ko", "EC")):
    """Read v2 eggNOG output and retain no-hit and hit-without-function queries.

    Unknown means no value in any declared field (KO or EC by default). Seed
    hits, GO terms, descriptions and all raw eggNOG fields remain available.
    Unknown query IDs, duplicate rows and malformed tables fail explicitly.
    Return normalized rows in query order and the exact unknown query roster.
    """
    if not unknown_fields or len(set(unknown_fields)) != len(unknown_fields):
        raise ValueError("Choose one or more distinct annotation fields for unknown selection")
    header, annotated = None, {}
    with Path(path).open() as stream:
        for line in stream:
            if line.startswith("#query\t"):
                if header is not None:
                    raise ValueError("Duplicate eggNOG annotation header")
                header = line[1:].rstrip("\r\n").split("\t")
                if len(header) != len(set(header)) or not set(unknown_fields) <= set(header):
                    raise ValueError("eggNOG header lacks distinct requested annotation fields")
                continue
            if line.startswith("#") or not line.strip():
                continue
            if header is None:
                raise ValueError("eggNOG data precedes its #query header")
            values = line.rstrip("\r\n").split("\t")
            if len(values) != len(header):
                raise ValueError("Malformed eggNOG row width")
            row = dict(zip(header, values))
            query = row["query"]
            if query not in queries or query in annotated:
                raise ValueError(f"Foreign or duplicate eggNOG query: {query}")
            annotated[query] = row
    if header is None:
        raise ValueError("Missing eggNOG #query header")
    rows, unknown = [], []
    for query in sorted(queries):
        row = {field: "" for field in header}
        row.update(annotated.get(query, {}))
        row["query"] = query
        assigned = any(row[field].strip() not in MISSING for field in unknown_fields)
        row["annotation_status"] = (
            "assigned"
            if assigned
            else "hit_without_requested_terms"
            if query in annotated
            else "no_hit"
        )
        row["unknown_under_declared_rule"] = not assigned
        rows.append(row)
        if not assigned:
            unknown.append(query)
    return rows, unknown


def annotate_fasta(
    fasta,
    output,
    *,
    emapper,
    data_dir,
    threads=2,
    emapper_python=None,
    sensitivity="sensitive",
    tax_scope="auto",
    max_proteins=100000,
    unknown_fields=("KEGG_ko", "EC"),
):
    """Run eggNOG-mapper 2.1.12 on exact unique selected sequences.

    The annotation database stays on disk (no --dbmem). DIAMOND uses block size
    0.5; threads are explicit. Automatic taxonomic scope retains mixed-domain
    inputs instead of silently restricting the whole microbiome to bacteria.
    Database files, tool, queries and results are hashed. A failed command keeps
    its logs and never writes COMPLETE.json. External databases are not downloaded.
    """
    if type(threads) is not int or threads < 1:
        raise ValueError("Annotation threads must be a positive integer")
    if sensitivity not in {"sensitive", "more-sensitive", "very-sensitive"}:
        raise ValueError("Choose sensitive, more-sensitive or very-sensitive")
    if not re.fullmatch(r"auto|[1-9][0-9]*(?:,[1-9][0-9]*)*", tax_scope):
        raise ValueError("Taxonomic scope must be auto or a comma-separated list of taxon IDs")
    out = Path(output)
    if out.exists() or out.is_symlink():
        raise FileExistsError(f"Output exists: {out}; choose a new directory")
    emapper, data = Path(emapper).resolve(strict=True), Path(data_dir).resolve(strict=True)
    command = [str(local_interpreter(emapper_python))] if emapper_python else []
    command.append(str(emapper))
    try:
        version = subprocess.check_output(
            command + ["--data_dir", str(data), "--version"], text=True, stderr=subprocess.STDOUT
        )
    except subprocess.CalledProcessError as error:
        raise RuntimeError(
            "eggNOG version check failed: " + (error.output or "")[-1500:]
        ) from error
    if not re.search(r"\b2\.1\.12\b", version):
        raise ValueError(f"This adapter requires eggNOG-mapper 2.1.12; found {version.strip()}")
    if not re.search(r"Installed eggNOG DB version:\s*5\.0\.2\b", version):
        raise ValueError("This adapter requires the installed eggNOG database version 5.0.2")
    database_files = [
        data / name
        for name in (
            "eggnog.db",
            "eggnog.taxa.db",
            "eggnog.taxa.db.traverse.pkl",
            "eggnog_proteins.dmnd",
        )
    ]
    for path in database_files:
        if not path.is_file() or not path.stat().st_size:
            raise ValueError(f"Missing or empty eggNOG database file: {path}")
    identities = {str(p): sha256(p) for p in database_files + [emapper]}
    sequences, fasta_hash = prepare_annotation_targets(fasta, out, max_proteins=max_proteins)
    query_hash = sha256(out / "queries.fasta")
    argv = command + [
        "-i",
        str(out.resolve() / "queries.fasta"),
        "--itype",
        "proteins",
        "-o",
        "annotation",
        "--output_dir",
        str(out.resolve()),
        "--data_dir",
        str(data),
        "-m",
        "diamond",
        "--cpu",
        str(threads),
        "--sensmode",
        sensitivity,
        "--block_size",
        "0.5",
        "--tax_scope",
        tax_scope,
        "--target_orthologs",
        "all",
        "--seed_ortholog_evalue",
        "0.001",
        "--seed_ortholog_score",
        "60",
    ]
    provenance = {
        "fasta_sha256": fasta_hash,
        "query_sha256": query_hash,
        "source_files": identities,
        "version_output": version,
        "argv": argv,
        "unknown_fields": list(unknown_fields),
        "unique_sequences": len(sequences),
        "annotation_database_in_memory": False,
    }
    (out / "PROVENANCE.json").write_text(json.dumps(provenance, indent=2) + "\n")
    with (out / "eggnog.stdout").open("x") as stdout, (out / "eggnog.stderr").open("x") as stderr:
        resources = run_recorded(argv, stdout=stdout, stderr=stderr)
    (out / "RESOURCES.json").write_text(json.dumps(resources, indent=2) + "\n")
    if resources["exit_code"]:
        raise RuntimeError(f"eggNOG failed; inspect {out / 'eggnog.stderr'}")
    annotation = out / "annotation.emapper.annotations"
    rows, unknown = read_eggnog_annotations(annotation, sequences, unknown_fields=unknown_fields)
    write_table(out / "annotations.tsv", list(rows[0]), rows)
    with (out / "unknown.fasta").open("x") as stream:
        for query in unknown:
            stream.write(f">{query}\n{sequences[query]}\n")
    if (
        sha256(fasta) != fasta_hash
        or sha256(out / "queries.fasta") != query_hash
        or any(sha256(p) != h for p, h in identities.items())
    ):
        raise ValueError("Annotation source changed during execution")
    result = {
        "status": "PASS",
        "unique_sequences": len(sequences),
        "unknown_sequences": len(unknown),
        "unknown_fields": list(unknown_fields),
        "resources": resources,
        "files": {p.name: sha256(p) for p in out.iterdir() if p.is_file()},
    }
    (out / "COMPLETE.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


def read_annotation_settings(path):
    """Validate optional annotation settings; resolve local paths beside the JSON file.

    The configuration is deliberately separate from search settings. It contains
    local resources and annotation choices, never URLs to download during a run.
    """
    source = Path(path).resolve(strict=True)
    settings = json.loads(source.read_text())
    allowed = {
        "emapper",
        "emapper_python",
        "data_dir",
        "threads",
        "sensitivity",
        "tax_scope",
        "max_proteins",
        "unknown_fields",
        "generate_embeddings",
        "model",
        "device",
        "max_residues",
        "checkpoint",
    }
    if not isinstance(settings, dict) or set(settings) - allowed:
        raise ValueError("Annotation settings contain unknown keys; see ANNOTATION.md")
    if not {"emapper", "data_dir"} <= set(settings):
        raise ValueError("Annotation settings require emapper and data_dir local paths")
    for key in ("emapper", "emapper_python", "data_dir", "checkpoint"):
        if key in settings:
            if not isinstance(settings[key], str) or not settings[key]:
                raise ValueError(f"Annotation {key} must name a local path")
            value = Path(settings[key]).expanduser()
            value = source.parent / value if not value.is_absolute() else value
            value = (
                local_interpreter(value) if key == "emapper_python" else value.resolve(strict=True)
            )
            if (key == "data_dir" and not value.is_dir()) or (
                key != "data_dir" and not value.is_file()
            ):
                raise ValueError(f"Invalid local annotation path: {key}")
            settings[key] = value
    for key, maximum in (("threads", None), ("max_proteins", None), ("max_residues", 2046)):
        if key in settings and (
            type(settings[key]) is not int
            or settings[key] < 1
            or (maximum and settings[key] > maximum)
        ):
            raise ValueError(f"Invalid positive integer annotation setting: {key}")
    if type(settings.get("generate_embeddings", False)) is not bool:
        raise ValueError("generate_embeddings must be a JSON boolean")
    if settings.get("sensitivity", "sensitive") not in {
        "sensitive",
        "more-sensitive",
        "very-sensitive",
    }:
        raise ValueError("Unsupported annotation sensitivity")
    if not isinstance(settings.get("tax_scope", "auto"), str) or not re.fullmatch(
        r"auto|[1-9][0-9]*(?:,[1-9][0-9]*)*", settings.get("tax_scope", "auto")
    ):
        raise ValueError("Taxonomic scope must be auto or comma-separated taxon IDs")
    fields = settings.get("unknown_fields", ["KEGG_ko", "EC"])
    allowed_fields = {"KEGG_ko", "EC", "GOs", "PFAMs", "COG_category", "Description"}
    if not isinstance(fields, list) or not fields or any(not isinstance(x, str) for x in fields):
        raise ValueError("unknown_fields must be a nonempty list of annotation field names")
    if len(set(fields)) != len(fields) or not set(fields) <= allowed_fields:
        raise ValueError("unknown_fields must contain distinct supported annotation fields")
    if settings.get("model", "esmc_300m") not in {"esmc_300m", "esmc_600m"}:
        raise ValueError("Choose esmc_300m or esmc_600m")
    if settings.get("device", "cpu") not in {"cpu", "cuda"}:
        raise ValueError("Embedding device must be cpu or cuda")
    if settings.get("generate_embeddings", False):
        import importlib.metadata

        from fasta_lake.embeddings import resolve_checkpoint

        resolve_checkpoint(settings.get("checkpoint"), settings.get("model", "esmc_300m"))
        if importlib.metadata.version("esm") != "3.2.3":
            raise ValueError("Embedding settings require the optional esm==3.2.3 environment")
    return settings


def annotate_study(groups, settings_path, output):
    """Annotate the exact representatives exported by a completed study grouping.

    Representative annotations are not automatically propagated to all possible
    group members. Embedding is optional and uses the same recorded unknown rule.
    """
    import csv

    from fasta_lake.embeddings import embed_unknown

    source = Path(groups).resolve(strict=True)
    settings_hash = sha256(settings_path)
    settings = read_annotation_settings(settings_path)
    summary_path = source / "summary.json"
    summary_hash = sha256(summary_path)
    summary = json.loads(summary_path.read_text())
    expected = summary.get("annotation_targets", {})
    fasta = source / "representatives.fasta"
    if expected.get("fasta") != fasta.name or expected.get("sha256") != sha256(fasta):
        raise ValueError(
            "Study representatives are missing or changed; rerun grouping into a new folder"
        )
    group_path = source / "study_groups.tsv"
    group_hash = sha256(group_path)
    with group_path.open() as stream:
        groups = [r["study_group"] for r in csv.DictReader(stream, delimiter="\t")]
    records = read_fasta(fasta, max_records=settings.get("max_proteins", 100000))
    if len(groups) != len(set(groups)) or set(groups) != set(records):
        raise ValueError("Annotation representatives do not match the study dictionary")
    generate = settings.pop("generate_embeddings", False)
    embedding_keys = {"model", "device", "max_residues", "checkpoint"}
    embedding = {key: settings.pop(key) for key in embedding_keys if key in settings}
    result = annotate_fasta(fasta, output, **settings)
    if generate:
        result["embeddings"] = embed_unknown(
            output, Path(output) / "embeddings", threads=settings.get("threads", 2), **embedding
        )
    if (
        sha256(settings_path) != settings_hash
        or sha256(summary_path) != summary_hash
        or sha256(group_path) != group_hash
    ):
        raise ValueError("Study annotation inputs changed during annotation")
    result["study"] = {
        "summary_sha256": summary_hash,
        "groups_sha256": group_hash,
        "settings_sha256": settings_hash,
        "scope": "Reporting representatives; shared origin remains unresolved",
    }
    (Path(output) / "STUDY_COMPLETE.json").write_text(json.dumps(result, indent=2) + "\n")
    return result
