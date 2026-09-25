"""Commands for portable peptide/protein evidence and optional CKG views."""

import json
from pathlib import Path

import click


@click.group("network")
def network():
    """Export peptide memberships and render experimental CKG network views."""


@network.command("export")
@click.option(
    "--groups", required=True, type=click.Path(exists=True, file_okay=False, path_type=Path)
)
@click.option("--out", required=True, type=click.Path(path_type=Path))
@click.option("--max-rows", default=1_000_000, show_default=True, type=click.IntRange(min=1))
@click.option("--max-edges", default=1_000_000, show_default=True, type=click.IntRange(min=1))
def export_cmd(groups, out, max_rows, max_edges):
    """Export complete accepted memberships from a group-study result."""
    from fasta_lake.networks import export_evidence_network

    try:
        result = export_evidence_network(groups, out, max_rows=max_rows, max_edges=max_edges)
    except (OSError, ValueError, RuntimeError) as error:
        raise click.ClickException(str(error)) from error
    click.echo(json.dumps(result, indent=2))


@network.command("render")
@click.option(
    "--bundle", required=True, type=click.Path(exists=True, file_okay=False, path_type=Path)
)
@click.option("--out", required=True, type=click.Path(path_type=Path))
@click.option(
    "--annotations",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    help="Checked network eggnog output to show in the evidence view.",
)
@click.option(
    "--ckg-python",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Python executable from a compatible CKG environment.",
)
@click.option(
    "--ckg-source",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    help="Optional CKG source checkout to use with that environment.",
)
@click.option(
    "--max-nodes",
    default=150,
    show_default=True,
    type=click.IntRange(min=2),
    help="Limit the view; complete exported evidence is retained.",
)
@click.option("--max-edges", default=300, show_default=True, type=click.IntRange(min=1))
def render_cmd(bundle, out, annotations, ckg_python, ckg_source, max_nodes, max_edges):
    """Create an experimental CKG report from a completed export."""
    from fasta_lake.networks import render_ckg_network

    try:
        result = render_ckg_network(
            bundle,
            out,
            ckg_python=ckg_python,
            ckg_source=ckg_source,
            max_nodes=max_nodes,
            max_edges=max_edges,
            annotations=annotations,
        )
    except (OSError, ValueError, RuntimeError) as error:
        raise click.ClickException(str(error)) from error
    click.echo(json.dumps(result, indent=2))


@network.command("eggnog")
@click.option(
    "--bundle", required=True, type=click.Path(exists=True, file_okay=False, path_type=Path)
)
@click.option(
    "--annotation",
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help="Completed FastaLake annotation folder or native eggNOG output.",
)
@click.option(
    "--fasta",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Study representatives.fasta or a recorded search FASTA.",
)
@click.option(
    "--query-fasta",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Original annotation queries; required only for native eggNOG output.",
)
@click.option("--out", required=True, type=click.Path(path_type=Path))
@click.option(
    "--field",
    "fields",
    multiple=True,
    default=("KEGG_ko", "EC", "GOs", "PFAMs"),
    show_default=True,
    help="Annotation namespace; repeat for multiple namespaces.",
)
@click.option("--max-records", default=100_000, type=click.IntRange(min=1), show_default=True)
def eggnog_cmd(bundle, annotation, fasta, query_fasta, out, fields, max_records):
    """Attach eggNOG annotations by exact protein sequence; retain unknowns."""
    from fasta_lake.network_annotations import annotate_evidence_network

    try:
        result = annotate_evidence_network(
            bundle,
            annotation,
            fasta,
            out,
            query_fasta=query_fasta,
            fields=fields,
            max_records=max_records,
        )
    except (OSError, ValueError, RuntimeError) as error:
        raise click.ClickException(str(error)) from error
    click.echo(json.dumps(result, indent=2))
