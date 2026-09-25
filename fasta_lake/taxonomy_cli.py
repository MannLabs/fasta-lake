"""Use Unipept query files and native taxonomy output with a checked study."""

import json
from pathlib import Path

import click


@click.group("taxonomy")
def taxonomy():
    """Export Unipept queries and plot its saved taxonomic assignments locally."""


@taxonomy.command("prepare")
@click.option(
    "--groups", required=True, type=click.Path(exists=True, file_okay=False, path_type=Path)
)
@click.option("--out", required=True, type=click.Path(path_type=Path))
def prepare(groups, out):
    """Write canonical peptides for Unipept pept2lca --all --equate."""
    from fasta_lake.taxonomy import prepare_unipept_queries

    try:
        result = prepare_unipept_queries(groups, out)
    except (ValueError, OSError) as error:
        raise click.ClickException(str(error)) from error
    click.echo(json.dumps(result, indent=2))


@taxonomy.command("report")
@click.option(
    "--groups", required=True, type=click.Path(exists=True, file_okay=False, path_type=Path)
)
@click.option(
    "--unipept-output", required=True, type=click.Path(exists=True, dir_okay=False, path_type=Path)
)
@click.option("--out", required=True, type=click.Path(path_type=Path))
@click.option(
    "--database", required=True, help="Reference/version reported by the selected Unipept server."
)
@click.option(
    "--equate-il", is_flag=True, help="Declare that Unipept was run with --equate (required)."
)
@click.option("--host-taxon-id", default=9606, show_default=True, type=click.IntRange(min=1))
@click.option(
    "--contaminant-peptides",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Optional explicit reference: one unmodified peptide per line.",
)
def report(groups, unipept_output, out, database, equate_il, host_taxon_id, contaminant_peptides):
    """Render upstream sunburst/treemap/tree with explicit counts, signal and citations."""
    from fasta_lake.taxonomy import taxonomy_report

    try:
        result = taxonomy_report(
            groups,
            unipept_output,
            out,
            database=database,
            equate_il=equate_il,
            host_taxon_id=host_taxon_id,
            contaminant_peptides=contaminant_peptides,
        )
    except (ValueError, OSError) as error:
        raise click.ClickException(str(error)) from error
    click.echo(json.dumps(result, indent=2))


@taxonomy.command("prepare-reference")
@click.option(
    "--fasta", required=True, type=click.Path(exists=True, dir_okay=False, path_type=Path)
)
@click.option("--out", required=True, type=click.Path(path_type=Path))
@click.option(
    "--taxonomy-map",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Explicit accession/ncbi_taxon_id TSV; unmapped proteins stay unresolved.",
)
@click.option("--max-proteins", default=100_000, show_default=True, type=click.IntRange(min=1))
def prepare_reference(fasta, out, taxonomy_map, max_proteins):
    """Experimental: export a selected FASTA for Unipept's native index builder."""
    from fasta_lake.taxonomy_reference import prepare_unipept_reference

    try:
        result = prepare_unipept_reference(
            fasta, out, taxonomy=taxonomy_map, max_proteins=max_proteins
        )
    except (ValueError, OSError) as error:
        raise click.ClickException(str(error)) from error
    click.echo(json.dumps(result, indent=2))
