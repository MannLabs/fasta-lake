"""Installed commands for specimen-matched molecular evidence."""

import json

import click


def checked(function, *args, **kwargs):
    """Translate expected input and execution failures into readable Click errors."""
    try:
        return function(*args, **kwargs)
    except (OSError, ValueError, RuntimeError, KeyError) as error:
        raise click.ClickException(str(error)) from error


@click.group("molecular")
def molecular():
    """Prepare TPM inputs, select exact proteins and compare matched searches."""


@molecular.command("prepare-source")
@click.option("--out", required=True, type=click.Path())
@click.option("--specimen", required=True)
@click.option("--reference-id", required=True)
@click.option(
    "--reference-scope",
    type=click.Choice(["complete", "available"]),
    default="complete",
    show_default=True,
    help="Available explicitly reports genes without sequences and excludes them.",
)
@click.option("--tpm", required=True, type=click.Path(exists=True, dir_okay=False))
@click.option("--proteins", required=True, type=click.Path(exists=True, dir_okay=False))
@click.option("--provenance", required=True, type=click.Path(exists=True, dir_okay=False))
@click.option("--alternate-proteins", type=click.Path(exists=True, dir_okay=False))
@click.option("--annotations", type=click.Path(exists=True, dir_okay=False))
@click.option("--metabolites", type=click.Path(exists=True, dir_okay=False))
@click.option("--compound-links", type=click.Path(exists=True, dir_okay=False))
def prepare_source(out, **kwargs):
    """Check and copy specimen inputs into a new portable source directory."""
    from .sources import prepare_source as prepare

    click.echo(json.dumps(checked(prepare, out, **kwargs), indent=2))


@molecular.command("validate-source")
@click.option("--source", required=True, type=click.Path(exists=True, dir_okay=False))
def validate_source(source):
    """Check source hashes, gene agreement, annotation links and declared scope."""
    from .sources import load_reference

    reference = checked(load_reference, source)
    click.echo(
        json.dumps(
            dict(
                status="PASS",
                specimen=reference.manifest["specimen"],
                reference_id=reference.manifest["reference_id"],
                genes=len(reference.available),
                exact_protein_sequences=len(
                    {reference.proteins[g][2] for g in reference.available}
                ),
                coverage=reference.coverage,
                scope="File integrity and gene agreement; provider supplies provenance.",
            ),
            indent=2,
        )
    )


def selection_options(default_dna):
    """Build a Click decorator for independent DNA/RNA selection thresholds."""
    def decorate(function):
        """Attach source, specimen, output and molecular-selection options in order."""
        options = [
            click.option("--source", required=True, type=click.Path(exists=True, dir_okay=False)),
            click.option(
                "--specimen", required=True, help="Specimen matching the source manifest."
            ),
            click.option("--out", required=True, type=click.Path()),
            click.option(
                "--background",
                type=click.Path(exists=True, dir_okay=False),
                help="Identical target-only host/contaminant background for compared arms.",
            ),
            click.option(
                "--dna-percent",
                type=float,
                default=default_dna,
                show_default=True,
                help="Top percent of positive available DNA genes; zero disables.",
            ),
            click.option(
                "--rna-percent",
                type=float,
                default=0,
                show_default=True,
                help="Top percent of positive available RNA genes; zero disables.",
            ),
            click.option(
                "--dna-absolute",
                type=float,
                help="Strict DNA TPM threshold; requires dna-percent=0.",
            ),
            click.option(
                "--rna-absolute",
                type=float,
                help="Strict RNA TPM threshold; requires rna-percent=0.",
            ),
            click.option(
                "--rank-rounding",
                type=click.Choice(["floor", "ceil"]),
                default="floor",
                show_default=True,
            ),
            click.option(
                "--metabolites",
                is_flag=True,
                help="Also select through positive compounds and direct reaction/KO links.",
            ),
        ]
        for option in reversed(options):
            function = option(function)
        return function

    return decorate


@molecular.command("add")
@click.option("--baseline", required=True, type=click.Path(exists=True, dir_okay=False))
@selection_options(1.0)
def add(baseline, source, specimen, out, rank_rounding, **kwargs):
    """Keep a final de novo FASTA and append DNA/RNA/metabolite-supported proteins."""
    from .selection import build_selection

    result = checked(
        build_selection, source, specimen, out, baseline=baseline, rounding=rank_rounding, **kwargs
    )
    click.echo(json.dumps(result, indent=2))


@molecular.command("select")
@selection_options(0.0)
def select(source, specimen, out, rank_rounding, **kwargs):
    """Build a standalone molecular FASTA with an optional common background."""
    from .selection import build_selection

    result = checked(build_selection, source, specimen, out, rounding=rank_rounding, **kwargs)
    click.echo(json.dumps(result, indent=2))


@molecular.command("compare")
@click.option("--left", required=True, type=click.Path(exists=True, file_okay=False))
@click.option("--right", required=True, type=click.Path(exists=True, file_okay=False))
@click.option("--out", required=True, type=click.Path())
@click.option("--left-label", default="metaG", show_default=True)
@click.option("--right-label", default="de_novo", show_default=True)
@click.option("--peptide-q", type=float, default=0.01, show_default=True)
@click.option("--protein-q", type=float)
@click.option("--compress-tables", is_flag=True, help="Write overlap tables as gzip TSV.")
def compare_searches(
    left, right, out, left_label, right_label, peptide_q, protein_q, compress_tables
):
    """Compare exact proteins, accepted peptides and LFQ in two matched Sage searches."""
    from .complementarity import compare

    try:
        result = compare(
            left,
            right,
            out,
            q=peptide_q,
            protein_q=protein_q,
            labels=(left_label, right_label),
            compress_tables=compress_tables,
        )
    except (OSError, ValueError, RuntimeError, KeyError) as error:
        raise click.ClickException(str(error)) from error
    click.echo(json.dumps(result, indent=2))


@molecular.command("increments")
@click.option("--left", required=True, type=click.Path(exists=True, file_okay=False))
@click.option("--right", required=True, type=click.Path(exists=True, file_okay=False))
@click.option(
    "--reference",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help="Full declared target reference including the shared host background.",
)
@click.option("--host-reference", type=click.Path(exists=True, dir_okay=False))
@click.option("--out", required=True, type=click.Path())
@click.option("--left-label", default="de_novo", show_default=True)
@click.option("--right-label", default="with_molecular", show_default=True)
@click.option("--peptide-q", type=float, default=0.01, show_default=True)
@click.option("--protein-q", type=float)
def increments(
    left, right, reference, host_reference, out, left_label, right_label, peptide_q, protein_q
):
    """List gained/lost peptide and protein support against one complete reference."""
    from .increments import increment_report

    try:
        result = increment_report(
            left,
            right,
            reference,
            out,
            host_reference=host_reference,
            q=peptide_q,
            protein_q=protein_q,
            labels=(left_label, right_label),
        )
    except (OSError, ValueError, RuntimeError, KeyError) as error:
        raise click.ClickException(str(error)) from error
    click.echo(json.dumps(result, indent=2))
