"""Installed optional experimental functional-profile workflow."""

import json
from pathlib import Path

import click


@click.command("explore-study")
@click.option(
    "--matrix", required=True, type=click.Path(exists=True, dir_okay=False, path_type=Path)
)
@click.option(
    "--annotations", required=True, type=click.Path(exists=True, dir_okay=False, path_type=Path)
)
@click.option("--out", required=True, type=click.Path(path_type=Path))
@click.option(
    "--level", required=True, help="One declared annotation level, such as KO, EC or genus."
)
@click.option("--metadata", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@click.option(
    "--group-column", help="Explicitly request one-way PERMANOVA for this metadata column."
)
@click.option(
    "--block-column", help="Restrict label permutations to donor/repeated-measure blocks."
)
@click.option("--permutations", type=click.IntRange(min=1), default=999, show_default=True)
@click.option("--seed", type=click.IntRange(min=0), default=0, show_default=True)
@click.option(
    "--backend", type=click.Choice(["numpy", "numba"]), default="numpy", show_default=True
)
def explore_cmd(
    matrix,
    annotations,
    out,
    level,
    metadata,
    group_column,
    block_column,
    permutations,
    seed,
    backend,
):
    """EXPERIMENTAL: inspect annotated signal, coverage, diversity and CLR PCA."""
    from fasta_lake.exploration import explore_study

    try:
        result = explore_study(
            matrix,
            annotations,
            out,
            level=level,
            metadata=metadata,
            group_column=group_column,
            block_column=block_column,
            permutations=permutations,
            seed=seed,
            backend=backend,
        )
    except (OSError, ValueError, RuntimeError, ImportError) as error:
        raise click.ClickException(str(error)) from error
    click.echo(
        json.dumps(
            {"status": result["status"], "experimental": True, "outputs": result["outputs"]},
            indent=2,
        )
    )
