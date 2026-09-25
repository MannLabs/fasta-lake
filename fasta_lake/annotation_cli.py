"""Installed entry points for optional protein annotation and embeddings."""

import json
from pathlib import Path

import click


@click.group("annotate")
def annotation():
    """Annotate selected proteins and inspect explicit annotation gaps."""


@annotation.command("eggnog")
@click.option(
    "--fasta", type=click.Path(exists=True, dir_okay=False, path_type=Path), required=True
)
@click.option("--out", type=click.Path(path_type=Path), required=True)
@click.option(
    "--emapper", type=click.Path(exists=True, dir_okay=False, path_type=Path), required=True
)
@click.option("--emapper-python", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@click.option(
    "--data-dir", type=click.Path(exists=True, file_okay=False, path_type=Path), required=True
)
@click.option("--threads", type=click.IntRange(min=1), default=2, show_default=True)
@click.option(
    "--sensitivity",
    type=click.Choice(["sensitive", "more-sensitive", "very-sensitive"]),
    default="sensitive",
    show_default=True,
)
@click.option("--tax-scope", default="auto", show_default=True)
@click.option("--max-proteins", type=click.IntRange(min=1), default=100000, show_default=True)
@click.option(
    "--unknown-field",
    "unknown_fields",
    multiple=True,
    default=["KEGG_ko", "EC"],
    type=click.Choice(["KEGG_ko", "EC", "GOs", "PFAMs", "COG_category", "Description"]),
    help="Unknown means no value in any selected field; repeat to choose fields.",
)
@click.option(
    "--embed-unknown",
    "generate_embeddings",
    is_flag=True,
    help="After annotation, generate local ESM-C vectors for its unknown proteins.",
)
@click.option("--model", type=click.Choice(["esmc_300m", "esmc_600m"]), default="esmc_300m")
@click.option(
    "--checkpoint",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Explicit local model weights; resources download installs them separately.",
)
@click.option("--device", type=click.Choice(["cpu", "cuda"]), default="cpu")
@click.option("--max-residues", type=click.IntRange(min=1, max=2046), default=1024)
def eggnog_cmd(
    fasta,
    out,
    emapper,
    emapper_python,
    data_dir,
    threads,
    sensitivity,
    tax_scope,
    max_proteins,
    unknown_fields,
    generate_embeddings,
    model,
    device,
    max_residues,
    checkpoint,
):
    """Run eggNOG-mapper 2.1.12 with versioned inputs and complete coverage records."""
    from fasta_lake.annotation import annotate_fasta
    from fasta_lake.embeddings import embed_unknown, resolve_checkpoint

    try:
        if generate_embeddings:
            resolve_checkpoint(checkpoint, model)
            import importlib.metadata

            if importlib.metadata.version("esm") != "3.2.3":
                raise ValueError("Automatic embedding requires the esm==3.2.3 optional environment")
        result = annotate_fasta(
            fasta,
            out,
            emapper=emapper,
            emapper_python=emapper_python,
            data_dir=data_dir,
            threads=threads,
            sensitivity=sensitivity,
            tax_scope=tax_scope,
            max_proteins=max_proteins,
            unknown_fields=unknown_fields,
        )
        if generate_embeddings:
            result["embeddings"] = embed_unknown(
                out,
                out / "embeddings",
                model=model,
                device=device,
                threads=threads,
                max_residues=max_residues,
                checkpoint=checkpoint,
            )
    except (OSError, ValueError, RuntimeError, ImportError) as error:
        raise click.ClickException(str(error)) from error
    click.echo(json.dumps(result, indent=2))


@annotation.command("embed-unknown")
@click.option(
    "--annotation",
    "source",
    required=True,
    type=click.Path(exists=True, file_okay=False, path_type=Path),
)
@click.option("--out", type=click.Path(path_type=Path), required=True)
@click.option(
    "--model", type=click.Choice(["esmc_300m", "esmc_600m"]), default="esmc_300m", show_default=True
)
@click.option("--device", type=click.Choice(["cpu", "cuda"]), default="cpu", show_default=True)
@click.option("--threads", type=click.IntRange(min=1), default=2, show_default=True)
@click.option(
    "--max-residues", type=click.IntRange(min=1, max=2046), default=1024, show_default=True
)
@click.option(
    "--checkpoint",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Explicit local model weight file; analysis never downloads weights.",
)
def embed_cmd(source, out, model, device, threads, max_residues, checkpoint):
    """Embed the recorded unknown roster; preserve every residue through windowing."""
    from fasta_lake.embeddings import embed_unknown

    try:
        result = embed_unknown(
            source,
            out,
            model=model,
            device=device,
            threads=threads,
            max_residues=max_residues,
            checkpoint=checkpoint,
        )
    except (OSError, ValueError, RuntimeError, ImportError) as error:
        raise click.ClickException(str(error)) from error
    click.echo(json.dumps(result, indent=2))


@annotation.command("study")
@click.option(
    "--groups", required=True, type=click.Path(exists=True, file_okay=False, path_type=Path)
)
@click.option(
    "--settings", required=True, type=click.Path(exists=True, dir_okay=False, path_type=Path)
)
@click.option("--out", required=True, type=click.Path(path_type=Path))
def study_cmd(groups, settings, out):
    """Annotate a study's verified reporting representatives, with optional ESM-C."""
    from fasta_lake.annotation import annotate_study

    try:
        result = annotate_study(groups, settings, out)
    except (OSError, ValueError, RuntimeError, ImportError) as error:
        raise click.ClickException(str(error)) from error
    click.echo(json.dumps(result, indent=2))
