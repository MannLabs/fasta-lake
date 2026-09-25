"""Human-readable resource discovery, explicit download and bounded unpacking."""

import json
import lzma
import math
import tarfile
from pathlib import Path
from urllib.error import URLError

import click

from fasta_lake.downloads import (
    catalogue,
    download_resource,
    plan_download,
    unpack_resource,
)


def size_label(size):
    """Format bytes as GiB and exact bytes; preserve an unknown estimate as text."""
    return "unknown" if size is None else f"{size / 1024**3:.3f} GiB ({size:,} bytes)"


@click.group("resources")
def resources():
    """Choose, size and download public reference databases and model weights."""


@resources.command("list")
@click.option("--json", "as_json", is_flag=True)
def list_cmd(as_json):
    """Show installed choices and dated sizes; this command is offline."""
    data = catalogue()
    if as_json:
        click.echo(json.dumps(data, indent=2))
        return
    click.echo("Sizes are dated catalogue snapshots; plan checks current transfer metadata.\n")
    for key, spec in data.items():
        sizes = [f["bytes"] for f in spec["files"]]
        total = sum(sizes) if all(x is not None for x in sizes) else None
        approximate = "; ".join(
            x["approximate_download_size"]
            for x in spec["files"]
            if "approximate_download_size" in x
        )
        click.echo(f"{key:24s} {spec['title']}\n  Download: {size_label(total)} {approximate}")


def show_plan(plan):
    """Print the selected resource, destinations, size estimates and checksums."""
    click.echo(f"{plan['title']}\nVersion: {plan['version']}")
    click.echo(f"Destination: {Path(plan['destination']) / plan['id']}")
    click.echo(f"Download: {size_label(plan['transfer_bytes'])}")
    click.echo(f"Expanded estimate: {size_label(plan.get('expanded_bytes_estimate'))}")
    click.echo(f"Disk currently free: {size_label(plan['free_disk_bytes'])}")
    click.echo(f"Source: {plan['source']}\n{plan['notes']}")
    for item in plan["files"]:
        click.echo(f"  {item['name']}: {size_label(item['remote']['bytes'])}")
        click.echo(f"  {item['url']}")
        click.echo(
            f"  Upstream checksum: {item.get('checksum') or 'unavailable; local SHA256 recorded'}"
        )


@resources.command("plan")
@click.argument("name")
@click.option("--directory", type=click.Path(path_type=Path), required=True)
@click.option("--json", "as_json", is_flag=True)
def plan_cmd(name, directory, as_json):
    """Show release, download bytes and disk needs without downloading data."""
    try:
        plan = plan_download(name, directory)
    except (OSError, ValueError, URLError) as error:
        raise click.ClickException(str(error)) from error
    if as_json:
        click.echo(json.dumps(plan, indent=2))
    else:
        show_plan(plan)


@resources.command("download")
@click.argument("name")
@click.option("--directory", type=click.Path(path_type=Path), required=True)
@click.option("--resume", is_flag=True, help="Resume an unchanged interrupted resource.")
def download_cmd(name, directory, resume):
    """Show the plan, then explicitly download one chosen resource; do not unpack."""
    try:
        plan = plan_download(name, directory)
        show_plan(plan)
        result = download_resource(
            plan,
            resume=resume,
            progress=lambda name, done, total: click.echo(
                f"{name}: {size_label(done)} / {size_label(total)}"
            ),
        )
    except (OSError, ValueError, URLError) as error:
        raise click.ClickException(str(error)) from error
    click.echo(f"Verified download: {result['directory']}/DOWNLOAD_COMPLETE.json")


@resources.command("unpack")
@click.argument(
    "archives",
    nargs=-1,
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option("--directory", type=click.Path(path_type=Path), required=True)
@click.option(
    "--max-gib",
    type=click.FloatRange(min=0, min_open=True),
    required=True,
    help="Maximum expanded disk bytes in GiB; original archives also occupy disk.",
)
def unpack_cmd(archives, directory, max_gib):
    """Unpack archives into a new directory within an explicit disk budget."""
    try:
        if not math.isfinite(max_gib):
            raise ValueError("Choose a finite expanded disk budget")
        result = unpack_resource(archives, directory, max_bytes=int(max_gib * 1024**3))
    except (OSError, ValueError, EOFError, tarfile.TarError, lzma.LZMAError) as error:
        raise click.ClickException(str(error)) from error
    click.echo(f"Expanded {size_label(result['expanded_bytes'])} into {directory}")
