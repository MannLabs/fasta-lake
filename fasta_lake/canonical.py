"""
Canonical database management for FastaLake baseline comparisons.

Downloads and manages canonical SwissProt FASTAs for the baseline search
that FastaLake is benchmarked against.
"""

from __future__ import annotations

import logging
import os
import shutil
import urllib.request
from pathlib import Path

from fasta_lake.helpers import count_fasta_headers

logger = logging.getLogger(__name__)

# UniProt asks callers to identify themselves, and this header used to carry a named
# author's institutional e-mail address. That is one person's personal data, hard-coded
# into a package, sent to a third-party server on every download, by every user who
# installs it -- none of whom is that person. It was also the only piece of site-specific
# configuration in the core package.
#
# Identify the SOFTWARE by default. A contact is a property of the deployment, not of the
# source, so it comes from the environment and is simply absent when unset.
CONTACT_ENV = "FASTALAKE_CONTACT"
USER_AGENT_BASE = "FastaLake (+https://github.com/mannlabs/FastaLake)"


def _user_agent() -> str:
    """Build the request identity, appending the optional configured contact."""
    contact = os.environ.get(CONTACT_ENV, "").strip()
    return f"{USER_AGENT_BASE[:-1]}; {contact})" if contact else USER_AGENT_BASE


# UniProt download URLs
UNIPROT_URLS = {
    "human": "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28%28organism_id%3A9606%29+AND+%28reviewed%3Atrue%29%29",
    "mouse": "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28%28organism_id%3A10090%29+AND+%28reviewed%3Atrue%29%29",
    "rat": "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28%28organism_id%3A10116%29+AND+%28reviewed%3Atrue%29%29",
}


def download_canonical(
    organism: str = "human",
    output_dir: str | Path = ".",
    force: bool = False,
) -> Path:
    """Download the latest canonical SwissProt FASTA for an organism.

    Parameters
    ----------
    organism : str
        'human', 'mouse', or 'rat'.
    output_dir : path
        Where to save the FASTA.
    force : bool
        Re-download even if file exists.

    Returns
    -------
    Path
        Path to the downloaded FASTA.
    """
    if organism not in UNIPROT_URLS:
        raise ValueError(f"Unknown organism: {organism}. Available: {list(UNIPROT_URLS)}")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / f"swissprot_{organism}_canonical.fasta"

    if output_path.exists() and not force:
        n = count_fasta_headers(output_path)
        logger.info("Canonical already exists: %s (%d entries)", output_path, n)
        return output_path

    url = UNIPROT_URLS[organism]
    logger.info("Downloading canonical %s from UniProt...", organism)

    try:
        req = urllib.request.Request(url)
        req.add_header("User-Agent", _user_agent())

        with urllib.request.urlopen(req) as response:
            content = response.read()

        output_path.write_bytes(content)
        n = count_fasta_headers(output_path)
        logger.info("Downloaded: %s (%d entries)", output_path.name, n)

    except Exception as e:
        logger.error("Download failed: %s", e)
        # Fall back to local copy if available
        local = _find_local_canonical(organism)
        if local:
            logger.info("Using local copy: %s", local)
            shutil.copy2(local, output_path)
        else:
            raise

    return output_path


def _find_local_canonical(organism: str) -> Path | None:
    """Try to find a local canonical FASTA on the cluster."""
    candidates = [
        # FastaLake databases
        Path(__file__).parent.parent / "databases" / "isoforms" / organism / "uniprot",
    ]

    species_name = organism.replace("human", "homo_sapiens").replace("mouse", "mus_musculus")
    name_patterns = [
        f"swissprot_{organism}_reviewed",
        f"swissprot_{species_name}_swissprot",
    ]

    for cand_dir in candidates:
        if not cand_dir.exists():
            continue
        for f in cand_dir.iterdir():
            if f.suffix == ".fasta" and any(p in f.name.lower() for p in name_patterns):
                return f
            # Also check for generic swissprot files
            if (
                "swissprot" in f.name.lower()
                and "reviewed" in f.name.lower()
                and f.suffix == ".fasta"
            ):
                return f

    return None


def get_canonical(organism: str = "human") -> Path | None:
    """Get path to canonical FASTA, preferring local copy.

    Does NOT download — returns None if not found locally.
    Use download_canonical() to fetch from UniProt.
    """
    local = _find_local_canonical(organism)
    if local:
        return local

    # Check common cluster paths
    import os

    peter = os.environ.get("PETER")
    if not peter:
        return None
    cluster_paths = [
        Path(peter) / "FastaLake" / "databases" / "isoforms" / organism / "uniprot",
    ]

    for p in cluster_paths:
        if not p.exists():
            continue
        for f in sorted(p.glob("swissprot_*reviewed*.fasta"), reverse=True):
            return f

    return None
