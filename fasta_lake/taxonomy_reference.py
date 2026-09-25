"""Experimental sequence-preserving adapter for the upstream Unipept index builder.

Native index format: https://github.com/unipept/unipept-index
(commit 718a81e8285d9889839710129389575efe55014a in the small lookup pilot).
This is a FastaLake format/identity adapter; upstream builds and searches the
index. Taxonomic mapping and a complete custom LCA service are separate.
"""

import csv
import json
from pathlib import Path

from fasta_lake.multi_omics.exact import read_fasta, write_table
from fasta_lake.networks import _sha


def prepare_unipept_reference(fasta, output, *, taxonomy=None, max_proteins=100_000):
    """Write sa-builder's four-field TSV from a selected protein FASTA.

    The optional map has accession and ncbi_taxon_id columns. IDs are supplied
    by the caller, not guessed from accession names or GTDB species strings.
    Unmapped sequences remain in the index with taxon 0 (unresolved). This
    prepares peptide-to-protein lookup; a full Unipept LCA service additionally
    needs a compatible taxonomy/database backend. It is not created here.
    """
    if type(max_proteins) is not int or max_proteins < 1:
        raise ValueError("max_proteins must be a positive integer")
    out = Path(output)
    if out.exists() or out.is_symlink():
        raise FileExistsError(f"Output exists: {out}; choose a new directory")
    paths = [Path(fasta).resolve(strict=True)]
    if taxonomy:
        paths.append(Path(taxonomy).resolve(strict=True))
    before = {str(p): _sha(p) for p in paths}
    proteins = read_fasta(paths[0], max_records=max_proteins)
    if not proteins or any(k.startswith("rev_") for k in proteins):
        raise ValueError("Supply a nonempty target-only FASTA")
    taxa = {}
    if taxonomy:
        with paths[1].open() as f:
            rows = csv.DictReader(f, delimiter="\t")
            if rows.fieldnames != ["accession", "ncbi_taxon_id"]:
                raise ValueError("Taxonomy map requires accession and ncbi_taxon_id columns")
            for row in rows:
                if None in row or any(v is None for v in row.values()):
                    raise ValueError("Malformed taxonomy mapping row")
                accession = row["accession"]
                if accession in taxa or accession not in proteins:
                    raise ValueError("Duplicate or foreign accession in taxonomy map")
                value = row["ncbi_taxon_id"]
                if not value.isascii() or not value.isdecimal() or not 1 <= int(value) <= 2**32 - 1:
                    raise ValueError("Mapped NCBI taxon IDs must be positive uint32 integers")
                taxa[accession] = int(value)
    if before != {str(p): _sha(p) for p in paths}:
        raise ValueError("Reference inputs changed while being read")
    out.mkdir(parents=True)
    with (out / "proteins.tsv").open("w") as f:
        for accession, (_, sequence, _) in proteins.items():
            # No header: this is the exact upstream sa-builder input format.
            # Empty annotations remain empty; taxonomic mapping does not invent function.
            f.write(f"{accession}\t{taxa.get(accession, 0)}\t{sequence}\t\n")
    write_table(
        out / "protein_identity.tsv",
        ["accession", "sequence_sha256", "ncbi_taxon_id", "taxonomy_status"],
        [
            dict(
                accession=k,
                sequence_sha256=v[2],
                ncbi_taxon_id=taxa.get(k, 0),
                taxonomy_status="caller_supplied_mapping" if k in taxa else "unresolved",
            )
            for k, v in proteins.items()
        ],
    )
    receipt = dict(
        status="PASS_REFERENCE_EXPORT",
        experimental=True,
        proteins=len(proteins),
        mapped_proteins=len(taxa),
        unresolved_proteins=len(proteins) - len(taxa),
        inputs=before,
        format="accession, uint32 taxon ID, exact sequence, empty annotations; TSV without header",
        taxonomy_scope="Caller-supplied NCBI mapping; GTDB names are not converted automatically. "
        "Taxon 0 is unresolved.",
        service_scope="Upstream peptide-to-protein index input; "
        "a full LCA/API service is a separate integration.",
        files={p.name: _sha(p) for p in out.iterdir() if p.is_file()},
    )
    (out / "COMPLETE.json").write_text(json.dumps(receipt, indent=2) + "\n")
    return receipt
