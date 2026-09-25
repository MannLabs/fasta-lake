"""Portable specimen-qualified DNA/RNA and optional metabolite inputs."""

from __future__ import annotations

import csv
import gzip
import json
import re
import shutil
from dataclasses import dataclass
from pathlib import Path

from .exact import load_source, numeric, read_fasta, read_tpm, sha256, write_table


def table_rows(path, required):
    """Yield string-valued rows from plain/gzipped TSV with a complete unique header.

    ``required`` names mandatory columns. Malformed rows raise ValueError; numeric
    conversion and missing-value interpretation belong to the calling parser.
    """
    opener = gzip.open if Path(path).suffix == ".gz" else open
    with opener(path, "rt", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        fields = reader.fieldnames or []
        if len(fields) != len(set(fields)) or not set(required) <= set(fields):
            raise ValueError(f"Expected unique TSV columns: {', '.join(required)}")
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise ValueError("Malformed TSV row")
            yield row


@dataclass
class MolecularSource:
    """Verified specimen inputs and their gene-level relationships.

    ``genes`` retains parsed DNA/RNA values and original TPM text. ``proteins``
    maps gene IDs to (header, normalized sequence, exact-sequence SHA256).
    ``available`` is their intersection; coverage retains the excluded genes.
    Optional annotations, compounds and direct reaction links remain separate.
    """

    manifest: dict
    inputs: dict
    genes: dict
    proteins: dict
    annotations: dict
    compounds: dict
    links: list

    @property
    def available(self):
        """Return sorted genes with both a TPM row and a protein sequence."""
        return sorted(set(self.genes) & set(self.proteins))

    @property
    def coverage(self):
        """Summarize available sequences, unavailable genes and missing molecular values."""
        available = set(self.available)
        missing = set(self.genes) - available
        return {
            "reference_scope": self.manifest.get("reference_scope", "complete"),
            "raw_tpm_genes": len(self.genes),
            "available_genes": len(available),
            "unavailable_genes": len(missing),
            "proteins_without_tpm": len(set(self.proteins) - set(self.genes)),
            "unavailable_positive_dna": sum(
                self.genes[g]["dna"] is not None and self.genes[g]["dna"] > 0 for g in missing
            ),
            "unavailable_positive_rna": sum(
                self.genes[g]["rna"] is not None and self.genes[g]["rna"] > 0 for g in missing
            ),
            "missing_dna_in_available": sum(self.genes[g]["dna"] is None for g in available),
            "missing_rna_in_available": sum(self.genes[g]["rna"] is None for g in available),
            "available_genes_with_ko": sum(bool(self.annotations.get(g)) for g in available),
            "positive_compounds": sum(v is not None and v > 0 for v in self.compounds.values()),
        }


def _inputs(obj, base):
    """Resolve source-manifest paths and verify required inputs and pinned hashes."""
    inputs = {}
    for key in (
        "tpm",
        "proteins",
        "provenance",
        "alternate_proteins",
        "annotations",
        "metabolites",
        "compound_links",
    ):
        if key not in obj:
            if key in ("tpm", "proteins", "provenance"):
                raise ValueError(f"Missing source input: {key}")
            continue
        record = obj[key]
        path = Path(record["path"])
        if not path.is_absolute():
            path = base / path
        path = path.resolve(strict=True)
        if sha256(path) != record["sha256"]:
            raise ValueError(f"{key} checksum mismatch")
        inputs[key] = path
    return inputs


def _load_v2(obj, base):
    """Validate a specimen-qualified v2 source, its sequence coverage and annotations."""
    for key in ("specimen", "reference_id"):
        if not isinstance(obj.get(key), str) or not obj[key].strip():
            raise ValueError(f"Missing {key}")
    if obj.get("reference_scope") not in ("complete", "available"):
        raise ValueError("Reference scope must be complete or explicitly available")
    inputs = _inputs(obj, base)
    genes = read_tpm(inputs["tpm"], nan_is_missing=True)
    proteins = read_fasta(inputs["proteins"])
    common = set(genes) & set(proteins)
    if not common:
        raise ValueError("No shared TPM/protein genes")
    if obj["reference_scope"] == "complete" and set(genes) != set(proteins):
        raise ValueError(
            "All-gene reference mismatch; use available scope only with an exclusion ledger"
        )
    if "alternate_proteins" in inputs:
        alternate = read_fasta(inputs["alternate_proteins"])
        shared = set(proteins) & set(alternate)
        if not shared:
            raise ValueError("Alternate protein reference has no shared genes")
        if any(proteins[g][2] != alternate[g][2] for g in shared):
            raise ValueError("Normalized protein sequences disagree between source copies")
    annotations = {}
    if "annotations" in inputs:
        for row in table_rows(inputs["annotations"], ("prodigal_protein", "KEGG_ko")):
            gene = row["prodigal_protein"]
            if gene not in genes or gene in annotations:
                raise ValueError("Foreign or duplicate annotation gene")
            tokens = re.split(r"[,;\s]+", row["KEGG_ko"].strip())
            kos = set()
            for token in tokens:
                if token.lower() in ("", "-", "na", "nan", "none"):
                    continue
                token = token.removeprefix("ko:")
                if not re.fullmatch(r"K\d{5}", token):
                    raise ValueError(f"Invalid KO identifier: {token}")
                kos.add(token)
            annotations[gene] = kos
            # The benchmark annotation export stores TPM to six decimal places.
            for column, layer in (("metaG_tpm", "dna"), ("metaT_tpm", "rna")):
                if column in row:
                    observed = numeric(row[column], nan_is_missing=True)
                    raw = genes[gene][layer]
                    if (raw is None) != (observed is None) or (
                        raw is not None and abs(raw - observed) > 5.1e-7
                    ):
                        raise ValueError("Raw/annotated shared TPM values disagree")
    compounds = {}
    if "metabolites" in inputs:
        for row in table_rows(inputs["metabolites"], ("specimen", "compound", "value")):
            if row["specimen"] != obj["specimen"]:
                raise ValueError("Metabolite/source specimen mismatch")
            compound = row["compound"]
            if not re.fullmatch(r"C\d{5}", compound) or compound in compounds:
                raise ValueError("Invalid or duplicate compound identifier")
            compounds[compound] = numeric(row["value"], nan_is_missing=True)
    links, seen = [], set()
    if "compound_links" in inputs:
        for row in table_rows(inputs["compound_links"], ("compound", "reaction", "KO")):
            edge = tuple(row[k] for k in ("compound", "reaction", "KO"))
            if not all(
                re.fullmatch(prefix + r"\d{5}", value)
                for prefix, value in zip(("C", "R", "K"), edge)
            ):
                raise ValueError("Invalid compound/reaction/KO link")
            if edge in seen:
                raise ValueError("Duplicate compound/reaction/KO link")
            seen.add(edge)
            links.append(dict(zip(("compound", "reaction", "KO"), edge)))
    return MolecularSource(obj, inputs, genes, proteins, annotations, compounds, links)


def load_reference(manifest_path):
    """Load a v1/v2 source.json and verify its pinned files and reference contract.

    Relative input paths resolve beside the manifest. Complete scope requires
    matching TPM/protein gene populations; explicitly available scope retains
    exclusions for coverage reporting. Return a MolecularSource; reject changed
    files, inconsistent sequences or unsupported schemas before selection.
    Provider records supply biological specimen provenance.
    """
    path = Path(manifest_path).resolve()
    obj = json.loads(path.read_text())
    if obj.get("schema") == "fastalake.molecular-source.v1":
        obj, inputs, genes, proteins = load_source(path)
        return MolecularSource(obj, inputs, genes, proteins, {}, {}, [])
    if obj.get("schema") != "fastalake.molecular-source.v2":
        raise ValueError("Unsupported molecular source schema")
    return _load_v2(obj, path.parent)


def write_coverage(source, output):
    """Write missing-gene and measurement coverage to a caller-owned output directory.

    Preserve original TPM values for genes without sequences, and separately
    report proteins without TPM rows. This function does not select candidates.
    Callers create a fresh directory before writing these three report files.
    """
    output = Path(output)
    missing = sorted(set(source.genes) - set(source.proteins))
    write_table(
        output / "unavailable_genes.tsv",
        ["gene_id", "metaG_tpm", "metaT_tpm"],
        (
            dict(
                gene_id=g,
                metaG_tpm=source.genes[g]["raw_dna"],
                metaT_tpm=source.genes[g]["raw_rna"],
            )
            for g in missing
        ),
    )
    write_table(
        output / "proteins_without_tpm.tsv",
        ["gene_id", "protein_sha256"],
        (
            dict(gene_id=g, protein_sha256=source.proteins[g][2])
            for g in sorted(set(source.proteins) - set(source.genes))
        ),
    )
    (output / "REFERENCE_COVERAGE.json").write_text(json.dumps(source.coverage, indent=2) + "\n")


def prepare_source(
    output,
    *,
    specimen,
    reference_id,
    tpm,
    proteins,
    provenance,
    reference_scope="complete",
    alternate_proteins=None,
    annotations=None,
    metabolites=None,
    compound_links=None,
):
    """Package checked specimen-matched molecular inputs for portable reuse.

    Parameters
    ----------
    output : path-like
        New source directory; existing output is rejected.
    specimen, reference_id : str
        Explicit biological specimen and shared gene-reference identifiers.
    tpm, proteins, provenance : path-like
        Gene TPM table, encoded-protein FASTA and provider linkage record.
    reference_scope : {"complete", "available"}
        Complete requires matching TPM/protein gene sets. Available retains
        exclusions for a deliberately incomplete sequence reference.
    alternate_proteins : path-like or None
        Optional second FASTA for validating shared gene sequences.
    annotations, metabolites, compound_links : path-like or None
        Optional gene annotations, measured compounds and reaction/KO links.

    Returns
    -------
    dict
        Path to the written source manifest and sequence/measurement coverage counts.

    Notes
    -----
    Input consistency is checked before packaging. Provider records establish
    specimen provenance; filenames alone do not. Missing TPM is retained as
    missing, separately from genes whose protein sequence is unavailable.
    """
    output = Path(output)
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    obj = dict(
        schema="fastalake.molecular-source.v2",
        specimen=specimen,
        reference_id=reference_id,
        reference_scope=reference_scope,
    )
    supplied = dict(
        tpm=tpm,
        proteins=proteins,
        provenance=provenance,
        alternate_proteins=alternate_proteins,
        annotations=annotations,
        metabolites=metabolites,
        compound_links=compound_links,
    )
    for key, value in supplied.items():
        if value is not None:
            path = Path(value).resolve(strict=True)
            obj[key] = dict(path=str(path), sha256=sha256(path))
    source = _load_v2(obj, Path.cwd())
    output.mkdir(parents=True, exist_ok=False)
    for key, path in source.inputs.items():
        suffix = (
            ".faa"
            if key in ("proteins", "alternate_proteins")
            else (".json" if key == "provenance" else ".tsv")
        )
        name = key + suffix + (".gz" if path.suffix == ".gz" else "")
        shutil.copyfile(path, output / name)
        if sha256(output / name) != obj[key]["sha256"]:
            raise RuntimeError(f"Source changed while copying: {key}")
        obj[key]["path"] = name
    (output / "source.json").write_text(json.dumps(obj, indent=2) + "\n")
    write_coverage(source, output)
    write_reference(source, output / "reference.fasta")
    return dict(source_manifest=str(output / "source.json"), **source.coverage)


def write_reference(source, path):
    """Write every available exact sequence, including genes with zero/missing TPM."""
    path = Path(path)
    sequences = {source.proteins[g][2]: source.proteins[g][1] for g in source.available}
    with path.open("x") as stream:
        for digest, sequence in sorted(sequences.items()):
            stream.write(">FMO_SHA256_" + digest + "\n" + sequence + "\n")
    return {"sequences": len(sequences), "sha256": sha256(path)}
