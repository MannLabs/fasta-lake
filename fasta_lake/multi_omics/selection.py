"""Gene-level evidence rules followed by an auditable exact-sequence FASTA union."""

from __future__ import annotations

import gzip
import json
import math
from pathlib import Path

from fasta_lake.hashing import HASH_VERSION

from .exact import read_fasta, sha256, top_cutoff, write_table
from .sources import load_reference, write_coverage


def validate_rules(
    dna_percent=0, rna_percent=0, dna_absolute=None, rna_absolute=None, rounding="floor"
):
    """Validate independent DNA/RNA selection modes without loading any inputs.

    A layer uses a percentage of positive genes or a finite nonnegative absolute
    TPM threshold. Set its percentage to zero when using an absolute threshold.
    Invalid ranges, rounding rules and conflicting modes raise ValueError.
    """
    for layer, percent, absolute in (
        ("DNA", dna_percent, dna_absolute),
        ("RNA", rna_percent, rna_absolute),
    ):
        top_cutoff([], percent, rounding)
        if absolute is not None:
            if percent != 0:
                raise ValueError(f"Choose {layer} percent or absolute mode; set percent=0")
            if not math.isfinite(absolute) or absolute < 0:
                raise ValueError("Absolute TPM thresholds must be finite and nonnegative")


def _targets(path, *, reject_duplicates):
    """Read target-only proteins and index exact sequences by SHA-256.

    Reject decoy-prefixed IDs and, when requested, duplicate baseline sequences.
    Return the accession records and one accession for each sequence hash.
    """
    records = read_fasta(path)
    hashes = {}
    for accession, (_, _, digest) in records.items():
        if accession.startswith(("rev_", "decoy_")):
            raise ValueError("Inputs must contain targets only; Sage generates new decoys")
        if digest in hashes and reject_duplicates:
            raise ValueError("Baseline contains duplicate exact protein sequences")
        hashes.setdefault(digest, accession)
    return records, hashes


def build_selection(
    source_manifest,
    specimen,
    output,
    *,
    baseline=None,
    background=None,
    dna_percent=0,
    rna_percent=0,
    dna_absolute=None,
    rna_absolute=None,
    rounding="floor",
    metabolites=False,
):
    """Select gene-supported sequences and write a traceable target FASTA union.

    ``source_manifest`` must match ``specimen``. Selection operates on available
    genes before exact-sequence deduplication. DNA/RNA percentages rank positive
    values independently; absolute thresholds use strict greater-than tests.
    Missing values never supply support. Optional metabolites require the source's
    compound, reaction/KO and gene annotation tables.

    Supply the de novo ``baseline`` for the standard additive workflow. Its bytes,
    sequences and headers are retained; new targets are appended deterministically.
    Omitting it creates a standalone molecular comparator. ``output`` must be new.
    Return the same counts, rules and hashes written to COMPLETE.json, alongside
    search.fasta, gene evidence and reference-coverage ledgers. A completed union
    is a candidate-selection result, not a peptide identification or FDR estimate.
    """
    output = Path(output)
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    validate_rules(dna_percent, rna_percent, dna_absolute, rna_absolute, rounding)
    source = load_reference(source_manifest)
    if source.manifest["specimen"] != specimen:
        raise ValueError("Acquisition/source specimen mismatch")
    if metabolites and not {"annotations", "metabolites", "compound_links"} <= set(source.inputs):
        raise ValueError(
            "Metabolite selection requires annotations, metabolites and compound_links"
        )
    originals, baseline_hashes = (
        _targets(baseline, reject_duplicates=True) if baseline else ({}, {})
    )
    background_records, background_hashes = (
        _targets(background, reject_duplicates=False) if background else ({}, {})
    )
    for accession, record in background_records.items():
        if accession in originals and record[2] != originals[accession][2]:
            raise ValueError("Background accession collides with a different baseline sequence")
    available = source.available
    cutoffs = {
        layer: top_cutoff((source.genes[g][layer] for g in available), percent, rounding)
        for layer, percent in (("dna", dna_percent), ("rna", rna_percent))
    }
    active_links = [
        r
        for r in source.links
        if source.compounds.get(r["compound"]) is not None and source.compounds[r["compound"]] > 0
    ]
    active_kos = {r["KO"] for r in active_links} if metabolites else set()
    retained_hashes = dict(background_hashes)
    retained_hashes.update(baseline_hashes)
    retained_accessions = {**background_records, **originals}
    selected, ledger = {}, []
    for gene in available:
        row = source.genes[gene]
        header, sequence, digest = source.proteins[gene]
        if gene.startswith(("rev_", "decoy_")):
            raise ValueError("Reference must contain targets only")
        passes = {}
        for layer, absolute in (("dna", dna_absolute), ("rna", rna_absolute)):
            value, cutoff = row[layer], cutoffs[layer]
            passes[layer] = value is not None and (
                value > absolute
                if absolute is not None
                else cutoff is not None and value > 0 and value >= cutoff
            )
        matched_kos = source.annotations.get(gene, set()) & active_kos
        metab_pass = bool(matched_kos)
        if (passes["dna"] or passes["rna"] or metab_pass) and digest not in retained_hashes:
            accession = "FMO_SHA256_" + digest
            if accession in retained_accessions and retained_accessions[accession][2] != digest:
                raise ValueError("Molecular accession collides with a different retained sequence")
            selected[digest] = (accession, sequence)
        ledger.append(
            dict(
                specimen=specimen,
                reference_id=source.manifest["reference_id"],
                gene_id=gene,
                original_header=header,
                protein_sha256=digest,
                metaG_tpm=row["raw_dna"],
                metaT_tpm=row["raw_rna"],
                dna_pass=int(passes["dna"]),
                rna_pass=int(passes["rna"]),
                metabolite_pass=int(metab_pass),
                selected_kos=";".join(sorted(matched_kos)),
                in_baseline=int(digest in baseline_hashes),
                in_background=int(digest in background_hashes),
            )
        )
    for row in ledger:
        digest = row["protein_sha256"]
        row["included"] = int(digest in retained_hashes or digest in selected)
        row["final_accession"] = (
            retained_hashes.get(digest, "FMO_SHA256_" + digest) if row["included"] else ""
        )
    expected = set(retained_hashes) | set(selected)
    if not expected:
        raise ValueError(
            "Selection contains no target sequences; inspect cutoffs or supply background"
        )
    output.mkdir(parents=True, exist_ok=False)
    target = output / "search.fasta"
    if baseline:
        baseline = Path(baseline)
        if baseline.suffix == ".gz":
            with gzip.open(baseline, "rb") as stream:
                raw = stream.read()
        else:
            raw = baseline.read_bytes()
    else:
        raw = b""
    appended = []
    for digest in sorted(set(background_hashes) - set(baseline_hashes)):
        accession = background_hashes[digest]
        header, sequence, _ = background_records[accession]
        appended.append((header, sequence))
    appended.extend(selected[d] for d in sorted(selected))
    with target.open("wb") as stream:
        stream.write(raw)
        if raw and appended and not raw.endswith(b"\n"):
            stream.write(b"\n")
        for header, sequence in appended:
            stream.write((">" + header + "\n" + sequence + "\n").encode("utf-8"))
    final = read_fasta(target)
    if {r[2] for r in final.values()} != expected or len(final) != len(expected):
        raise RuntimeError("Actual FASTA does not equal the complete exact-sequence union")
    if any(final[a] != record for a, record in originals.items()):
        raise RuntimeError("Union changed a baseline sequence/header")
    if baseline and not appended and target.read_bytes() != raw:
        raise RuntimeError("No-addition FASTA differs from baseline bytes")
    write_table(output / "gene_evidence.tsv", list(ledger[0]), ledger)
    write_coverage(source, output)
    if metabolites:
        write_table(output / "metabolite_links.tsv", ["compound", "reaction", "KO"], active_links)
        write_table(
            output / "compound_evidence.tsv",
            ["specimen", "compound", "value", "positive", "KOs"],
            (
                dict(
                    specimen=specimen,
                    compound=c,
                    value="" if value is None else value,
                    positive=int(value is not None and value > 0),
                    KOs=";".join(sorted({r["KO"] for r in active_links if r["compound"] == c})),
                )
                for c, value in sorted(source.compounds.items())
            ),
        )
    receipt = dict(
        schema="fastalake.molecular-union.v2",
        hash_version=HASH_VERSION,
        specimen=specimen,
        reference_id=source.manifest["reference_id"],
        reference_coverage=source.coverage,
        baseline_sequences=len(baseline_hashes),
        background_sequences=len(background_hashes),
        additions=len(expected - set(baseline_hashes)),
        molecular_additions=len(selected),
        final_sequences=len(final),
        genes=len(available),
        dna_measured_genes=sum(source.genes[g]["dna"] is not None for g in available),
        rna_measured_genes=sum(source.genes[g]["rna"] is not None for g in available),
        selected_genes=sum(
            bool(r["dna_pass"] or r["rna_pass"] or r["metabolite_pass"]) for r in ledger
        ),
        dna_cutoff=cutoffs["dna"],
        rna_cutoff=cutoffs["rna"],
        dna_percent=dna_percent,
        rna_percent=rna_percent,
        dna_absolute=dna_absolute,
        rna_absolute=rna_absolute,
        rank_rounding=rounding,
        metabolites_enabled=metabolites,
        active_metabolite_kos=len(active_kos),
        selection_code_sha256=sha256(__file__),
        source_code_sha256=sha256(Path(__file__).with_name("sources.py")),
        input_parser_code_sha256=sha256(Path(__file__).with_name("exact.py")),
        source_manifest_sha256=sha256(source_manifest),
        baseline_sha256=sha256(baseline) if baseline else None,
        background_sha256=sha256(background) if background else None,
        fasta_sha256=sha256(target),
        gene_evidence_sha256=sha256(output / "gene_evidence.tsv"),
        input_sha256={k: sha256(v) for k, v in source.inputs.items()},
        output_sha256={p.name: sha256(p) for p in sorted(output.iterdir()) if p.is_file()},
    )
    # No success receipt is emitted if any pinned source changed during selection.
    if any(receipt["input_sha256"][key] != source.manifest[key]["sha256"] for key in source.inputs):
        raise RuntimeError("Molecular inputs changed during selection")
    (output / "COMPLETE.json").write_text(json.dumps(receipt, indent=2, allow_nan=False) + "\n")
    return receipt
