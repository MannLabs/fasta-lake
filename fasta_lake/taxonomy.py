"""Join Unipept's own peptide taxonomy output to measured FastaLake evidence.

Taxonomy assignment stays with Unipept. This adapter preserves its lineage,
cutoffs and reference identity, and reuses its open-source visualisations.


Upstream taxonomy: https://github.com/unipept/unipept-cli; Mesuere et al.
(2016), https://doi.org/10.1093/bioinformatics/btw039.
Upstream visuals: https://github.com/unipept/unipept-visualizations;
Verschaffelt et al. (2022), https://doi.org/10.1093/bioinformatics/btab590.
The original JS bundle is retained under assets/ with its MIT/ISC notices.
FastaLake performs the checked evidence join and preserves internal-rank mass.
"""

from __future__ import annotations

import base64
import csv
import gzip
import html
import json
import math
import shutil
from collections import defaultdict
from pathlib import Path

from fasta_lake.networks import _identity, _sha, read_evidence_network

RANKS = (
    "domain",
    "realm",
    "kingdom",
    "subkingdom",
    "superphylum",
    "phylum",
    "subphylum",
    "superclass",
    "class",
    "subclass",
    "infraclass",
    "superorder",
    "order",
    "suborder",
    "infraorder",
    "parvorder",
    "superfamily",
    "family",
    "subfamily",
    "tribe",
    "subtribe",
    "genus",
    "subgenus",
    "species_group",
    "species_subgroup",
    "species",
    "subspecies",
    "strain",
    "varietas",
    "forma",
)
UNKNOWN = {
    "not_returned": (-1, "Not returned by Unipept"),
    "lookup_cutoff": (-2, "Lookup cutoff reached"),
    "cutoff_not_reported": (-3, "Lookup completeness unreported"),
    "unresolved": (-4, "Unresolved taxonomy"),
}
BALANCE = ("Host", "Bacteria", "Archaea", "Other assigned taxa", "Unresolved")


def canonical_peptide(value):
    """Use the study's I/L-equivalent peptide identity, without removing mass gaps."""
    if (
        not isinstance(value, str)
        or not value
        or any(c not in "ACDEFGHIKLMNPQRSTVWY" for c in value)
    ):
        raise ValueError(f"Expected an unmodified canonical amino-acid peptide: {value!r}")
    return value.replace("I", "L")


def read_unipept_output(path, peptides, *, max_rows=1_000_000, max_bytes=256 * 1024**2):
    """Read native pept2lca --all JSON/CSV; keep absent and truncated lookups visible.

    Caller must have used --equate for FastaLake's I/L-equivalent peptide keys.
    Rank assignments are never made more specific than Unipept's returned LCA.
    """
    source = Path(path)
    if source.stat().st_size > max_bytes:
        raise ValueError("Unipept output exceeds the explicit file-size limit")
    with source.open(encoding="utf-8-sig", newline="") as handle:
        if source.suffix.lower() == ".json":
            rows = json.load(handle)
            if not isinstance(rows, list):
                raise ValueError("Unipept JSON must contain a list of peptide results")
        else:
            reader = csv.DictReader(handle)
            if not reader.fieldnames or len(reader.fieldnames) != len(set(reader.fieldnames)):
                raise ValueError("Unipept CSV needs distinct column names")
            rows = []
            for raw in reader:
                if len(rows) >= max_rows:
                    raise ValueError("Unipept output exceeds the explicit row limit")
                rows.append(raw)
    if len(rows) > max_rows:
        raise ValueError("Unipept output exceeds the explicit row limit")
    records = {}
    for raw in rows:
        if not isinstance(raw, dict) or None in raw:
            raise ValueError("Malformed Unipept output row")
        needed = {
            "peptide",
            "taxon_id",
            "taxon_name",
            "taxon_rank",
            "genus_id",
            "species_id",
            "phylum_id",
        }
        if not needed <= raw.keys():
            raise ValueError("Unipept output needs pept2lca --all with full taxonomy lineage")
        peptide = canonical_peptide(raw["peptide"])
        if peptide not in peptides or peptide in records:
            raise ValueError("Foreign or duplicate I/L-equivalent peptide in Unipept output")
        cutoff = raw.get("cutoff_used", None)
        if cutoff in (True, 1, "1", "true", "True"):
            status = "lookup_cutoff"
        elif cutoff in (False, 0, "0", "false", "False", ""):
            status = "assigned"
        elif cutoff is None:
            status = "cutoff_not_reported"
        else:
            raise ValueError("Invalid Unipept cutoff_used flag")
        lineage = []
        lca_id = raw["taxon_id"]
        if lca_id in (None, "", 0, "0"):
            lca_id = None
            if status == "assigned":
                status = "unresolved"
        else:
            try:
                if isinstance(lca_id, (bool, float)):
                    raise ValueError("Taxonomy IDs must be integer identifiers")
                lca_id = int(lca_id)
            except (TypeError, ValueError) as error:
                raise ValueError("Invalid Unipept taxonomy ID") from error
            if lca_id < 1:
                raise ValueError("Invalid Unipept taxonomy ID")
        lca_seen = False
        lca_rank = "domain" if raw["taxon_rank"] == "superkingdom" else raw["taxon_rank"]
        for rank in RANKS:
            key = "superkingdom" if rank == "domain" and "domain_id" not in raw else rank
            taxid, name = raw.get(key + "_id"), raw.get(key + "_name")
            if taxid in (None, "", 0, "0"):
                continue
            if isinstance(taxid, (bool, float)):
                raise ValueError("Taxonomy IDs must be integer identifiers")
            taxid = int(taxid)
            if taxid < 1 or not isinstance(name, str) or not name.strip():
                raise ValueError("Unipept lineage IDs require corresponding nonempty names")
            _identity(name, "taxon name")
            if lca_seen or (lca_rank in RANKS and RANKS.index(rank) > RANKS.index(lca_rank)):
                raise ValueError("Unipept output contains descendants below its assigned LCA")
            if rank == lca_rank and taxid != lca_id:
                raise ValueError("Unipept LCA disagrees with its lineage")
            lineage.append(dict(id=taxid, name=name, rank=rank))
            lca_seen = taxid == lca_id
        if lca_id and lca_id != 1 and not lca_seen:
            name = _identity(raw["taxon_name"], "LCA name")
            rank = _identity(raw["taxon_rank"], "LCA rank")
            lineage.append(dict(id=lca_id, name=name, rank=rank))
        if lca_id == 1 and status == "assigned":
            status = "unresolved"
        records[peptide] = dict(peptide=peptide, status=status, lineage=lineage, raw=raw)
    for peptide in sorted(peptides - records.keys()):
        records[peptide] = dict(peptide=peptide, status="not_returned", lineage=[], raw={})
    return records


def build_taxonomy_tree(records, weights):
    """Make Unipept's count/selfCount tree with each peptide's weight counted once."""
    root = dict(id=1, name="All peptide evidence", count=0.0, selfCount=0.0, children=[])
    nodes = {(): root}
    for peptide, weight in sorted(weights.items()):
        if not math.isfinite(weight) or weight < 0:
            raise ValueError("Taxonomy weights must be finite and nonnegative")
        if not weight:
            continue
        record = records[peptide]
        lineage = record["lineage"] if record["status"] == "assigned" else []
        if not lineage:
            ident, name = UNKNOWN[record["status"] if record["status"] in UNKNOWN else "unresolved"]
            lineage = [dict(id=ident, name=name, rank="unresolved")]
        root["count"] += weight
        parent, key = root, ()
        for taxon in lineage:
            key = (*key, taxon["id"])
            if key not in nodes:
                nodes[key] = dict(**taxon, count=0.0, selfCount=0.0, children=[])
                parent["children"].append(nodes[key])
            node = nodes[key]
            if node["name"] != taxon["name"]:
                raise ValueError("Conflicting names for the same taxonomy lineage")
            node["count"] += weight
            parent = node
        parent["selfCount"] += weight
    for node in nodes.values():
        node["children"].sort(key=lambda child: (-child["count"], child["name"]))
        if not math.isclose(
            node["count"],
            node["selfCount"] + math.fsum(child["count"] for child in node["children"]),
            rel_tol=1e-10,
            abs_tol=1e-8,
        ):
            raise ValueError("Taxonomy tree does not conserve peptide weights")
    return root


def _balance(record, host_taxon_id):
    """Classify an assigned lineage as host, bacteria, archaea or other; retain unresolved."""
    if record["status"] != "assigned":
        return "Unresolved"
    ancestors = {node["id"] for node in record["lineage"]}
    if host_taxon_id in ancestors:
        return "Host"
    if 2 in ancestors:
        return "Bacteria"
    if 2157 in ancestors:
        return "Archaea"
    return "Other assigned taxa"


def prepare_unipept_queries(groups, output):
    """Export the checked study peptide roster for the upstream Unipept CLI."""
    out = Path(output)
    if out.exists():
        raise FileExistsError(f"Output exists: {out}")
    network = read_evidence_network(groups)
    peptides = sorted(n["label"] for n in network["nodes"] if n["node_type"] == "peptide")
    out.mkdir(parents=True)
    (out / "peptides.txt").write_text("".join(p + "\n" for p in peptides))
    record = dict(
        status="PASS",
        graph_id=network["graph_id"],
        peptides=len(peptides),
        peptide_sha256=_sha(out / "peptides.txt"),
        required_settings="pept2lca --all --equate",
        note="Running the Unipept CLI contacts its selected server; default is public UniProt.",
    )
    (out / "COMPLETE.json").write_text(json.dumps(record, indent=2) + "\n")
    return record


def taxonomy_report(
    groups,
    unipept_output,
    output,
    *,
    database,
    equate_il,
    host_taxon_id=9606,
    contaminant_peptides=None,
):
    """Create local taxonomy/holobiont plots from Unipept output, without web requests.

    Counts are distinct accepted canonical peptides per acquisition. Signal is
    summed assigned LFQ-feature intensity, independent of directLFQ reporting-
    group estimates. Host/taxon fractions are signal shares, not cell biomass.
    Contaminant reference overlap is reported separately from host identity.
    """
    if not equate_il:
        raise ValueError(
            "Use Unipept --equate and declare equate_il=True for study peptide identities"
        )
    database = _identity(database, "Unipept database/version label")
    if type(host_taxon_id) is not int or host_taxon_id < 1:
        raise ValueError("Declare a positive host NCBI taxonomy ID")
    out, source = Path(output), Path(groups).resolve(strict=True)
    if out.exists() or out.is_symlink():
        raise FileExistsError(f"Output exists: {out}; choose a new directory")
    network = read_evidence_network(source)
    paths = [source / n for n in network["source_files"]]
    paths += [source / "feature_assignments.tsv.gz", Path(unipept_output).resolve(strict=True)]
    if contaminant_peptides:
        paths.append(Path(contaminant_peptides).resolve(strict=True))
    identities = {str(p): _sha(p) for p in paths}
    if any(
        identities[str(source / name)] != digest for name, digest in network["source_files"].items()
    ):
        raise ValueError("A grouping input changed while preparing taxonomy")
    peptides = {n["label"] for n in network["nodes"] if n["node_type"] == "peptide"}
    assigned = {
        n["label"]: n["assigned_group"] for n in network["nodes"] if n["node_type"] == "peptide"
    }
    records = read_unipept_output(unipept_output, peptides)
    samples = network["source_summary"]["acquisition_names"]
    observed = {s: {} for s in samples}
    with gzip.open(source / "peptide_evidence.tsv.gz", "rt") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            observed[row["sample"]][row["canonical_peptide"]] = 1.0
    quantities = defaultdict(list)
    seen = set()
    with gzip.open(source / "feature_assignments.tsv.gz", "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"sample", "input_line", "canonical_peptide", "intensity", "study_group"}
        if not required <= set(reader.fieldnames or []):
            raise ValueError("LFQ feature ledger lacks the required identity/intensity fields")
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise ValueError("Malformed LFQ feature ledger row")
            sample, peptide = row["sample"], row["canonical_peptide"]
            key = (sample, row["input_line"])
            if sample not in observed or peptide not in observed[sample] or key in seen:
                raise ValueError("LFQ feature ledger has foreign or duplicate evidence")
            if row["study_group"] != assigned[peptide]:
                raise ValueError("LFQ feature assignment disagrees with the study dictionary")
            seen.add(key)
            value = float(row["intensity"])
            if not math.isfinite(value) or value < 0:
                raise ValueError("LFQ feature ledger has invalid intensity")
            quantities[(sample, peptide)].append(value)
    if len(seen) != network["source_summary"]["lfq_features"]:
        raise ValueError("LFQ ledger feature count disagrees with the grouping summary")
    total = math.fsum(v for values in quantities.values() for v in values)
    if not math.isclose(
        total, network["source_summary"]["intensity_out"], rel_tol=1e-9, abs_tol=1e-6
    ):
        raise ValueError("LFQ signal does not reconcile to the grouping summary")
    contaminants = set()
    if contaminant_peptides:
        contaminants = {
            canonical_peptide(line.strip())
            for line in Path(contaminant_peptides).read_text().splitlines()
            if line.strip()
        }
    reports, table, balances = [], [], []
    for sample in samples:
        weights = {p: math.fsum(quantities[(sample, p)]) for p in observed[sample]}
        totals = dict(peptides=len(observed[sample]), signal=math.fsum(weights.values()))
        trees = dict(
            peptides=build_taxonomy_tree(records, observed[sample]),
            signal=build_taxonomy_tree(records, weights),
        )
        for peptide in observed[sample]:
            record = records[peptide]
            table.append(
                dict(
                    sample=sample,
                    peptide=peptide,
                    status=record["status"],
                    taxon_id=record["raw"].get("taxon_id", ""),
                    taxon_name=record["raw"].get("taxon_name", ""),
                    taxon_rank=record["raw"].get("taxon_rank", ""),
                    lfq_sum_intensity=weights[peptide],
                    contaminant_reference_overlap=peptide in contaminants,
                )
            )
        balance = {}
        for category in BALANCE:
            keys = [p for p in observed[sample] if _balance(records[p], host_taxon_id) == category]
            count, signal = len(keys), math.fsum(weights[p] for p in keys)
            balance[category] = dict(peptides=count, signal=signal)
            balances.append(
                dict(
                    sample=sample,
                    category=category,
                    peptides=count,
                    signal=signal,
                    peptide_fraction=count / totals["peptides"] if totals["peptides"] else "",
                    signal_fraction=signal / totals["signal"] if totals["signal"] else "",
                )
            )
        overlap = dict(
            peptides=len(set(observed[sample]) & contaminants),
            signal=math.fsum(weights[p] for p in weights if p in contaminants),
        )
        reports.append(
            dict(sample=sample, totals=totals, trees=trees, balance=balance, overlap=overlap)
        )
    if identities != {str(p): _sha(p) for p in paths}:
        raise ValueError("Taxonomy report inputs changed during processing")
    out.mkdir(parents=True)
    shutil.copyfile(unipept_output, out / ("unipept_original" + Path(unipept_output).suffix))
    from fasta_lake.multi_omics.exact import write_table

    write_table(
        out / "peptide_taxonomy.tsv",
        [
            "sample",
            "peptide",
            "status",
            "taxon_id",
            "taxon_name",
            "taxon_rank",
            "lfq_sum_intensity",
            "contaminant_reference_overlap",
        ],
        table,
    )
    write_table(
        out / "holobiont_balance.tsv",
        ["sample", "category", "peptides", "signal", "peptide_fraction", "signal_fraction"],
        balances,
    )
    payload = dict(
        database=database,
        samples=reports,
        host_taxon_id=host_taxon_id,
        contaminant_reference_supplied=bool(contaminant_peptides),
        metrics={
            "peptides": "Distinct accepted peptides",
            "signal": "Summed LFQ-feature intensity (a.u.)",
        },
    )
    from fasta_lake.community_plots import plot_community_context

    static_plots = plot_community_context(payload, out)
    payload["static_plots"] = static_plots
    (out / "taxonomy.json").write_text(json.dumps(payload, indent=2) + "\n")
    _write_report(out, payload)
    receipt = dict(
        status="PASS",
        graph_id=network["graph_id"],
        database_label=database,
        database_label_source="User-declared; native Unipept output preserved",
        equate_il=True,
        host_taxon_id=host_taxon_id,
        inputs=identities,
        source_implementation=_sha(__file__),
        plotting_implementation=_sha(Path(__file__).with_name("community_plots.py")),
        template_implementation=_sha(Path(__file__).parent / "assets/taxonomy_report.html"),
        scope="Unipept reference-dependent LCA; peptide counts and summed LFQ signal; "
        "not organism biomass or directLFQ taxon abundance",
        visualizations="Unipept Visualizations 3.0.0 (MIT)",
        static_plots=static_plots,
        files={p.name: _sha(p) for p in out.iterdir() if p.is_file()},
    )
    (out / "COMPLETE.json").write_text(json.dumps(receipt, indent=2) + "\n")
    return receipt


def _write_report(out, payload):
    """Embed upstream Unipept visualisations and data for an offline HTML report."""
    assets = Path(__file__).with_name("assets")
    module = base64.b64encode((assets / "unipept-visualizations.js").read_bytes()).decode()
    for name in ("UNIPEPT_LICENSE.md", "D3_LICENSE.md", "UNIPEPT_SOURCE.json"):
        shutil.copyfile(assets / name, out / name)
    data = json.dumps(payload).replace("<", "\\u003c").replace("&", "\\u0026")
    template = Path(__file__).with_name("assets") / "taxonomy_report.html"
    document = (
        template.read_text()
        .replace("__UNIPEPT_MODULE__", module)
        .replace("__TAXONOMY_DATA__", data)
    )
    document = document.replace("__DATABASE_LABEL__", html.escape(payload["database"]))
    (out / "index.html").write_text(document)
