"""Render a bounded peptide/protein view through an explicitly selected CKG install.

This worker runs in CKG's Python environment. It reads the portable FastaLake
export and uses CKG analytics/viz functions; no Neo4j connection is needed.


CKG source: https://github.com/MannLabs/CKG; Santos et al. (2022),
https://doi.org/10.1038/s41587-021-01145-6. Calls its analytics and viz modules
with NetworkX (https://github.com/networkx/networkx), Plotly
(https://github.com/plotly/plotly.py) and PyVis (https://github.com/WestHealth/pyvis).
FastaLake supplies checked graph data and the bounded component layout.
The optional compatibility changes are recorded in integrations/ckg/.
"""

from __future__ import annotations

import argparse
import html
import importlib.metadata
import inspect
import json
import math
import re
import sys
from collections import Counter
from pathlib import Path

if __package__:
    from .networks import (
        _sha,
        add_function_annotations,
        load_evidence_network,
        load_network_annotations,
        select_network_view,
    )
else:
    from networks import (
        _sha,
        add_function_annotations,
        load_evidence_network,
        load_network_annotations,
        select_network_view,
    )


def remove_unused_bootstrap(document):
    """Remove PyVis's unused remote Bootstrap assets from its inline template."""
    document = re.sub(
        r'<link\b[^>]*href="https://cdn\.jsdelivr\.net/npm/bootstrap@[^\"]*"[^>]*>',
        "",
        document,
        flags=re.S,
    )
    document = re.sub(
        r'<script\b[^>]*src="https://cdn\.jsdelivr\.net/npm/bootstrap@[^\"]*"[^>]*>\s*</script>',
        "",
        document,
        flags=re.S,
    )
    if re.search(r"<(?:script|link)\b[^>]*(?:src|href)=[\"\'](?:https?:)?//", document, re.I):
        raise ValueError("CKG/PyVis output still loads external scripts or styles")
    return document


def _json_for_script(value):
    """Serialize data safely for a JavaScript string literal in an HTML script."""
    return (
        json.dumps(value, ensure_ascii=False)
        .replace("<", "\\u003c")
        .replace(">", "\\u003e")
        .replace("&", "\\u0026")
        .replace("\u2028", "\\u2028")
        .replace("\u2029", "\\u2029")
    )


def component_layout(graph):
    """Keep disconnected membership components readable in a deterministic grid."""
    import networkx as nx

    components = sorted(nx.connected_components(graph), key=lambda c: (-len(c), sorted(c)))
    radii = [max(0.45, math.sqrt(len(c)) * 0.23) for c in components]
    row_width = max(3, math.sqrt(sum((2 * r + 0.5) ** 2 for r in radii)) * 1.3)
    positions = {}
    x_offset, y_offset, row_height = 0, 0, 0
    for component, radius in zip(components, radii):
        width = 2 * radius + 0.5
        if x_offset and x_offset + width > row_width:
            x_offset = 0
            y_offset += row_height
            row_height = 0
        part = graph.subgraph(sorted(component)).copy()
        local = nx.spring_layout(part, seed=42, weight=None, scale=radius)
        for node, (x, y) in local.items():
            positions[node] = (float(x_offset + radius + x), float(-y_offset - radius + y))
        x_offset += width
        row_height = max(row_height, width)
    return positions


def render(bundle, output, ckg_source=None, max_nodes=150, max_edges=300, annotations=None):
    """Plot reported membership and topological communities without re-inferring proteins."""
    out = Path(output)
    if out.exists() or out.is_symlink():
        raise FileExistsError(f"CKG output exists: {out}; choose a new directory")
    if ckg_source:
        source = Path(ckg_source).resolve(strict=True)
        if not (source / "ckg/analytics_core/viz/viz.py").is_file():
            raise ValueError("CKG source must contain ckg/analytics_core/viz/viz.py")
        sys.path.insert(0, str(source))
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import networkx as nx
    import pandas as pd
    import plotly.graph_objects as go
    import plotly.io as pio
    from ckg.analytics_core.analytics import analytics
    from ckg.analytics_core.viz import viz

    network = load_evidence_network(bundle)
    layer = load_network_annotations(annotations, network["graph_id"]) if annotations else None
    reserve = (
        min(20, max_nodes // 3, sum(n["namespace"] in {"KEGG_ko", "EC"} for n in layer["nodes"]))
        if layer
        else 0
    )
    view = select_network_view(
        network,
        max_nodes=max(2, max_nodes - reserve),
        max_edges=max(1, max_edges - reserve),
        preferred_proteins={edge["source"] for edge in layer["edges"]} if layer else (),
    )
    if layer:
        view = add_function_annotations(view, layer, max_nodes=max_nodes, max_edges=max_edges)
    graph = nx.Graph()
    for node in view["nodes"]:
        peptide = node["node_type"] == "peptide"
        function = node["node_type"] == "annotated_function"
        # Tooltips interpret HTML. Escape labels; canonical IDs never contain raw labels.
        title = html.escape(node["label"]) + "<br>" + html.escape(node["node_type"])
        if not function:
            title += f"<br>Reported memberships in full export: {node['reported_degree']}"
        if "annotation_status" in node:
            title += "<br>eggNOG: " + html.escape(node["annotation_status"])
        if function:
            title += "<br>Namespace: " + html.escape(node["namespace"])
        graph.add_node(
            node["id"],
            label=node["label"],
            title=title,
            node_type=node["node_type"],
            color="#7562a7" if function else "#246b82" if peptide else "#b66b37",
            shape="diamond" if function else "dot" if peptide else "square",
            size=9 if peptide else 14,
        )
    for edge in view["edges"]:
        annotation_edge = edge.get("relation") == "orthology_annotation"
        assigned = edge.get("assigned_representative", False)
        graph.add_edge(
            edge["source"],
            edge["target"],
            width=2.0 if assigned else 1.0,
            dashes=not assigned,
            color="#7562a7" if annotation_edge else "#64748b" if assigned else "#bac3cc",
            relation=edge.get("relation", "reported_membership"),
            assigned_representative=assigned,
            observed_acquisitions=edge.get("observed_acquisitions", 0),
            title="eggNOG orthology annotation; not measured activity"
            if annotation_edge
            else f"Reported in {edge['observed_acquisitions']} acquisition(s); "
            + ("assigned representative" if assigned else "other reported member"),
        )
    communities = (
        analytics.get_network_communities(
            graph, {"communities_algorithm": "greedy_modularity", "values": None}
        )
        if graph.number_of_edges()
        else {key: n for n, key in enumerate(graph)}
    )
    if set(communities) != set(graph):
        raise ValueError("CKG communities lost or added graph nodes")
    nx.set_node_attributes(graph, communities, "community")
    nodes, edges = viz.network_to_tables(graph, "source", "target")
    if len(nodes) != len(view["nodes"]) or len(edges) != len(view["edges"]):
        raise ValueError("CKG graph/table round-trip counts disagree")
    out.mkdir(parents=True, exist_ok=False)
    nodes.to_csv(out / "view_nodes.tsv", sep="\t", index=False)
    edges.to_csv(out / "view_edges.tsv", sep="\t", index=False)
    (out / "VIEW.json").write_text(json.dumps(view, indent=2) + "\n")

    counts = Counter(n["node_type"] for n in network["nodes"])
    bars = pd.DataFrame(
        {
            "Evidence": ["Canonical peptides", "Reported protein accessions", "Membership edges"],
            "Count": [counts["peptide"], counts["reported_protein"], len(network["edges"])],
        }
    )
    chart = viz.get_barplot(bars, "full-evidence-counts", {"x": "Evidence", "y": "Count"})
    figure = go.Figure(chart.figure)
    figure.update_layout(title="Full exported evidence", template="plotly_white")
    pio.write_html(figure, out / "evidence_counts.html", include_plotlyjs=True, auto_open=False)
    if graph:
        pyvis = viz.get_notebook_network_pyvis(
            graph, {"width": "100%", "height": "720px", "cdn_resources": "in_line"}
        )
        # A deterministic initial layout gives an immediate readable view.
        positions = component_layout(graph)
        for node in pyvis.nodes:
            x, y = positions[node["id"]]
            # Canvas y increases downward; preserve the static plot's component order.
            node.update(x=float(x * 140), y=float(-y * 140))
        # CKG enables the PyVis configuration editor for notebook exploration.
        # A report starts with a calm overview; complete labels remain in hover
        # details and tables, with an explicit label toggle for closer inspection.
        labels = {}
        for node in pyvis.nodes:
            full_label = node["label"]
            labels[node["id"]] = (
                full_label if len(full_label) <= 26 else full_label[:15] + "…" + full_label[-7:]
            )
            if graph.nodes[node["id"]]["node_type"] != "annotated_function":
                node["label"] = ""
        pyvis.conf = False  # PyVis 0.3.2 template's notebook configuration panel.
        pyvis.set_options(
            json.dumps(
                {
                    "configure": {"enabled": False},
                    "physics": {"enabled": False},
                    "interaction": {"hover": True, "tooltipDelay": 100},
                    "nodes": {"font": {"size": 18, "face": "Arial", "color": "#213547"}},
                    "edges": {"smooth": False},
                }
            )
        )
        document = remove_unused_bootstrap(pyvis.generate_html())
        original_labels = {node["id"]: node["label"] for node in pyvis.nodes}
        controls = """<style>
body{font:16px/1.5 system-ui;color:#213547;margin:24px auto;max-width:1280px;padding:0 20px}
h1{font-size:24px;margin-bottom:8px}.fl-controls{padding:12px 0}
button{background:white;color:#246b82;border:1px solid #b6c8d0;border-radius:5px;
padding:8px 14px;font:inherit;cursor:pointer}.fl-legend{color:#475569;margin:8px 0}
#mynetwork{border:1px solid #dde5e9!important;border-radius:8px;float:none!important}
</style><h1>Peptide–protein evidence</h1>
<p>Hover for full identities; drag or zoom to explore. Labels can be shown when needed.</p>
<p class="fl-legend">Blue circles: peptides · Orange squares: reported proteins ·
Purple diamonds: orthology annotations. Solid edges mark assigned representatives;
dashed edges show other reported memberships. These links are not protein interactions.</p>
<div class="fl-controls"><button id="fl-labels" aria-pressed="false">Show labels</button>
<button id="fl-fit">Fit network</button></div>"""
        document = document.replace("<body>", "<body>" + controls, 1)
        controls_js = (
            "<script>\nconst flLabels = "
            + _json_for_script(labels)
            + ", flOverview = "
            + _json_for_script(original_labels)
            + ";\n"
            + """
const flButton = document.getElementById('fl-labels');
flButton.addEventListener('click', () => {
 const show = flButton.getAttribute('aria-pressed') !== 'true';
 const labels = show ? flLabels : flOverview;
 nodes.update(Object.entries(labels).map(([id, label]) => ({id, label})));
 flButton.setAttribute('aria-pressed', String(show));
 flButton.textContent = show ? 'Hide labels' : 'Show labels';
});
document.getElementById('fl-fit').addEventListener('click', () => network.fit());
</script>"""
        )
        document = document.replace("</body>", controls_js + "</body>", 1)
    else:
        positions = {}
        document = (
            '<!doctype html><meta charset="utf-8">'
            "<p>No accepted peptide memberships are available.</p>"
        )
    # User labels remain data inside JS; prevent a label from closing a script tag.
    for node in view["nodes"]:
        if "</script" in node["label"].lower():
            raise ValueError("CKG view cannot safely render this accession label")
    credits = (
        '<p style="font:13px system-ui;color:#475569;padding:12px">'
        "Clinical Knowledge Graph: Santos et al. (2022), "
        '<a href="https://doi.org/10.1038/s41587-021-01145-6">'
        "doi:10.1038/s41587-021-01145-6</a>. "
        "Evidence integration: FastaLake."
    )
    if layer:
        credits += (
            " eggNOG-mapper v2: Cantalapiedra et al. (2021), "
            '<a href="https://doi.org/10.1093/molbev/msab293">'
            "doi:10.1093/molbev/msab293</a>."
        )
    credits += "</p>"
    document = (
        document.replace("</body>", credits + "</body>")
        if "</body>" in document
        else (document + credits)
    )
    (out / "network.html").write_text(document)
    fig, ax = plt.subplots(figsize=(12, 8))
    if graph:
        for kind, color, shape in [
            ("peptide", "#246b82", "o"),
            ("reported_protein", "#b66b37", "s"),
            ("annotated_function", "#7562a7", "D"),
        ]:
            selected = [n for n, d in graph.nodes(data=True) if d["node_type"] == kind]
            nx.draw_networkx_nodes(
                graph,
                positions,
                nodelist=selected,
                node_color=color,
                node_shape=shape,
                node_size=28 if kind == "peptide" else 70,
                ax=ax,
                label={
                    "peptide": "Peptides",
                    "reported_protein": "Reported proteins",
                    "annotated_function": "eggNOG annotations",
                }[kind],
            )
        function_labels = {
            n: d["label"]
            for n, d in graph.nodes(data=True)
            if d["node_type"] == "annotated_function"
        }
        nx.draw_networkx_labels(
            graph,
            positions,
            function_labels,
            font_size=7,
            ax=ax,
            bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "none"},
        )
        for assigned in (False, True):
            pairs = [
                (a, b)
                for a, b, d in graph.edges(data=True)
                if d["assigned_representative"] == assigned
                and d["relation"] != "orthology_annotation"
            ]
            nx.draw_networkx_edges(
                graph,
                positions,
                edgelist=pairs,
                style="solid" if assigned else "dashed",
                width=1 if assigned else 0.65,
                alpha=0.65,
                edge_color="#64748b",
                ax=ax,
            )
        nx.draw_networkx_edges(
            graph,
            positions,
            ax=ax,
            edge_color="#7562a7",
            style="dotted",
            edgelist=[
                (a, b)
                for a, b, d in graph.edges(data=True)
                if d["relation"] == "orthology_annotation"
            ],
            alpha=0.65,
        )
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.02), ncol=3, frameon=False)
        ax.margins(0.1)
    else:
        ax.text(0.5, 0.5, "No accepted peptide memberships", ha="center", transform=ax.transAxes)
    ax.set_title(
        f"Peptide–protein evidence: {len(view['nodes'])}/{view['total_nodes']} nodes, "
        f"{len(view['edges'])}/{view['total_edges']} edges",
        loc="left",
        pad=20,
    )
    ax.set_axis_off()
    fig.text(
        0.06,
        0.035,
        "Solid: assigned reporting representative. Dashed: other reported membership.\n"
        "Purple dotted: orthology annotation. Full evidence remains in the exported tables.",
        fontsize=10,
        color="#475569",
    )
    citation = "CKG: Santos et al. (2022), doi:10.1038/s41587-021-01145-6."
    if layer:
        citation += "\neggNOG-mapper: Cantalapiedra et al. (2021), doi:10.1093/molbev/msab293."
    fig.text(0.06, -0.02, citation, fontsize=8, color="#475569")
    for ext in ("png", "pdf", "svg"):
        fig.savefig(out / f"network.{ext}", dpi=180, bbox_inches="tight")
    plt.close(fig)
    scope = html.escape(network["scope"])
    annotation_note = ""
    if layer:
        import shutil

        for name in ["protein_annotations.tsv", "function_nodes.tsv", "function_edges.tsv"]:
            shutil.copyfile(Path(annotations) / name, out / name)
        annotation_note = (
            "<h2>Functional annotations</h2><p>"
            + html.escape(layer["scope"])
            + "</p><p>Purple diamonds and dotted edges show eggNOG KO/EC assignments. "
            "Additional namespaces remain in the complete annotation export. "
            "No-hit and unmatched proteins remain in the coverage table.</p>"
            + '<p><a href="protein_annotations.tsv">Protein annotation coverage</a> · '
            '<a href="function_edges.tsv">All protein–function links</a></p>'
        )
    (out / "index.html").write_text(f"""<!doctype html><html lang="en"><meta charset="utf-8">
<title>FastaLake peptide–protein evidence</title><style>body{{font:17px/1.6 system-ui;
max-width:1080px;
margin:48px auto;
padding:0 24px;
color:#213547}}h1{{line-height:1.2}}a{{color:#246b82}}img,iframe{{width:100%;
border:0}}aside{{background:#f2f6f8;
padding:20px}}small{{color:#475569}}</style>
<h1>Peptide–protein evidence</h1><p>{scope}.</p>
<aside>Showing {len(view["nodes"])} of {view["total_nodes"]} nodes and {len(view["edges"])} of
{view["total_edges"]} edges.
Selection: {html.escape(view["selection"])}. Solid edges mark the assigned reporting representative;
 dashed grey edges show other reported memberships. Purple edges show optional annotations.
CKG finds {len(set(communities.values()))} topological communities within this view;
 these are not functional pathways.</aside>
<p><a href="network.html">Open interactive network</a> · <a href="evidence_counts.html">Explore
full evidence counts</a> · <a href="view_nodes.tsv">View nodes</a> · <a
href="view_edges.tsv">View edges</a></p>
{annotation_note}
<img src="network.png" alt="Bounded peptide and reported-protein membership network">
<iframe src="network.html" height="850" title="Interactive evidence network"></iframe>
<p><small>Network analysis and interactive plotting use Clinical Knowledge Graph (CKG),
NetworkX, Plotly and PyVis. Complete source identities and rendering settings are in
RENDER_COMPLETE.json.</small></p>{credits}</html>""")
    files = [p.name for p in out.iterdir() if p.is_file()]
    source_hashes = {
        name: {
            "path": str(Path(inspect.getfile(module)).resolve()),
            "sha256": _sha(inspect.getfile(module)),
        }
        for name, module in [("analytics", analytics), ("viz", viz)]
    }
    receipt = dict(
        status="PASS",
        graph_id=network["graph_id"],
        full_nodes=view["total_nodes"],
        full_edges=view["total_edges"],
        visible_nodes=len(view["nodes"]),
        visible_edges=len(view["edges"]),
        hidden_nodes=view["hidden_nodes"],
        hidden_edges=view["hidden_edges"],
        view_components=nx.number_connected_components(graph) if graph else 0,
        view_communities=len(set(communities.values())),
        community_scope="Unweighted bounded-view topology",
        annotation_layer_sha256=_sha(Path(annotations) / "annotations.json")
        if annotations
        else None,
        ckg_sources=source_hashes,
        versions={
            n: importlib.metadata.version(n) for n in ["networkx", "pandas", "plotly", "pyvis"]
        },
        files={name: _sha(out / name) for name in files},
        scope=network["scope"],
    )
    (out / "RENDER_COMPLETE.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(receipt, indent=2))
    return receipt


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--ckg-source", type=Path)
    parser.add_argument("--annotations", type=Path)
    parser.add_argument("--max-nodes", type=int, default=150)
    parser.add_argument("--max-edges", type=int, default=300)
    args = parser.parse_args()
    render(
        args.bundle, args.output, args.ckg_source, args.max_nodes, args.max_edges, args.annotations
    )
