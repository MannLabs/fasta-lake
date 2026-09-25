# Optional CKG rendering adapter

FastaLake calls [MannLabs/CKG](https://github.com/MannLabs/CKG) through a separate
Python environment. The tested installation is a locally modernised checkout
based on `6b5b625`, with Python 3.12, NetworkX 3.6, Plotly 6.5, Dash 3.3 and
PyVis 0.3.2. This is a compatibility-tested local source, not a claim that a fresh
upstream installation has the same dependencies or behavior.

The supplied `compatibility.patch` proposes four focused corrections:

- Render PyVis through its public HTML API.
- Preserve caller options and graph data, and pass width/height by name.
- Use pandas' supported `records` orientation.
- Return an edge, rather than a node, for Girvan–Newman's edge callback.

These changes were tested in an isolated copy. The original working CKG source
and configuration were preserved. Inspect and apply the patch only to a
compatible copy; `git apply --check compatibility.patch` checks whether its
context matches. It does not install CKG or reproduce all local modernisation.
The standalone compatibility regression is `test_compatibility.py`.

The adapter imports `ckg.analytics_core.analytics.analytics` and
`ckg.analytics_core.viz.viz`. It needs no Neo4j server or database credentials.
The [network tutorial](../../docs/how-to/evidence-networks.md) shows the commands and
output meanings. The default Docker image exports networks and joins eggNOG
annotations; rendering through CKG is an optional separate environment.

Upstream source and licence remain CKG's. `LICENSE.rst` retains its licence;
the patch is an explicitly identified FastaLake compatibility contribution.
No upstream issue or pull request has been submitted by this workflow.

Santos A et al. *A knowledge graph to interpret clinical proteomics data.*
Nature Biotechnology 40:692–702 (2022).
[doi:10.1038/s41587-021-01145-6](https://doi.org/10.1038/s41587-021-01145-6).
