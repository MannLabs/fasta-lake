# Community context and evidence links

Run this after the [measured CAMPI example](../campi/README.md):

```bash
python examples/community/run_example.py --groups campi_run/study --out community
```

The complete Docker demo already runs it. Open `community/taxonomy/index.html`
for Unipept's interactive sunburst, treemap and taxonomy tree. The companion
holobiont PDF/SVG/PNG shows both peptide counts and summed LFQ-feature intensity.

The example reads **recorded public Unipept results**, checks the exact query
roster and preserves counts and intensity. It makes no network request. Native
Unipept 4.2.2 output used UniProt 2026.02 with `--all --equate`; the timestamp,
command and hashes are in [input provenance](inputs/PROVENANCE.json).

It also exports the peptide–protein graph and attaches one real eggNOG-mapper
enolase annotation by full protein sequence. The remaining proteins are
deliberately unannotated in this teaching fixture. This is partial coverage,
not a complete functional analysis. The portable eggNOG table retains the native
header and exact row bytes; only execution-comment lines were removed, with the
untouched original hash recorded. The query FASTA contains the one matched query.

[Use your own Unipept results](../../docs/how-to/taxonomy.md) ·
[Add annotations or render with CKG](../../docs/how-to/evidence-networks.md).
Both guides give interpretation limits and tool citations; citations also
appear in the taxonomy report and exported figures.
