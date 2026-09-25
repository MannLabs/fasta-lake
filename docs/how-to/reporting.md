# Report and share a study

This page lists what to report with a completed study so that others can repeat
it, and what a public release needs. You need the completed run directory and
its provenance records.

## Report the provenance

Supply the exact software version or source manifest, binary versions and hashes, acquisition manifest, catalogue provenance, prediction provenance, complete search configuration, exact FASTAs searched and grouping dictionary. Explain any host or contaminant additions and every exclusion rule. Report the counting unit beside each depth or overlap statistic.

## Describe the design

For repeatability, state which acquisitions share biological material and which variation the design can estimate. For comparisons, retain both arms for each paired acquisition, including low-yield outliers. Separate the effect of reference content from the effect of sample-specific selection when the scientific claim concerns one of those factors.

## Release publicly

For a public release, deposit shareable data and software in an appropriate versioned archive and use the assigned persistent identifiers in the final manuscript. The local release bundle and this guide do not establish that a GitHub release, data deposition or DOI has already been created. An independent colleague should run the example in a clean environment before the package is described as broadly portable.

## Documentation sources

The documentation structure takes inspiration from [SAGE's algorithm and result documentation](https://sage-docs.vercel.app/docs/how_it_works) and [AlphaPeptDeep's tutorials and API documentation](https://alphapeptdeep.readthedocs.io/en/latest/). The FastaLake explanations and examples here describe this implementation in their own wording.
