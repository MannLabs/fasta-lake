# Interpret core biological functions

A functional profile helps ask what a microbial community can do, and which of
those capabilities have protein evidence. A useful starting panel separates
**cellular housekeeping**, **functions common in the sampled ecosystem**, and
**specialised functions**. These are interpretation categories, not a universal
list of genes required in every organism.

## A starting panel for gut metaproteomics

| Panel | Examples | Interpretation |
|---|---|---|
| Information processing | Translation/ribosomes, transcription, DNA replication and repair | Cellular maintenance and growth machinery |
| Central metabolism and energy | Central carbon reactions, ATP-generating machinery | Conserved biochemical processes; routes differ between organisms |
| Building blocks and cofactors | Amino-acid, nucleotide, vitamin and cofactor metabolism | Biosynthesis and salvage; auxotrophy and cross-feeding matter |
| Protein and cell maintenance | Protein folding/turnover, membrane and envelope processes | Cellular upkeep; not every component is universally essential |
| Gut ecosystem functions | Glycan utilisation, fermentation and short-chain-fatty-acid pathways | Functions distributed across community members; a different organism can provide the same reaction |
| Specialised functions | Selected xenobiotic transformations or other niche-specific pathways | Show a stated reference and biological rationale; do not label every non-core annotation as novel |

This is a proposed reporting panel, not a completed automated KO-to-core mapping.
The initial labels follow the distinction between the minimal gut bacterial genome
and the community metagenome in [Qin et al. (2010)](https://www.nature.com/articles/nature08821).
The expanded Human Microbiome Project distinguishes broadly distributed pathways
from human-associated and body-site-enriched functions; its reported core depends
on a specified population and detection rule.
[Lloyd-Price et al. (2017)](https://www.nature.com/articles/nature23889).

## Keep three questions separate

1. **Reference classification:** is this a recognised housekeeping or ecosystem function?
2. **Observed prevalence:** in how many eligible specimens was there evidence here?
3. **Molecular support:** was the candidate selected by de novo, metaG, metaT or another declared source, and was it subsequently identified by MS?

!!! warning
    Do not transfer a DNA prevalence threshold directly to proteomic detection.
    Undetected protein evidence does not establish biological absence. A KO hit does
    not establish a complete pathway or metabolic flux. Pathway/module completion
    requires the relevant alternative and required reactions to be evaluated.
    See the [KEGG module definitions](https://www.kegg.jp/kegg/module.html).

For paired multiomics, compare the same specimen and reference definitions, retain
unannotated signal, and keep DNA potential, RNA evidence and measured protein
quantities distinct. A newly selected FASTA member is not automatically a newly
identified function. Demonstrating gains needs the paired baseline search too.
