# Measure what DNA and RNA add

This page measures what a specimen's DNA and RNA add to a de novo search, at
the level of exact protein sequences and peptides. You need a completed de novo
search and a validated molecular source. Keep the final de novo FASTA, append exact sequences supported by a specimen's
DNA/RNA measurements, then search the same spectra again. Compare both gained
and lost evidence: keeping every baseline protein does not guarantee retaining
every accepted peptide.

## Build explicit alternatives

First validate the specimen-qualified, explicitly scoped quantification reference using
the [molecular input contract](molecular-inputs.md). DNA and RNA share the same
gene-to-sequence table. Missing values are unavailable, not zero. The current
prepare-source interface supports complete references and explicit available-gene
subsets; available scope exports missing genes instead of silently dropping them.

These commands use all positive TPM genes as explicit comparison arms.
They are not universal optimized cutoffs; the existing add-command default
remains a separate configured setting.

~~~bash
fasta-lake molecular validate-source --source source.json
fasta-lake molecular select --source source.json --specimen specimen_1 \
  --dna-percent 100 --background human.fasta --out dna_alone
fasta-lake molecular select --source source.json --specimen specimen_1 \
  --dna-percent 0 --rna-absolute 0 --background human.fasta --out rna_alone
fasta-lake molecular add --baseline denovo.fasta --source source.json \
  --specimen specimen_1 --dna-percent 100 --out dna_union
fasta-lake molecular add --baseline denovo.fasta --source source.json \
  --specimen specimen_1 --dna-percent 0 --rna-absolute 0 --out rna_union
fasta-lake molecular add --baseline denovo.fasta --source source.json \
  --specimen specimen_1 --dna-percent 100 --rna-absolute 0 --out dna_rna_union
~~~

Use the resulting search.fasta in a copy of the acquisition's Sage configuration.
Preserve spectra, search settings, decoy generation, confidence thresholds and
the same host background. Keep every comparison, including losses. The add
command exports the per-gene TPM and exact-protein provenance needed below.

## Measure sequence-level support

~~~bash
fasta-lake molecular increments \
  --left searches/denovo --right searches/denovo_dna \
  --reference full_declared_reference.fasta --host-reference human.fasta \
  --left-label de_novo --right-label de_novo_DNA --out increments/dna
~~~

The declared reference must include every exact protein sequence in both search
FASTAs and all available alternatives relevant to specificity, including the
human background. Do not pass just the selected union if a wider reference is
available. The software checks sequence containment; the provider must justify
the declared reference's scope and biological provenance.

Outputs:

- `protein_support.tsv.gz`: every reference sequence compatible with accepted
  peptides in either arm; exact hash, reference accessions, candidate membership,
  reported PSM support, full supporting peptide lists, gained/lost peptides and
  reference-specific support. Join protein_sha256 to gene_evidence.tsv from
  molecular add or another explicitly qualified annotation dictionary.
- `peptide_changes.tsv.gz`: every retained, gained and lost canonical peptide,
  all compatible reference hashes, and host/nonhost ambiguity if a host reference
  was supplied. “Nonhost” does not independently establish microbial taxonomy.
- `summary.json`: peptide identity retention in both directions, positive LFQ
  feature counts, gains/losses and shared-feature agreement; compatible sequence
  counts at one peptide, two peptides, or two including one reference-specific
  peptide. Added candidates are separate from additional peptide evidence.
- `COMPLETE.json`: exact input, implementation and output hashes. Existing
  output folders are never overwritten.

Protein identity uses exact normalized sequence hashes. Peptide compatibility
uses I/L equivalence and ignores bracketed modifications; LFQ alignment retains
the exact peptidoform and charge. A peptide specific to one sequence in the
declared reference is not universally unique. Compatible-sequence counts and
the two-peptide threshold are not protein-level FDR, a protein-grouping method,
unique organism attribution, or measured enzyme activity. Original same-run LFQ
agreement is distinct from joint LFQ, MBR and independent quantitative accuracy.

## A simple benchmark figure

Use paired acquisitions, a fixed donor roster and donor-level summaries.
Show (1) yield for de novo, molecular alone and their union, (2) the fraction of
molecular-search peptide identities also recovered by de novo, and (3) added
and lost protein support. Include candidate database sizes and shared-peptide
ambiguity. Equal counts are not identical proteomes.

Two reference scopes answer different questions:

- **Public-reservoir de novo versus metagenome-derived searches:** can recovery
  approach the metagenome comparator without that specimen's sequencing?
- **Same specimen reference, different selectors:** how efficiently does de novo
  choose among already available genes, and what do measured DNA/RNA add?

!!! warning
    Do not use the second experiment to prove the first. An available-gene research
    subset with documented missing genes must not be labelled a verified original
    quantification reference to bypass molecular validate-source.
