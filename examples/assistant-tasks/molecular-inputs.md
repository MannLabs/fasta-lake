# Add matched DNA or RNA evidence

Help me add specimen-matched molecular evidence to a FastaLake study.

My inputs:
- Existing proteomic study inputs/manifest: [paths]
- Specimen-to-acquisition mapping from the provider: [path]
- Shared gene protein FASTA and DNA/RNA TPM table per specimen: [paths]
- Provider's reference/provenance notes: [path]
- DNA selection rule: [chosen percentage or absolute threshold]
- RNA selection rule: [chosen percentage or absolute threshold]
- New output folder: [path]

Read AGENTS.md, docs/how-to/molecular-inputs.md and
docs/how-to/molecular-docker.md. Inspect examples/multiomics/README.md;
its molecular values are synthetic teaching data.

Prepare and validate each source before adding it to the acquisition manifest.
Require matching gene IDs and verified specimen linkage. Keep missing
measurements distinct from zero and report sequence-coverage gaps. Ask me to
resolve absent source information or unspecified scientific thresholds.

Keep every de novo baseline sequence and header. Add molecular candidates
after selection parsimony and before Sage. Inspect gene_evidence.tsv and the
final FASTA, then distinguish added candidates from accepted search peptides.
For a comparison, retain both gained and lost identifications.
