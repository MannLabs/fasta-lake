# Stages in the supported sharing workflow

This page resolves the stage-reference link retained in the component parameter notes.
The sharing release uses the portable manifest runner described in the
[user guide](../../docs/how-to/acquisition-manifest.md). The historical laboratory scheduler chain
is outside this distribution. The older component notes' reference to an exhaustive
site configuration inventory therefore does not describe this page.

| Stage | Executable or script | Inputs and resulting object |
|---|---|---|
| Reference construction | `lake_builder` | Source FASTAs become exact normalized sequence identities and header mappings |
| Pooled evidence selection | `fasta_extractor` | Reference lake and pooled eligible predictions become an evidence lake |
| Acquisition-specific drawing | `fasta_extractor` | Evidence lake and one acquisition's ranked predictions become a drawn FASTA |
| Razor selection | `parsimony_engine` | Drawn FASTA and eligible inference peptides become a selected search FASTA |
| Search and LFQ | SAGE 0.14.6 | One mzML and its selected FASTA become search and quantitative observations |
| Filtered aggregation | `aggregation_engine` | Accepted target observations become peptide-feature and member-set matrices |
| Study reporting | `tools/group_study.py` | Pooled accepted evidence becomes a study dictionary and dual-key table |

The supported workflow uses no sequence clustering. Construction uses the top-30%
rank rule, 9–50-residue evidence and I/L matching; inference uses eligible mapped
predictions without that rank cutoff, with `razor` and `hash-acc` specified explicitly.
The [Methods](../../docs/concepts/methods.md) records the full validated search settings
and stage-specific counting rules. The [CAMPI walkthrough](../../examples/campi/README.md)
runs every stage on included measured inputs.

Use each executable's `--help` for its complete current CLI parameter list.
The portable runner writes actual commands, parameters, input hashes and executable
checksums into each new output directory. These execution records establish the
settings of a particular run; a default listed in an older component note does not.
