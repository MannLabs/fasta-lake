# Troubleshoot a run

This page lists the messages a run can show, their likely causes and the next
step. You need the run's logs and QC summary at hand.

## Messages you may see first

| Message or observation | Next step |
|---|---|
| Output already exists | Choose a new `--out`; keep the earlier logs for comparison |
| Missing executable | Use Docker, or supply the binary directory/explicit Sage path |
| Specimen mismatch | Correct the acquisition-to-source mapping before searching |
| Missing sequences in molecular inputs | Inspect the coverage ledger; choose incomplete-reference mode only for a declared available-reference analysis |
| PCA skipped | Read `qc/summary.json`; inspect coverage and the complete varying group count |
| Process killed for memory | Keep the failed logs, inspect the failing stage, and use the resource guide to size the next pilot |

A selected protein is a search candidate. Read its peptide evidence and
quantification before making a biological claim.

For existing DIANovo CSVs, read the [sequence and probability input contract](de-novo-inputs.md).

## Symptoms, causes and actions

| Symptom | Likely explanation | Action |
|---|---|---|
| Output exists | A previous or partial run occupies that path | Use a new output directory and retain the old logs |
| Executable not found | Python installed without Rust builds or PATH configuration | Build required crates and inspect binary discovery |
| No matched proteins | Input format, score direction, length filter or reference coverage | Inspect eligible predictions and test a known peptide containment example |
| Missing results table | Search failed or an incomplete directory was supplied | Read SAGE logs and job termination status |
| Several LFQ intensity columns | Several acquisitions were searched in one invocation | Use the supported one-acquisition search layout or implement an explicit wide-table adapter |
| Many UNCLASSIFIED sources | Sequence hashes carry no catalogue semantics | Supply the source mapping sidecar |
| Very large search memory | Candidate or metagenome database exceeds the pilot allocation | Profile and rerun in a new directory with adequate memory |
| More member-set rows than expected proteins | Different accession combinations form different keys | Report member sets honestly or construct a declared study dictionary |
| Low overlap between technical acquisitions | Stochastic evidence and representative selection contribute | Compare canonical peptides, exact sequences and study groups separately |

## Limits of the validated workflow

!!! note
    Clustering is absent from the validated workflow. Experimental Python inference, taxonomy budgets, mass-k-mer rescue, raw conversion and alternative engine exporters remain separate interfaces with different assumptions. They should not be substituted into the tutorial while retaining the same Methods description. A new species, acquisition mode, score convention or reference design requires its own pilot and validation.
