# Work with an agent or LLM

Choose the task closest to what you want to do. Fill in the short input block,
then give the file to an assistant that can read this repository. The tasks work
as starting instructions for a coding agent or as prompts for a chat assistant;
without terminal access, you run the commands and return the output.

| I want to… | Task file | Files to start from |
|---|---|---|
| Install and see a result | [Install and demo](install-and-demo.md) | [Demo guide](../../docs/get-started/demo.md) |
| Run my own spectra | [First study](first-study.md) | [Working acquisition manifest](../campi/inputs/samples.tsv) |
| Add matched DNA/RNA | [Molecular inputs](molecular-inputs.md) | [Multi-omics example](../multiomics/README.md) |
| Plan a cluster run | [Cluster pilot](cluster-pilot.md) | [Site configuration template](../site/README.md) |
| Fix a failed run | [Troubleshooting](troubleshoot.md) | Your command and stage logs |
| Understand the results | [Read results](read-results.md) | Your study matrix and QC |

[AGENTS.md](../../AGENTS.md) gives assistants the repository-wide starting
instructions. [The assistant guide](../../docs/how-to/assistants.md) explains
how to use these tasks.

These files are prompts, not executable configuration. Replace bracketed
placeholders with real inputs. Keep credentials out of prompts. The molecular
example uses synthetic DNA/RNA values; it is a format example, not your study.

For a runnable sample manifest, use the bundled CAMPI TSV above rather than
copying visually spaced text into a table. For a configurable run, see
[the configuration guide](../../docs/reference/configuration.md) and its
generated example files.
