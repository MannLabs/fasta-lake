# Names used in FastaLake

A biological sample, a database entry and a reported protein group describe
different things. These definitions follow the workflow from input to result.

| Term | Meaning | Where you see it |
|---|---|---|
| **Specimen** | Biological material linking proteomic and matched DNA/RNA measurements. One specimen can have several acquisitions. | Molecular source manifest |
| **Acquisition** | One mass-spectrometry measurement, supplied as one mzML file with its predictions. The manifest's `sample` column names acquisitions. | `samples.tsv` |
| **Protein reservoir** | The prepared collection of complete protein sequences from which candidates are retrieved. Also called the **reference lake** in commands. | `--lake` |
| **De novo query** | An eligible predicted peptide retained for exact lookup. It is selection evidence, before database-search acceptance. | Prediction CSV |
| **Study-supported sequence space** | Reservoir proteins containing a retained query from the declared study. Also called the **evidence lake** in output files. | `evidence.fasta` |
| **Acquisition-specific candidate set** | Proteins retrieved from the study-supported space by one acquisition's retained queries. | `sample_lake.fasta` |
| **Selection parsimony** | Deterministic razor assignment of mapped predicted peptides to candidates; proteins receiving an assignment are retained in full. | Rust `parsimony_engine` |
| **Selected acquisition-specific database** | Complete proteins retained after selection parsimony. This is the de novo baseline. | `<sample>_razor.fasta` |
| **Final target FASTA** | Proteins supplied to the search, including any requested molecular additions. Sage generates decoys separately. | Search `config.json` |
| **Accepted peptide** | A target peptide passing the declared search-confidence threshold. Canonical sequence summaries remove modification notation and treat I/L as equivalent. | Peptide evidence tables |
| **Study-level protein group** | A reporting representative and the accepted peptides assigned to it across acquisitions. Shared evidence can leave biological origin ambiguous. | `study_groups.tsv` |
| **Protein-group quantity** | The sum or directLFQ estimate for a study group in an acquisition. An unquantified cell remains missing. | `study_group_matrix.tsv` |

## Two kinds of sequence identity

Peptide lookup treats isoleucine and leucine as equivalent. The prepared lake's
protein hashes preserve the original I/L distinction: one matching peptide can
occur in several distinct protein sequences. See
[identity rules](algorithms.md#four-identity-and-counting-rules).

## Selection and reporting

Before searching, selection parsimony uses predicted peptides to choose a
compact FASTA. After searching, study grouping uses accepted evidence to build
a common reporting dictionary. They use different evidence and tie rules.

The [workflow](pipeline-stages.md) shows both operations. The
[function map](../reference/function-map.md) points to their implementations.

## Where molecular additions enter

Appending specimen-matched DNA/RNA-supported sequences **after selection**
differs from pooling assembly proteins into the reservoir **before lookup**.
Report where each source entered; “uses metagenomics” alone does not specify
the method.
