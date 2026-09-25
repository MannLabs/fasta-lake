# Configuration file

`fasta-lake init --template <name>` writes an example `fastalake.json`; `fasta-lake validate` checks one. Three fields are required: `database.sources`, `peptide_evidence.predictions_dir` and `output.directory`. Everything else has a default. The templates below are generated from the code; the field semantics are documented in the [`fasta_lake.config`](api/config.md) module.

=== "minimal"

    ```json
    --8<-- "docs/reference/generated/fastalake_minimal.json"
    ```

=== "metaproteomics"

    ```json
    --8<-- "docs/reference/generated/fastalake_metaproteomics.json"
    ```

=== "clinical"

    ```json
    --8<-- "docs/reference/generated/fastalake_clinical.json"
    ```

=== "full"

    ```json
    --8<-- "docs/reference/generated/fastalake_full.json"
    ```
