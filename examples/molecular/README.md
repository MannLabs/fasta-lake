# Synthetic molecular-evidence example

After installing this source version of FastaLake, run from the repository root:

```bash
python examples/molecular/run_example.py --out molecular_example_01
```

No spectra, network, search engine or Rust binary is needed for this example.
It creates four synthetic gene records, a TPM table, a checksummed source manifest,
and a one-protein de novo baseline. Gene B passes the DNA cutoff; gene C passes
the RNA cutoff. Genes A and B encode the same sequence, so their original gene
records remain separate while only one protein is added. The expected final FASTA
contains three exact sequences, with two additions and the baseline header intact.

These invented values test the file contract and selection behavior only. They
provide no evidence of biological performance or an optimal abundance threshold.
See [the comparison and sweep documentation](../../docs/concepts/molecular-complementarity.md)
for measured-search requirements and the complete input contract.
