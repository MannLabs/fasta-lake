# Add molecular evidence to the de novo workflow

Run after the bundled CAMPI example:

```bash
python examples/multiomics/run_example.py --campi campi_run_01 --out multiomics_run_01
```

This executes real Sage searches on the bundled CAMPI spectra. Molecular TPM,
KO and compound values are synthetic test inputs over real candidate sequences;
this example is an execution check, not a biological benchmark. Every original
de novo sequence/header must remain, three molecular candidates must be added,
and Sage must search the union before study grouping/directLFQ.
