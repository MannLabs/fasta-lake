# Chunked workflow acceptance

Run after completing the bundled public CAMPI example:

```bash
python examples/chunked/run_example.py --campi campi_run_01 \
  --out chunked_check_01 --laptop
```

The test repeats the actual runner with reference pieces, then compares the
merged evidence and per-acquisition drawn/razor FASTAs byte for byte against
the unchunked runs. It also checks the study matrix axes/values and automatic QC.
`--laptop` also tests resource detection, automatic threads and the saved plan.
Docker acceptance runs this test and a chunked additive metaG/metaT workflow
under a 4 GiB container limit and two CPU slots.
These are cropped measured inputs; see the [resource guide](../../docs/how-to/compute.md)
for the separate full-reference/laptop-capacity question.
