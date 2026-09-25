# Prepare and run my first study

Help me run my own acquisitions through FastaLake.

My inputs:
- Repository and execution environment: [paths/platform]
- Prepared protein FASTA: [path]
- Prediction CSVs and matching mzMLs: [paths]
- Acquisition manifest: [path, or help needed]
- New output folder and resource limits: [values]
- Quantification: [sum, or directLFQ if deliberately selected]

Read AGENTS.md, docs/get-started/first-study.md,
docs/how-to/acquisition-manifest.md and docs/how-to/compute.md.
Use one manifest row per acquisition. The working example is
examples/campi/inputs/samples.tsv; paths are relative to the manifest.
Confirm ambiguous pairings with me. Do not infer biological specimen identity
from a filename.

Validate with tools/run_manifest.py --validate-only and the intended settings.
Profile one representative full acquisition before sizing the full study.
Keep every input unchanged and use a fresh output directory for each attempt.

After running, inspect COMPLETE.json, commands.json, the study matrix and QC.
Explain missingness and PCA eligibility. Record the actual configuration,
revision, input identities, measurements and remaining questions.
