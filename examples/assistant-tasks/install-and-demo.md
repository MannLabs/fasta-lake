# Install FastaLake and inspect the demo

Help me get a checked first result with FastaLake.

My setup:
- Repository: [local checkout, or authorized GitHub clone]
- Computer: [operating system and architecture]
- New result folder: [path]
- Limits: [available disk/RAM or known cluster rules]

Read AGENTS.md, docs/get-started/install.md and docs/get-started/demo.md.
Inspect the actual environment before choosing Docker or the source route.
Use the commands from this checkout, and preserve existing environments.

Run the bundled examples in a new directory. On a Docker-capable machine,
the repository command is:
`bash tools/check_local_docker.sh docker_check_01`

The first build needs downloads; the examples run offline afterward.
Check exit status, the final PASS and each example's VALIDATION.json.
Then show me the CAMPI QC image and study matrix and explain their rows,
columns and missing values. Give me the revision, result paths and any
unavailable checks. If you cannot run commands, guide me through them and
inspect the output I return.
