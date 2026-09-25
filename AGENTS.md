# FastaLake instructions for coding and analysis assistants

Help the user reach a checked result with understandable commands and clear
evidence. Follow the user's task and existing authorization; ask for missing
scientific inputs or access when needed, without repeatedly reconfirming agreed work.

## Orient before running

- Read [README.md](README.md), then [the assistant guide](docs/how-to/assistants.md).
  For code changes, also read [.github/CONTRIBUTING.md](.github/CONTRIBUTING.md).
- The repository root is the supported source; run Python/Rust development commands
  and `bash tools/check_local_docker.sh docker_check_01` from there.
- Record the repository revision and working-tree state. Inspect the actual OS,
  architecture, available RAM/CPU/disk, Docker or existing source environment,
  and scheduler access. State which machine your terminal can actually reach.
- FastaLake was developed in a cluster environment, and we optimise primarily
  for that setting. Use allocated compute for substantial work; follow the
  site's scheduler rules. Inspect existing jobs and completion receipts before
  submitting another run. Do not guess a partition, account or GPU allocation.

## Installation and first result

- Use [the Docker guide](docs/get-started/docker.md) for the bundled environment or
  [source installation](docs/get-started/install.md#source-installation) where
  Docker is unavailable. Read the checked-out commands; do not invent a PyPI
  installation, registry image tag or executable name.
- Use an isolated Python environment and the recorded dependency/tool versions.
  A Python installation alone does not supply the separate Rust and Sage binaries.
- On a Docker-capable machine, the repository-root check above builds the image
  and runs the complete offline demo. It requires a new output directory. Inspect
  exit status and the example validation records before reporting success.
- In a chat without terminal access, supply commands for the user to run and
  inspect their returned output. Do not claim to have executed them.

## A study with the user's data

1. Confirm prepared protein FASTA, supported de novo predictions, mzML files and
   an acquisition manifest. RAW conversion and de novo model inference are upstream.
   Resolve paths relative to the manifest; do not infer specimen pairing from names.
2. Keep source data unchanged and choose a new result directory. Prefer read-only
   input mounts in Docker. Read the relevant tutorial and the current command's
   `--help`; `python tools/run_manifest.py` is the source-run entry point inside
   the repository root, and `/opt/fastalake/tools/run_manifest.py` is its image path.
3. Run `--validate-only` with the intended settings before execution. Preflight is
   not a completed analysis or a guarantee of sufficient RAM. Use a representative
   full acquisition to size the main study; consult the measured resource guide.
4. Preserve de novo razor targets when adding specimen-matched metaG/metaT or
   metabolite-supported candidates. Follow [molecular inputs](docs/how-to/molecular-inputs.md)
   for shared gene references, independent layer cutoffs and missing-sequence policy.
   Do not silently enable incomplete-reference mode or substitute missing TPM with zero.
5. Keep default sum quantification unless another method is selected. directLFQ
   and basic QC/PCA are supported; PCA eligibility and exclusions must be reported.
   Annotation, ESM-C and experimental downstream analysis are explicit extensions.
6. Download only the chosen resources into the chosen local destination. Report
   transfer size and disk requirements, including unknown estimates. Retain the
   source/version/checksum records. Do not start every catalogue/model download.
7. Execute within the agreed resources. Threads share process data; independent
   jobs can duplicate it. `--laptop` plans capacity but does not enforce a RAM cap.
   Keep failed logs; inspect the failing stage before changing the resource request.

## Verify and hand back

- Check process exit codes, scheduler state and completion/validation receipts.
  Read `commands.json`, the exact search configuration and `qc/summary.json`.
  A queued job, a generated plot or an existing directory is not a workflow pass.
- Report revision/image identity, commands, inputs/settings, output paths,
  verification performed, skipped checks and remaining limits. Link the matrix,
  QC and peptide/molecular evidence. Separate example checks from full-study claims.
- Preserve source data, existing outputs and unfavourable results. The manuscript
  and its review records live outside this repository and are not part of a run.
- Selection support is not peptide identification. Shared peptides can leave
  protein/species/function attribution ambiguous. Search q values do not establish
  whole-workflow FDR; independent entrapment remains unresolved. Do not change the
  declared method to match expected numbers or treat embeddings as function labels.
- Keep credentials and unapproved research data out of chat and external posts.
  The code licence does not grant redistribution rights to every bundled input.

## Changes to the software

Follow [functions in execution order](docs/reference/function-map.md). Keep scientific
contracts explicit: exact sequence identity, acquisition membership, cutoff ties,
missingness, deterministic ordering/sums and complete provenance. Add meaningful
regressions for changed behavior. Run the affected Python/Rust checks and the
installed/Docker route when workflow or packaging changes. Use CI's build steps
for cross-language tests, and report unavailable-data skips explicitly.
