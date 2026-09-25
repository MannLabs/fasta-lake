# Diagnose a failed or incomplete run

Help me understand this FastaLake failure.

- Revision and environment: [values]
- Command or scheduler job ID: [value]
- Output folder: [path]
- Error or unexpected result: [description]

Read AGENTS.md and the relevant guide. Inspect the actual process/scheduler
state, commands.json, stage logs and completion records. Identify the first
failing stage and explain the evidence for its cause.

Preserve the failed output and input files. Do not change scientific cutoffs
to obtain a passing result. Fix the demonstrated problem and validate the fix;
use a new output directory if a rerun is needed. Distinguish execution failure
from an expected exclusion, such as too few acquisitions for PCA.

Report what failed, what changed, which checks passed and where the new
results and original logs are.
