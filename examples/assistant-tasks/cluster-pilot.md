# Size a study on my cluster

Help me prepare a small, informative FastaLake pilot.

My setup:
- Repository/environment: [path]
- Input manifest and prepared lake: [paths]
- Scheduler, account and permitted partition: [site details]
- Per-job RAM/CPU/time and concurrency limits: [values]
- New output directory: [path]

Read AGENTS.md and docs/how-to/compute.md. Inspect current jobs and completion
records before submitting anything. Use the site's known scheduler settings;
ask for settings I have not supplied.

Validate inputs, then run one representative full acquisition within the
stated allocation. Record memory, time and disk by stage, including selection
parsimony and Sage. Reference chunking and a planning budget do not enforce
a total process-memory cap. Threads and independent jobs have different memory costs.

Use the pilot to propose and execute the remaining study within the agreed
resources. Keep failed logs and use a new directory for any retry. Finish with
scheduler status, output checks, resource measurements and exact run identities.
