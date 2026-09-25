# Security policy

## Reporting a vulnerability

Use GitHub's private vulnerability reporting for this repository
(**Security → Report a vulnerability**). Do not open a public issue for a security problem.
Include the version (`fasta-lake --version`), the command or input that triggers the problem,
and what an attacker could achieve. You will get a first response within ten working days.

## Scope

FastaLake reads FASTA and prediction tables that you supply, runs local binaries, and downloads
public reference resources over HTTPS from the catalogued hosts in `fasta_lake/resources.py`.
It has no network service, no authentication and stores no credentials. Reports about the
following are in scope:

- a crafted input file (FASTA, CSV, TSV, mzML, JSON configuration) that makes a command read or
  write outside its output directory, execute code, or exhaust memory without a clear error
- a download path that accepts a non-HTTPS URL, skips the size or checksum checks in
  `fasta_lake/downloads.py`, or writes outside the target directory
- the container image running with more privilege than `USER 10001`, or a CI workflow that
  exposes a token to untrusted input

## What is checked before every release

- every GitHub Action is pinned to a commit hash and runs with `contents: read`; checkouts do
  not persist credentials
- Rust dependencies pass `cargo deny` with `rust/deny.toml` (advisories, licences, sources)
- Python dependencies are audited with `pip-audit` on the release pins; the container base
  images are pinned by digest
- `tools/check_no_site_paths.py` fails CI if any tracked text file carries an internal path

## Supported versions

Only the latest release on `main` receives fixes.
