# Environment and site configuration

The study runner takes the reference FASTA, the acquisition manifest, the output directory and the Sage executable as explicit arguments. Paths inside the manifest resolve relative to the manifest. No laboratory filesystem layout is assumed, and no tracked file may contain one: a scan in CI refuses site-specific paths, and `examples/site/fastalake.site.example.sh` shows where such values belong on a shared installation.

## Environment variables

| Variable | Purpose |
|---|---|
| `FASTALAKE_BIN_DIR` | Directory containing the Rust executables. Explicit executable arguments take precedence; otherwise the package looks here, then in the repository build directories, then on `PATH`. |
| `FASTALAKE_CONVERSION_SCRIPTS` | Location of independently configured raw-conversion and AlphaNovo scripts for the experimental `convert` command. Those tools and their model weights are not bundled and are outside the real-data acceptance workflow. |
| `FASTALAKE_TEST_DATA`, `FL_JANKO_DIR`, `FASTALAKE_TEST_*` | Optional historical test fixtures. Each variable is named at the top of the test that uses it; a missing fixture produces a declared skip, and the bundled synthetic tests and tutorials never need them. |

The Python wheel contains configuration templates, not Rust executables or a Sage installation; see [Install](install.md) for how the executables are built or obtained.

## Site-specific values

`tools/check_no_site_paths.py` scans the distributed text files for hardcoded laboratory paths and runs in CI. When a script needs a site value, it reads it from the untracked `fastalake.site.sh`, never from a literal in a tracked file. Source manifests (`release/SOURCE_MANIFEST.json`) identify the exact files a release contains. Internal study provenance stays in the private review record and is not rewritten to simulate portable input locations.
