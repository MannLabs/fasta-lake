# Contributing to FastaLake

Thank you for helping make metaproteomics easier to inspect and reproduce.
Useful contributions include a clearer error message, a smaller reproducer,
a tutorial correction, a measured speedup or a new analysis.

## Understand the code in a few minutes

1. Read [the principle and diagram](../concepts/principle.md).
2. Follow [the functions in execution order](../reference/function-map.md).
3. Open the relevant module beside its tests. The whole repository is the supported source; the manuscript and its review
   records are kept privately, outside this repository.

The [development review](development.md) distinguishes engineering
work from scientific questions. An open question is welcome when its scope and
evidence are clear.

## Set up and check a change

From the repository root, use Python 3.12 on Linux x86_64 for the locked analysis
environment (the hash-pinned set the container and CI use), then add the
development tools on top:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --require-hashes -r requirements/container.txt
python -m pip install --no-deps -e '.[analysis,directlfq]'
python -m pip install -e '.[dev]'
ruff check fasta_lake tests
python -m pytest tests -q -rs
python tools/source_manifest.py --check   # on a release branch only: --write, and commit release/SOURCE_MANIFEST.json
```

Rust uses the toolchain in `rust-toolchain.toml`. For a changed crate:

```bash
cargo fmt --manifest-path rust/parsimony_engine/Cargo.toml --all --check
cargo clippy --locked --manifest-path rust/parsimony_engine/Cargo.toml --all-targets -- -D warnings
cargo test --locked --manifest-path rust/parsimony_engine/Cargo.toml
```

Replace the crate name as needed. Cross-language tests also need release
executables; [CI](https://github.com/MannLabs/fasta-lake/blob/main/.github/workflows/ci.yml) shows the exact build commands and
`FASTALAKE_BIN_DIR` setup. Tests requiring unavailable external data report a
skip. Inspect `-rs`: skipped data checks do not establish empirical acceptance.

For workflow or packaging changes, run the [complete Docker check](../get-started/install.md#verify-the-installation).
It executes the included examples offline and checks outputs. A small fixture
tests integration; representative full acquisitions are needed for resource or
scientific performance claims.

## Make the change reviewable

- Name a function for the operation it performs. Separate numeric transformations
  from file handling and orchestration where that helps reading.
- Document non-obvious inputs, units, missing values, identity rules and outputs.
  A short example often explains more than a long docstring.
- Add a regression for the actual failure or scientific invariant. Check exact
  sequences, sample membership, missingness and provenance as well as counts.
- Keep existing outputs and adverse results. Use a new destination for a replay;
  record source/input identities and explain intended numerical changes.
- In a pull request, describe the trigger, resulting behavior and validation.
  A minimal synthetic reproducer is preferable to an inaccessible private path.

## Report an issue or suggest an upstream improvement

Use [FastaLake issues](https://github.com/MannLabs/fasta-lake/issues) for this
package. Include the command, revision/image identity, platform, relevant log,
expected result and smallest shareable input. For resource issues, add the
allocation, stage and measured peak memory. Remove credentials and private data.

Before approaching another project, check its current source, releases and
existing issues. Reproduce the problem without FastaLake where possible. Credit
an existing fix; propose a focused change with a regression when one is needed.

FastaLake code uses the [MIT licence](https://github.com/MannLabs/fasta-lake/blob/main/LICENSE). Preserve attribution
and the separate terms of external dependencies, references and example data.

## Credit code where it is used

Place the upstream GitHub link and relevant paper beside each scientific adapter
or reused implementation, in its module/function docstring or a nearby comment.
Name the upstream method and distinguish our input checks or adaptations. Preserve
licence notices and record versions/changes for bundled or modified source. Carry
relevant citations into generated reports and exported figures; a central references
list alone is not sufficient.

## AI-assisted contributions

AI coding tools are part of FastaLake's development history. Describe their
role when it matters to a contribution, review the resulting code and retain
its tests and provenance. The contributor remains responsible for correctness,
scientific assumptions and upstream attribution.
