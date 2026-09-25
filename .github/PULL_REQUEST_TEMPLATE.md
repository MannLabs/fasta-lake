## What this changes

<!-- The trigger (issue, failure, or measurement), the resulting behaviour, and what stays unchanged. -->

## How it was validated

<!-- Commands run and their outcome. A minimal synthetic reproducer beats a private path. -->

## Checklist

- [ ] `ruff check fasta_lake tests` and `python -m pytest tests -q -rs` pass; skipped data checks are listed, not hidden
- [ ] Changed Rust crates pass `cargo fmt --check`, `cargo clippy -D warnings` and `cargo test`
- [ ] `python tools/render_reference.py --write` was run and its output is committed; `python tools/source_manifest.py --write` too if this is a release branch
- [ ] A regression test covers the actual failure or scientific invariant
- [ ] Documentation and docstrings describe the new behaviour; `mkdocs build --strict` passes
- [ ] Existing outputs, expected results and unfavourable results are preserved, not edited to fit
