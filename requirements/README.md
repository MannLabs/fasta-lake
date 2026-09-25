# Python dependency pins

- `container.txt`: hash-locked Python 3.12 / Linux x86_64 runtime, shared by
  Docker and the source examples. Install from the repository root with
  `python -m pip install --require-hashes -r requirements/container.txt`.
- `release.txt`: core dependency constraints used to compile that lock and by
  the release dependency audit.

The generation command is recorded at the top of `container.txt`. Update the
runtime as a set; do not change an individual version without regenerating its
dependency closure and hashes. See [contributing](../docs/project/contributing.md).
