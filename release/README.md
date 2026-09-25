# Release records

`SOURCE_MANIFEST.json` records the files and hashes captured for a release.
`RELEASE_PROVENANCE.json` preserves its build/source provenance. The existing
records belong to **1.1.0rc2**; their bytes and historical paths are preserved.
They do not claim that the current development tree matches that release.

When preparing the next release, stage the intended source tree, then run from
the repository root:

```bash
python tools/source_manifest.py --write
git add release/SOURCE_MANIFEST.json
python tools/source_manifest.py --check
```

The manifest uses paths relative to the repository root and excludes itself.
Update release provenance as part of the release process. CI enforces manifest
consistency on release branches and version tags; ordinary development branches
retain the preceding release snapshot until a new one is prepared.
