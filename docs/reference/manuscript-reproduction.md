# Software runs and manuscript reproduction

A software example, a figure rebuild and a replay of the original scientific
analysis answer different questions. Choose the level you need.

| Goal | Starting point | What to retain |
|---|---|---|
| Learn or test the supported workflow | [CAMPI example](../get-started/first-study.md) | Example validation and command logs |
| Analyse a new study | [Acquisition manifest](../how-to/acquisition-manifest.md) and [search settings](../how-to/search-settings.md) | Exact inputs, settings, executable versions and output provenance |
| Reproduce a manuscript result | The manuscript version's source-data package and producing analysis code | Manuscript identity, original executable/configuration identities, input and output hashes |

The companion [paper repository](https://github.com/MannLabs/fasta-lake-paper)
is maintained separately and may require collaborator access. Its
`manuscript/CURRENT.json` identifies its packaged manuscript; verify that it
matches the document you are reading. A newer manuscript circulated separately
does not automatically update that repository. Historical plotting profiles
reproduce their declared snapshots.

## Distinguish the peptide populations

For the September manuscript workflow:

| Operation | Convention | Manifest-runner option |
|---|---|---|
| Exact lookup | Eligible predictions of 9–50 residues; top 30% by score, including cutoff ties; I/L equivalent | `--min-length 9 --max-length 50 --top-fraction 0.30` |
| Selection parsimony | All eligible acquisition predictions mapped onto candidates; fixed evidence counts; accession-hash ties | Runner uses `hash-acc` |
| Main search | Search peptides of 5–50 residues; subsequent target peptide q ≤ 0.01 | `--search-min-length 5 --search-max-length 50` |
| Seven-residue sensitivity search | A separate search with confidence re-estimated | `--search-min-length 7`, with the remaining archived settings |
| Nested 9–50 analysis | Subset of accepted peptides from the original 5–50 search, retaining that search's q values | A downstream subset, not a new minimum-9 search |

These are selected conventions, not a complete reproduction command. Use the
exact archived reservoir, prediction exports, acquisition roster and search
configuration for each comparison. Support for another predictor's input format
does not establish equivalent scientific performance.

## Preserve the producing software

Ordinary identification searches and corrected joint-LFQ analyses have distinct
Sage provenance. Use the matching source, patch and build records for the
particular result; the bundled example executable does not stand in for every
historical producing build. Experimental bag/substitution analyses also have
their own producing code and are separate from the supported exact-lookup route.

In a supplied paper checkout, `python paper.py verify` checks packaged identities
and declared producer coverage. Follow that version's instructions for numerical
recomputation or figure generation, using a new output directory. Byte integrity,
saved-result arithmetic and full raw-data replay are distinct checks.

The [validation overview](../concepts/validation.md) explains how these checks
relate to scientific interpretation.
