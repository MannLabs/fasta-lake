# Choose and download your local resources

Start with the proteins relevant to your study. You do not need every catalogue.
FastaLake downloads a resource only when you run `resources download`; construction,
searching and annotation use the local files you configure.

## First download: a small reference you can inspect

From an installed FastaLake environment:

```bash
fasta-lake resources list
fasta-lake resources plan uniprot-swissprot --directory ./resources
fasta-lake resources download uniprot-swissprot --directory ./resources
fasta-lake resources unpack \
  ./resources/uniprot-swissprot/uniprot_sprot.fasta.gz \
  --directory ./resources/swissprot-proteins --max-gib 1
```

The plan prints the version, official URL, destination, transfer bytes, available
disk space and expanded size when known. Downloading keeps the original archive.
Unpacking uses a new directory and an explicit **disk budget**, not a RAM allocation.
The expanded file is `resources/swissprot-proteins/uniprot_sprot.fasta`.

### If you use Docker

With the image from the [installation page](../get-started/install.md), mount the chosen
local resource folder and run the same helper inside the container:

```bash
mkdir -p resources
docker run --rm --platform linux/amd64 \
  --user "$(id -u):$(id -g)" \
  --mount "type=bind,source=$PWD/resources,target=/resources" \
  fastalake:docker fasta-lake resources plan uniprot-swissprot \
  --directory /resources
```

Replace `plan` with `download` when ready to transfer that resource. These two
commands need network access. To unpack, replace the command after
`fastalake:docker` with:

```bash
fasta-lake resources unpack /resources/uniprot-swissprot/uniprot_sprot.fasta.gz \
  --directory /resources/swissprot-proteins --max-gib 1
```

The files remain in your computer's `resources/` folder after the container exits.

A real check of release 2026_03 transferred 93,801,562 bytes and expanded
288,201,213 bytes, containing 575,748 sequences. UniProt's MD5 matched;
FastaLake also recorded SHA256. The catalogue uses UniProt’s official EBI mirror
(the original FTP host timed out during later endpoint checks). These numbers describe that release, not future releases.

## Which resource belongs in your lake?

| Choice | Purpose | Download snapshot, 15 September 2026 |
|---|---|---|
| `gmsc-90`, `gmsc-100` | [GMSC](https://gmsc.big-data-biology.org/): small-ORF proteins, at 90% family identity or exact sequence identity | 9.33 / 20.60 GB compressed |
| `gmgc-90nr`, `gmgc-95nr`, `gmgc-100nr` | [GMGC](https://gmgc.embl.de/download.cgi): global microbial protein sequences; choose a clustering level | Upstream rounded listing: 24G / 31G / 82G; exact streaming length may be unavailable |
| `gmgc-95nr.complete`, `gmgc-95nr.no-rare` | GMGC subsets of complete genes or genes excluding rare entries | Upstream rounded listing: 13G / 16G |
| `uhgp-50`, `uhgp-90`, `uhgp-95`, `uhgp-100` | [UHGP proteins associated with UHGG human-gut genomes](https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v2.0.2/protein_catalogue/) | 5.16 / 10.57 / 8.88 / 28.34 GB archives, v2.0.2 |
| `uniprot-swissprot`, `uniprot-trembl` | [UniProt](https://www.uniprot.org/help/downloads): reviewed or unreviewed protein records | 0.094 / 40.63 GB compressed; current release resolved at download |
| `uniprot-human` | [Human reference proteome UP000005640](https://www.uniprot.org/proteomes/UP000005640) for host proteins | 0.048 GB compressed snapshot; current release resolved at download |
| `eggnog-5.0.2` | Local orthology annotation with eggNOG-mapper 2.1.12 | 12.06 GB transfer; about 51 GB expanded, plus archives |
| `esmc-300m`, `esmc-600m` | Local ESM-C embeddings for explicitly selected annotation gaps | 1.33 / 2.30 GB weights; no unpacking needed |

GB in this table means 10⁹ bytes. The command reports GiB (2³⁰ bytes) and exact
byte counts where available. An archive can include cluster membership tables as
well as its representative FASTA; that is why archive sizes need not decrease
monotonically with clustering threshold.

GMSC 90AA/100AA are identity thresholds, not protein lengths. GMGC's 95nr
catalogue was clustered at 95% nucleotide identity; its 90nr catalogue is at
90% amino-acid identity. UHGG supplies genomes; UHGP is the protein route
for FastaLake. Downloading all underlying genomes is unnecessary for protein search.

These choices do not replace the paper's frozen reference or deduplicate a combined
lake automatically. Keep source identities, choose compatible headers and follow
[database preparation](#prepare-a-reference-lake) and [chunked construction](compute.md).
Annotation/quality tables are separate resources; a larger reference is not by
itself stronger identification evidence.

## Large catalogue workflow

```bash
fasta-lake resources plan uhgp-95 --directory /your/external-disk/references
fasta-lake resources download uhgp-95 --directory /your/external-disk/references
```

Choose an unpacking budget after checking free disk space. When the upstream
expanded size is unavailable, the plan says `unknown`; it does not estimate it
from the final four bytes of a gzip file, which overflow for large FASTAs.

```bash
# Example ceiling: 80 GiB of additional disk space, if that much is available.
fasta-lake resources unpack \
  /your/external-disk/references/uhgp-95/uhgp-95.tar.gz \
  --directory /your/external-disk/references/uhgp95-unpacked --max-gib 80
```

Inspect the unpacked file roster in `UNPACK_COMPLETE.json` and use the protein
FASTA as your reference input. Storage for the archive, expanded reference,
evidence and search outputs is separate from the [RAM planning budget](compute.md).

## Interruptions and integrity

An incomplete transfer stays in `.partial`; it never becomes a finished resource.
Repeat the same command with `--resume`. Resume requires an unchanged strong ETag
and a server that honours byte ranges. If the server cannot support this, use a
new destination for a fresh transfer; the partial file is retained for inspection.
GMGC's streaming endpoint may not support resume or expose an exact transfer size.

A killed process may leave an empty `.lock` directory. Check that the old download
has stopped before removing that lock and resuming. Existing resources are never
overwritten. `--resume` verifies already completed files against their receipts.

Each file has a `.download.json` receipt, and the resource has a saved plan and
`DOWNLOAD_COMPLETE.json`. SHA256 is always recorded. Where an upstream checksum
exists (including UniProt Swiss-Prot/TrEMBL and ESM-C weights), it must match. A
local checksum alone records the downloaded bytes; it is labelled accordingly.
Archive expansion rejects paths outside the destination, links, duplicate files
and data exceeding the chosen disk budget.

## Annotation and model setup

The supported eggNOG 5.0.2 files are available from the official legacy
[eggNOG endpoint](http://eggnog5.embl.de/download/emapperdb-5.0.2/). Its newer HTTPS
replacement returned 403 for these files at verification; the catalogue records
the working HTTP source explicitly. No upstream digest was available there.
The mapper checks database compatibility, and annotation records full local hashes.

```bash
fasta-lake resources download eggnog-5.0.2 --directory ./resources
fasta-lake resources unpack ./resources/eggnog-5.0.2/*.gz \
  --directory ./resources/eggnog-data --max-gib 55
fasta-lake resources download esmc-300m --directory ./resources
```

Point annotation to `eggnog-data` and the local checkpoint
`resources/esmc-300m/esmc_300m_2024_12_v0.pth`. The model is pinned to a published
revision and SHA256 from the [ESM-C model repository](https://huggingface.co/Biohub/esmc-300m-2024-12).
Tools and model dependencies are installed separately; see [annotation](annotation.md).

## Prepare a reference lake

For a first, small reference, convert the downloaded Swiss-Prot FASTA to a lake
with sequence-derived identifiers and a source-header sidecar:

```bash
mkdir curated_reference_01
lake_builder --source SWISSPROT:1:resources/swissprot-proteins/uniprot_sprot.fasta \
  --uniform-tag FL --track-all-sources \
  --output curated_reference_01/lake.fasta \
  --output-headers curated_reference_01/source_headers.tsv \
  --output-stats curated_reference_01/build_stats.txt
```

`lake_builder` is included in Docker and the Rust build. Repeat `--source` to
combine chosen sources; each source has a tag, numeric priority and local FASTA
path. Lower numbers have higher priority when the same sequence occurs in
multiple sources. Keep the source-header sidecar and download receipts with the
lake. Hash accessions describe sequence identity; recover biological source names
from the sidecar instead of parsing those hashes as taxonomic labels.

Use `--lake curated_reference_01/lake.fasta` in the manifest runner. The
[worked synthetic construction](https://github.com/MannLabs/fasta-lake/blob/main/demo/run_demo.sh) shows every stage. Initial
global catalogue deduplication is separate from chunked evidence extraction and
has not been shown to fit every worldwide catalogue within the laptop budget.
