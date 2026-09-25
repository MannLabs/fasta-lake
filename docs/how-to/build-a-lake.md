# Build a reference lake

A reference lake defines the sequence universe from which FastaLake can select. Use catalogues appropriate to the study and retain their versions, download dates, file checksums and licences. Include expected host and contaminant sequences deliberately if the scientific design requires them. FastaLake cannot recover a missing reference sequence merely because a spectrum supports it.

## Build one source

One source can be built as follows:

```bash
mkdir lake_build_01
lake_builder --source REF:1:reference.fasta --uniform-tag REF \
  --track-all-sources --output lake_build_01/lake.fasta \
  --output-headers lake_build_01/joint_headers.tsv \
  --output-stats lake_build_01/lake.stats.txt
```

## Source arguments and deduplication

The source argument has the form `TAG:PRIORITY:PATH`. Use separate `--source` arguments for multiple FASTAs. The priority and source-header policy affect retained descriptive provenance. Exact sequence deduplication uses the normalized sequence hash. `--track-all-sources` retains the relation to all observed source headers, which matters when the same sequence occurs in several catalogues. Choose source tags that are meaningful and consistent within a study.

## Sequence normalization

Normalization keeps ASCII letters, converts lowercase to uppercase, and removes other characters before hashing. Ambiguity letters remain literal. This permissive identity contract also removes internal stop markers. A sequence containing a biological stop or a malformed translation therefore needs source-quality review before lake construction. Normalization is reproducible preprocessing, not validation that the translated protein is biologically correct.

## Memory and duplicate accessions

The builder starts with a modest index allocation that grows with the input. A tiny FASTA does not trigger a reference-scale preallocation. Large lakes still require substantial memory. Duplicate original accessions are acceptable across source catalogues when the lake builder replaces them with content identifiers and retains source mappings. Direct parsimony input with duplicate accessions is rejected because an accession-indexed map cannot represent two different sequence records safely.
