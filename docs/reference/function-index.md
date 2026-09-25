# Python function index

For a first reading, use [functions in execution order](function-map.md).
This generated index contains **494 functions and methods**, including private
and nested helpers, from the Python package and the two main workflow tools.
Sections follow workflow purpose; modules and names are alphabetical within them.
Links open the exact definition. The [Python API](api/index.md) renders full
docstrings. Legacy and experimental functions remain listed; their presence
here does not make them part of the default runner.

Regenerate with `python tools/render_function_index.py --write`.
CI checks the index against the source.

## 1. Prepare inputs and plan resources

### fasta_lake.canonical

| Function | Purpose |
|---|---|
| [`_find_local_canonical`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/canonical.py#L106) | Try to find a local canonical FASTA on the cluster. |
| [`_user_agent`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/canonical.py#L32) | Build the request identity, appending the optional configured contact. |
| [`download_canonical`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/canonical.py#L46) | Download the latest canonical SwissProt FASTA for an organism. |
| [`get_canonical`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/canonical.py#L136) | Get path to canonical FASTA, preferring local copy. |

### fasta_lake.capacity

| Function | Purpose |
|---|---|
| [`_cgroup_limits`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/capacity.py#L66) | Read v1/v2 constraints at every ancestor up to the mounted hierarchy root. |
| [`_host_memory`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/capacity.py#L33) | Return total and available host memory in bytes from macOS or procfs. |
| [`_integer`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/capacity.py#L16) | Parse a positive capacity below 2**60, or return None for an unusable limit. |
| [`_read`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/capacity.py#L25) | Read a stripped system-limit file, returning an empty string on I/O failure. |
| [`_unescape_mount`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/capacity.py#L61) | Decode octal escapes in a Linux mount path. |
| [`detect_capacity`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/capacity.py#L128) | Observe capacity available to this process, not the cluster node's full RAM. |
| [`plan_laptop`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/capacity.py#L170) | Choose shared-process threads and reference pieces; never infer jobs per GB. |

### fasta_lake.complexity

| Function | Purpose |
|---|---|
| [`_recommend_strategy`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/complexity.py#L197) | Auto-recommend inference strategy based on complexity metrics. |
| [`ComplexityReport.summary`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/complexity.py#L84) | Return summary dict for JSON stats output. |
| [`estimate_complexity`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/complexity.py#L106) | Estimate sample complexity by mapping peptides to a reference database. |

### fasta_lake.contaminants

| Function | Purpose |
|---|---|
| [`_extract_accession`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/contaminants.py#L142) | Extract UniProt-style accession from various header formats. |
| [`ContaminantDB.__post_init__`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/contaminants.py#L100) | Populate the default contaminant accessions when none were supplied. |
| [`ContaminantDB.is_contaminant`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/contaminants.py#L105) | True if this protein ID corresponds to a contaminant. |
| [`ContaminantDB.load`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/contaminants.py#L118) | Load additional accessions from a text file (one per line). |
| [`ContaminantDB.load_from_fasta`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/contaminants.py#L129) | Extract accessions from a FASTA file. |
| [`exclude_contaminants`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/contaminants.py#L200) | Filter out PSM rows whose protein is a contaminant. |
| [`is_contaminant_id`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/contaminants.py#L193) | Check if a protein_id is a contaminant. Cheap wrapper for FDR filters. |
| [`tag_contaminants`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/contaminants.py#L171) | Rewrite protein IDs with ``CON_`` prefix where they match contaminants. |

### fasta_lake.dianovo

| Function | Purpose |
|---|---|
| [`parse_dianovo_prediction`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/dianovo.py#L16) | Return a complete uppercase peptide and mean confidence, including zeros. |

### fasta_lake.downloads

| Function | Purpose |
|---|---|
| [`_hashes`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/downloads.py#L153) | Stream SHA-256 and an optional expected MD5/SHA-256; reject mismatches. |
| [`_lock`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/downloads.py#L175) | Hold a directory lock for a download and release it on context exit. |
| [`_read_metadata`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/downloads.py#L40) | Fetch at most 2 MiB of resource metadata, raising for a larger response. |
| [`_safe_name`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/downloads.py#L147) | Reject empty names, directory traversal and embedded path separators. |
| [`catalogue`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/downloads.py#L33) | Return the installed, dated resource catalogue without making network calls. |
| [`download_file`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/downloads.py#L191) | Stream one planned file; publish only after size/checksum validation. |
| [`download_resource`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/downloads.py#L289) | Download a resolved resource into its named local folder, with per-file receipts. |
| [`free_bytes`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/downloads.py#L121) | Read free disk space from the closest existing destination parent. |
| [`inspect_remote`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/downloads.py#L84) | Read transfer metadata without downloading the resource body. |
| [`plan_download`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/downloads.py#L131) | Resolve files, transfer bytes, disk availability and expanded-size uncertainty. |
| [`resolve_resource`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/downloads.py#L49) | Resolve a catalogue choice, including live release identity for UniProt. |
| [`unpack_resource`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/downloads.py#L353) | Expand gzip/xz/tar archives into a fresh folder with a hard output-byte budget. |
| [`unpack_resource.copy_stream`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/downloads.py#L374) | Extract one safe archive member within the byte budget and record its hash. |

### fasta_lake.hashing

| Function | Purpose |
|---|---|
| [`normalize_sequence`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/hashing.py#L64) | Normalize an amino-acid sequence per the frozen hashing contract. |
| [`sequence_hash`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/hashing.py#L84) | Return the v5.3 SHA256 hex digest of a normalized amino-acid sequence. |

### fasta_lake.headers.__init__

| Function | Purpose |
|---|---|
| [`_build_gene_map`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/headers/__init__.py#L946) | Build accession → gene and description → gene maps from canonical FASTA. |
| [`_parse_disease_variant`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/headers/__init__.py#L609) | Parse disease variant format: sp\|P02649_v21K\|APOE_VAR description. |
| [`_parse_ens_variant`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/headers/__init__.py#L409) | Parse personal variant format: ens\|ENST_Mut\|GENE description |
| [`_parse_ensembl_pep`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/headers/__init__.py#L503) | Parse Ensembl pep format: ENSP... pep chromosome:... gene:ENSG... gene_symbol:XXX |
| [`_parse_gencode`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/headers/__init__.py#L385) | Parse GENCODE format: ENSP\|ENST\|ENSG\|OTTHUMG\|OTTHUMT\|transcript-NNN\|GENE\|length |
| [`_parse_gnomad`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/headers/__init__.py#L329) | Parse gnomAD variant format: var-POP\|ACC\|GENE_Mut\|pos [sp\|...\|... desc GN=GENE] |
| [`_parse_lake_hybrid`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/headers/__init__.py#L580) | Parse lake_builder output: PROJECT [SWISSPROT] sp\|ACC\|NAME desc GN=GENE md5=... |
| [`_parse_mgyg`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/headers/__init__.py#L473) | Parse UHGP/MGYG format: MGYG000000000_00000 description |
| [`_parse_ncbi_microbe`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/headers/__init__.py#L438) | Parse NCBI microbe format: SPECIES\|WP_ACC\|LOCUS [gene=X] description OS=... |
| [`_parse_prodigal`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/headers/__init__.py#L535) | Parse Prodigal format: 1_1 # start # end # strand # ID=1_1;partial=... |
| [`_parse_sample_ref`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/headers/__init__.py#L551) | Parse SAMPLE\|REF\| prefix format from Stage 2 output. |
| [`_parse_simple_microbe`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/headers/__init__.py#L491) | Parse simple microbe formats: GMGC, SMORF, Assembly. |
| [`_parse_uniprot`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/headers/__init__.py#L297) | Parse UniProt sp\| or tr\| format. |
| [`_resolve_gene`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/headers/__init__.py#L973) | Try to resolve a missing gene name using canonical mapping. |
| [`audit_fasta`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/headers/__init__.py#L1003) | Audit a FASTA file's headers without modifying it. |
| [`generate_decoys`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/headers/__init__.py#L781) | Generate reversed decoy sequences and append to FASTA. |
| [`heal_fasta`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/headers/__init__.py#L664) | Heal a FASTA file for a specific search engine. |
| [`parse_header`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/headers/__init__.py#L146) | Parse a FASTA header line (without leading '>') into a ParsedHeader. |
| [`ParsedHeader.format_for_diann`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/headers/__init__.py#L63) | Format header for DIA-NN compatibility. |
| [`ParsedHeader.format_for_msfragger`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/headers/__init__.py#L120) | Format header for MSFragger/Philosopher compatibility. |
| [`ParsedHeader.format_for_sage`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/headers/__init__.py#L104) | Format header for SAGE compatibility. |
| [`prepare_fasta_for_engine`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/headers/__init__.py#L844) | Prepare a FASTA file for a specific search engine. |
| [`print_audit`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/headers/__init__.py#L1042) | Pretty-print audit results. |

### fasta_lake.headers.__main__

| Function | Purpose |
|---|---|
| [`main`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/headers/__main__.py#L8) | Parse and dispatch the FASTA header healing, audit and preparation commands. |

### fasta_lake.headers.fasta_io

| Function | Purpose |
|---|---|
| [`load_fasta_full`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/headers/fasta_io.py#L21) | Load FASTA file preserving full headers. |
| [`load_fasta_full._flush`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/headers/fasta_io.py#L41) | Store the pending FASTA record, rejecting empty sequences and duplicate IDs. |
| [`load_fasta_proteins`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/headers/fasta_io.py#L80) | Load FASTA file returning protein_id → sequence (backward compatible). |
| [`write_fasta`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/headers/fasta_io.py#L90) | Write FASTA file with full headers. |
| [`write_fasta_simple`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/headers/fasta_io.py#L122) | Write a subset of proteins to FASTA with full headers. |

### fasta_lake.headers.validation

| Function | Purpose |
|---|---|
| [`detect_experiment`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/headers/validation.py#L905) | Auto-detect experiment type from de novo CSVs and mzML filenames. |
| [`print_validation`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/headers/validation.py#L854) | Pretty-print validation results. |
| [`validate_denovo_csv`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/headers/validation.py#L139) | Validate a de novo prediction CSV before processing. |
| [`validate_evidence_lake`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/headers/validation.py#L595) | Validate evidence lake output from fasta_extractor. |
| [`validate_inference_output`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/headers/validation.py#L745) | Validate Stage 3 inference output. |
| [`validate_mzml`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/headers/validation.py#L54) | Validate an mzML file before processing. |
| [`validate_persample_fasta`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/headers/validation.py#L672) | Validate per-sample FASTA from Stage 2. |
| [`validate_pipeline`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/headers/validation.py#L817) | Run all validation checks on available pipeline outputs. |
| [`validate_source_fasta`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/headers/validation.py#L422) | Validate a source FASTA file before lake building. |

### fasta_lake.helpers

| Function | Purpose |
|---|---|
| [`_accession_from_header`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/helpers.py#L184) | Read the first header token, reporting the source line for an empty header. |
| [`_store_record`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/helpers.py#L192) | Refuse duplicate accessions and empty sequences instead of losing records. |
| [`clean_sequence`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/helpers.py#L47) | Clean a protein sequence by uppercasing and removing stop codons. |
| [`collect_peptides_with_scores`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/helpers.py#L238) | Collect peptides with their best AlphaNovo scores from a CSV file. |
| [`count_fasta_headers`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/helpers.py#L136) | Number of '>' records in a FASTA file, read with the file closed afterwards. |
| [`get_db_type`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/helpers.py#L64) | Classify a protein by its database source based on ID prefix. |
| [`get_species_id`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/helpers.py#L93) | Extract a species/genome identifier from a protein ID, or ``None``. |
| [`load_fasta_proteins`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/helpers.py#L142) | Load proteins from a FASTA file into a dict. |
| [`load_fasta_proteins._flush`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/helpers.py#L163) | Pass the pending FASTA record to the shared record validator. |
| [`map_peptides_to_proteins`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/helpers.py#L316) | Map peptides to proteins using Aho-Corasick multi-pattern matching. |
| [`normalize_il`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/helpers.py#L30) | Normalize isoleucine (I) to leucine (L) for mass spec equivalence. |
| [`require_columns`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/helpers.py#L210) | Raise ValueError when a table lacks the columns a reader depends on. |

### fasta_lake.lake_builder

| Function | Purpose |
|---|---|
| [`_parse_gencode_header`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/lake_builder.py#L100) | GENCODE: >ENSP00000269305.4\|ENST...\|ENSG...\|gene_name\|desc |
| [`_parse_generic_header`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/lake_builder.py#L113) | Preserve the first header token as accession and the remaining description. |
| [`_parse_uhgp_header`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/lake_builder.py#L88) | UHGP: >MGYP000123 description |
| [`_parse_uniprot_header`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/lake_builder.py#L55) | Parse UniProt metadata into a source record, or return None for other headers. |
| [`build_lake`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/lake_builder.py#L174) | Build a merged lake FASTA + provenance sidecar from multiple sources. |
| [`iter_fasta`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/lake_builder.py#L147) | Yield (SourceRecord, sequence) pairs from a FASTA, tagging every record |
| [`parse_header`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/lake_builder.py#L120) | Parse a FASTA header line (without leading ``>``). |

### fasta_lake.lake_merge

| Function | Purpose |
|---|---|
| [`_merge_description`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/lake_merge.py#L110) | Prefer the longest *informative* description; fall back to longest. |
| [`_merge_gene`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/lake_merge.py#L102) | First non-empty gene in priority order. |
| [`_merge_organism`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/lake_merge.py#L124) | If all sources agree → that organism. Else → 'ambiguous', list conflicts. |
| [`_merge_pe`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/lake_merge.py#L139) | UniProt PE: LOWER is stronger evidence. Return min across sources. |
| [`_merge_sv`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/lake_merge.py#L145) | Return the highest recorded sequence version, or None when all are absent. |
| [`_merge_tax_ids`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/lake_merge.py#L134) | Collect distinct recorded taxonomy IDs in numeric order. |
| [`_primary_source`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/lake_merge.py#L210) | Return the SourceRecord that won primary status, or None if not found. |
| [`_priority_index`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/lake_merge.py#L94) | Return a source's priority position, sorting unknown tags after known ones. |
| [`_record_to_jsonable`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/lake_merge.py#L280) | Convert MergedRecord to a JSON-serialisable dict (flattens dataclasses). |
| [`format_fasta_header`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/lake_merge.py#L218) | Format the merged FASTA header line (without leading ``>``). |
| [`merge_records`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/lake_merge.py#L151) | Merge N source records that share an identical sequence into one. |
| [`priority_win`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/lake_merge.py#L349) | Legacy behaviour: pick one winner by priority, discard other sources. |
| [`read_provenance_sidecar`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/lake_merge.py#L319) | Read a provenance sidecar JSONL (or JSONL.gz) back into MergedRecords. |
| [`write_provenance_sidecar`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/lake_merge.py#L301) | Write merged records as JSONL (or JSONL.gz if compress=True). |

## 2. Select acquisition databases

### fasta_lake.blosum

| Function | Purpose |
|---|---|
| [`_parse_matrix`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/blosum.py#L100) | Parse a whitespace-delimited substitution matrix into residue-pair scores. |
| [`generate_all_variants`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/blosum.py#L192) | Generate all single-substitution variants for a sequence. |
| [`generate_variants`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/blosum.py#L166) | Generate single-AA variants at a position. |
| [`score`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/blosum.py#L116) | Substitution score for aa1 -> aa2. |
| [`substitutions`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/blosum.py#L130) | Possible substitutions for an amino acid, sorted by score descending. |

### fasta_lake.chunked

| Function | Purpose |
|---|---|
| [`_pieces`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/chunked.py#L53) | Yield complete-record pieces without buffering a whole FASTA record. |
| [`_predictions`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/chunked.py#L37) | Mirror Rust's recursive, case-sensitive CSV roster, including symlink files. |
| [`_sha`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/chunked.py#L22) | Stream a file into a lowercase SHA-256 hex digest. |
| [`_stamp`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/chunked.py#L31) | Record a file's byte count and SHA-256 without loading it into memory. |
| [`extract_chunked`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/chunked.py#L102) | Extract and merge without altering the prepared input lake's identities. |

### fasta_lake.ics

| Function | Purpose |
|---|---|
| [`add_ics`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/ics.py#L88) | Add ICS column to search results. Auto-detects engine from columns. |
| [`add_ics_alphapept`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/ics.py#L77) | Add ICS column to AlphaPept results. |
| [`add_ics_sage`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/ics.py#L42) | Add ICS column to SAGE results. |
| [`clean_sequence`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/ics.py#L26) | Strip modifications, keep uppercase AA only. |
| [`ics`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/ics.py#L35) | Core ICS formula. Returns 0-1. |
| [`summarize`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/ics.py#L124) | ICS summary statistics. |

### fasta_lake.inference.runner

| Function | Purpose |
|---|---|
| [`_rescue_failed_peptides`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/inference/runner.py#L46) | Run mass k-mer anchoring on peptides that failed exact matching. |
| [`run_inference`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/inference/runner.py#L178) | Run protein inference on a single sample. |

### fasta_lake.inference.strategies

| Function | Purpose |
|---|---|
| [`infer_bayesian`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/inference/strategies.py#L541) | Bayesian protein probability: iterative evidence redistribution. |
| [`infer_info_score`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/inference/strategies.py#L416) | Information-content parsimony: keep proteins above an info score threshold. |
| [`infer_razor`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/inference/strategies.py#L209) | Razor parsimony: assign shared peptides to protein with most unique peptides. |
| [`infer_razor._protein_sort_key`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/inference/strategies.py#L241) | Rank the Python razor candidates by support, with lexical accession ties. |
| [`infer_species_budget`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/inference/strategies.py#L25) | Species-budget parsimony: cross-species aware shared peptide resolution. |
| [`infer_species_budget._protein_sort_key`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/inference/strategies.py#L103) | Rank by optional unique-score support, then evidence counts and accession. |
| [`infer_uniform`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/inference/strategies.py#L297) | Uniform parsimony: keep ALL proteins with any peptide evidence. |
| [`infer_uniform_2pep`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/inference/strategies.py#L353) | Uniform 2-peptide: require at least ``min_peptides`` per protein. |
| [`score_lake_bayesian`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/inference/strategies.py#L766) | Score every protein in the evidence lake using Bayesian EM across ALL samples. |

### fasta_lake.refine

| Function | Purpose |
|---|---|
| [`generate_variant_fasta`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/refine.py#L107) | Write a proposal FASTA: one best-scoring substitution per candidate protein. |
| [`identify_refinement_candidates`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/refine.py#L38) | Find low-ICS PSMs that are candidates for variant refinement. |
| [`run_refinement`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/refine.py#L229) | Not implemented: a validated refinement pass does not exist yet. |

## 3. Add molecular evidence

### fasta_lake.multi_omics.__main__

| Function | Purpose |
|---|---|
| [`_add_gate_args`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/multi_omics/__main__.py#L20) | Add options for the historical clustered FL_MO evidence-gate command. |
| [`main`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/multi_omics/__main__.py#L51) | Dispatch the historical FL_MO gate or frozen-median reproduction command. |

### fasta_lake.multi_omics.cli

| Function | Purpose |
|---|---|
| [`add`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/multi_omics/cli.py#L132) | Keep a final de novo FASTA and append DNA/RNA/metabolite-supported proteins. |
| [`checked`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/multi_omics/cli.py#L8) | Translate expected input and execution failures into readable Click errors. |
| [`compare_searches`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/multi_omics/cli.py#L161) | Compare exact proteins, accepted peptides and LFQ in two matched Sage searches. |
| [`increments`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/multi_omics/cli.py#L197) | List gained/lost peptide and protein support against one complete reference. |
| [`molecular`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/multi_omics/cli.py#L17) | Prepare TPM inputs, select exact proteins and compare matched searches. |
| [`prepare_source`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/multi_omics/cli.py#L39) | Check and copy specimen inputs into a new portable source directory. |
| [`select`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/multi_omics/cli.py#L144) | Build a standalone molecular FASTA with an optional common background. |
| [`selection_options`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/multi_omics/cli.py#L71) | Build a Click decorator for independent DNA/RNA selection thresholds. |
| [`selection_options.decorate`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/multi_omics/cli.py#L73) | Attach source, specimen, output and molecular-selection options in order. |
| [`validate_source`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/multi_omics/cli.py#L48) | Check source hashes, gene agreement, annotation links and declared scope. |

### fasta_lake.multi_omics.complementarity

| Function | Purpose |
|---|---|
| [`compare`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/multi_omics/complementarity.py#L232) | Compare two searches of the same acquisition with identical non-file settings. |
| [`correlation`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/multi_omics/complementarity.py#L220) | Return Pearson correlation, or None for fewer than two or constant values. |
| [`finite`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/multi_omics/complementarity.py#L31) | Convert a value to float and reject NaN or infinity with a field label. |
| [`load_search`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/multi_omics/complementarity.py#L64) | Read accepted target evidence and positive LFQ features for one acquisition. |
| [`peptide`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/multi_omics/complementarity.py#L23) | Validate a Sage peptide and remove modification syntax before I/L collapse. |
| [`rank`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/multi_omics/complementarity.py#L205) | Return one-based ranks in input order, averaging the ranks of tied values. |
| [`search_inputs`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/multi_omics/complementarity.py#L39) | Resolve a single-acquisition search and separate settings from file paths. |

### fasta_lake.multi_omics.evidence_gate

| Function | Purpose |
|---|---|
| [`build_evidence_lookup`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/multi_omics/evidence_gate.py#L165) | Build the per-sample evidence attribution table (the FL_MO core). |
| [`EvidenceRow.has_tpm`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/multi_omics/evidence_gate.py#L98) | Report whether either historical DNA or RNA support flag is set. |
| [`EvidenceRow.in_gate`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/multi_omics/evidence_gate.py#L102) | Whether this representative is admitted under the given gate. |
| [`EvidenceRow.to_tsv`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/multi_omics/evidence_gate.py#L112) | Serialize this legacy gate row, formatting TPM to four decimal places. |
| [`gate_members`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/multi_omics/evidence_gate.py#L281) | Return the set of cluster_rep_ids admitted under `gate`. |
| [`passes_tpm`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/multi_omics/evidence_gate.py#L156) | Return (passes_metaG_top1, passes_metaT_1) for a single assembly protein. |
| [`top1_threshold`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/multi_omics/evidence_gate.py#L128) | Per-sample metaG cutoff: smallest value v* such that values >= v* are the |

### fasta_lake.multi_omics.exact

| Function | Purpose |
|---|---|
| [`build_union`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/multi_omics/exact.py#L183) | Keep every baseline sequence and append independently selected molecular proteins. |
| [`load_source`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/multi_omics/exact.py#L144) | Load one explicitly specimen-qualified quantification reference. |
| [`numeric`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/multi_omics/exact.py#L66) | Parse a nonnegative finite TPM while preserving explicit missing tokens. |
| [`read_fasta`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/multi_omics/exact.py#L29) | Return accession -> (full header, normalized sequence, SHA-256). |
| [`read_fasta.save`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/multi_omics/exact.py#L37) | Store a normalized record with its exact-sequence hash, enforcing ID limits. |
| [`read_tpm`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/multi_omics/exact.py#L84) | Read plain or gzipped gene TPM rows without converting missing values to zero. |
| [`sha256`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/multi_omics/exact.py#L20) | Hash file bytes in 8 MiB blocks and return the SHA-256 hex digest. |
| [`top_cutoff`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/multi_omics/exact.py#L114) | Top percentage among positive GENES; zero disables, ties all included. |
| [`write_table`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/multi_omics/exact.py#L130) | Write dictionary rows as a TSV, using deterministic gzip for a .gz path. |

### fasta_lake.multi_omics.increments

| Function | Purpose |
|---|---|
| [`evidence_tables`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/multi_omics/increments.py#L77) | Return exhaustive support rows and bidirectional counts. |
| [`increment_report`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/multi_omics/increments.py#L212) | Validate matched searches and export complete reference-wide evidence. |
| [`reference_map`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/multi_omics/increments.py#L53) | Map I/L-canonical peptides to distinct exact protein sequences. |

### fasta_lake.multi_omics.reproduce

| Function | Purpose |
|---|---|
| [`per_sample_median`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/multi_omics/reproduce.py#L48) | Return (median, n_samples) of an integer depth column in a per-sample TSV. |
| [`reproduce`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/multi_omics/reproduce.py#L56) | Recompute and check the published FL_MO razor median. |

### fasta_lake.multi_omics.runner

| Function | Purpose |
|---|---|
| [`_denovo_rep_ids`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/multi_omics/runner.py#L124) | First-token header ids from the de-novo extracted FASTA. |
| [`_f`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/multi_omics/runner.py#L46) | Parse a legacy TPM value, treating blank or malformed text as zero. |
| [`_pull_rep_seqs`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/multi_omics/runner.py#L134) | Stream the rep FASTA, returning {rep_id: raw_sequence_bytes} for needed. |
| [`index_fasta`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/multi_omics/runner.py#L83) | Return {first_whitespace_token_of_header: sequence_str} for a FASTA. |
| [`load_hash_to_rep`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/multi_omics/runner.py#L102) | Load the FL_MO hash_to_rep table: {sha256_hex: cluster_rep_id}. |
| [`load_tpm_tsv`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/multi_omics/runner.py#L58) | Yield (assembly_protein_id, metaG_tpm, metaT_tpm) from a TPM TSV. |
| [`run_fl_mo_sample`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/multi_omics/runner.py#L154) | Run the FL_MO gate for one sample and write the augmented FASTA + lookup. |

### fasta_lake.multi_omics.selection

| Function | Purpose |
|---|---|
| [`_targets`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/multi_omics/selection.py#L37) | Read target-only proteins and index exact sequences by SHA-256. |
| [`build_selection`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/multi_omics/selection.py#L54) | Select gene-supported sequences and write a traceable target FASTA union. |
| [`validate_rules`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/multi_omics/selection.py#L16) | Validate independent DNA/RNA selection modes without loading any inputs. |

### fasta_lake.multi_omics.sources

| Function | Purpose |
|---|---|
| [`_inputs`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/multi_omics/sources.py#L81) | Resolve source-manifest paths and verify required inputs and pinned hashes. |
| [`_load_v2`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/multi_omics/sources.py#L108) | Validate a specimen-qualified v2 source, its sequence coverage and annotations. |
| [`load_reference`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/multi_omics/sources.py#L182) | Load a v1/v2 source.json and verify its pinned files and reference contract. |
| [`MolecularSource.available`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/multi_omics/sources.py#L53) | Return sorted genes with both a TPM row and a protein sequence. |
| [`MolecularSource.coverage`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/multi_omics/sources.py#L58) | Summarize available sequences, unavailable genes and missing molecular values. |
| [`prepare_source`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/multi_omics/sources.py#L233) | Package checked specimen-matched molecular inputs for portable reuse. |
| [`table_rows`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/multi_omics/sources.py#L16) | Yield string-valued rows from plain/gzipped TSV with a complete unique header. |
| [`write_coverage`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/multi_omics/sources.py#L201) | Write missing-gene and measurement coverage to a caller-owned output directory. |
| [`write_reference`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/multi_omics/sources.py#L317) | Write every available exact sequence, including genes with zero/missing TPM. |

## 4. Search, group and quantify

### fasta_lake.aggregate

| Function | Purpose |
|---|---|
| [`build_quant_matrix`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/aggregate.py#L94) | Build a quantification matrix across samples. |
| [`compute_cross_sample_stats`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/aggregate.py#L138) | Compute cross-sample statistics. |
| [`load_diann_sample`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/aggregate.py#L47) | Load a single-sample DIA-NN report. |
| [`print_cohort_summary`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/aggregate.py#L271) | Print a summary table of all samples. |
| [`run_rust_aggregation`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/aggregate.py#L196) | Run the Rust aggregation_engine binary. |

### fasta_lake.quantification

| Function | Purpose |
|---|---|
| [`_baseline`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/quantification.py#L53) | Validate the sum matrix without requiring numerical packages for sum mode. |
| [`_directlfq`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/quantification.py#L186) | Run directLFQ on ions assigned to the fixed study dictionary. |
| [`_peptide_support`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/quantification.py#L356) | Write observed peptide support before/after estimation, one acquisition at a time. |
| [`_read_frame`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/quantification.py#L179) | Read a TSV with pandas while preserving identifier-like NA strings. |
| [`quantify_study`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/quantification.py#L73) | Quantify an existing study dictionary without assigning peptides again. |
| [`require_directlfq`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/quantification.py#L39) | Fail before running a workflow if the optional backend is unavailable. |
| [`stamp`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/quantification.py#L30) | Return basename, byte size and a streamed SHA256 for a pathlib.Path input. |

### fasta_lake.search

| Function | Purpose |
|---|---|
| [`_dda_config`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/search.py#L13) | The explicitly versioned DDA settings used in the release acceptance runs. |
| [`sage_config`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/search.py#L60) | Build the complete configuration for one Sage 0.14.6 acquisition. |

### fasta_lake.study

| Function | Purpose |
|---|---|
| [`_members`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/study.py#L205) | Collect distinct nonempty accessions, excluding the rev_ decoy prefix. |
| [`_passing`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/study.py#L121) | Accept only a finite q value between zero and the inclusive threshold. |
| [`_peptide_evidence`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/study.py#L533) | Retain accepted identifications, including peptides without positive LFQ features. |
| [`_quantify`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/study.py#L398) | Assign positive finite LFQ features to the fixed peptide-to-group dictionary. |
| [`_read_fasta`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/study.py#L165) | Verify accession identity before merging acquisition evidence. |
| [`_read_fasta.retain`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/study.py#L170) | Merge a nonempty FASTA record after checking local IDs and sequence identity. |
| [`_rows`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/study.py#L135) | Hash the actual bytes read, then parse strict, tab-delimited records. |
| [`_rows.lines`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/study.py#L140) | Decode TSV lines while hashing exactly the bytes consumed by the parser. |
| [`_sage_filename`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/study.py#L130) | Sage emits the basename including mzML extension; accept recorded full paths. |
| [`_write_table`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/study.py#L210) | Create a TSV exclusively, preserving the caller's column and row order. |
| [`canonical_peptide`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/study.py#L26) | Remove Sage modifications and use the declared I/L equivalence. |
| [`greedy_set_cover`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/study.py#L40) | Choose study representatives by maximum uncovered peptide gain. |
| [`group_searches`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/study.py#L218) | Build study groups and sum quantities from completed acquisition searches. |
| [`member_set_key`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/study.py#L35) | Preserve the full distinct accession membership, without source-tag stripping. |
| [`study_dictionary`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/study.py#L76) | Assign accepted peptides to a fixed set of study representatives. |

## 5. Inspect and interpret results

### fasta_lake.annotation

| Function | Purpose |
|---|---|
| [`annotate_fasta`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/annotation.py#L116) | Run eggNOG-mapper 2.1.12 on exact unique selected sequences. |
| [`annotate_study`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/annotation.py#L325) | Annotate the exact representatives exported by a completed study grouping. |
| [`prepare_annotation_targets`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/annotation.py#L27) | Deduplicate a selected FASTA by exact sequence; retain every input alias. |
| [`read_annotation_settings`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/annotation.py#L244) | Validate optional annotation settings; resolve local paths beside the JSON file. |
| [`read_eggnog_annotations`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/annotation.py#L62) | Read v2 eggNOG output and retain no-hit and hit-without-function queries. |

### fasta_lake.community_plots

| Function | Purpose |
|---|---|
| [`plot_community_context`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/community_plots.py#L21) | Plot both count and summed-signal fractions; keep empty acquisitions visible. |

### fasta_lake.downstream

| Function | Purpose |
|---|---|
| [`_names`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/downstream.py#L51) | Read unique nonblank IDs and require them to belong to the expected roster. |
| [`_plot`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/downstream.py#L374) | Save descriptive QC panels and the eligible PCA result using local plot styles. |
| [`_plot.draw_pca`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/downstream.py#L411) | Draw available finite PCA scores, or show why no PCA coordinates were saved. |
| [`analyze_study`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/downstream.py#L63) | Summarize study coverage and plot eligible PCA without imputing quantities. |
| [`require_analysis`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/downstream.py#L24) | Load the pinned QC backend and return AlphaPeptTools, AnnData, NumPy and pandas. |

### fasta_lake.embeddings

| Function | Purpose |
|---|---|
| [`embed_unknown`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/embeddings.py#L65) | Embed every unknown from a completed, unchanged FastaLake annotation run. |
| [`residue_embedding_sum`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/embeddings.py#L29) | Sum residue embeddings in float64, excluding the ESM-C BOS/EOS tokens. |
| [`resolve_checkpoint`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/embeddings.py#L45) | Validate an explicit local checkpoint against its published upstream SHA256. |
| [`sequence_windows`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/embeddings.py#L19) | Yield nonoverlapping full-coverage windows; never silently truncate a protein. |

### fasta_lake.exploration

| Function | Purpose |
|---|---|
| [`_plot`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/exploration.py#L228) | Render coverage beside annotated-feature diversity; label experimental scope. |
| [`explore_study`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/exploration.py#L37) | Export functional signal, coverage, diversity, Bray–Curtis and CLR PCA. |
| [`read_group_terms`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/exploration.py#L21) | Read a long table with study_group and term columns; deduplicate exact pairs. |

### fasta_lake.metaproteomics

| Function | Purpose |
|---|---|
| [`_matrix`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/metaproteomics.py#L24) | Validate a nonempty feature-by-acquisition frame and return numeric values. |
| [`_within_group_ss`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/metaproteomics.py#L169) | Sum within-group squared distances, dividing each pair by its group size. |
| [`aggregate_terms`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/metaproteomics.py#L44) | Split observed group signal equally across distinct declared functional terms. |
| [`bray_curtis`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/metaproteomics.py#L131) | Return a labelled Bray–Curtis distance matrix and explicitly excluded samples. |
| [`clr_complete`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/metaproteomics.py#L152) | Return sample × feature natural-log CLR on the complete positive core. |
| [`diversity_summary`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/metaproteomics.py#L107) | Return observed richness, Shannon entropy and effective Shannon diversity. |
| [`permanova`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/metaproteomics.py#L182) | One-way experimental PERMANOVA with optional within-block label permutations. |
| [`permanova.align`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/metaproteomics.py#L213) | Align complete acquisition labels to the distance matrix and encode groups. |
| [`permanova.statistic`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/metaproteomics.py#L256) | Compute the PERMANOVA pseudo-F statistic for the supplied group labels. |
| [`permanova.within`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/metaproteomics.py#L245) | Compute within-group dispersion for one observed or permuted label vector. |
| [`permanova.within`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/metaproteomics.py#L250) | Compute within-group dispersion for one observed or permuted label vector. |
| [`relative_signal`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/metaproteomics.py#L96) | Close observed signal to one per nonempty acquisition; preserve empty columns. |

### fasta_lake.network_annotations

| Function | Purpose |
|---|---|
| [`annotate_evidence_network`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/network_annotations.py#L145) | Attach annotations to exact sequences from a recorded study FASTA. |
| [`function_terms`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/network_annotations.py#L120) | Split one annotation namespace without creating terms for missing values. |
| [`load_eggnog_annotations`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/network_annotations.py#L43) | Read a completed FastaLake annotation folder or native eggNOG v2 output. |

### fasta_lake.networks

| Function | Purpose |
|---|---|
| [`_identity`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/networks.py#L53) | Validate a nonblank identifier without padding or control characters. |
| [`_integer`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/networks.py#L42) | Parse an integer at least minimum, rejecting noncanonical numeric text. |
| [`_rows`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/networks.py#L27) | Yield exact-schema TSV rows, enforcing the caller's row limit. |
| [`_sha`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/networks.py#L18) | Return a streamed SHA-256 digest of the evidence file. |
| [`add_function_annotations`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/networks.py#L395) | Extend a bounded view with typed orthology edges without changing assignments. |
| [`export_evidence_network`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/networks.py#L254) | Write complete typed tables and a checked portable graph JSON in a new folder. |
| [`load_evidence_network`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/networks.py#L310) | Read an unchanged, completed network export; reject partial or modified bundles. |
| [`load_network_annotations`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/networks.py#L378) | Read an unchanged functional layer belonging to the requested evidence graph. |
| [`read_evidence_network`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/networks.py#L65) | Read and reconcile a schema-1 study's complete reported-membership graph. |
| [`read_evidence_network.node_id`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/networks.py#L193) | Derive a kind-prefixed node ID from the graph identity and node label. |
| [`render_ckg_network`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/networks.py#L442) | Render a portable export with the user's separate compatible CKG environment. |
| [`select_network_view`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/networks.py#L330) | Select a deterministic bounded view, prioritizing shared peptide memberships. |

### fasta_lake.proteomics.cascade

| Function | Purpose |
|---|---|
| [`_canonical_peptide_set`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/proteomics/cascade.py#L22) | Pool tryptic peptides from the canonical FASTA under the supplied digest limits. |
| [`build_smart_tier`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/proteomics/cascade.py#L30) | canonical backbone + smart isoforms + Ig + variants, gene-healed. |

### fasta_lake.proteomics.grouping

| Function | Purpose |
|---|---|
| [`_norm_name`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/proteomics/grouping.py#L145) | Lowercase a protein name and remove fragment and isoform qualifiers. |
| [`_set_gene`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/proteomics/grouping.py#L243) | Replace the first GN= token, or append it when absent. |
| [`assign_gene`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/proteomics/grouping.py#L152) | Return (gene, method). |
| [`build_gene_index`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/proteomics/grouping.py#L124) | Build the discriminative peptide->gene index from canonical SwissProt. |
| [`heal_fasta_genes`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/proteomics/grouping.py#L196) | Rewrite GN= in input_fasta with peptide-voted real gene symbols. Returns stats. |
| [`is_accession_gene`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/proteomics/grouping.py#L83) | True if GN= is actually a UniProt accession (i.e. a fake gene). |
| [`is_ig`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/proteomics/grouping.py#L44) | Test whether the header matches the configured immunoglobulin pattern. |
| [`iter_fasta`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/proteomics/grouping.py#L49) | Yield (header_without_gt, sequence). |
| [`parse_gene`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/proteomics/grouping.py#L64) | Return the first GN= gene token, or None when the header has none. |
| [`parse_protein_name`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/proteomics/grouping.py#L70) | Protein description between the id token and the first KEY= tag. |
| [`tryptic_peptides`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/proteomics/grouping.py#L88) | Cut after K/R (no proline rule — matches the run's --cut 'K*,R*'). I/L folded. |

### fasta_lake.qc_plots

| Function | Purpose |
|---|---|
| [`intensity_distribution_summary`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/qc_plots.py#L11) | Summarise measured positive reporting-group intensities per acquisition. |
| [`plot_intensity_distributions`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/qc_plots.py#L47) | Save paginated violins, missingness bars and their exact summary table. |

### fasta_lake.quality

| Function | Purpose |
|---|---|
| [`compute_sequence_entropy`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/quality.py#L118) | Compute normalized Shannon entropy of a protein sequence. |
| [`filter_sequences_by_quality`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/quality.py#L266) | Apply sequence quality filters to a protein dictionary. |
| [`has_poly_run`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/quality.py#L76) | Check if a sequence contains a homopolymeric run. |
| [`is_fragment`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/quality.py#L168) | Check if a protein is likely a fragment. |
| [`noncanonical_fraction`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/quality.py#L68) | Fraction of residues outside the 22 genetically encoded amino acids. |
| [`normalise_residues`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/quality.py#L44) | Clean a protein sequence: strip trailing stop symbols, replace non-canonical |
| [`QualityReport.summary`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/quality.py#L244) | Return summary dict for JSON stats output. |

### fasta_lake.taxonomy

| Function | Purpose |
|---|---|
| [`_balance`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/taxonomy.py#L218) | Classify an assigned lineage as host, bacteria, archaea or other; retain unresolved. |
| [`_write_report`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/taxonomy.py#L441) | Embed upstream Unipept visualisations and data for an offline HTML report. |
| [`build_taxonomy_tree`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/taxonomy.py#L179) | Make Unipept's count/selfCount tree with each peptide's weight counted once. |
| [`canonical_peptide`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/taxonomy.py#L70) | Use the study's I/L-equivalent peptide identity, without removing mass gaps. |
| [`prepare_unipept_queries`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/taxonomy.py#L232) | Export the checked study peptide roster for the upstream Unipept CLI. |
| [`read_unipept_output`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/taxonomy.py#L81) | Read native pept2lca --all JSON/CSV; keep absent and truncated lookups visible. |
| [`taxonomy_report`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/taxonomy.py#L253) | Create local taxonomy/holobiont plots from Unipept output, without web requests. |

### fasta_lake.taxonomy_reference

| Function | Purpose |
|---|---|
| [`prepare_unipept_reference`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/taxonomy_reference.py#L17) | Write sa-builder's four-field TSV from a selected protein FASTA. |

## 6. Run the workflow and configure tools

### fasta_lake.annotation_cli

| Function | Purpose |
|---|---|
| [`annotation`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/annotation_cli.py#L10) | Annotate selected proteins and inspect explicit annotation gaps. |
| [`eggnog_cmd`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/annotation_cli.py#L57) | Run eggNOG-mapper 2.1.12 with versioned inputs and complete coverage records. |
| [`embed_cmd`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/annotation_cli.py#L133) | Embed the recorded unknown roster; preserve every residue through windowing. |
| [`study_cmd`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/annotation_cli.py#L160) | Annotate a study's verified reporting representatives, with optional ESM-C. |

### fasta_lake.cli

| Function | Purpose |
|---|---|
| [`_maybe_inject_entrapment`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/cli.py#L1451) | Apply a config's entrapment section to a database, returning what to search. |
| [`_setup_logging`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/cli.py#L35) | Configure logging for CLI usage. |
| [`analyze_study_cmd`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/cli.py#L1890) | Export AnnData and run AlphaPeptTools QC, log2 PCA and optional replicate CV. |
| [`cli`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/cli.py#L48) | FASTA Lake: Sample-specific protein database construction for metaproteomics. |
| [`contam_list`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/cli.py#L133) | List the bundled contaminant accessions. |
| [`contam_tag`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/cli.py#L166) | Scan a FASTA, prefix CON_ on headers matching the contaminant list. |
| [`contaminants_grp`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/cli.py#L122) | Inspect or apply the contaminant database (v5.1). |
| [`convert`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/cli.py#L1196) | Convert raw mass spec data to de novo predictions. |
| [`diagnose_cmd`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/cli.py#L895) | Produce FDR calibration curve + outlier sample list. |
| [`entrapment_classes_cmd`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/cli.py#L256) | List entrapment classes with their difficulty and interpretation. |
| [`entrapment_fdp`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/cli.py#L624) | Compute the FDP from search results, with n, median, IQR and bootstrap CI. |
| [`entrapment_generate`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/cli.py#L384) | Produce an entrapment FASTA, deterministically. |
| [`entrapment_grp`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/cli.py#L211) | Generate, inject and measure entrapment sequences. |
| [`entrapment_inject`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/cli.py#L490) | Concatenate entrapment into a database and write a manifest. |
| [`entrapment_prefixes`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/cli.py#L296) | Report the entrapment prefixes actually present in a FASTA. |
| [`entrapment_prepare`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/cli.py#L764) | Generate + inject entrapment as configured by a fastalake.json. |
| [`evidence_lake`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/cli.py#L1352) | Stage 1: Create evidence lake from peptide predictions. |
| [`export`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/cli.py#L1129) | Prepare a FASTA for a specific search engine. |
| [`extract_chunked_cmd`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/cli.py#L1944) | Build evidence in bounded reference pieces and merge before inference. |
| [`extract_sample`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/cli.py#L1408) | Stage 2: Extract per-sample protein database. |
| [`fl_mo_gate`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/cli.py#L1776) | Run the FL_MO gate for one sample at Stage 3. |
| [`fl_mo_grp`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/cli.py#L1722) | Legacy clustered FL_MO reproduction (archived Methods §M5.5). |
| [`fl_mo_reproduce`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/cli.py#L1817) | Reproduce the published FL_MO razor per-sample median (13,776). |
| [`infer`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/cli.py#L1027) | Experimental Python inference, distinct from the validated Rust workflow. |
| [`init`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/cli.py#L1518) | Generate an example fastalake.json configuration. |
| [`lake_build`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/cli.py#L1326) | Stage 0: Build reference protein lake from databases. |
| [`lake_inspect_cmd`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/cli.py#L1642) | Query a lake provenance sidecar produced by lake-merge. |
| [`lake_merge_cmd`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/cli.py#L1563) | Merge multiple source FASTAs into one lake with field-level annotation merge (v5.2). |
| [`plan_resources_cmd`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/cli.py#L1989) | Show laptop planning choices without starting a workflow. |
| [`preset_grp`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/cli.py#L58) | Inspect or export validated parameter presets (v5). |
| [`preset_list`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/cli.py#L67) | List registered parameter presets and their search thresholds. |
| [`preset_sage`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/cli.py#L96) | Emit a SAGE config with a preset applied. |
| [`preset_show`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/cli.py#L78) | Show all parameters for a preset. |
| [`quantify_study_cmd`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/cli.py#L1842) | Quantify completed study groups and retain peptide/ion support diagnostics. |
| [`tune_cmd`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/cli.py#L821) | Recommend SAGE parameters by post-hoc sweep over existing entrapment results. |
| [`validate`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/cli.py#L1497) | Validate a FASTA Lake JSON configuration file. |

### fasta_lake.config

| Function | Purpose |
|---|---|
| [`ClusteringConfig.__post_init__`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/config.py#L413) | Require identity and coverage fractions in the interval (0, 1]. |
| [`ClusteringConfig.min_seq_id`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/config.py#L430) | The value to pass as MMseqs2 ``--min-seq-id``. |
| [`ClusteringConfig.mmseqs_args`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/config.py#L434) | The MMseqs2 ``easy-cluster`` argument list driven by this config. |
| [`ClusteringConfig.skip_clustering`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/config.py#L421) | True when clustering is a no-op for this config. |
| [`create_example_config`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/config.py#L674) | Generate an example JSON configuration string. |
| [`DatabaseConfig.__post_init__`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/config.py#L106) | Convert source dictionaries into validated DatabaseSource objects. |
| [`DatabaseSource.__post_init__`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/config.py#L78) | Normalize the source path and reject unsupported source categories. |
| [`EntrapmentConfig.__post_init__`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/config.py#L372) | Normalize the injection stage and validate seed, count and FDR settings. |
| [`FastaLakeConfig.__post_init__`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/config.py#L573) | Normalize and validate acquisition and search-engine names. |
| [`FastaLakeConfig.get_sources_by_type`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/config.py#L620) | Return database sources filtered by type. |
| [`FastaLakeConfig.has_contaminants`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/config.py#L632) | Check if contaminant databases are configured. |
| [`FastaLakeConfig.has_entrapment`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/config.py#L624) | Check if entrapment is configured (via config or source tags). |
| [`FastaLakeConfig.has_variants`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/config.py#L628) | Check if variant databases are configured. |
| [`FastaLakeConfig.validate`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/config.py#L611) | Validate that required fields are set and paths are sensible. |
| [`InferenceConfig.__post_init__`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/config.py#L289) | Reject an inference strategy absent from the supported strategy registry. |
| [`load_config`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/config.py#L642) | Load a FASTA Lake configuration from a JSON file. |
| [`OutputConfig.__post_init__`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/config.py#L488) | Normalize the configured output path to a string. |

### fasta_lake.download_cli

| Function | Purpose |
|---|---|
| [`download_cmd`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/download_cli.py#L86) | Show the plan, then explicitly download one chosen resource; do not unpack. |
| [`list_cmd`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/download_cli.py#L32) | Show installed choices and dated sizes; this command is offline. |
| [`plan_cmd`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/download_cli.py#L70) | Show release, download bytes and disk needs without downloading data. |
| [`resources`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/download_cli.py#L26) | Choose, size and download public reference databases and model weights. |
| [`show_plan`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/download_cli.py#L50) | Print the selected resource, destinations, size estimates and checksums. |
| [`size_label`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/download_cli.py#L20) | Format bytes as GiB and exact bytes; preserve an unknown estimate as text. |
| [`unpack_cmd`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/download_cli.py#L117) | Unpack archives into a new directory within an explicit disk budget. |

### fasta_lake.exploration_cli

| Function | Purpose |
|---|---|
| [`explore_cmd`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/exploration_cli.py#L32) | EXPERIMENTAL: inspect annotated signal, coverage, diversity and CLR PCA. |

### fasta_lake.gzip_io

| Function | Purpose |
|---|---|
| [`open_gzip_text`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/gzip_io.py#L22) | Open ``path`` for writing gzip-compressed text with a zero MTIME header. |

### fasta_lake.network_cli

| Function | Purpose |
|---|---|
| [`eggnog_cmd`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/network_cli.py#L111) | Attach eggNOG annotations by exact protein sequence; retain unknowns. |
| [`export_cmd`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/network_cli.py#L21) | Export complete accepted memberships from a group-study result. |
| [`network`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/network_cli.py#L10) | Export peptide memberships and render experimental CKG network views. |
| [`render_cmd`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/network_cli.py#L61) | Create an experimental CKG report from a completed export. |

### fasta_lake.presets

| Function | Purpose |
|---|---|
| [`apply_preset`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/presets.py#L153) | Apply preset fields to an existing RunConfig-like object. |
| [`describe_preset`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/presets.py#L178) | Human-readable description of what a preset does. |
| [`get_preset`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/presets.py#L145) | Return the named preset or raise ValueError listing the available names. |
| [`Preset.as_dict`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/presets.py#L50) | Return a dataclass dictionary of this preset's settings. |
| [`preset_sage_config`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/presets.py#L164) | Apply preset to a SAGE JSON config dict. Returns modified copy. |

### fasta_lake.provenance

| Function | Purpose |
|---|---|
| [`summarise_provenance`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/provenance.py#L113) | Aggregate a provenance JSONL into summary counts. |
| [`write_provenance`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/provenance.py#L38) | Write per-protein provenance sidecar. |

### fasta_lake.resources

| Function | Purpose |
|---|---|
| [`local_interpreter`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/resources.py#L10) | Keep a Python environment's launch path, including its executable symlink. |
| [`run_recorded`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/resources.py#L23) | Record POSIX child usage when available. |

### fasta_lake.run_config

| Function | Purpose |
|---|---|
| [`RunConfig.capture_environment`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/run_config.py#L150) | Auto-populate environment fields (git commit, binary paths, etc.). |
| [`RunConfig.load`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/run_config.py#L212) | Load config from JSON/YAML. |
| [`RunConfig.max_length`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/run_config.py#L117) | Read the legacy length alias, or set both de novo and search maxima. |
| [`RunConfig.max_length`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/run_config.py#L122) | Read the legacy length alias, or set both de novo and search maxima. |
| [`RunConfig.min_length`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/run_config.py#L105) | Read the legacy length alias, or set both de novo and search minima. |
| [`RunConfig.min_length`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/run_config.py#L110) | Read the legacy length alias, or set both de novo and search minima. |
| [`RunConfig.save`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/run_config.py#L202) | Save config to YAML. |
| [`RunConfig.summary`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/run_config.py#L241) | Human-readable summary. |
| [`RunConfig.update_cluster`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/run_config.py#L226) | Update config with clustering results. |
| [`RunConfig.update_evidence`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/run_config.py#L219) | Update config with evidence lake results. |
| [`RunConfig.update_sample`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/run_config.py#L230) | Update config with per-sample extraction results. |
| [`RunConfig.update_search`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/run_config.py#L235) | Update config with search results for one sample. |

### fasta_lake.rust.binaries

| Function | Purpose |
|---|---|
| [`find_binary`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/rust/binaries.py#L29) | Locate a compiled Rust binary. |
| [`run_fasta_extractor`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/rust/binaries.py#L89) | Run the Rust ``fasta_extractor`` binary. |
| [`run_lake_builder`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/rust/binaries.py#L198) | Run the Rust ``lake_builder`` binary. |

### fasta_lake.taxonomy_cli

| Function | Purpose |
|---|---|
| [`prepare`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/taxonomy_cli.py#L19) | Write canonical peptides for Unipept pept2lca --all --equate. |
| [`prepare_reference`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/taxonomy_cli.py#L80) | Experimental: export a selected FASTA for Unipept's native index builder. |
| [`report`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/taxonomy_cli.py#L50) | Render upstream sunburst/treemap/tree with explicit counts, signal and citations. |
| [`taxonomy`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/taxonomy_cli.py#L10) | Export Unipept queries and plot its saved taxonomic assignments locally. |

### tools.group_study

| Function | Purpose |
|---|---|
| [`main`](https://github.com/MannLabs/fasta-lake/blob/main/tools/group_study.py#L14) | Parse grouping options, validate optional backends, then group and report the study. |

### tools.run_manifest

| Function | Purpose |
|---|---|
| [`load_manifest`](https://github.com/MannLabs/fasta-lake/blob/main/tools/run_manifest.py#L41) | Validate acquisitions and resolve inputs relative to the manifest. |
| [`main`](https://github.com/MannLabs/fasta-lake/blob/main/tools/run_manifest.py#L123) | Validate and execute the acquisition-manifest workflow. |
| [`main.execute`](https://github.com/MannLabs/fasta-lake/blob/main/tools/run_manifest.py#L370) | Run one command with resource measurements and append its result to commands.json. |
| [`main.extract`](https://github.com/MannLabs/fasta-lake/blob/main/tools/run_manifest.py#L422) | Run whole-reference or chunked exact lookup with the same declared evidence settings. |
| [`sha256`](https://github.com/MannLabs/fasta-lake/blob/main/tools/run_manifest.py#L32) | Stream file bytes into a SHA-256 digest for the run's provenance. |

## 7. Diagnostics and comparisons

### fasta_lake._ckg_worker

| Function | Purpose |
|---|---|
| [`_json_for_script`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/_ckg_worker.py#L65) | Serialize data safely for a JavaScript string literal in an HTML script. |
| [`component_layout`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/_ckg_worker.py#L77) | Keep disconnected membership components readable in a deterministic grid. |
| [`remove_unused_bootstrap`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/_ckg_worker.py#L46) | Remove PyVis's unused remote Bootstrap assets from its inline template. |
| [`render`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/_ckg_worker.py#L101) | Plot reported membership and topological communities without re-inferring proteins. |

### fasta_lake.bags

| Function | Purpose |
|---|---|
| [`bag_plan`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/bags.py#L141) | Validate recorded spans and choose an anchor for all-hit streaming lookup. |
| [`compile_queries`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/bags.py#L212) | Compile native AlphaNovo alternatives without turning variants into observations. |
| [`compile_queries.write`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/bags.py#L261) | Emit one route-tagged pattern, failing before truncating a query. |
| [`control_form`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/bags.py#L195) | Shuffle residues reproducibly while preserving bag spans and modifications. |
| [`dump`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/bags.py#L48) | Atomically replace a JSON receipt, rejecting non-finite numbers. |
| [`evidence_for`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/bags.py#L329) | Return original-query support only when the arm admits this protein. |
| [`fnv1a`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/bags.py#L321) | Compute the production unsigned 64-bit FNV-1a accession tie-breaker. |
| [`load_queries`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/bags.py#L89) | Keep original row ranking/ties; collapse selection evidence by original sequence. |
| [`low_complexity`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/bags.py#L78) | Apply the archived composition and alternating-repeat exclusion rules. |
| [`metadata`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/bags.py#L42) | Record the absolute path, byte size and SHA-256 of an input or output. |
| [`norm`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/bags.py#L28) | Canonicalize residue case and collapse indistinguishable I/L spellings. |
| [`select_databases`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/bags.py#L336) | Apply FastaLake's fixed-count razor rule to original-query/protein edges. |
| [`sha`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/bags.py#L33) | Hash a complete file while bounding read-buffer memory. |
| [`tokens`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/bags.py#L165) | Separate residue and terminal tokens and check their unmodified spelling. |
| [`upstream`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/bags.py#L56) | Load the actual pinned implementation; never select another installed copy. |

### fasta_lake.benchmark

| Function | Purpose |
|---|---|
| [`_is_entrapment_accession`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/benchmark.py#L177) | One rule for both loaders: an entrapment record carries ENTRAP or SHUFFLED |
| [`_pearson`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/benchmark.py#L450) | Pearson correlation coefficient. |
| [`build_gene_map`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/benchmark.py#L183) | Build accession → gene name map from a FASTA file. |
| [`compare`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/benchmark.py#L303) | Compare FastaLake search results against a baseline. |
| [`is_real_gene_symbol`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/benchmark.py#L30) | Check if a gene name is a real HGNC/MGI symbol, not an accession fallback. |
| [`load_diann_report`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/benchmark.py#L122) | Load a DIA-NN report.tsv and extract peptides, genes, protein groups. |
| [`load_sage_report`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/benchmark.py#L211) | Load a SAGE results.sage.tsv and extract peptides, genes, protein groups. |
| [`print_benchmark`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/benchmark.py#L400) | Pretty-print benchmark results. |

### fasta_lake.canonical_comparison

| Function | Purpose |
|---|---|
| [`_load_diann_gene_peptides`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/canonical_comparison.py#L152) | Load DIA-NN report: gene → {peptides, intensity}. |
| [`classify_losses`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/canonical_comparison.py#L48) | Classify every gene lost from canonical → FastaLake. |
| [`print_loss_report`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/canonical_comparison.py#L185) | Print a detailed loss classification report. |

### fasta_lake.diagnose

| Function | Purpose |
|---|---|
| [`calibration_curve`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/diagnose.py#L97) | Compute empirical FDP vs nominal q-value curve across samples. |
| [`calibration_curve.target_db_size_fn`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/diagnose.py#L108) | Supply the historical 75,000-protein fallback when no size callback exists. |
| [`outlier_samples`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/diagnose.py#L177) | Identify samples with FDP > mean + z*std (outliers). |
| [`per_sample_fdp_at_q`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/diagnose.py#L31) | Compute per-sample PSM- and protein-level FDP at a given q threshold. |
| [`write_report`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/diagnose.py#L223) | Write calibration JSON/CSV and outlier JSON into the requested directory. |

### fasta_lake.entrapment

| Function | Purpose |
|---|---|
| [`_validate_counts`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/entrapment.py#L40) | Reject negative counts, non-integers and booleans before FDP accounting. |
| [`estimate_fdp`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/entrapment.py#L225) | Estimate False Discovery Proportion using entrapment hits. |
| [`estimate_fdp_paired`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/entrapment.py#L324) | Wen 2025 Eq 4 paired FDP estimator — a *tighter* valid upper bound. |
| [`generate_noble_paired`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/entrapment.py#L180) | Noble paired entrapment: split proteins into target and entrapment sets. |
| [`generate_shuffled_entrapment`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/entrapment.py#L47) | Generate shuffled entrapment proteins from source sequences. |
| [`load_entrapment_fasta`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/entrapment.py#L137) | Load user-supplied entrapment proteins from a FASTA file. |
| [`prepare_entrapment`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/entrapment.py#L512) | Prepare entrapment proteins based on configuration. |
| [`validate_entrapment_fdr`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/entrapment.py#L433) | Validate FDR calibration using entrapment analysis. |

### fasta_lake.entrapment_classes

| Function | Purpose |
|---|---|
| [`class_names`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/entrapment_classes.py#L154) | Registered class names, preserving the historical display order. |
| [`EntrapmentClass.as_metadata`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/entrapment_classes.py#L44) | Serialisable metadata for embedding in an FDP report. |
| [`get_class`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/entrapment_classes.py#L144) | Look up an entrapment class by name. |

### fasta_lake.entrapment_stages

| Function | Purpose |
|---|---|
| [`EntrapmentStage.as_metadata`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/entrapment_stages.py#L93) | Serialisable metadata for embedding in a manifest or FDP record. |
| [`get_stage`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/entrapment_stages.py#L202) | Look up a stage, accepting aliases. |
| [`normalize_stage`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/entrapment_stages.py#L186) | Map an accepted spelling onto a canonical stage name. |
| [`stages_comparable`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/entrapment_stages.py#L207) | True only when two stages are the same injection point. |

### fasta_lake.entrapment_workflow

| Function | Purpose |
|---|---|
| [`_config_class_name`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/entrapment_workflow.py#L1796) | Map an :class:`EntrapmentConfig` onto a registered class name. |
| [`_find_search_db`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/entrapment_workflow.py#L1232) | Return the first sorted matching search database, or None if absent. |
| [`_reconcile`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/entrapment_workflow.py#L1646) | Prefer an explicit value, but never let it silently contradict a manifest. |
| [`_sha256_file`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/entrapment_workflow.py#L334) | SHA-256 of a file, streamed. |
| [`bootstrap_median_ci`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/entrapment_workflow.py#L1104) | Seeded percentile bootstrap CI on the median. |
| [`classify_psms`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/entrapment_workflow.py#L725) | Classify accepted PSMs into target / entrapment / shared. |
| [`compute_fdp`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/entrapment_workflow.py#L1240) | Compute one FDP record for a (dataset, class, stage) triple. |
| [`compute_fdp._governing_manifest`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/entrapment_workflow.py#L1359) | Find the nearest governing injection manifest, with a single-manifest fallback. |
| [`conserved_peptide_scan`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/entrapment_workflow.py#L854) | Test entrapment-attributed peptides for occurrence in the target space. |
| [`conserved_peptide_scan.norm`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/entrapment_workflow.py#L882) | Apply I/L equivalence only when requested for this diagnostic scan. |
| [`ConservedScan.control_rate`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/entrapment_workflow.py#L849) | Fraction of reversed entrapment peptides found -- should be ~0. |
| [`ConservedScan.rate`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/entrapment_workflow.py#L844) | Fraction of entrapment peptides also present in the target space. |
| [`count_fasta_headers`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/entrapment_workflow.py#L200) | Count proteins in a FASTA by counting **header lines only**. |
| [`detect_prefixes`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/entrapment_workflow.py#L247) | Report the entrapment-like prefixes actually present in a FASTA. |
| [`expected_results_subdir`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/entrapment_workflow.py#L1165) | The per-sample subdirectory the run templates write for an entrapment class. |
| [`fdp_estimates`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/entrapment_workflow.py#L966) | Compute every FDP estimator, returning ``None`` where undefined. |
| [`find_sample_results`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/entrapment_workflow.py#L1174) | Find ``(sample_id, results_tsv)`` under ``results_dir``. |
| [`generate_entrapment_fasta`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/entrapment_workflow.py#L400) | Write an entrapment FASTA for one entrapment class. |
| [`inject_entrapment`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/entrapment_workflow.py#L537) | Concatenate entrapment into a database and record what was done. |
| [`iter_fasta`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/entrapment_workflow.py#L178) | Stream ``(header, sequence)`` pairs from a FASTA file. |
| [`load_manifests`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/entrapment_workflow.py#L642) | Load entrapment manifests, refusing to mix injection stages or classes. |
| [`partition_counts`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/entrapment_workflow.py#L219) | Return ``(n_target, n_entrapment)`` protein counts for a FASTA. |
| [`prepare_from_config`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/entrapment_workflow.py#L1748) | Honour ``entrapment.enabled`` in a FastaLake JSON config. |
| [`shuffle_sequence`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/entrapment_workflow.py#L351) | Shuffle a sequence, optionally holding K/R residues in place. |
| [`stage_measures`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/entrapment_workflow.py#L168) | One-line statement of what an FDP at this injection stage estimates. |
| [`strip_modifications`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/entrapment_workflow.py#L698) | Reduce a modified peptide string to bare uppercase residues. |
| [`summarize`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/entrapment_workflow.py#L1130) | n, median, IQR and a seeded bootstrap CI for a list of per-sample values. |
| [`verify_prefix_matches`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/entrapment_workflow.py#L301) | Assert the configured prefix finds the entrapment partition that is declared. |
| [`write_conserved_report`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/entrapment_workflow.py#L947) | Write the conserved peptides and their target proteins to TSV. |
| [`write_fdp_reports`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/entrapment_workflow.py#L1656) | Write the per-sample TSV, the (dataset, class, stage) record, and conserved peptides. |

### fasta_lake.tune

| Function | Purpose |
|---|---|
| [`apply_filter`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/tune.py#L57) | Count target / decoy / entrapment PSMs under a filter. |
| [`compute_fdp_combined`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/tune.py#L102) | Noble 2025 combined-method FDP as percentage; ``None`` when undefined. |
| [`default_filter_grid`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/tune.py#L46) | Standard grid. len × q × score cartesian product. |
| [`find_sample_results`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/tune.py#L114) | Yield (sample_id, entrap_type, results.sage.tsv path) for each sample. |
| [`parse_tsv`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/tune.py#L76) | Parse required Sage fields for the historical spectrum-q tuning diagnostic. |
| [`tune`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/tune.py#L132) | Tune parameters by post-hoc filter sweep. |
| [`tune.target_db_size_fn`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/tune.py#L168) | Supply the historical 75,000-protein fallback when no size callback exists. |
| [`write_reports`](https://github.com/MannLabs/fasta-lake/blob/main/fasta_lake/tune.py#L240) | Write tuning CSV, recommended settings and supporting reports. |
