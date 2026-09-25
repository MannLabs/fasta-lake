# Rust function index

Use [functions in execution order](function-map.md) to find the engine for a step.
This index links **170 function declarations**, sorted by file and name.
It omits test directories and conventional in-file test modules.
It is a source navigation aid, not a list of public or enabled APIs.
Repeated names belong to different files or implementations.

Regenerate with `python tools/render_function_index.py --write`.

## rust/aggregation_engine/build.rs

| Function | Definition |
|---|---|
| `main` | [Line 4](https://github.com/MannLabs/fasta-lake/blob/main/rust/aggregation_engine/build.rs#L4) |

## rust/aggregation_engine/src/main.rs

| Function | Definition |
|---|---|
| `apply_bh_correction` | [Line 924](https://github.com/MannLabs/fasta-lake/blob/main/rust/aggregation_engine/src/main.rs#L924) |
| `build_functional_group_matrix` | [Line 1484](https://github.com/MannLabs/fasta-lake/blob/main/rust/aggregation_engine/src/main.rs#L1484) |
| `build_functional_peptide_matrix` | [Line 1552](https://github.com/MannLabs/fasta-lake/blob/main/rust/aggregation_engine/src/main.rs#L1552) |
| `build_peptide_matrix` | [Line 774](https://github.com/MannLabs/fasta-lake/blob/main/rust/aggregation_engine/src/main.rs#L774) |
| `build_protein_group_matrix` | [Line 842](https://github.com/MannLabs/fasta-lake/blob/main/rust/aggregation_engine/src/main.rs#L842) |
| `canonical_for_group` | [Line 1401](https://github.com/MannLabs/fasta-lake/blob/main/rust/aggregation_engine/src/main.rs#L1401) |
| `classify_group` | [Line 331](https://github.com/MannLabs/fasta-lake/blob/main/rust/aggregation_engine/src/main.rs#L331) |
| `classify_group_with_hosts` | [Line 1147](https://github.com/MannLabs/fasta-lake/blob/main/rust/aggregation_engine/src/main.rs#L1147) |
| `find_sample_dirs` | [Line 446](https://github.com/MannLabs/fasta-lake/blob/main/rust/aggregation_engine/src/main.rs#L446) |
| `group_value_votes` | [Line 1362](https://github.com/MannLabs/fasta-lake/blob/main/rust/aggregation_engine/src/main.rs#L1362) |
| `is_classifiable_accession` | [Line 322](https://github.com/MannLabs/fasta-lake/blob/main/rust/aggregation_engine/src/main.rs#L322) |
| `is_host_member` | [Line 368](https://github.com/MannLabs/fasta-lake/blob/main/rust/aggregation_engine/src/main.rs#L368) |
| `load_functional_lookup` | [Line 1257](https://github.com/MannLabs/fasta-lake/blob/main/rust/aggregation_engine/src/main.rs#L1257) |
| `load_host_accessions` | [Line 1112](https://github.com/MannLabs/fasta-lake/blob/main/rust/aggregation_engine/src/main.rs#L1112) |
| `load_search_q` | [Line 470](https://github.com/MannLabs/fasta-lake/blob/main/rust/aggregation_engine/src/main.rs#L470) |
| `main` | [Line 1801](https://github.com/MannLabs/fasta-lake/blob/main/rust/aggregation_engine/src/main.rs#L1801) |
| `normalize_protein_group` | [Line 293](https://github.com/MannLabs/fasta-lake/blob/main/rust/aggregation_engine/src/main.rs#L293) |
| `parse_intensity` | [Line 437](https://github.com/MannLabs/fasta-lake/blob/main/rust/aggregation_engine/src/main.rs#L437) |
| `parse_lfq_sample` | [Line 539](https://github.com/MannLabs/fasta-lake/blob/main/rust/aggregation_engine/src/main.rs#L539) |
| `parse_psm_sample` | [Line 681](https://github.com/MannLabs/fasta-lake/blob/main/rust/aggregation_engine/src/main.rs#L681) |
| `parse_q` | [Line 425](https://github.com/MannLabs/fasta-lake/blob/main/rust/aggregation_engine/src/main.rs#L425) |
| `parse_q_thresholds` | [Line 157](https://github.com/MannLabs/fasta-lake/blob/main/rust/aggregation_engine/src/main.rs#L157) |
| `representative_protein` | [Line 406](https://github.com/MannLabs/fasta-lake/blob/main/rust/aggregation_engine/src/main.rs#L406) |
| `resolve_functional_levels` | [Line 1787](https://github.com/MannLabs/fasta-lake/blob/main/rust/aggregation_engine/src/main.rs#L1787) |
| `resolve_level` | [Line 1448](https://github.com/MannLabs/fasta-lake/blob/main/rust/aggregation_engine/src/main.rs#L1448) |
| `sanitize_level` | [Line 1616](https://github.com/MannLabs/fasta-lake/blob/main/rust/aggregation_engine/src/main.rs#L1616) |
| `short` | [Line 229](https://github.com/MannLabs/fasta-lake/blob/main/rust/aggregation_engine/src/main.rs#L229) |
| `split_functional_value` | [Line 1328](https://github.com/MannLabs/fasta-lake/blob/main/rust/aggregation_engine/src/main.rs#L1328) |
| `strip_source_tags` | [Line 263](https://github.com/MannLabs/fasta-lake/blob/main/rust/aggregation_engine/src/main.rs#L263) |
| `validate_search_schema` | [Line 525](https://github.com/MannLabs/fasta-lake/blob/main/rust/aggregation_engine/src/main.rs#L525) |
| `write_functional_group_matrix` | [Line 1629](https://github.com/MannLabs/fasta-lake/blob/main/rust/aggregation_engine/src/main.rs#L1629) |
| `write_functional_peptide_matrix` | [Line 1686](https://github.com/MannLabs/fasta-lake/blob/main/rust/aggregation_engine/src/main.rs#L1686) |
| `write_group_metadata` | [Line 1170](https://github.com/MannLabs/fasta-lake/blob/main/rust/aggregation_engine/src/main.rs#L1170) |
| `write_group_to_level_map` | [Line 1746](https://github.com/MannLabs/fasta-lake/blob/main/rust/aggregation_engine/src/main.rs#L1746) |
| `write_peptide_matrix` | [Line 999](https://github.com/MannLabs/fasta-lake/blob/main/rust/aggregation_engine/src/main.rs#L999) |
| `write_protein_group_matrix` | [Line 1054](https://github.com/MannLabs/fasta-lake/blob/main/rust/aggregation_engine/src/main.rs#L1054) |

## rust/common/dianovo.rs

| Function | Definition |
|---|---|
| `parse_prediction` | [Line 11](https://github.com/MannLabs/fasta-lake/blob/main/rust/common/dianovo.rs#L11) |

## rust/fasta_export/build.rs

| Function | Definition |
|---|---|
| `main` | [Line 4](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_export/build.rs#L4) |

## rust/fasta_export/src/main.rs

| Function | Definition |
|---|---|
| `compute_fdp_combined` | [Line 447](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_export/src/main.rs#L447) |
| `compute_fdp_lower` | [Line 465](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_export/src/main.rs#L465) |
| `count_entrapment_hits` | [Line 475](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_export/src/main.rs#L475) |
| `export_fasta` | [Line 380](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_export/src/main.rs#L380) |
| `extract_entry_gene` | [Line 219](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_export/src/main.rs#L219) |
| `extract_gmgc_gene` | [Line 244](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_export/src/main.rs#L244) |
| `extract_gn` | [Line 191](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_export/src/main.rs#L191) |
| `extract_gnomad_gene` | [Line 231](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_export/src/main.rs#L231) |
| `extract_os` | [Line 206](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_export/src/main.rs#L206) |
| `generate_decoys` | [Line 318](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_export/src/main.rs#L318) |
| `heal_fasta` | [Line 410](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_export/src/main.rs#L410) |
| `heal_for_diann` | [Line 270](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_export/src/main.rs#L270) |
| `heal_for_msfragger` | [Line 297](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_export/src/main.rs#L297) |
| `heal_for_sage` | [Line 260](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_export/src/main.rs#L260) |
| `heal_header` | [Line 303](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_export/src/main.rs#L303) |
| `main` | [Line 535](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_export/src/main.rs#L535) |
| `parse_header` | [Line 89](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_export/src/main.rs#L89) |

## rust/fasta_extractor_v2/build.rs

| Function | Definition |
|---|---|
| `main` | [Line 4](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/build.rs#L4) |

## rust/fasta_extractor_v2/src/bin/fasta_sort_by_id.rs

| Function | Definition |
|---|---|
| `main` | [Line 58](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/src/bin/fasta_sort_by_id.rs#L58) |

## rust/fasta_extractor_v2/src/bin/hash_fasta.rs

| Function | Definition |
|---|---|
| `emit` | [Line 55](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/src/bin/hash_fasta.rs#L55) |
| `hash_hex` | [Line 50](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/src/bin/hash_fasta.rs#L50) |
| `main` | [Line 62](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/src/bin/hash_fasta.rs#L62) |
| `normalize` | [Line 38](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/src/bin/hash_fasta.rs#L38) |

## rust/fasta_extractor_v2/src/bin/lake_builder.rs

| Function | Definition |
|---|---|
| `extract_clean_sequence` | [Line 329](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/src/bin/lake_builder.rs#L329) |
| `fmt_num` | [Line 193](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/src/bin/lake_builder.rs#L193) |
| `main` | [Line 363](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/src/bin/lake_builder.rs#L363) |
| `parse_fasta_regions` | [Line 297](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/src/bin/lake_builder.rs#L297) |
| `parse_source_spec` | [Line 239](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/src/bin/lake_builder.rs#L239) |
| `progress_bar` | [Line 341](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/src/bin/lake_builder.rs#L341) |
| `sha256_hex` | [Line 205](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/src/bin/lake_builder.rs#L205) |
| `sort_joint_headers` | [Line 219](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/src/bin/lake_builder.rs#L219) |
| `sort_joint_headers_source_blind` | [Line 234](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/src/bin/lake_builder.rs#L234) |

## rust/fasta_extractor_v2/src/bin/mass_corasick.rs

| Function | Definition |
|---|---|
| `aa_mass` | [Line 34](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/src/bin/mass_corasick.rs#L34) |
| `build` | [Line 84](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/src/bin/mass_corasick.rs#L84) |
| `load_predictions` | [Line 197](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/src/bin/mass_corasick.rs#L197) |
| `main` | [Line 378](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/src/bin/mass_corasick.rs#L378) |
| `mass_to_bin` | [Line 65](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/src/bin/mass_corasick.rs#L65) |
| `sweep_protein` | [Line 128](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/src/bin/mass_corasick.rs#L128) |

## rust/fasta_extractor_v2/src/bin/mass_corasick_v2.rs

| Function | Definition |
|---|---|
| `aa_mass` | [Line 23](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/src/bin/mass_corasick_v2.rs#L23) |
| `build` | [Line 102](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/src/bin/mass_corasick_v2.rs#L102) |
| `get_bin` | [Line 203](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/src/bin/mass_corasick_v2.rs#L203) |
| `load_predictions` | [Line 266](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/src/bin/mass_corasick_v2.rs#L266) |
| `main` | [Line 426](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/src/bin/mass_corasick_v2.rs#L426) |
| `mass_to_bin` | [Line 71](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/src/bin/mass_corasick_v2.rs#L71) |
| `sweep_protein` | [Line 213](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/src/bin/mass_corasick_v2.rs#L213) |

## rust/fasta_extractor_v2/src/bin/mass_corasick_v3.rs

| Function | Definition |
|---|---|
| `aa_mass` | [Line 28](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/src/bin/mass_corasick_v3.rs#L28) |
| `build` | [Line 106](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/src/bin/mass_corasick_v3.rs#L106) |
| `compute_kmer_bins` | [Line 71](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/src/bin/mass_corasick_v3.rs#L71) |
| `find_csv_files` | [Line 272](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/src/bin/mass_corasick_v3.rs#L272) |
| `get_bin` | [Line 165](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/src/bin/mass_corasick_v3.rs#L165) |
| `is_valid_peptide` | [Line 53](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/src/bin/mass_corasick_v3.rs#L53) |
| `load_csv_peptides` | [Line 222](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/src/bin/mass_corasick_v3.rs#L222) |
| `load_proteins` | [Line 295](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/src/bin/mass_corasick_v3.rs#L295) |
| `main` | [Line 892](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/src/bin/mass_corasick_v3.rs#L892) |
| `mass_to_bin` | [Line 48](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/src/bin/mass_corasick_v3.rs#L48) |
| `run_merge` | [Line 690](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/src/bin/mass_corasick_v3.rs#L690) |
| `run_prep` | [Line 342](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/src/bin/mass_corasick_v3.rs#L342) |
| `run_rescue` | [Line 505](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/src/bin/mass_corasick_v3.rs#L505) |
| `sweep_protein` | [Line 174](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/src/bin/mass_corasick_v3.rs#L174) |

## rust/fasta_extractor_v2/src/bin/mass_kmer_anchor.rs

| Function | Definition |
|---|---|
| `aa_mass` | [Line 32](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/src/bin/mass_kmer_anchor.rs#L32) |
| `build` | [Line 152](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/src/bin/mass_kmer_anchor.rs#L152) |
| `exact_match` | [Line 486](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/src/bin/mass_kmer_anchor.rs#L486) |
| `kmer_mass` | [Line 91](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/src/bin/mass_kmer_anchor.rs#L91) |
| `load_fasta` | [Line 292](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/src/bin/mass_kmer_anchor.rs#L292) |
| `load_predictions` | [Line 338](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/src/bin/mass_kmer_anchor.rs#L338) |
| `load_predictions_from_file` | [Line 392](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/src/bin/mass_kmer_anchor.rs#L392) |
| `main` | [Line 571](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/src/bin/mass_kmer_anchor.rs#L571) |
| `mass_to_bin` | [Line 64](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/src/bin/mass_kmer_anchor.rs#L64) |
| `new` | [Line 76](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/src/bin/mass_kmer_anchor.rs#L76) |
| `new` | [Line 117](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/src/bin/mass_kmer_anchor.rs#L117) |
| `peptide_mass` | [Line 95](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/src/bin/mass_kmer_anchor.rs#L95) |
| `query` | [Line 192](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/src/bin/mass_kmer_anchor.rs#L192) |
| `theoretical_mass` | [Line 101](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/src/bin/mass_kmer_anchor.rs#L101) |

## rust/fasta_extractor_v2/src/bin/parsimony.rs

| Function | Definition |
|---|---|
| `collect_csvs` | [Line 694](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/src/bin/parsimony.rs#L694) |
| `load_fasta` | [Line 656](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/src/bin/parsimony.rs#L656) |
| `main` | [Line 105](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/src/bin/parsimony.rs#L105) |
| `normalize` | [Line 649](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/src/bin/parsimony.rs#L649) |
| `parse_prediction_csv` | [Line 713](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/src/bin/parsimony.rs#L713) |
| `run_bayesian_em` | [Line 454](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/src/bin/parsimony.rs#L454) |
| `run_info_score` | [Line 565](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/src/bin/parsimony.rs#L565) |
| `run_razor` | [Line 601](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/src/bin/parsimony.rs#L601) |
| `walk` | [Line 696](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/src/bin/parsimony.rs#L696) |

## rust/fasta_extractor_v2/src/main.rs

| Function | Definition |
|---|---|
| `add` | [Line 374](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/src/main.rs#L374) |
| `collect_csvs` | [Line 388](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/src/main.rs#L388) |
| `extract_sequence` | [Line 294](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/src/main.rs#L294) |
| `is_low_complexity` | [Line 311](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/src/main.rs#L311) |
| `load_peptides` | [Line 494](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/src/main.rs#L494) |
| `load_peptides_single_csv` | [Line 541](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/src/main.rs#L541) |
| `load_peptides_top_percent` | [Line 413](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/src/main.rs#L413) |
| `load_peptides_top_percent_single_csv` | [Line 573](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/src/main.rs#L573) |
| `load_peptides_txt` | [Line 630](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/src/main.rs#L630) |
| `main` | [Line 673](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/src/main.rs#L673) |
| `new` | [Line 364](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/src/main.rs#L364) |
| `normalize_il_seq` | [Line 305](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/src/main.rs#L305) |
| `parse_fasta_regions` | [Line 259](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/src/main.rs#L259) |
| `parse_prediction_csv` | [Line 141](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/src/main.rs#L141) |
| `walk` | [Line 390](https://github.com/MannLabs/fasta-lake/blob/main/rust/fasta_extractor_v2/src/main.rs#L390) |

## rust/lake_stats/build.rs

| Function | Definition |
|---|---|
| `main` | [Line 4](https://github.com/MannLabs/fasta-lake/blob/main/rust/lake_stats/build.rs#L4) |

## rust/lake_stats/src/main.rs

| Function | Definition |
|---|---|
| `is_clean` | [Line 42](https://github.com/MannLabs/fasta-lake/blob/main/rust/lake_stats/src/main.rs#L42) |
| `main` | [Line 152](https://github.com/MannLabs/fasta-lake/blob/main/rust/lake_stats/src/main.rs#L152) |
| `new` | [Line 95](https://github.com/MannLabs/fasta-lake/blob/main/rust/lake_stats/src/main.rs#L95) |
| `next` | [Line 105](https://github.com/MannLabs/fasta-lake/blob/main/rust/lake_stats/src/main.rs#L105) |
| `trypsin_into` | [Line 55](https://github.com/MannLabs/fasta-lake/blob/main/rust/lake_stats/src/main.rs#L55) |

## rust/parsimony_engine/build.rs

| Function | Definition |
|---|---|
| `main` | [Line 4](https://github.com/MannLabs/fasta-lake/blob/main/rust/parsimony_engine/build.rs#L4) |

## rust/parsimony_engine/src/main.rs

| Function | Definition |
|---|---|
| `classify_by_prefix` | [Line 388](https://github.com/MannLabs/fasta-lake/blob/main/rust/parsimony_engine/src/main.rs#L388) |
| `credit_order` | [Line 960](https://github.com/MannLabs/fasta-lake/blob/main/rust/parsimony_engine/src/main.rs#L960) |
| `extract_mgyg_species` | [Line 439](https://github.com/MannLabs/fasta-lake/blob/main/rust/parsimony_engine/src/main.rs#L439) |
| `fnv1a64` | [Line 924](https://github.com/MannLabs/fasta-lake/blob/main/rust/parsimony_engine/src/main.rs#L924) |
| `fnv1a_pep_acc` | [Line 913](https://github.com/MannLabs/fasta-lake/blob/main/rust/parsimony_engine/src/main.rs#L913) |
| `get_db_type` | [Line 372](https://github.com/MannLabs/fasta-lake/blob/main/rust/parsimony_engine/src/main.rs#L372) |
| `is_low_complexity` | [Line 448](https://github.com/MannLabs/fasta-lake/blob/main/rust/parsimony_engine/src/main.rs#L448) |
| `load_fasta` | [Line 635](https://github.com/MannLabs/fasta-lake/blob/main/rust/parsimony_engine/src/main.rs#L635) |
| `load_peptides` | [Line 484](https://github.com/MannLabs/fasta-lake/blob/main/rust/parsimony_engine/src/main.rs#L484) |
| `load_source_map` | [Line 418](https://github.com/MannLabs/fasta-lake/blob/main/rust/parsimony_engine/src/main.rs#L418) |
| `load_weights` | [Line 279](https://github.com/MannLabs/fasta-lake/blob/main/rust/parsimony_engine/src/main.rs#L279) |
| `main` | [Line 1475](https://github.com/MannLabs/fasta-lake/blob/main/rust/parsimony_engine/src/main.rs#L1475) |
| `map_peptides_parallel` | [Line 692](https://github.com/MannLabs/fasta-lake/blob/main/rust/parsimony_engine/src/main.rs#L692) |
| `normalize_il` | [Line 273](https://github.com/MannLabs/fasta-lake/blob/main/rust/parsimony_engine/src/main.rs#L273) |
| `run_bayesian` | [Line 1289](https://github.com/MannLabs/fasta-lake/blob/main/rust/parsimony_engine/src/main.rs#L1289) |
| `run_info_score` | [Line 1242](https://github.com/MannLabs/fasta-lake/blob/main/rust/parsimony_engine/src/main.rs#L1242) |
| `run_razor` | [Line 791](https://github.com/MannLabs/fasta-lake/blob/main/rust/parsimony_engine/src/main.rs#L791) |
| `run_species_budget` | [Line 985](https://github.com/MannLabs/fasta-lake/blob/main/rust/parsimony_engine/src/main.rs#L985) |
| `run_uniform` | [Line 1122](https://github.com/MannLabs/fasta-lake/blob/main/rust/parsimony_engine/src/main.rs#L1122) |
| `run_uniform_2pep` | [Line 1215](https://github.com/MannLabs/fasta-lake/blob/main/rust/parsimony_engine/src/main.rs#L1215) |
| `strip_source_tags` | [Line 333](https://github.com/MannLabs/fasta-lake/blob/main/rust/parsimony_engine/src/main.rs#L333) |
| `write_fasta` | [Line 1444](https://github.com/MannLabs/fasta-lake/blob/main/rust/parsimony_engine/src/main.rs#L1444) |

## rust/shared_peptide_classifier/build.rs

| Function | Definition |
|---|---|
| `main` | [Line 4](https://github.com/MannLabs/fasta-lake/blob/main/rust/shared_peptide_classifier/build.rs#L4) |

## rust/shared_peptide_classifier/src/main.rs

| Function | Definition |
|---|---|
| `load_cluster_subset` | [Line 115](https://github.com/MannLabs/fasta-lake/blob/main/rust/shared_peptide_classifier/src/main.rs#L115) |
| `load_lake_subset` | [Line 155](https://github.com/MannLabs/fasta-lake/blob/main/rust/shared_peptide_classifier/src/main.rs#L155) |
| `load_rescued` | [Line 76](https://github.com/MannLabs/fasta-lake/blob/main/rust/shared_peptide_classifier/src/main.rs#L76) |
| `main` | [Line 203](https://github.com/MannLabs/fasta-lake/blob/main/rust/shared_peptide_classifier/src/main.rs#L203) |
