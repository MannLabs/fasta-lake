# Methods: sample-specific database construction and release validation

!!! note "Current grouping contract"

    The active utility now uses [direct peptide assignment](output-contract.md). Original dual-key descriptions and prebuilt reading copies below record the previous method; they are not the current output contract.

## Workflow and scope

FastaLake constructs an acquisition-specific protein sequence database by matching de novo peptide evidence to a pooled reference catalogue and resolving redundant sequence explanations before database searching. The workflow described here used exact sequence deduplication, I/L-insensitive peptide containment, Rust razor inference with an accession-hash tie rule, and SAGE database searching. Sequence clustering, mass-k-mer rescue and taxonomy-based inference budgets were not applied.

The release validation used ten HSTOOL technical acquisitions and ten MICROB-PREDICT acquisitions with matched metagenome FASTAs. Existing prediction CSVs and mzML files were used as inputs. Raw conversion, de novo model inference, metagenome assembly and reference-catalogue generation were not repeated. Those upstream operations therefore require their original provenance for a complete study-level Methods description. This page describes the software operations performed in the release validation and provides a parameterized basis for the final manuscript.

## Reference sequence identity and provenance

Reference protein sequences were represented using the FastaLake sequence-hashing contract. Normalization retained ASCII letters, converted lowercase letters to uppercase, and removed nonletter characters. Ambiguity codes were retained literally. SHA-256 was computed on the resulting bytes, providing a deterministic accession for each normalized sequence. Isoleucine and leucine remained distinct at this identity layer. Source-header sidecars retained the relationship between normalized sequences and their original catalogue records.

The validation reused the existing cohort reference lakes. The HSTOOL lake contained 233,497,702 sequence records and occupied 61,260,701,006 bytes. The PREDICT lake contained 272,919,212 sequence records and occupied 71,668,046,270 bytes. These counts describe the reference files used in this validation. They do not claim that every catalogue sequence represented a protein expressed in the acquisitions. The acquisition manifests and retained run provenance identify the inputs and downstream file checksums. Full-file SHA-256 digests of both reference lakes are supplied in the release source data.

## De novo evidence selection

Prediction CSVs supplied a peptide sequence and a scalar confidence score for each prediction row. Rows were eligible for construction when the cleaned sequence contained 9–50 residues and the score was positive and finite. Low-complexity sequences were excluded when one residue accounted for more than 60% of sequence length, when a sequence of at least 12 residues contained at most two distinct residue types, or when a sequence of at least 10 residues comprised a repeating dipeptide.

Eligible rows were ranked within each acquisition by decreasing score. The score at rank \(\lceil 0.30 \cdot n \rceil\), where \(n\) denotes the number of eligible rows, defined the threshold. All rows at or above that threshold were retained, including score ties. Repeated prediction rows contributed to the rank distribution before the retained peptide collection was deduplicated. Consequently the retained fraction could exceed 30% when scores tied at the cutoff. This rank filter was not interpreted as a calibrated peptide error probability.

Isoleucine was replaced by leucine for peptide-to-protein matching. The same transformation was applied to the protein sequence during containment testing, while the stored protein sequence and sequence-derived accession retained the original residues. Matching therefore admitted I/L-equivalent evidence without collapsing exact reference sequence identity.

## Evidence lake and acquisition-specific drawing

An Aho-Corasick automaton was constructed from the pooled retained peptide patterns and applied to the reference lake. A protein sequence entered the cohort evidence lake when it contained at least one pattern as an exact contiguous substring under the matching normalization. The resulting evidence lakes contained 8,785,676 HSTOOL sequences and 14,998,182 PREDICT sequences, supported by 52,993 and 114,923 distinct pooled construction peptides, respectively.

For each acquisition, its own retained peptide patterns were then matched against the cohort evidence lake. The resulting drawn database included every evidence-lake sequence containing at least one acquisition-specific pattern. No sequence clustering intervened between the pooled and acquisition-specific stages. Construction matching did not require an enzymatic peptide boundary. Enzyme specificity was imposed during the later database search.

## Pre-search razor inference

The Rust parsimony engine constructed a bipartite graph connecting distinct eligible prediction peptides to proteins in the acquisition-specific drawn database. Peptides in the 9–50-residue interval were supplied without a top-fraction argument at this stage. The inference evidence collection therefore included eligible lower-ranked predictions that were not used to draw the database.

Each peptide was assigned to one explaining protein. Proteins were preferred in descending order of the number of unique mapped peptides, followed by the total number of distinct mapped peptides. A unique peptide was defined relative to the current drawn database as a peptide mapping to exactly one protein. Residual ties were resolved by the smallest 64-bit FNV-1a hash of the accession bytes, with lexical accession order resolving a hash collision. Peptides were traversed in descending mapping degeneracy and then lexical sequence order. Under the accession-only tie rule, assignment priorities did not depend on previous assignments.

Every protein receiving at least one peptide assignment was retained in the search FASTA. This procedure preserves an explanation for every mapped inference peptide but does not guarantee a minimum-cardinality protein set. Proteins supported only by shared peptides can survive. The retained sequences were treated as search hypotheses, and no claim of protein identification was made from construction support alone.

## Database searching

Each acquisition was searched separately with SAGE 0.14.6 against its selected FASTA. The ten PREDICT acquisitions were additionally searched against their matched metagenome FASTAs using the same mzML inputs, SAGE executable and search settings. Only the supplied FASTA and corresponding output location differed between paired search configurations.

| Parameter | Value |
|---|---|
| Enzyme cleavage | C-terminal to K or R, restricted before P |
| Missed cleavages | 2 |
| Search peptide length | 5–50 residues |
| Peptide neutral mass | 500–5,000 Da |
| Fragment m/z interval | 200–2,000 |
| Fragment ion series | b and y |
| Minimum ion index | 2 |
| Precursor tolerance | −20 to +20 ppm |
| Fragment tolerance | −20 to +20 ppm |
| Isotope error range | −1 to +3 |
| Fixed modification | C, +57.0215 Da |
| Variable modification | M, +15.9949 Da |
| Maximum variable modifications | 2 |
| Minimum / maximum retained peaks | 15 / 150 |
| Minimum matched peaks | 4 |
| Maximum fragment charge | 2 |
| Reported PSMs per spectrum | 1 |
| Decoy generation | Enabled, prefix rev_ |
| Retention-time prediction | Enabled |
| Deisotoping, chimeric and wide-window modes | Disabled |
| Fragment index bucket size | 32,768 |
| Label-free quantification | Enabled |
| LFQ peak scoring / integration | Hybrid / Sum |
| LFQ spectral-angle setting | 0.7 |

The five-residue search minimum differed from the nine-residue construction minimum. This allowed the search engine to evaluate shorter enzymatically admissible peptides from proteins already included in the candidate database. The complete JSON configuration and actual FASTA path were retained for each search. SAGE was invoked with PIN-file output enabled. The implementation and statistical procedures of SAGE are described by [Lazear, Journal of Proteome Research (2023)](https://doi.org/10.1021/acs.jproteome.3c00486).

## Identification and intensity accounting

Canonical peptide identification counts used target rows with label equal to 1 and finite peptide_q values between zero and 0.01 inclusive. Bracketed modification annotations were removed and I was replaced by L before sequence deduplication. PSM counts used the spectrum-level q field when that counting unit was reported. These filters were applied at their stated level rather than treating the different confidence columns as interchangeable.

The core LFQ aggregator required the corresponding search-results table and one acquisition intensity column per search directory. Peptide-feature matrices were gated on search-level peptide_q, and protein-member-set matrices on search-level protein_q, each at 0.01. The MS1-specific LFQ q_value column did not supply these gates. Observations were filtered before their intensities were added to the relevant key. Total intensity could therefore differ between the two matrix types because the identification gates differed. Nonfinite and negative intensity values did not contribute accepted signal. A protein-member-set key comprised the sorted, deduplicated accessions reported for the feature. Member-set rows were not interpreted as uniquely resolved biological proteins.

For the canonical-peptide repeatability analysis, positive finite LFQ feature intensities were summed for target canonical peptides accepted in the corresponding acquisition. No additional filter on the MS1-integration-specific LFQ q_value field was imposed in that analysis. The identification gate and quantification definition were therefore reported separately. No missing-value imputation was applied.

## Study-level reporting groups (previous dual-key method)

This section records the previous method; the current rule is direct peptide assignment, described in [the output contract](output-contract.md).

The study reporting dictionary was restricted to the union of accessions in the exact FASTAs searched. Target peptide observations passing the stated peptide-q filter defined a pooled canonical peptide-to-protein graph. Proteins with no passing observed peptide remained in the searched universe but were not assigned an evidence-derived study group.

Study representatives were selected by greedy set cover. At each iteration, the protein explaining the largest number of uncovered peptides was selected, with lexical accession order resolving equal gains. Each peptide was assigned to the first selected representative containing it. Nonselected proteins were mapped to the representative explaining the largest number of their peptide assignments, with selected order resolving ties. The number of contributing representative groups and the fraction of peptide assignments supporting the chosen mapping were retained as ambiguity measures.

For quantitative reporting, each accepted LFQ feature retained a sample key based on the lexical first accession in its member set. Its protein members were mapped through the study dictionary, and the earliest selected available study representative supplied the study key. Distinct canonical peptide counts and LFQ feature counts were retained separately. Each feature contributed intensity once, and total included intensity was checked for conservation before and after assignment. This operation harmonized reporting under a declared cohort dictionary and did not recalculate study-wide FDR.

## Empirical acceptance design and statistics

The HSTOOL subset comprised the ten acquisitions labelled TechnicalReplicates 1–10 in the retained cohort metadata. The PREDICT subset comprised the first ten declared eligible acquisitions with matched metagenome files, corresponding to MicrobPred 10–19. Selection was specified before comparing search outcomes. Both low-yield PREDICT acquisitions remained in the analysis.

Pairwise overlap between HSTOOL acquisitions was measured with Jaccard similarity, defined as intersection size divided by union size. Exact selected sequence-accession sets and identified canonical peptide sets were evaluated separately. The 45 pairwise similarities shared acquisitions and were treated as descriptive, dependent comparisons.

For peptide intensity repeatability, the coefficient of variation was calculated as the sample standard deviation divided by the mean across ten acquisitions. This analysis used the 12,057 canonical peptides with positive accepted intensity in all ten runs. An additional sensitivity analysis divided each acquisition's intensities by the median ratio to the corresponding peptide's across-acquisition median, implemented on the log2 scale. Raw and normalized results were retained separately. Neither complete-case selection nor normalization was interpreted as recovering missing measurements or establishing biological reproducibility.

PREDICT yield was compared within acquisition as the ratio of identified canonical peptide counts from the FastaLake and matched metagenome searches. The median of the ten paired ratios summarized the comparison. No significance test was performed on this software-acceptance subset. Reference content and sample-specific selection differed simultaneously between the two arms, so their isolated effects could not be estimated from this comparison alone.

## Error calibration and diagnostic scope

The acceptance experiment did not introduce a new independent entrapment or holdout design. Search q-values were interpreted within the search procedure and supplied candidate set. Construction evidence, search evidence and post-search reporting groups were distinguished throughout.

Entrapment diagnostics require a control satisfying the absence and exchangeability assumptions at the level being counted. Empty controls were treated as undefined, and a random split of a real reference was not treated as an absent control. Construction-stage shuffled matches were reported descriptively without FDR certification. A combined entrapment estimator of the form \(N_E(1 + 1/r)/(N_T + N_E)\) requires a declared effective entrapment-to-target ratio \(r\) and appropriate assumptions. Its point estimate is not a confidence bound for an individual experiment. The methodological basis and distinctions among entrapment designs are discussed by [Wen et al., Nature Methods (2025)](https://doi.org/10.1038/s41592-025-02719-x).

## Software and reproducibility

The reviewed release candidate was designated 1.1.0rc1. Validation used Linux compute nodes, Python 3.12.4 and the Rust 1.91.1 toolchain. Rust dependencies were resolved from retained lockfiles. Input manifests, executable identities, commands, logs and comparison tables were preserved in separate output directories. Failed initial attempts and resource-related retries were retained alongside successful runs.

The original and hardened aggregation engines were applied to the same thirty completed search outputs. At the 0.01 LFQ threshold, all six resulting peptide-feature and protein-member-set matrix comparisons were numerically identical. Candidate extraction and razor selection were also compared against the original twenty FastaLake acquisitions, with exact FASTA comparisons retained. These checks establish the observed behavior for the tested data and supported configuration. They do not prove correctness for every possible input.

## Bundled reproducibility example

A small executable example was prepared from the CAMPI SIHUMIx mock-community acquisitions S01 and S02 deposited under PRIDE PXD023217, associated with [Van Den Bossche et al., Nature Communications 12, 7305 (2021)](https://doi.org/10.1038/s41467-021-27542-8). The retention-time interval 45 ≤ RT < 50 minutes was specified before inspecting the new search outcomes. All MS1 and MS2 spectra within the interval were retained, together with any referenced precursor MS1 spectrum. Measured binary-array payloads and native scan identifiers were preserved, and XML indices were rebuilt. S01 contained 381 MS1 and 5,689 MS2 spectra; S02 contained 375 MS1 and 5,615 MS2 spectra. Each retained MS2 spectrum was paired to one previously generated AlphaNovo prediction by native scan number, with precursor m/z agreement within 0.059 ppm. The prediction provenance and source-file checksums accompany the example.

All 29,673 target entries from the deposited SIHUMI_DB1UNIPROT.faa were retained, while its 29,673 DECOY_ entries were excluded before SAGE decoy generation. Exact normalization and deduplication produced 29,635 reference sequences. Construction, inference, search and grouping used the settings described above, with the construction rank cutoff recomputed within each cropped acquisition. The reference execution produced 538 pooled evidence sequences, selected databases of 348 and 287 sequences, 548 and 426 canonical target peptides at peptide_q ≤ 0.01, and 448 post-search study groups. Database identities were checked by exact SHA-256 comparison. Search counts, peptide-set concordance and conservation of included intensity were checked by the accompanying verification script against its declared reference and tolerances.

This example provides a distributable execution check in addition to the separate 30-search acceptance experiment. Cropping changes the evidence ranking, calibration population and chromatographic boundaries. Its results were therefore not interpreted as full-acquisition yields, full-run quantitative measurements or independent error-rate validation. Neural-network training, raw conversion and model inference were not repeated in this exercise.

## Information required for the final study manuscript

The final manuscript must link these implementation methods to the original prediction-generation provenance, including the actual AlphaNovo model and checkpoint, preprocessing, score definition and training-data relationship to the evaluation spectra. It must also supply reference-catalogue releases, upstream assembly methods, confirmed biological pairing and replicate definitions, and deposited data identifiers. These details cannot be reconstructed reliably from a prediction filename or inferred from a current source checkout. The release guide deliberately does not invent them.
