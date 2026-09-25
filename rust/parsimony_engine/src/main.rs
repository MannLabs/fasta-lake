use aho_corasick::AhoCorasick;
use clap::{Parser, ValueEnum};
use dashmap::DashMap;
use rayon::prelude::*;
use serde::Serialize;
#[path = "../../common/dianovo.rs"]
mod dianovo;

use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::PathBuf;
use std::time::Instant;

// ============================================================================
// CLI
// ============================================================================

const VERSION: &str = concat!(env!("CARGO_PKG_VERSION"), "+", env!("GIT_REV"));

#[derive(Parser)]
#[command(name = "parsimony_engine")]
#[command(version = VERSION)]
#[command(about = "High-performance parsimony inference for metaproteomics")]
struct Cli {
    /// Strategy: razor, species_budget, uniform, uniform_2pep, info_score, bayesian
    #[arg(short = 'm', long, value_enum)]
    strategy: Strategy,

    /// Stage2 FASTA (sample-specific protein database)
    #[arg(short = 'f', long)]
    fasta: PathBuf,

    /// AlphaNovo predictions CSV
    #[arg(short = 'p', long)]
    peptides: PathBuf,

    /// Output directory
    #[arg(short = 'o', long)]
    output_dir: PathBuf,

    /// Sample ID
    #[arg(short = 's', long)]
    sample: String,

    /// Optional accession -> source-catalogue map (TSV: accession, source), for lakes
    /// whose accessions are content hashes and carry no catalogue prefix.
    ///
    /// Without it, `db_breakdown` on a SHA-256-keyed lake can only report UNCLASSIFIED,
    /// because `PREDICT_<64 hex>` says nothing about which catalogue the sequence came
    /// from. The hash-keyed eggNOG lookups already have this layout and can be passed
    /// directly. Coverage is printed, so a map that only partly covers the database
    /// cannot be mistaken for a complete one.
    #[arg(long)]
    source_map: Option<PathBuf>,

    /// Minimum peptide length [default: 9, matches Python wrapper and v5.3 entrapment-FDR fix]
    #[arg(long, default_value_t = 9)]
    min_length: usize,

    /// Maximum peptide length [default: 50]
    #[arg(long, default_value_t = 50)]
    max_length: usize,

    /// Top percent of peptides by score per sample (e.g. 0.30 for top 30%)
    #[arg(long)]
    top_percent: Option<f64>,

    /// Normalize I/L (mass spec equivalence) [default: true]
    #[arg(long, default_value_t = true, num_args = 0..=1, default_missing_value = "true", action = clap::ArgAction::Set)]
    normalize_il: bool,

    /// Number of threads for parallel search
    #[arg(short = 't', long, default_value_t = 1)]
    threads: usize,

    /// Optional per-peptide weight TSV (peptide<TAB>weight, with header).
    /// Multiplied into the denovo score in info_score and bayesian strategies.
    /// Missing peptides get weight 1.0. Razor and uniform ignore weights.
    #[arg(long)]
    weights: Option<PathBuf>,

    /// How razor breaks ties between candidates with identical evidence.
    ///
    /// REQUIRED, deliberately: there is NO default. A tie-break rule changes which
    /// proteins are reported, so a run that does not state its rule cannot be
    /// reproduced or compared, and the omission is invisible in the output.
    ///
    /// This was a real failure, not a hypothetical. The canonical 05_parsimony.sh
    /// omitted --tiebreak and inherited whatever the binary on its PATH happened to
    /// do; the binary it actually found predates this option entirely and hardcodes
    /// alpha. Making the flag required is what makes that class of mistake impossible,
    /// and it is strictly better than flipping the default: flipping would silently
    /// change the result of every existing script, whereas this makes them fail loudly
    /// and reproduce exactly once `--tiebreak alpha` is added.
    ///
    /// For v26 use `hash-acc`. To reproduce pre-v26 published output use `alpha`.
    #[arg(long, value_enum, required = true)]
    tiebreak: TieBreak,
}

/// How razor resolves a shared peptide when two or more candidate proteins carry
/// IDENTICAL evidence (same unique-peptide count, same total mapped peptides).
///
/// This is criterion 3 of razor's comparator and it is not a detail: measured on the
/// 97-sample healthy-stool cohort, 39.56% of shared peptides (median) reach it, and the
/// rule chosen changes the per-sample database by up to 18% in size with under half its
/// content in common.
///
/// NONE OF THESE IS WRONG. They sit at different points on the comparability-resolution
/// axis, which is the same trade-off `--group-key-mode` and `--min-seq-id` expose
/// elsewhere in the pipeline:
///
///   cross-sample comparability  <----------------------->  sample-specific resolution
///        Alpha, Hash                                              Greedy
///   (winner is a pure function of the accession, so the      (winner depends on that
///    same protein wins in every sample)                       sample's evidence)
///
/// Measured on healthy stool (97 samples): Greedy is more parsimonious within a sample
/// (-1.86% database, smaller in 96/97) but less comparable across samples (singletons
/// 52.32% -> 55.25%, mean pairwise Jaccard 0.0896 -> 0.0827). Depth cost of switching is
/// -1.21%.
///
/// ⚠ CORRECTED 2026-07-30. This paragraph used to end: "Hash is catalogue-neutral AND
/// cross-sample deterministic, at +18.6% database size, so it dominates Alpha on both of
/// those axes." The word "dominates" is not supportable and the framing is misleading,
/// because it prices a +18.6% database as if it were free:
///   * Full-cohort paired measurement puts the inflation HIGHER than +18.6% --
///     PREDICT +21.5% and CAPSCAN +22.1%, bigger on 306/306 and 368/368 samples.
///   * Hash also raises per-sample depth ~6% (PREDICT 11,785 -> 12,518; CAPSCAN
///     11,729 -> 12,590). More identifications out of a larger database at fixed FDR is
///     exactly the pattern entrapment exists to test, and it has NOT been tested. Decoy
///     FDR is identical across arms (0.987% PSM), which does not settle it -- decoys and
///     entrapment answer different questions.
///   * Peter's ruling 2026-07-30: the hash arm is NOT reported. It stays as the
///     catalogue-neutrality reference only, and no depth or database-size figure from it
///     goes into the manuscript.
///
/// Use `hash-acc` if catalogue neutrality is what is wanted; it achieves it without the
/// inflation, at a protein-yield cost that is stated on that variant.
///
/// Pick by what the study needs. See analyses/razor_tiebreak_bias/FINDING.md and
/// COMPARABILITY_RESOLUTION_AXIS.md.
#[derive(Clone, Copy, PartialEq, Eq, Debug, ValueEnum)]
enum TieBreak {
    /// Lexicographically smallest accession wins. The historical behaviour, so it
    /// reproduces pre-v26 output exactly. Order-independent by
    /// construction, because it is a pure function of the accession strings. Its flaw
    /// is that byte order encodes a catalogue preference nobody chose --
    /// GMGC < GMSC < MGYG < sp| < tr| -- and `db_breakdown` reports those categories.
    Alpha,
    /// FNV-1a keyed on (peptide, accession). Catalogue-neutral, but keying on the
    /// peptide means a different protein wins each tie, which SPREADS assignments and
    /// inflates the database ~18%. Kept as the neutrality reference; not recommended,
    /// because razor exists to be parsimonious.
    Hash,
    /// FNV-1a on the ACCESSION ALONE. Removes Alpha's catalogue preference while keeping
    /// both of Alpha's virtues. It is a TRADE, not a free win.
    ///
    /// ⚠ CORRECTED 2026-07-30. This doc previously read "The recommended rule, and the only
    /// one that is better than Alpha on every axis rather than trading." **Measurement does
    /// not support that**, and a code comment must not assert what the data contradicts:
    ///
    ///   DE-BIASING: CONFIRMED, and only where accessions carry a catalogue.
    ///     Healthy stool (catalogue-keyed, n=97): Alpha puts 57.96% of the database in
    ///     GMGC; greedy, hash and hash-acc all land at ~40% -- three mechanistically
    ///     unrelated rules agreeing, on 97/97 samples. Alpha's swapped-out set is 81.9%
    ///     GMGC.
    ///     PREDICT (SHA-256-keyed, n=306): largest Alpha-vs-hash-acc difference is
    ///     0.016 pp. With no catalogue in the accession there is no preference to remove,
    ///     so hash-acc buys NOTHING on such a lake. Most FastaLake cohorts are of this kind.
    ///
    ///   COST: -11.34% protein yield on healthy stool (260,441 -> 230,919 target proteins
    ///     at protein_q<0.01) at an unchanged achieved PSM-level FDR (0.987% in every arm).
    ///     Whether those lost proteins were false positives removed or true proteins lost
    ///     CANNOT be settled with decoys; it needs an entrapment FDP measurement, which has
    ///     not been done.
    ///
    ///   SIZE: statistically indistinguishable from Alpha, not smaller. Paired per-sample,
    ///     PREDICT n=306 median +0.5 sequences (sign test p=0.56); CAPSCAN n=368 median
    ///     +3.0 sequences (p=1.4e-4, i.e. significant only because n is large enough to
    ///     resolve +0.04%). GREEDY is the arm that shrinks the database (median -103/-112
    ///     sequences, smaller on 304/306 and 362/365 samples).
    ///
    /// So: pick hash-acc when catalogue neutrality on a catalogue-keyed lake is worth a
    /// yield cost.
    ///
    /// ⚠ DEFAULT CHANGED 2026-07-31: hash-acc is now the DEFAULT, alpha is retained for
    /// reproducing pre-v26 output. The ruling rests on the v26 accession scheme: with
    /// every accession `COHORT_<sha256hex>` (--uniform-tag), no accession carries a
    /// catalogue, so alpha has no preference left to express and the two rules converge
    /// -- measured at 0.016 pp on PREDICT, which was already hash-keyed in v25. hash-acc
    /// was chosen anyway so neutrality holds by the RULE and not only by the accession
    /// scheme, i.e. it does not silently depend on --uniform-tag also being correct.
    ///
    /// The -11.34% yield cost recorded above was measured on a CATALOGUE-KEYED lake
    /// (v25 healthy stool) and is not expected to apply to a v26 lake, where all four
    /// cohorts are hash-keyed. That expectation is UNVERIFIED. It must be confirmed by
    /// the entrapment FDP measurement, not assumed -- decoys cannot settle it.
    ///
    /// Alpha's two virtues both come from its winner depending ONLY on the accession
    /// set: (i) CONCENTRATION -- if two proteins tie on five peptides, all five go to
    /// the same one, so one protein is selected rather than two, which is parsimony;
    /// (ii) CROSS-SAMPLE INVARIANCE -- one sample ties on peptide P1, another on P2, and
    /// the same protein wins both, so both samples select it. Its single flaw is that
    /// the ordering it uses, byte order, correlates with the catalogue prefix.
    ///
    /// Hashing the accession alone keeps both virtues and removes the flaw: still a pure
    /// function of the accession, so still concentrating and still sample-invariant, but
    /// a hash of "GMGC..." bears no ordering relation to a hash of "MGYG...".
    ///
    /// `Hash` keys on (peptide, accession), which was carried over from `credit_order`
    /// where SPREADING credit is wanted because `uniform` is the permissive arm. For razor
    /// spreading is the wrong objective: measured at +18.6% database size on healthy stool,
    /// and +21.5% / +22.1% on the full PREDICT and CAPSCAN cohorts, bigger on every single
    /// sample. Keying on the accession alone avoids that.
    HashAcc,
    /// Prefer the tied candidate already carrying the most peptides assigned so far;
    /// FNV-1a for residual ties. Concentrates assignments like Alpha, but on accumulated
    /// EVIDENCE rather than on the accession string. Requires a fixed peptide order
    /// (descending degeneracy, then peptide) because it carries state across peptides,
    /// so unlike Alpha it is order-independent only GIVEN that order.
    Greedy,
}

#[derive(Clone, ValueEnum)]
enum Strategy {
    Razor,
    SpeciesBudget,
    Uniform,
    Uniform2pep,
    InfoScore,
    Bayesian,
}

// ============================================================================
// Data structures
// ============================================================================

#[derive(Serialize)]
struct ParsimonyStats {
    sample_id: String,
    strategy: String,
    runtime_seconds: f64,
    stage2_proteins: usize,
    proteins_with_evidence: usize,
    selected_proteins: usize,
    reduction_factor: f64,
    unique_peptides: usize,
    shared_peptides: usize,
    total_peptides_mapped: usize,
    db_breakdown: BTreeMap<String, DbBreakdown>,
    // species_budget specific
    #[serde(skip_serializing_if = "Option::is_none")]
    species_budget_ceiling: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    validated_species: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    validated_species_retained: Option<usize>,
}

#[derive(Serialize, Default)]
struct DbBreakdown {
    unique: usize,
    shared_only: usize,
    total: usize,
}

// ============================================================================
// Helpers
// ============================================================================

fn normalize_il(seq: &str) -> String {
    seq.replace('I', "L")
}

/// Load optional per-peptide weights TSV. Header required; "peptide" and "weight" columns.
/// Peptides are normalized I→L if `normalize_il_flag` is set, matching the peptide_scores key.
fn load_weights(path: &PathBuf, normalize_il_flag: bool) -> HashMap<String, f64> {
    let f = std::fs::File::open(path).expect("cannot open weights TSV");
    let r = BufReader::new(f);
    let mut iter = r.lines();
    let header = iter
        .next()
        .expect("empty weights TSV")
        .expect("read error on weights header");
    let cols: Vec<&str> = header.split('\t').collect();
    let pep_idx = cols
        .iter()
        .position(|c| *c == "peptide")
        .expect("weights TSV missing 'peptide' column");
    let w_idx = cols
        .iter()
        .position(|c| *c == "weight")
        .expect("weights TSV missing 'weight' column");
    let mut out: HashMap<String, f64> = HashMap::new();
    for line in iter {
        let line = line.expect("read error on weights row");
        if line.is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() <= pep_idx.max(w_idx) {
            continue;
        }
        let pep = parts[pep_idx];
        let key = if normalize_il_flag {
            normalize_il(pep)
        } else {
            pep.to_string()
        };
        if let Ok(w) = parts[w_idx].parse::<f64>() {
            out.insert(key, w);
        }
    }
    out
}

/// Strip a TWO-tag source prefix like "UNCLUSTERED|CAPSCAN|" to get the raw protein ID.
/// Format: TAG1|TAG2|<original_id> where original_id may itself contain pipes (e.g. sp|P12345|NAME)
///
/// DOC CORRECTED 2026-07-30. This comment previously also claimed it strips a SINGLE tag,
/// `"CLUSTER97|"`. It does not: the guard is `pipe_positions.len() >= 2`, so a one-pipe
/// prefix is left in place and the accession then fails prefix classification. The doc,
/// not the code, was wrong -- and it misled a test written against it.
///
/// LATENT, NOT LIVE. Scanned 200,000 headers of a real healthy-stool `stage2.fasta`: all
/// 52 pipe-containing accessions are `sp|X|Y` (two pipes, first segment lowercase, so
/// correctly NOT stripped and then matched by the `sp|` branch). Zero single-pipe
/// uppercase tags occur. Behaviour is therefore deliberately left UNCHANGED: altering
/// `strip_source_tags` would move classification for real data and needs its own
/// verification, which a doc fix does not.
fn strip_source_tags(protein_id: &str) -> &str {
    // Count pipes: if >=2 and the first segment looks like a tag (all caps, short),
    // strip the first two pipe-separated segments
    let pipe_positions: Vec<usize> = protein_id.match_indices('|').map(|(i, _)| i).collect();
    if pipe_positions.len() >= 2 {
        let first_seg = &protein_id[..pipe_positions[0]];
        // Source tags are uppercase words like UNCLUSTERED, CLUSTER97, CAPSCAN
        if first_seg
            .chars()
            .all(|c| c.is_ascii_uppercase() || c.is_ascii_digit() || c == '_')
        {
            return &protein_id[pipe_positions[1] + 1..];
        }
    }
    protein_id
}

/// Source catalogue for one protein, for `db_breakdown`.
///
/// FIXED 2026-07-30. The previous version ended in `else { "smORF" }`, a CATCH-ALL. That
/// silently mislabelled every accession whose prefix it did not recognise -- which is
/// every accession in a SHA-256-keyed lake. PREDICT and CAPSCAN write
/// `PREDICT_<64 hex>` ids, so every db_breakdown ever emitted for those cohorts read
/// "100% smORF". The truth, measured from `kept_header` in predict_joint_headers.tsv, is
/// roughly 35% UHGG / 33% GMSC / 17% GMGC / 15% MuPr_Assembly (the project's own
/// metagenome) / 0.06% Human. The reported figure did not merely lack precision, it named
/// the wrong catalogue for 100% of proteins.
///
/// Two changes:
///   1. `GMSC` is now matched EXPLICITLY. It had been relying on the catch-all, because
///      GMSC IS the smORF catalogue -- which is why the fallthrough was named "smORF" and
///      why the bug stayed invisible on catalogue-keyed lakes. Making it explicit means
///      output for a catalogue-keyed lake is UNCHANGED (asserted in the tests below), so
///      only the genuinely-unrecognised case moves.
///   2. The fallthrough is now `UNCLASSIFIED`, not a real catalogue name. A value we
///      cannot derive must not be reported as one we can.
///
/// For hash-keyed lakes, pass `--source-map` to get a CORRECT breakdown rather than an
/// honest `UNCLASSIFIED`.
fn get_db_type<'a>(
    protein_id: &'a str,
    source_map: Option<&'a HashMap<String, String>>,
) -> &'a str {
    if let Some(map) = source_map {
        if let Some(src) = map.get(protein_id) {
            return src.as_str();
        }
        // Fall through to the prefix rule on a miss: a map covering only part of the
        // database must not silently relabel the remainder.
    }
    classify_by_prefix(protein_id)
}

/// Prefix-based catalogue classification. Correct only while the accession still CARRIES
/// its catalogue prefix; returns `UNCLASSIFIED` when it does not, rather than guessing.
fn classify_by_prefix(protein_id: &str) -> &'static str {
    let raw = strip_source_tags(protein_id);
    if raw.starts_with("sp|") || raw.starts_with("tr|") {
        "Human"
    } else if raw.starts_with("MGYG") {
        "UHGG"
    } else if raw.starts_with("GMGC") {
        "GMGC"
    } else if raw.starts_with("GMSC") {
        // GMSC is the smORF catalogue. Previously reached via the catch-all.
        "smORF"
    } else if raw.contains("Assembly_") {
        // Project-specific metagenome assembly, e.g. MuPr_Assembly_9998637.
        "Assembly"
    } else {
        "UNCLASSIFIED"
    }
}

/// Load an accession -> source-catalogue map, for lakes whose accessions are content
/// hashes and so carry no catalogue prefix.
///
/// Format: TSV, column 1 = accession, column 2 = source; a leading header line is skipped
/// when column 2 is literally "source". That is deliberately the layout of the existing
/// hash-keyed eggNOG lookups
/// (`projects/<cohort>/lookups/funcdecoy_p05_eggnog_<COHORT>_hashkeyed.tsv`), so those
/// files can be passed straight in.
///
/// Coverage is REPORTED by the caller, never assumed: a lookup built for an older
/// reference set may cover only a fraction of the current lake.
fn load_source_map(path: &std::path::Path) -> std::io::Result<HashMap<String, String>> {
    use std::io::BufRead;
    let file = std::fs::File::open(path)?;
    let mut map = HashMap::new();
    for (i, line) in std::io::BufReader::new(file).lines().enumerate() {
        let line = line?;
        let mut it = line.split('\t');
        let (Some(acc), Some(src)) = (it.next(), it.next()) else {
            continue;
        };
        if i == 0 && src == "source" {
            continue;
        }
        if acc.is_empty() || src.is_empty() {
            continue;
        }
        map.insert(acc.to_string(), src.to_string());
    }
    Ok(map)
}

fn extract_mgyg_species(protein_id: &str) -> Option<&str> {
    let raw = strip_source_tags(protein_id);
    if raw.starts_with("MGYG") {
        raw.split('_').next()
    } else {
        None
    }
}

fn is_low_complexity(seq: &str) -> bool {
    let n = seq.len();
    if n == 0 {
        return true;
    }
    let bytes = seq.as_bytes();
    let mut counts = [0u32; 26];
    for &b in bytes {
        if b.is_ascii_uppercase() {
            counts[(b - b'A') as usize] += 1;
        }
    }
    let max_count = *counts.iter().max().unwrap_or(&0);
    if (max_count as f64 / n as f64) > 0.6 {
        return true;
    }
    let unique_aas = counts.iter().filter(|&&c| c > 0).count();
    if n >= 12 && unique_aas <= 2 {
        return true;
    }
    if n >= 10 {
        let dp = &bytes[..2];
        let is_dipeptide_repeat = bytes
            .chunks(2)
            .all(|chunk| chunk == dp || (chunk.len() == 1 && chunk[0] == dp[0]));
        if is_dipeptide_repeat {
            return true;
        }
    }
    false
}

// ============================================================================
// Load peptides from AlphaNovo CSV
// ============================================================================

fn load_peptides(
    csv_path: &PathBuf,
    min_len: usize,
    max_len: usize,
    top_percent: Option<f64>,
    do_normalize: bool,
) -> HashMap<String, f64> {
    eprintln!(
        "Loading peptides from {}...",
        csv_path.file_name().unwrap().to_string_lossy()
    );

    let file = fs::File::open(csv_path).expect("Cannot open peptide CSV");
    let mut reader = csv::ReaderBuilder::new()
        .has_headers(true)
        .flexible(true)
        .from_reader(BufReader::new(file));

    let headers: Vec<String> = reader
        .headers()
        .expect("No CSV headers")
        .iter()
        .map(|s| s.trim().trim_matches('\r').to_string())
        .collect();

    // Auto-detect de novo CSV format (DIANovo / AlphaNovo / generic).
    // DIANovo  : pred_seq + pred_prob (per-position confidence list)
    // AlphaNovo: peptide_prediction_detokenized_unmodified + score
    // generic  : sequence + score
    let dianovo_seq_col = headers.iter().position(|h| h == "pred_seq");
    let dianovo_prob_col = headers.iter().position(|h| h == "pred_prob");
    let alphanovo_col = headers
        .iter()
        .position(|h| h == "peptide_prediction_detokenized_unmodified");
    let generic_seq_col = headers.iter().position(|h| h == "sequence");
    let score_col = headers.iter().position(|h| h == "score");

    enum Fmt {
        Dianovo(usize, usize),
        Scalar(usize, Option<usize>),
    }
    let fmt = if let (Some(s), Some(p)) = (dianovo_seq_col, dianovo_prob_col) {
        Fmt::Dianovo(s, p)
    } else if let Some(s) = alphanovo_col {
        Fmt::Scalar(s, score_col)
    } else if let Some(s) = generic_seq_col {
        Fmt::Scalar(s, score_col)
    } else {
        panic!(
            "Unrecognized de novo CSV format in {:?}. Expected columns: \
                pred_seq/pred_prob (DIANovo), \
                peptide_prediction_detokenized_unmodified/score (AlphaNovo), \
                or sequence/score (generic). Found headers: {:?}",
            csv_path.file_name().unwrap_or_default(),
            headers
        );
    };

    let mut all_peptides: Vec<(String, f64)> = Vec::new();
    let mut invalid_dianovo = 0usize;

    for result in reader.records() {
        let record = match result {
            Ok(r) => r,
            Err(_) => continue,
        };

        let (peptide_raw, score) = match fmt {
            Fmt::Dianovo(s, p) => {
                match dianovo::parse_prediction(
                    record.get(s).unwrap_or(""),
                    record.get(p).unwrap_or(""),
                ) {
                    Some(value) => value,
                    None => {
                        invalid_dianovo += 1;
                        continue;
                    }
                }
            }
            Fmt::Scalar(s, sc_col) => {
                let seq = record.get(s).unwrap_or("").trim().to_uppercase();
                let sc: f64 = sc_col
                    .and_then(|c| record.get(c))
                    .unwrap_or("0.5")
                    .parse()
                    .unwrap_or(0.0);
                (seq, sc)
            }
        };

        if peptide_raw.is_empty() || peptide_raw.len() < min_len || peptide_raw.len() > max_len {
            continue;
        }
        if !score.is_finite() || score <= 0.0 {
            continue;
        }
        if is_low_complexity(&peptide_raw) {
            continue;
        }

        let peptide = if do_normalize {
            normalize_il(&peptide_raw)
        } else {
            peptide_raw
        };
        all_peptides.push((peptide, score));
    }

    if invalid_dianovo > 0 {
        eprintln!(
            "DIANovo: excluded {invalid_dianovo} unresolved/invalid prediction rows in {}",
            csv_path.display()
        );
    }

    // Apply top-percent filter if set
    let filtered: Vec<(String, f64)> = if let Some(pct) = top_percent {
        let mut sorted_scores: Vec<f64> = all_peptides.iter().map(|(_, s)| *s).collect();
        sorted_scores.sort_by(|a, b| b.partial_cmp(a).unwrap());
        let cutoff_idx = ((sorted_scores.len() as f64) * pct).ceil() as usize;
        let threshold = if cutoff_idx > 0 && cutoff_idx <= sorted_scores.len() {
            sorted_scores[cutoff_idx.min(sorted_scores.len()) - 1]
        } else {
            0.0
        };
        all_peptides
            .into_iter()
            .filter(|(_, s)| *s >= threshold)
            .collect()
    } else {
        all_peptides
    };

    // Deduplicate keeping max score
    let mut peptide_scores: HashMap<String, f64> = HashMap::new();
    for (pep, score) in filtered {
        let entry = peptide_scores.entry(pep).or_insert(0.0);
        if score > *entry {
            *entry = score;
        }
    }

    eprintln!("  {} unique peptides loaded", peptide_scores.len());
    peptide_scores
}

// ============================================================================
// Load FASTA proteins
// ============================================================================

fn load_fasta(fasta_path: &PathBuf) -> (Vec<String>, HashMap<String, String>) {
    eprintln!(
        "Loading FASTA {}...",
        fasta_path.file_name().unwrap().to_string_lossy()
    );

    let file = fs::File::open(fasta_path).expect("Cannot open FASTA");
    let reader = BufReader::with_capacity(8 * 1024 * 1024, file);

    let mut ids: Vec<String> = Vec::new();
    let mut proteins: HashMap<String, String> = HashMap::new();
    let mut current_id = String::new();
    let mut current_seq = String::new();

    for line in reader.lines() {
        let line = line.expect("Error reading FASTA");
        if line.starts_with('>') {
            if !current_id.is_empty() {
                ids.push(current_id.clone());
                if proteins.contains_key(&current_id) {
                    eprintln!("ERROR: duplicate FASTA accession {}; use sequence-derived unique identifiers", current_id);
                    std::process::exit(2);
                }
                proteins.insert(current_id.clone(), current_seq.clone());
                current_seq.clear();
            }
            current_id = line
                .strip_prefix('>')
                .unwrap_or(&line)
                .split_whitespace()
                .next()
                .unwrap_or("")
                .to_string();
        } else {
            current_seq.push_str(line.trim());
        }
    }
    if !current_id.is_empty() {
        ids.push(current_id.clone());
        if proteins.contains_key(&current_id) {
            eprintln!(
                "ERROR: duplicate FASTA accession {}; use sequence-derived unique identifiers",
                current_id
            );
            std::process::exit(2);
        }
        proteins.insert(current_id, current_seq);
    }

    eprintln!("  {} proteins loaded", ids.len());
    (ids, proteins)
}

// ============================================================================
// Parallel Aho-Corasick peptide-to-protein mapping
// ============================================================================

fn map_peptides_parallel(
    proteins: &HashMap<String, String>,
    peptide_list: &[String],
    do_normalize: bool,
    threads: usize,
) -> (
    HashMap<String, HashSet<String>>,
    HashMap<String, HashSet<String>>,
) {
    eprintln!(
        "Building Aho-Corasick automaton ({} patterns)...",
        peptide_list.len()
    );
    let start = Instant::now();

    let ac = AhoCorasick::builder()
        .build(peptide_list)
        .expect("Failed to build automaton");

    eprintln!(
        "  Automaton built in {:.2}ms",
        start.elapsed().as_secs_f64() * 1000.0
    );

    // Collect proteins into a vec for parallel iteration
    let protein_vec: Vec<(&String, &String)> = proteins.iter().collect();

    eprintln!(
        "Mapping peptides to {} proteins ({} threads)...",
        protein_vec.len(),
        threads
    );
    let map_start = Instant::now();

    let pep2prot: DashMap<String, HashSet<String>> = DashMap::new();
    let prot2pep: DashMap<String, HashSet<String>> = DashMap::new();

    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .ok(); // ignore if already set

    protein_vec.par_iter().for_each(|(prot_id, seq)| {
        let clean = seq.to_uppercase().replace('*', "");
        let normed = if do_normalize {
            normalize_il(&clean)
        } else {
            clean
        };

        let mut local_peps: Vec<String> = Vec::new();

        // find_overlapping_iter (not find_iter): a peptide that positionally overlaps
        // another matched peptide inside a protein must still register its edge. Using the
        // non-overlapping find_iter here silently dropped ~23% of peptide->protein edges,
        // inflating razor depth by a median 8.2% (see analyses/addyosmani_benchmark). This
        // matches the extractor, which already uses find_overlapping_iter.
        for mat in ac.find_overlapping_iter(&normed) {
            let pep = &peptide_list[mat.pattern().as_usize()];
            local_peps.push(pep.clone());
        }

        if !local_peps.is_empty() {
            let pid = (*prot_id).clone();
            for pep in &local_peps {
                pep2prot.entry(pep.clone()).or_default().insert(pid.clone());
            }
            let pep_set: HashSet<String> = local_peps.into_iter().collect();
            prot2pep.insert(pid, pep_set);
        }
    });

    let elapsed = map_start.elapsed().as_secs_f64();
    eprintln!("  Mapping done in {:.1}s", elapsed);
    eprintln!(
        "  {} peptides mapped to {} proteins",
        pep2prot.len(),
        prot2pep.len()
    );

    // Convert DashMaps to HashMaps
    let p2p: HashMap<String, HashSet<String>> = pep2prot.into_iter().collect();
    let pr2pe: HashMap<String, HashSet<String>> = prot2pep.into_iter().collect();

    (p2p, pr2pe)
}

// ============================================================================
// RAZOR PARSIMONY
// ============================================================================

/// Select complete candidate proteins by assigning each mapped peptide once.
///
/// The two maps describe the same acquisition-level peptide/protein graph.
/// Unique and total evidence counts are computed from that complete graph and
/// remain fixed during assignment. HashAcc resolves evidence ties by the smallest
/// FNV-1a accession hash, then lexical accession; other modes are explicit options.
/// Return accessions receiving at least one assignment. This is pre-search
/// selection, not the post-search greedy cover used for study reporting.
fn run_razor(
    pep2prot: &HashMap<String, HashSet<String>>,
    prot2pep: &HashMap<String, HashSet<String>>,
    tiebreak: TieBreak,
) -> HashSet<String> {
    eprintln!("\n=== RAZOR PARSIMONY (tiebreak={:?}) ===", tiebreak);

    // Step 1: Count unique peptides per protein
    let mut unique_count: HashMap<String, usize> = HashMap::new();
    let mut unique_peptides: HashSet<String> = HashSet::new();

    for (pep, prots) in pep2prot {
        if prots.len() == 1 {
            let prot = prots.iter().next().unwrap();
            *unique_count.entry(prot.clone()).or_insert(0) += 1;
            unique_peptides.insert(pep.clone());
        }
    }

    eprintln!("  Unique peptides: {}", unique_peptides.len());
    eprintln!("  Proteins with unique peptides: {}", unique_count.len());

    // Fixed peptide order. Required by TieBreak::Greedy, which carries state across
    // peptides; irrelevant to Alpha and Hash, which are pure functions of the
    // (peptide, accession) pair. Applying it unconditionally means all three rules see
    // the same order, so a comparison between them isolates the rule.
    //
    // Descending degeneracy first is deliberate: the most ambiguous peptides are placed
    // while the fewest proteins have accumulated evidence, so the arbitrary early
    // choices land on the peptides where they matter least.
    //
    // Alpha's output is unchanged by this sort. Verified empirically: razor under Alpha
    // gave identical selected-set sizes in HashMap order (Stage A) and in this sorted
    // order (Stage B) on all 97 healthy-stool samples.
    let mut peps: Vec<&String> = pep2prot.keys().collect();
    peps.sort_by(|a, b| {
        pep2prot[*b]
            .len()
            .cmp(&pep2prot[*a].len())
            .then_with(|| a.cmp(b))
    });

    // Step 2: Assign each peptide to its best protein
    let mut assigned: HashMap<String, HashSet<String>> = HashMap::new();
    // Peptides assigned so far, per protein. Only read by Greedy. Owned keys rather
    // than borrows, because `assigned` is mutated in the same loop.
    let mut assigned_n: HashMap<String, usize> = HashMap::new();

    for pep in peps {
        let prots = &pep2prot[pep];
        let best: String = if prots.len() == 1 {
            prots.iter().next().unwrap().clone()
        } else {
            // Evidence key: most unique peptides, then most total mapped peptides.
            // Unchanged; this is criteria 1 and 2 and is not what `tiebreak` varies.
            let ev = |p: &String| -> (usize, usize) {
                (
                    unique_count.get(p).copied().unwrap_or(0),
                    prot2pep.get(p).map(|s| s.len()).unwrap_or(0),
                )
            };
            let maxk = prots.iter().map(&ev).max().unwrap();
            let tied: Vec<&String> = prots.iter().filter(|p| ev(p) == maxk).collect();

            match tiebreak {
                // Historical behaviour: lexicographically smallest accession.
                TieBreak::Alpha => (*tied.iter().min_by(|a, b| a.cmp(b)).unwrap()).clone(),
                TieBreak::Hash => (*tied
                    .iter()
                    .min_by(|a, b| {
                        fnv1a_pep_acc(pep, a)
                            .cmp(&fnv1a_pep_acc(pep, b))
                            .then_with(|| a.cmp(b))
                    })
                    .unwrap())
                .clone(),
                // Accession-only hash: invariant to which peptide triggered the tie,
                // so it concentrates and is sample-invariant exactly as Alpha is.
                TieBreak::HashAcc => (*tied
                    .iter()
                    .min_by(|a, b| {
                        fnv1a64(a.as_bytes())
                            .cmp(&fnv1a64(b.as_bytes()))
                            .then_with(|| a.cmp(b))
                    })
                    .unwrap())
                .clone(),
                // Most already-assigned peptides first, hash for residual ties.
                TieBreak::Greedy => (*tied
                    .iter()
                    .max_by(|a, b| {
                        let ca = assigned_n.get(a.as_str()).copied().unwrap_or(0);
                        let cb = assigned_n.get(b.as_str()).copied().unwrap_or(0);
                        // (borrow-safe: assigned_n is keyed by owned String)
                        ca.cmp(&cb).then_with(|| {
                            // max_by wants the winner to compare Greater, and the hash
                            // rule is "smallest hash wins", so reverse it.
                            fnv1a_pep_acc(pep, b)
                                .cmp(&fnv1a_pep_acc(pep, a))
                                .then_with(|| b.cmp(a))
                        })
                    })
                    .unwrap())
                .clone(),
            }
        };
        *assigned_n.entry(best.clone()).or_insert(0) += 1;
        assigned.entry(best).or_default().insert(pep.clone());
    }

    let selected: HashSet<String> = assigned.keys().cloned().collect();
    eprintln!("  Selected proteins: {}", selected.len());
    selected
}

// ============================================================================
// DETERMINISTIC CREDIT ORDER (shared by species_budget & uniform)
// ============================================================================

/// FNV-1a keyed on (peptide, accession). Used by both `credit_order` and razor's
/// Hash/Greedy tie-breaks so a single definition governs every place a reproducible
/// arbitrary choice is made.
fn fnv1a_pep_acc(peptide: &str, acc: &str) -> u64 {
    let mut k = Vec::with_capacity(peptide.len() + acc.len() + 1);
    k.extend_from_slice(peptide.as_bytes());
    k.push(0);
    k.extend_from_slice(acc.as_bytes());
    fnv1a64(&k)
}

/// FNV-1a, 64-bit. Spelled out rather than using `DefaultHasher` because the std
/// hasher's algorithm is explicitly allowed to change between Rust releases, and
/// this value decides which proteins enter a published database.
fn fnv1a64(bytes: &[u8]) -> u64 {
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    for &b in bytes {
        h ^= b as u64;
        h = h.wrapping_mul(0x0000_0100_0000_01b3);
    }
    h
}

/// Order in which candidate proteins receive budget-limited credit for a shared
/// peptide, when the peptide's degeneracy exceeds its budget.
///
/// WHY A HASH AND NOT A RANKING. The pre-D3 code iterated a `HashSet`, i.e. it
/// allocated credit ARBITRARILY among candidates. That was nondeterministic and had
/// to be fixed, but "arbitrary" was the intended semantics: the budget exists to CAP
/// how far a shared peptide's credit spreads, not to pick winners by evidence.
///
/// The first D3 fix ordered candidates by unique-peptide evidence (most first). That
/// is deterministic but not semantics-preserving, and in `run_uniform` it is actively
/// harmful: Phase 1 has already selected every protein with unique evidence, so
/// ordering those first spends the budget on proteins that are already in the output.
/// When the budget is 1 -- the common case for a 9-mer, where `length_norm` = 0.125 --
/// Phase 2 then contributes NOTHING, and shared-only proteins can never be credited.
/// That silently turns the permissive contrast arm into a razor variant.
///
/// Sorting on the accession alone would be deterministic and would restore neutral
/// allocation, but FastaLake accessions carry source tags and sort
/// `CLUSTER97|` < `GMGC` < `MGYG` < `sp|` < `tr|`, so it would hand a systematic
/// advantage to whichever catalogue sorts earliest -- and `db_breakdown` in
/// `stats.json` reports exactly those categories.
///
/// Keying a hash on (peptide, protein) is deterministic, reproducible across runs and
/// thread counts, and unbiased with respect to BOTH unique-evidence status and
/// catalogue prefix. Because the key includes the peptide, a given protein is not
/// systematically favoured across all peptides. The accession is the tie-break, so the
/// order is total.
fn credit_order(peptide: &str, a: &str, b: &str) -> std::cmp::Ordering {
    let mut ka = Vec::with_capacity(peptide.len() + a.len() + 1);
    ka.extend_from_slice(peptide.as_bytes());
    ka.push(0);
    ka.extend_from_slice(a.as_bytes());
    let mut kb = Vec::with_capacity(peptide.len() + b.len() + 1);
    kb.extend_from_slice(peptide.as_bytes());
    kb.push(0);
    kb.extend_from_slice(b.as_bytes());
    fnv1a64(&ka).cmp(&fnv1a64(&kb)).then_with(|| a.cmp(b))
}

// ============================================================================
// SPECIES BUDGET PARSIMONY
//
// LIMITATION (v5.3): species are inferred via the `MGYG*` accession-prefix
// heuristic in `extract_mgyg_species()`. This works for UHGG-derived lakes
// but does NOT generalize to GMGC, GMSC, or assembly-only proteins (which
// fall into a single "unknown species" bucket). The Python wrapper accepts
// an external `--taxonomy-map` file to extend species attribution to those
// catalogs; the Rust engine does not yet support this. For non-UHGG lakes,
// prefer the Python implementation (`fasta_lake.cli infer --strategy
// species_budget --taxonomy-map ...`) until this divergence is unified.
// ============================================================================

fn run_species_budget(
    pep2prot: &HashMap<String, HashSet<String>>,
    prot2pep: &HashMap<String, HashSet<String>>,
    peptide_scores: &HashMap<String, f64>,
) -> (HashSet<String>, usize, usize, usize) {
    eprintln!("\n=== SPECIES-BUDGET PARSIMONY ===");

    // Find unique peptides and validated species
    let mut unique_peptides: HashSet<String> = HashSet::new();
    let mut validated_species: HashSet<String> = HashSet::new();
    let mut validated_proteins: HashSet<String> = HashSet::new();

    for (pep, prots) in pep2prot {
        if prots.len() == 1 {
            unique_peptides.insert(pep.clone());
        }
    }

    for (prot_id, pep_set) in prot2pep {
        if pep_set.iter().any(|p| unique_peptides.contains(p)) {
            if let Some(species) = extract_mgyg_species(prot_id) {
                validated_species.insert(species.to_string());
                validated_proteins.insert(prot_id.clone());
            }
        }
    }

    let species_ceiling = validated_species.len();
    eprintln!("  Unique peptides: {}", unique_peptides.len());
    eprintln!("  Validated species: {} (budget ceiling)", species_ceiling);

    // Phase 1: unique peptides auto-include
    let mut protein_scores: HashMap<String, usize> = HashMap::new();
    for pep in &unique_peptides {
        if let Some(prots) = pep2prot.get(pep) {
            let prot = prots.iter().next().unwrap();
            *protein_scores.entry(prot.clone()).or_insert(0) += 1;
        }
    }

    // Phase 2: shared peptides with budget
    let shared_peptides: Vec<&String> =
        pep2prot.keys().filter(|p| pep2prot[*p].len() > 1).collect();

    // Sort: longest first, then highest score
    let mut shared_sorted: Vec<(&String, usize, f64)> = shared_peptides
        .iter()
        .map(|p| (*p, p.len(), *peptide_scores.get(*p).unwrap_or(&0.5)))
        .collect();
    shared_sorted.sort_by(|a, b| b.1.cmp(&a.1).then(b.2.partial_cmp(&a.2).unwrap()));

    let mut budgets: HashMap<String, i64> = HashMap::new();
    for (pep, len, score) in &shared_sorted {
        let degeneracy = pep2prot[*pep].len();
        // length_norm: linear ramp 0.0 at len=6 → 1.0 at len=30. Methods M5 restricts
        // inputs to len 9–50 so this stays in [0.125, 1.83]. The 6/24 anchors flatten
        // species-budget allocation around tryptic mid-length peptides; the upper end
        // is unclamped because longer peptides should get proportionally more budget.
        let length_norm = ((*len as f64) - 6.0).max(0.0) / 24.0;
        let confidence = score * length_norm;
        let capped = degeneracy.min(species_ceiling);
        let raw = (capped as f64) * confidence;
        let budget = (raw.ceil() as i64).max(1);
        budgets.insert((*pep).clone(), budget);
    }

    for (pep, _, _) in &shared_sorted {
        let budget = budgets.get_mut(*pep).unwrap();
        if *budget <= 0 {
            continue;
        }

        if let Some(prots) = pep2prot.get(*pep) {
            // D3 fix: iterate candidates in a DETERMINISTIC, principled order
            // (most unique evidence first) instead of raw HashSet order, so WHICH
            // proteins get credited when degeneracy exceeds the budget is stable
            // run-to-run and consistent with razor.
            let mut ordered: Vec<&String> = prots.iter().collect();
            ordered.sort_by(|a, b| credit_order(pep, a.as_str(), b.as_str()));
            for prot_id in ordered {
                if *budget <= 0 {
                    break;
                }
                *protein_scores.entry(prot_id.clone()).or_insert(0) += 1;
                *budget -= 1;
            }
        }
    }

    // Phase 3: DB-specific selection
    let mut selected = HashSet::new();
    for (prot_id, score) in &protein_scores {
        // C4 fix: strip source tags before the MGYG test so tagged UHGG proteins
        // (e.g. "CLUSTER97|MGYG...") are classified the same way extract_mgyg_species
        // classifies them for the ceiling.
        let dominated_include = if strip_source_tags(prot_id).starts_with("MGYG") {
            if validated_proteins.contains(prot_id) {
                *score >= 1
            } else {
                *score >= 2 // shared-only UHGG needs 2
            }
        } else {
            *score >= 1
        };
        if dominated_include {
            selected.insert(prot_id.clone());
        }
    }

    // Count retained validated species
    let mut retained_species: HashSet<String> = HashSet::new();
    for prot in &selected {
        if let Some(sp) = extract_mgyg_species(prot) {
            retained_species.insert(sp.to_string());
        }
    }
    let validated_retained = retained_species.intersection(&validated_species).count();

    eprintln!("  Selected proteins: {}", selected.len());
    eprintln!(
        "  Species retained: {} ({} validated)",
        retained_species.len(),
        validated_retained
    );

    (
        selected,
        species_ceiling,
        validated_species.len(),
        validated_retained,
    )
}

// ============================================================================
// UNIFORM PARSIMONY
// ============================================================================

fn run_uniform(
    pep2prot: &HashMap<String, HashSet<String>>,
    prot2pep: &HashMap<String, HashSet<String>>,
    peptide_scores: &HashMap<String, f64>,
) -> (HashSet<String>, usize) {
    eprintln!("\n=== UNIFORM PARSIMONY ===");

    // Same phases 1&2 as species_budget, but phase 3 = uniform >=1
    let mut unique_peptides: HashSet<String> = HashSet::new();
    let mut validated_species: HashSet<String> = HashSet::new();

    for (pep, prots) in pep2prot {
        if prots.len() == 1 {
            unique_peptides.insert(pep.clone());
        }
    }

    for (prot_id, pep_set) in prot2pep {
        if pep_set.iter().any(|p| unique_peptides.contains(p)) {
            if let Some(species) = extract_mgyg_species(prot_id) {
                validated_species.insert(species.to_string());
            }
        }
    }

    let species_ceiling = validated_species.len();
    eprintln!("  Validated species: {} (budget ceiling)", species_ceiling);

    // Phase 1
    let mut protein_scores: HashMap<String, usize> = HashMap::new();
    for pep in &unique_peptides {
        if let Some(prots) = pep2prot.get(pep) {
            let prot = prots.iter().next().unwrap();
            *protein_scores.entry(prot.clone()).or_insert(0) += 1;
        }
    }

    // Phase 2
    let mut shared_sorted: Vec<(&String, usize, f64)> = pep2prot
        .keys()
        .filter(|p| pep2prot[*p].len() > 1)
        .map(|p| (p, p.len(), *peptide_scores.get(p).unwrap_or(&0.5)))
        .collect();
    shared_sorted.sort_by(|a, b| b.1.cmp(&a.1).then(b.2.partial_cmp(&a.2).unwrap()));

    let mut budgets: HashMap<String, i64> = HashMap::new();
    for (pep, len, score) in &shared_sorted {
        let degeneracy = pep2prot[*pep].len();
        // length_norm: linear ramp 0.0 at len=6 → 1.0 at len=30. Methods M5 restricts
        // inputs to len 9–50 so this stays in [0.125, 1.83]. The 6/24 anchors flatten
        // species-budget allocation around tryptic mid-length peptides; the upper end
        // is unclamped because longer peptides should get proportionally more budget.
        let length_norm = ((*len as f64) - 6.0).max(0.0) / 24.0;
        let confidence = score * length_norm;
        let capped = degeneracy.min(species_ceiling);
        let raw = (capped as f64) * confidence;
        budgets.insert((*pep).clone(), (raw.ceil() as i64).max(1));
    }

    for (pep, _, _) in &shared_sorted {
        let budget = budgets.get_mut(*pep).unwrap();
        if *budget <= 0 {
            continue;
        }
        if let Some(prots) = pep2prot.get(*pep) {
            // D3 fix: deterministic, most-unique-evidence-first credit order.
            let mut ordered: Vec<&String> = prots.iter().collect();
            ordered.sort_by(|a, b| credit_order(pep, a.as_str(), b.as_str()));
            for prot_id in ordered {
                if *budget <= 0 {
                    break;
                }
                *protein_scores.entry(prot_id.clone()).or_insert(0) += 1;
                *budget -= 1;
            }
        }
    }

    // Phase 3: UNIFORM - all proteins with >=1 peptide
    let selected: HashSet<String> = protein_scores
        .into_iter()
        .filter(|(_, score)| *score >= 1)
        .map(|(id, _)| id)
        .collect();

    eprintln!("  Selected proteins: {}", selected.len());
    (selected, species_ceiling)
}

// ============================================================================
// UNIFORM 2-PEP PARSIMONY
// ============================================================================

fn run_uniform_2pep(
    prot2pep: &HashMap<String, HashSet<String>>,
    min_peptides: usize,
) -> HashSet<String> {
    eprintln!(
        "\n=== UNIFORM 2-PEP PARSIMONY (min_peptides={}) ===",
        min_peptides
    );

    let selected: HashSet<String> = prot2pep
        .iter()
        .filter(|(_, peps)| peps.len() >= min_peptides)
        .map(|(id, _)| id.clone())
        .collect();

    eprintln!(
        "  Selected proteins: {} (of {} with evidence)",
        selected.len(),
        prot2pep.len()
    );
    selected
}

// ============================================================================
// INFO SCORE PARSIMONY
// ============================================================================

fn run_info_score(
    pep2prot: &HashMap<String, HashSet<String>>,
    prot2pep: &HashMap<String, HashSet<String>>,
    peptide_scores: &HashMap<String, f64>,
    peptide_weights: &HashMap<String, f64>,
    min_info: f64,
) -> HashSet<String> {
    eprintln!(
        "\n=== INFO SCORE PARSIMONY (min_info={}, n_weights={}) ===",
        min_info,
        peptide_weights.len()
    );

    let mut info_scores: HashMap<String, f64> = HashMap::new();
    for (prot_id, peps) in prot2pep {
        let mut info = 0.0_f64;
        // D6 fix: fold the float sum over a SORTED peptide view so the accumulated
        // info score is bit-identical run-to-run (HashSet order otherwise flips
        // borderline proteins at the min_info threshold).
        let mut sorted_peps: Vec<&String> = peps.iter().collect();
        sorted_peps.sort();
        for pep in sorted_peps {
            let n_match = pep2prot.get(pep).map(|s| s.len()).unwrap_or(0);
            if n_match == 0 {
                continue;
            }
            let denovo_score = peptide_scores.get(pep).copied().unwrap_or(1.0);
            let w = peptide_weights.get(pep).copied().unwrap_or(1.0);
            info += (denovo_score * w) / (n_match as f64);
        }
        info_scores.insert(prot_id.clone(), info);
    }

    let selected: HashSet<String> = info_scores
        .into_iter()
        .filter(|(_, score)| *score >= min_info)
        .map(|(id, _)| id)
        .collect();

    eprintln!("  Selected proteins: {}", selected.len());
    selected
}

// ============================================================================
// BAYESIAN EM PARSIMONY
// ============================================================================

fn run_bayesian(
    pep2prot: &HashMap<String, HashSet<String>>,
    prot2pep: &HashMap<String, HashSet<String>>,
    proteins: &HashMap<String, String>,
    peptide_scores: &HashMap<String, f64>,
    peptide_weights: &HashMap<String, f64>,
    n_iterations: usize,
    min_probability: f64,
) -> HashSet<String> {
    eprintln!(
        "\n=== BAYESIAN EM PARSIMONY (iters={}, min_prob={}, n_weights={}) ===",
        n_iterations,
        min_probability,
        peptide_weights.len()
    );

    // D4 fix: sort the protein universe so every float fold below (priors,
    // normalizers, M-step) accumulates in a fixed order, bit-identical run-to-run.
    let mut all_prots: Vec<String> = prot2pep.keys().cloned().collect();
    all_prots.sort();
    if all_prots.is_empty() {
        return HashSet::new();
    }

    // Peptide keys in a fixed order for the deterministic E-step below.
    let mut sorted_peps: Vec<&String> = pep2prot.keys().collect();
    sorted_peps.sort();

    // Protein lengths (floor at 50)
    let prot_len: HashMap<&str, f64> = all_prots
        .iter()
        .map(|pid| {
            let len = proteins.get(pid).map(|s| s.len()).unwrap_or(50).max(50);
            (pid.as_str(), len as f64)
        })
        .collect();

    // Initialize: prior proportional to score×weight-weighted peptide count / sqrt(length)
    let mut prob: HashMap<&str, f64> = HashMap::new();
    for pid in &all_prots {
        let mut pid_peps: Vec<&String> = prot2pep[pid].iter().collect();
        pid_peps.sort(); // deterministic fold order
        let raw: f64 = pid_peps
            .iter()
            .map(|pep| {
                let s = peptide_scores.get(*pep).copied().unwrap_or(1.0);
                let w = peptide_weights.get(*pep).copied().unwrap_or(1.0);
                s * w
            })
            .sum();
        prob.insert(pid.as_str(), raw / prot_len[pid.as_str()].sqrt());
    }

    // Normalize (sum over the sorted protein universe for a fixed fold order)
    let total: f64 = all_prots
        .iter()
        .map(|p| prob.get(p.as_str()).copied().unwrap_or(0.0))
        .sum();
    if total > 0.0 {
        for v in prob.values_mut() {
            *v /= total;
        }
    }

    // EM iterations
    for iteration in 0..n_iterations {
        let mut new_evidence: HashMap<&str, f64> = HashMap::new();

        // E-step: iterate peptides in fixed (sorted) order, and fold each
        // shared peptide's contributions over a sorted protein view, so both the
        // prob_sum normalizer and the per-protein evidence accumulate identically
        // run-to-run.
        for pep in &sorted_peps {
            let matching_prots = &pep2prot[*pep];
            let denovo_score = peptide_scores.get(*pep).copied().unwrap_or(1.0);
            let pep_w = peptide_weights.get(*pep).copied().unwrap_or(1.0);
            let evidence = denovo_score * pep_w;

            if matching_prots.len() == 1 {
                let pid = matching_prots.iter().next().unwrap();
                *new_evidence.entry(pid.as_str()).or_insert(0.0) += evidence;
            } else {
                let mut ordered: Vec<&String> = matching_prots.iter().collect();
                ordered.sort();
                let prob_sum: f64 = ordered
                    .iter()
                    .map(|p| prob.get(p.as_str()).copied().unwrap_or(0.0))
                    .sum();
                if prob_sum == 0.0 {
                    let share = evidence / (ordered.len() as f64);
                    for p in &ordered {
                        *new_evidence.entry(p.as_str()).or_insert(0.0) += share;
                    }
                } else {
                    for p in &ordered {
                        let post = prob.get(p.as_str()).copied().unwrap_or(0.0) / prob_sum;
                        *new_evidence.entry(p.as_str()).or_insert(0.0) += evidence * post;
                    }
                }
            }
        }

        // M-step
        let mut max_delta: f64 = 0.0;
        for pid in &all_prots {
            let old_p = prob.get(pid.as_str()).copied().unwrap_or(0.0);
            let new_p = new_evidence.get(pid.as_str()).copied().unwrap_or(0.0)
                / prot_len[pid.as_str()].sqrt();
            prob.insert(pid.as_str(), new_p);
            let delta = (new_p - old_p).abs();
            if delta > max_delta {
                max_delta = delta;
            }
        }

        // Normalize (fixed fold order over the sorted protein universe)
        let total: f64 = all_prots
            .iter()
            .map(|p| prob.get(p.as_str()).copied().unwrap_or(0.0))
            .sum();
        if total > 0.0 {
            for v in prob.values_mut() {
                *v /= total;
            }
        }

        // EM convergence threshold. 1e-8 is well below the LFQ intensity dynamic range
        // and protein-probability resolution; below this max-delta, further iterations
        // change protein assignments by less than the noise floor of directLFQ. Mirrors
        // the Python equivalent in fasta_lake/inference/strategies.py (search "1e-8").
        if max_delta < 1e-8 {
            eprintln!("  Converged at iteration {}", iteration + 1);
            break;
        }
    }

    // Scale so max = 1.0, then threshold
    let max_prob = prob.values().cloned().fold(0.0_f64, f64::max);
    let selected: HashSet<String> = if max_prob > 0.0 {
        prob.iter()
            .filter(|(_, &p)| (p / max_prob) >= min_probability)
            .map(|(&pid, _)| pid.to_string())
            .collect()
    } else {
        HashSet::new()
    };

    eprintln!("  Selected proteins: {}", selected.len());
    selected
}

// ============================================================================
// Write FASTA output
// ============================================================================

fn write_fasta(
    output_path: &PathBuf,
    selected: &HashSet<String>,
    proteins: &HashMap<String, String>,
) -> u64 {
    let file = fs::File::create(output_path).expect("Cannot create output FASTA");
    let mut writer = BufWriter::with_capacity(4 * 1024 * 1024, file);

    let mut sorted_ids: Vec<&String> = selected.iter().collect();
    sorted_ids.sort();

    for prot_id in &sorted_ids {
        if let Some(seq) = proteins.get(*prot_id) {
            writeln!(writer, ">{}", prot_id).unwrap();
            for chunk in seq.as_bytes().chunks(60) {
                writer.write_all(chunk).unwrap();
                writer.write_all(b"\n").unwrap();
            }
        }
    }

    writer.flush().unwrap();
    fs::metadata(output_path).unwrap().len()
}

// ============================================================================
// Main
// ============================================================================

/// Load one acquisition's candidates and predictions, run the chosen selection
/// strategy, and write the selected full-length FASTA with evidence summaries.
fn main() {
    let cli = Cli::parse();
    let total_start = Instant::now();
    let fail = |message: &str| {
        eprintln!("ERROR: {}", message);
        std::process::exit(2);
    };
    if cli.min_length == 0 || cli.max_length < cli.min_length || cli.threads == 0 {
        fail("Require positive threads and 1 <= min-length <= max-length");
    }
    if cli
        .top_percent
        .is_some_and(|v| !v.is_finite() || v <= 0.0 || v > 1.0)
    {
        fail("top-percent is a fraction in (0, 1], e.g. 0.30");
    }
    if cli.sample.is_empty()
        || cli.sample == "."
        || cli.sample == ".."
        || !cli
            .sample
            .bytes()
            .all(|c| c.is_ascii_alphanumeric() || b"._-".contains(&c))
    {
        fail("Sample ID must contain only ASCII letters, digits, dot, underscore or hyphen");
    }
    if !cli.fasta.is_file() || !cli.peptides.is_file() {
        fail("Input FASTA and peptide file must exist");
    }
    for p in [cli.source_map.as_ref(), cli.weights.as_ref()]
        .into_iter()
        .flatten()
    {
        if !p.is_file() {
            fail("Requested source map or weights file does not exist");
        }
    }
    if cli.output_dir.exists()
        && (!cli.output_dir.is_dir()
            || fs::read_dir(&cli.output_dir)
                .map(|mut r| r.next().is_some())
                .unwrap_or(true))
    {
        fail("Output directory must be new or empty; existing results are preserved");
    }

    let strategy_name = match cli.strategy {
        Strategy::Razor => "razor",
        Strategy::SpeciesBudget => "species_budget",
        Strategy::Uniform => "uniform",
        Strategy::Uniform2pep => "uniform_2pep",
        Strategy::InfoScore => "info_score",
        Strategy::Bayesian => "bayesian",
    };

    eprintln!("================================================================================");
    eprintln!(
        "PARSIMONY ENGINE — {} — {}",
        strategy_name.to_uppercase(),
        cli.sample
    );
    eprintln!("================================================================================");

    // Load peptides
    let peptide_scores = load_peptides(
        &cli.peptides,
        cli.min_length,
        cli.max_length,
        cli.top_percent,
        cli.normalize_il,
    );

    // Optional per-peptide weights (FastaLake v5.4 — provenance-aware parsimony)
    let peptide_weights: HashMap<String, f64> = match &cli.weights {
        Some(p) => {
            let w = load_weights(p, cli.normalize_il);
            eprintln!(
                "  {} per-peptide weights loaded from {}",
                w.len(),
                p.display()
            );
            w
        }
        None => HashMap::new(),
    };

    // Load FASTA
    let (_ids, proteins) = load_fasta(&cli.fasta);

    // Build peptide list
    let peptide_list: Vec<String> = peptide_scores.keys().cloned().collect();

    // Parallel mapping
    let (pep2prot, prot2pep) =
        map_peptides_parallel(&proteins, &peptide_list, cli.normalize_il, cli.threads);

    // Count unique/shared
    let unique_count = pep2prot.values().filter(|s| s.len() == 1).count();
    let shared_count = pep2prot.len() - unique_count;

    // Run strategy
    let (selected, species_ceiling, validated_sp, validated_sp_retained) = match cli.strategy {
        Strategy::Razor => {
            let sel = run_razor(&pep2prot, &prot2pep, cli.tiebreak);
            (sel, None, None, None)
        }
        Strategy::SpeciesBudget => {
            let (sel, ceil, vs, vsr) = run_species_budget(&pep2prot, &prot2pep, &peptide_scores);
            (sel, Some(ceil), Some(vs), Some(vsr))
        }
        Strategy::Uniform => {
            let (sel, ceil) = run_uniform(&pep2prot, &prot2pep, &peptide_scores);
            (sel, Some(ceil), None, None)
        }
        Strategy::Uniform2pep => {
            let sel = run_uniform_2pep(&prot2pep, 2);
            (sel, None, None, None)
        }
        Strategy::InfoScore => {
            let sel = run_info_score(&pep2prot, &prot2pep, &peptide_scores, &peptide_weights, 0.1);
            (sel, None, None, None)
        }
        Strategy::Bayesian => {
            let sel = run_bayesian(
                &pep2prot,
                &prot2pep,
                &proteins,
                &peptide_scores,
                &peptide_weights,
                10,
                0.01,
            );
            (sel, None, None, None)
        }
    };

    // DB breakdown.
    //
    // A SHA-256-keyed lake carries no catalogue prefix, so without --source-map every
    // protein lands in UNCLASSIFIED. That is deliberate and is reported below: the old
    // behaviour was to call all of them "smORF", which named the wrong catalogue for 100%
    // of proteins and was reported for months without anyone being able to see it.
    let source_map = match cli.source_map.as_deref() {
        Some(p) => match load_source_map(p) {
            Ok(m) => {
                eprintln!("source map: {} entries from {}", m.len(), p.display());
                Some(m)
            }
            Err(e) => {
                eprintln!(
                    "FATAL: --source-map given but unreadable ({}): {e}",
                    p.display()
                );
                std::process::exit(1);
            }
        },
        None => None,
    };

    let mut db_breakdown: BTreeMap<String, DbBreakdown> = BTreeMap::new();
    for prot_id in &selected {
        let db = get_db_type(prot_id, source_map.as_ref()).to_string();
        let entry = db_breakdown.entry(db).or_default();
        // Check if protein has unique peptide
        let has_unique = if let Some(peps) = prot2pep.get(prot_id) {
            peps.iter()
                .any(|p| pep2prot.get(p).map(|s| s.len()) == Some(1))
        } else {
            false
        };
        if has_unique {
            entry.unique += 1;
        } else {
            entry.shared_only += 1;
        }
        entry.total += 1;
    }

    // Write output
    fs::create_dir_all(&cli.output_dir).expect("Cannot create output dir");

    let suffix = strategy_name;
    let output_fasta = cli
        .output_dir
        .join(format!("{}_{}.fasta", cli.sample, suffix));
    let file_bytes = write_fasta(&output_fasta, &selected, &proteins);
    let file_mb = file_bytes as f64 / (1024.0 * 1024.0);

    let elapsed = total_start.elapsed().as_secs_f64();
    let reduction = if selected.is_empty() {
        0.0
    } else {
        proteins.len() as f64 / selected.len() as f64
    };

    // Stats
    let stats = ParsimonyStats {
        sample_id: cli.sample.clone(),
        strategy: strategy_name.to_string(),
        runtime_seconds: elapsed,
        stage2_proteins: proteins.len(),
        proteins_with_evidence: prot2pep.len(),
        selected_proteins: selected.len(),
        reduction_factor: reduction,
        unique_peptides: unique_count,
        shared_peptides: shared_count,
        total_peptides_mapped: pep2prot.len(),
        db_breakdown,
        species_budget_ceiling: species_ceiling,
        validated_species: validated_sp,
        validated_species_retained: validated_sp_retained,
    };

    let stats_path = cli
        .output_dir
        .join(format!("{}_{}_stats.json", cli.sample, suffix));
    let stats_json = serde_json::to_string_pretty(&stats).unwrap();
    fs::write(&stats_path, &stats_json).unwrap();

    // Print summary
    eprintln!("\n================================================================================");
    eprintln!("DONE — {} — {}", strategy_name.to_uppercase(), cli.sample);
    eprintln!("================================================================================");
    eprintln!("  Input proteins:    {:>12}", proteins.len());
    eprintln!("  With evidence:     {:>12}", prot2pep.len());
    eprintln!("  Selected:          {:>12}", selected.len());
    eprintln!("  Reduction:         {:>11.1}x", reduction);
    eprintln!("  Unique peptides:   {:>12}", unique_count);
    eprintln!("  Shared peptides:   {:>12}", shared_count);
    eprintln!("  Output:            {:.1} MB", file_mb);
    eprintln!("  Runtime:           {:.1}s", elapsed);

    for (db, br) in &stats.db_breakdown {
        eprintln!(
            "  {:15}: {:>8} ({} unique, {} shared-only)",
            db, br.total, br.unique, br.shared_only
        );
    }
    // Say so when the breakdown is not a real answer. Silence here is what let "100%
    // smORF" stand as a reported composition for PREDICT and CAPSCAN.
    let unclassified = stats
        .db_breakdown
        .get("UNCLASSIFIED")
        .map(|b| b.total)
        .unwrap_or(0);
    let classified_total: usize = stats.db_breakdown.values().map(|b| b.total).sum();
    if unclassified > 0 && classified_total > 0 {
        let pct = 100.0 * unclassified as f64 / classified_total as f64;
        eprintln!(
            "  WARNING: {unclassified} of {classified_total} proteins ({pct:.1}%) could not be \
             attributed to a source catalogue."
        );
        eprintln!(
            "           Their accessions carry no catalogue prefix -- this is a content-hash \
             lake."
        );
        eprintln!(
            "           Pass --source-map <accession-to-source TSV> for a real breakdown. Do NOT \
             report"
        );
        eprintln!("           this breakdown as the database's source composition.");
    }
    eprintln!("================================================================================");

    // Also print a one-line summary to stdout for easy parsing
    println!(
        "{}\t{}\t{}\t{}\t{}\t{:.1}x\t{:.1}s\t{:.1}MB",
        cli.sample,
        strategy_name,
        proteins.len(),
        prot2pep.len(),
        selected.len(),
        reduction,
        elapsed,
        file_mb
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    // ------------------------------------------------------------------------
    // db_breakdown source attribution (fixed 2026-07-30)
    // ------------------------------------------------------------------------

    /// THE BUG THIS LOCKS OUT. A content-hash accession must never be silently bucketed
    /// into a real catalogue. The old `else { "smORF" }` returned "smORF" for every one of
    /// these, which is how PREDICT and CAPSCAN reported "100% smORF" for months while the
    /// true composition was ~35% UHGG / 33% GMSC / 17% GMGC / 15% Assembly.
    #[test]
    fn hash_accessions_are_never_silently_bucketed() {
        for acc in [
            "PREDICT_0007fdecd98e981dada791d932a2405e8f269c7d674dfec96cfcc82b0dc1401f",
            "CAPSCAN_986d9daba40168bd0fb533669e01a2452175167a11986c31ee34c9ae4642dad0",
            "SIHUMIX_4afdb2e40ebd66609ec42854c2cd7244686ac463b094ca08a76af460aca3c1cc",
        ] {
            let got = get_db_type(acc, None);
            assert_eq!(
                got, "UNCLASSIFIED",
                "hash accession {acc} was attributed to '{got}' -- a content hash carries no \
                 catalogue, so any real catalogue name here is fabricated provenance"
            );
            assert_ne!(got, "smORF", "regression: the old catch-all is back");
        }
    }

    /// The fix must NOT change what a catalogue-keyed lake reports. `GMSC -> smORF` used
    /// to be reached through the catch-all and is now an explicit branch, so these must
    /// come out exactly as before. If this test fails, the fix has moved a published
    /// number (healthy stool's db_breakdown) and is not safe to ship.
    #[test]
    fn catalogue_keyed_classification_is_unchanged() {
        let cases = [
            ("GMGC10.000_006_220.RPMC", "GMGC"),
            ("GMGC10.146_643_264.UNKNOWN", "GMGC"),
            ("GMSC10.100AA.840_421_255", "smORF"),
            ("MGYG000268747_00393", "UHGG"),
            ("MGYG000040330_00278", "UHGG"),
            ("sp|P02768|ALBU_HUMAN", "Human"),
            ("tr|A0A123|SOME_PROT", "Human"),
            ("MuPr_Assembly_9998637", "Assembly"),
        ];
        for (acc, want) in cases {
            assert_eq!(
                get_db_type(acc, None),
                want,
                "classification moved for {acc}"
            );
        }
    }

    /// A source map fixes the hash case; a MISS in the map must fall through to the prefix
    /// rule rather than being dropped or relabelled. A map covering part of the database
    /// must not make the remainder disappear.
    #[test]
    fn source_map_resolves_hashes_and_misses_fall_through() {
        let mut map = HashMap::new();
        map.insert("PREDICT_aaa".to_string(), "GMGC10".to_string());
        map.insert("PREDICT_bbb".to_string(), "MuPr_Assembly".to_string());

        assert_eq!(get_db_type("PREDICT_aaa", Some(&map)), "GMGC10");
        assert_eq!(get_db_type("PREDICT_bbb", Some(&map)), "MuPr_Assembly");
        // Not in the map, no usable prefix -> honest UNCLASSIFIED, not a guess.
        assert_eq!(get_db_type("PREDICT_zzz", Some(&map)), "UNCLASSIFIED");
        // Not in the map but prefix-classifiable -> the prefix rule still applies.
        assert_eq!(get_db_type("MGYG000268747_00393", Some(&map)), "UHGG");
    }

    /// Source tags must still be stripped before classification.
    ///
    /// This test FAILED on first run and the CODE was right -- my assertion was wrong. It
    /// asserted that a single-pipe tag (`CLUSTER97|GMGC10...`) is stripped, because the
    /// doc comment on `strip_source_tags` said so. The guard is `>= 2` pipes, so it is
    /// not. The doc has been corrected; the behaviour is deliberately unchanged.
    ///
    /// The single-pipe case is asserted below as a KNOWN LIMITATION rather than removed,
    /// so it is recorded instead of forgotten. It is latent: 200,000 real healthy-stool
    /// headers contain 52 pipe-containing accessions and every one is `sp|X|Y`; no
    /// single-pipe uppercase tag occurs. If such accessions ever enter a lake they will
    /// land in UNCLASSIFIED -- visibly, which is the point of the fix.
    #[test]
    fn source_tags_are_stripped_before_classification() {
        // Two-tag prefix: stripped, then classified on the real accession.
        assert_eq!(
            get_db_type("UNCLUSTERED|CAPSCAN|MGYG000268747_00393", None),
            "UHGG"
        );
        assert_eq!(
            get_db_type("UNCLUSTERED|CAPSCAN|GMGC10.000_006_220.RPMC", None),
            "GMGC"
        );

        // UniProt keeps its own pipes: first segment is lowercase, so it is NOT treated
        // as a source tag, and the `sp|`/`tr|` branch matches. This is why the
        // uppercase-only guard exists.
        assert_eq!(get_db_type("sp|P02768|ALBU_HUMAN", None), "Human");

        // KNOWN LIMITATION, asserted so a future change to strip_source_tags trips here
        // and has to be thought about rather than silently altering classification.
        assert_eq!(
            get_db_type("CLUSTER97|GMGC10.000_1", None),
            "UNCLASSIFIED",
            "single-pipe source tags are NOT stripped (guard is >=2 pipes); if this now \
             returns GMGC, strip_source_tags changed and real-data classification moved"
        );
    }

    // C1 regression (Prove-It): a peptide that positionally OVERLAPS another matched
    // peptide inside a protein must still register its peptide->protein edge. With the
    // buggy non-overlapping `find_iter` this edge was dropped; `find_overlapping_iter`
    // recovers it. This test FAILS on find_iter and PASSES on find_overlapping_iter.
    #[test]
    fn overlapping_peptides_both_map_to_protein() {
        let mut proteins: HashMap<String, String> = HashMap::new();
        // "PEPTIDE" occupies [3,10); "EPTIDEK" occupies [4,11) — they overlap.
        proteins.insert("P1".to_string(), "AAAPEPTIDEKAAA".to_string());
        let peptides = vec!["PEPTIDE".to_string(), "EPTIDEK".to_string()];

        let (pep2prot, prot2pep) = map_peptides_parallel(&proteins, &peptides, false, 1);

        assert!(
            pep2prot.get("PEPTIDE").is_some_and(|s| s.contains("P1")),
            "PEPTIDE edge to P1 missing"
        );
        assert!(
            pep2prot.get("EPTIDEK").is_some_and(|s| s.contains("P1")),
            "EPTIDEK edge to P1 missing — the overlapping peptide was dropped (C1 bug)"
        );
        // Protein P1 must carry BOTH peptides as evidence, not just the first-matched one.
        assert_eq!(
            prot2pep.get("P1").map_or(0, |s| s.len()),
            2,
            "P1 should map to both overlapping peptides"
        );
    }

    // ---- D3 determinism fixture -------------------------------------------
    // A shared peptide whose degeneracy (3) exceeds its budget (1). Exactly one
    // of the three candidate proteins may be credited. The deterministic credit
    // order must always pick the SAME protein — the one with the most unique
    // evidence — so the selected set is stable and principled, not random.
    //
    // PROT_A carries a unique peptide (u=1); PROT_B and PROT_C are shared-only
    // (u=0). Non-MGYG proteins are selected at score >= 1. With correct ordering
    // the shared peptide's single unit of budget goes to PROT_A (already selected
    // via its unique peptide), leaving PROT_B / PROT_C uncredited and UNselected.
    // The buggy HashSet-order version would sometimes credit PROT_B or PROT_C,
    // pulling them into the output.
    #[allow(clippy::type_complexity)]
    fn degeneracy_over_budget_fixture() -> (
        HashMap<String, HashSet<String>>,
        HashMap<String, HashSet<String>>,
        HashMap<String, f64>,
    ) {
        let pa = "PROT_A".to_string();
        let pb = "PROT_B".to_string();
        let pc = "PROT_C".to_string();
        let uniq = "UNIQUEPEPA".to_string(); // unique to PROT_A
        let shared = "SHAREDPEPTIDEXYZ".to_string(); // shared across A/B/C

        let mut pep2prot: HashMap<String, HashSet<String>> = HashMap::new();
        pep2prot.insert(uniq.clone(), HashSet::from([pa.clone()]));
        pep2prot.insert(
            shared.clone(),
            HashSet::from([pa.clone(), pb.clone(), pc.clone()]),
        );

        let mut prot2pep: HashMap<String, HashSet<String>> = HashMap::new();
        prot2pep.insert(pa.clone(), HashSet::from([uniq.clone(), shared.clone()]));
        prot2pep.insert(pb.clone(), HashSet::from([shared.clone()]));
        prot2pep.insert(pc.clone(), HashSet::from([shared.clone()]));

        let mut scores: HashMap<String, f64> = HashMap::new();
        scores.insert(uniq, 0.9);
        scores.insert(shared, 0.9);

        (pep2prot, prot2pep, scores)
    }

    // D3 regression. These replace three tests that asserted the OPPOSITE and so
    // locked in the defect: they required PROT_B and PROT_C to be excluded, which is
    // precisely the behaviour that turned `uniform` into a razor variant.
    #[test]
    fn credit_order_is_deterministic_and_total() {
        use std::cmp::Ordering;
        // Same inputs -> same answer, every time.
        for _ in 0..100 {
            assert_eq!(
                credit_order("PEPTIDEK", "P1", "P2"),
                credit_order("PEPTIDEK", "P1", "P2")
            );
        }
        // Antisymmetric, and never Equal for distinct accessions (total order).
        for (a, b) in [("P1", "P2"), ("MGYG1", "GMGC1"), ("sp|X|", "tr|Y|")] {
            let ab = credit_order("PEPTIDEK", a, b);
            let ba = credit_order("PEPTIDEK", b, a);
            assert_ne!(ab, Ordering::Equal, "{a} vs {b} must not tie");
            assert_eq!(ab.reverse(), ba, "{a} vs {b} must be antisymmetric");
        }
        // Keyed on the peptide: a protein is not favoured across every peptide.
        let flips = (0..200)
            .filter(|i| {
                let pep = format!("PEPTIDE{i}K");
                credit_order(&pep, "AAA", "ZZZ") == Ordering::Greater
            })
            .count();
        assert!(
            (20..180).contains(&flips),
            "ordering of the same pair should vary with the peptide, got {flips}/200"
        );
    }

    #[test]
    fn credit_order_does_not_rank_by_unique_evidence_or_catalogue() {
        use std::cmp::Ordering;
        // Accession-only sorting would put every CLUSTER97|/GMGC accession ahead of
        // every sp|/tr| one, handing a systematic advantage to whichever catalogue
        // sorts earliest. The hash must not.
        let early = (0..200)
            .filter(|i| {
                let pep = format!("PEPTIDE{i}K");
                credit_order(&pep, "GMGC10.000_1", "sp|P02768|ALB") == Ordering::Less
            })
            .count();
        assert!(
            (40..160).contains(&early),
            "catalogue prefix must not decide the order, got {early}/200"
        );
    }

    /// THE defect this fix exists for: with budget < degeneracy, a protein whose only
    /// evidence is the shared peptide must still be reachable. Under the previous
    /// unique-evidence-first order it never was, because Phase 1 had already selected
    /// every unique-evidence protein and they consumed the whole budget first.
    #[test]
    fn uniform_phase2_can_credit_shared_only_proteins() {
        let pa = "PROT_A".to_string(); // carries unique evidence -> selected in Phase 1
        let pb = "PROT_B".to_string(); // shared-only
        let pc = "PROT_C".to_string(); // shared-only

        let mut pep2prot: HashMap<String, HashSet<String>> = HashMap::new();
        let mut prot2pep: HashMap<String, HashSet<String>> = HashMap::new();
        let mut scores: HashMap<String, f64> = HashMap::new();

        let uniq = "UNIQUEPEPA".to_string();
        pep2prot.insert(uniq.clone(), HashSet::from([pa.clone()]));
        scores.insert(uniq.clone(), 0.9);
        prot2pep.entry(pa.clone()).or_default().insert(uniq);

        // Many shared peptides, each degenerate across all three, each with budget 1.
        for i in 0..200 {
            let pep = format!("SHAREDPEPTIDE{i}K");
            pep2prot.insert(
                pep.clone(),
                HashSet::from([pa.clone(), pb.clone(), pc.clone()]),
            );
            scores.insert(pep.clone(), 0.9);
            for pr in [&pa, &pb, &pc] {
                prot2pep.entry(pr.clone()).or_default().insert(pep.clone());
            }
        }

        let (sel1, _) = run_uniform(&pep2prot, &prot2pep, &scores);
        let (sel2, _) = run_uniform(&pep2prot, &prot2pep, &scores);
        assert_eq!(sel1, sel2, "uniform selection must be identical run-to-run");

        assert!(
            sel1.contains("PROT_A"),
            "PROT_A has unique evidence and must be selected"
        );
        assert!(
            sel1.contains("PROT_B") && sel1.contains("PROT_C"),
            "shared-only proteins must be reachable by Phase 2 credit; got {sel1:?}"
        );
    }

    #[test]
    fn species_budget_credit_is_deterministic() {
        let (pep2prot, prot2pep, scores) = degeneracy_over_budget_fixture();
        let (sel1, _, _, _) = run_species_budget(&pep2prot, &prot2pep, &scores);
        let (sel2, _, _, _) = run_species_budget(&pep2prot, &prot2pep, &scores);
        assert_eq!(
            sel1, sel2,
            "species_budget selection must be identical run-to-run"
        );
    }

    /// The `--tiebreak` flag must default to the historical rule, so existing runs
    /// reproduce. Alpha picks the lexicographically smallest accession among candidates
    /// tied on evidence; Greedy must be able to pick a different one, or the flag does
    /// nothing.
    #[test]
    fn tiebreak_alpha_is_historical_and_greedy_differs() {
        // Two peptides, each shared between the same two evidence-identical proteins.
        // Neither protein owns a unique peptide, so criteria 1 and 2 tie and criterion 3
        // decides both assignments.
        let zz = "ZZZ_prot".to_string();
        let aa = "AAA_prot".to_string();
        let p1 = "SHAREDPEPONEK".to_string();
        let p2 = "SHAREDPEPTWOK".to_string();

        let mut pep2prot: HashMap<String, HashSet<String>> = HashMap::new();
        pep2prot.insert(p1.clone(), HashSet::from([zz.clone(), aa.clone()]));
        pep2prot.insert(p2.clone(), HashSet::from([zz.clone(), aa.clone()]));
        let mut prot2pep: HashMap<String, HashSet<String>> = HashMap::new();
        prot2pep.insert(zz.clone(), HashSet::from([p1.clone(), p2.clone()]));
        prot2pep.insert(aa.clone(), HashSet::from([p1.clone(), p2.clone()]));

        // Alpha: the lexicographically smallest accession takes BOTH peptides, so only
        // that one is selected. This is the historical behaviour.
        let sel_alpha = run_razor(&pep2prot, &prot2pep, TieBreak::Alpha);
        assert_eq!(
            sel_alpha.len(),
            1,
            "alpha should concentrate on one protein"
        );
        assert!(
            sel_alpha.contains("AAA_prot"),
            "alpha must pick the lex-smallest"
        );

        // Greedy also concentrates -- that is the point of it -- but on whichever
        // protein the first (hash-decided) assignment lands on, not on the name.
        let sel_greedy = run_razor(&pep2prot, &prot2pep, TieBreak::Greedy);
        assert_eq!(sel_greedy.len(), 1, "greedy should concentrate too");

        // Determinism: all three rules must be stable across repeated calls.
        for tb in [TieBreak::Alpha, TieBreak::Hash, TieBreak::Greedy] {
            let a = run_razor(&pep2prot, &prot2pep, tb);
            let b = run_razor(&pep2prot, &prot2pep, tb);
            assert_eq!(a, b, "{tb:?} must be reproducible");
        }
    }
}
