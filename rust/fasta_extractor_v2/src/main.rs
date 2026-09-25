//! High-performance parallel peptide-to-protein matcher with full sequence extraction
//!
//! Uses Aho-Corasick automaton + rayon parallelization + memory-mapped I/O
//! to efficiently search massive FASTA databases and extract matching proteins.

use aho_corasick::AhoCorasick;
use clap::Parser;
use dashmap::DashMap;
use indicatif::{ProgressBar, ProgressStyle};
use memmap2::MmapOptions;
use rayon::prelude::*;
use serde::Deserialize;
#[path = "../../common/dianovo.rs"]
mod dianovo;

use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::PathBuf;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Instant;

const VERSION: &str = concat!(env!("CARGO_PKG_VERSION"), "+", env!("GIT_REV"));

#[derive(Parser, Debug)]
#[command(name = "fasta_extractor")]
#[command(version = VERSION)]
#[command(about = "Extract proteins matching de novo peptides from massive FASTA databases")]
struct Args {
    /// Directory containing AlphaNovo prediction CSVs
    #[arg(short = 'p', long)]
    peptides_dir: PathBuf,

    /// FASTA database to search
    #[arg(short = 'd', long)]
    database: PathBuf,

    /// Output FASTA file for matched proteins
    #[arg(short = 'o', long)]
    output: PathBuf,

    /// Source tag for output headers (e.g., GMGC, GMSC, UHGP)
    #[arg(short = 's', long, default_value = "DB")]
    source_tag: String,

    /// Minimum peptide length. Default 9 (v5) — raised from 8 because SAGE digests
    /// FASTA to tryptic peptides ≥9; shorter de novo peptides admit proteins with
    /// low-quality evidence that inflate entrapment FDP. See DESIGN_V5.md.
    #[arg(long, default_value = "9")]
    min_length: usize,

    /// Maximum peptide length
    #[arg(long, default_value = "50")]
    max_length: usize,

    /// Minimum peptide score (flat threshold, ignored if --top-percent is set)
    #[arg(long, default_value = "0.0")]
    min_score: f64,

    /// Keep top X% of peptides by score PER SAMPLE (e.g., 0.30 for top 30%)
    #[arg(long)]
    top_percent: Option<f64>,

    /// Number of threads (default: all available)
    #[arg(short = 't', long)]
    threads: Option<usize>,

    /// Chunk size for parallel processing
    #[arg(long, default_value = "50000")]
    chunk_size: usize,

    /// Also output a TSV of peptide-protein hits
    #[arg(long)]
    output_hits: Option<PathBuf>,

    /// Normalize I/L (treat I and L as equivalent for mass spec) - enabled by default
    #[arg(long, default_value_t = true, action = clap::ArgAction::Set)]
    normalize_il: bool,

    /// Output per-protein evidence manifest TSV with confidence scores.
    /// Columns: protein_id, n_peptides, n_samples, mean_score, max_score, support_score
    #[arg(long)]
    output_evidence: Option<PathBuf>,

    /// Known-absent FASTA for an extraction-stage negative-control match diagnostic.
    /// These proteins are searched with the same automaton but excluded from
    /// the output FASTA. An entrapment report is written alongside the stats JSON.
    /// Example: Arabidopsis proteins in a human/gut experiment.
    #[arg(long)]
    entrapment: Option<PathBuf>,

    /// Entrapment method for auto-generation (use INSTEAD of --entrapment).
    /// "shuffled": shuffle sequences from the database preserving AA composition (Varela).
    /// "noble" is retained for CLI recognition but rejected: random splitting is not a null.
    /// When set, --entrapment-n controls how many proteins to generate.
    #[arg(long, value_parser = ["shuffled", "noble"])]
    entrapment_method: Option<String>,

    /// Number of entrapment proteins to generate (for --entrapment-method shuffled).
    /// Used for shuffled diagnostics. Default: 5000.
    #[arg(long, default_value = "5000")]
    entrapment_n: usize,

    /// Random seed for entrapment generation. Default: 42.
    #[arg(long, default_value = "42")]
    entrapment_seed: u64,

    /// Compatibility flag. Numeric annotations are removed by the CSV parser
    /// regardless of this flag. Named PTM notation is not supported.
    #[arg(long)]
    strip_mods: bool,

    /// Include FASTA (aka --exception / --hypothesis): always include these proteins
    /// in the output, regardless of peptide evidence in this sample. Two uses:
    /// (1) cross-sample evidence propagation — a "core proteome" of proteins detected
    /// across many samples put into every sample's database; (2) targeted inclusion —
    /// named candidate proteins (a hypothesis set, or a low-abundance organism's
    /// reference proteome) that bypass the de novo membership gate and enter the search
    /// database directly, so a low-abundance member missed by de novo can still be
    /// searched for. In both cases the search engine (SAGE/DIA-NN) still controls FDR;
    /// inclusion never asserts an identification, only the opportunity to find one.
    #[arg(long, alias = "exception", alias = "hypothesis")]
    include: Option<PathBuf>,
}

#[derive(Debug, Deserialize)]
struct AlphaNovoRow {
    #[serde(rename = "peptide_prediction_detokenized_unmodified")]
    sequence: Option<String>,
    #[serde(rename = "score")]
    score: Option<f64>,
}

#[derive(Debug, Deserialize)]
struct DIANovoRow {
    pred_seq: Option<String>,
    pred_prob: Option<String>,
}

/// Detect CSV format from headers and return (sequence, score) pairs
fn parse_prediction_csv(
    path: &std::path::Path,
    min_length: usize,
    max_length: usize,
    _do_normalize_il: bool,
) -> Vec<(String, f64)> {
    let file = match File::open(path) {
        Ok(f) => f,
        Err(_) => return Vec::new(),
    };

    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(true)
        .flexible(true)
        .from_reader(BufReader::new(file));

    let headers: Vec<String> = match rdr.headers() {
        Ok(h) => h
            .iter()
            .map(|s| s.trim().trim_matches('\r').to_string())
            .collect(),
        Err(_) => return Vec::new(),
    };

    let is_dianovo = headers.iter().any(|h| h == "pred_seq");
    let is_alphanova = headers
        .iter()
        .any(|h| h == "peptide_prediction_detokenized_unmodified");
    let is_casanovo =
        headers.iter().any(|h| h == "sequence") && headers.iter().any(|h| h == "score");

    let mut peptides: Vec<(String, f64)> = Vec::new();

    if is_dianovo {
        let mut invalid = 0usize;
        // Keep the explicit `if let Ok` so malformed CSV rows are skipped, not propagated;
        // `.flatten()` would silently discard deserialize errors. See CLAUDE.md.
        #[allow(clippy::manual_flatten)]
        for result in rdr.deserialize::<DIANovoRow>() {
            if let Ok(row) = result {
                if let (Some(seq), Some(prob_str)) = (row.pred_seq, row.pred_prob) {
                    if let Some((clean_seq, score)) = dianovo::parse_prediction(&seq, &prob_str) {
                        if clean_seq.len() >= min_length
                            && clean_seq.len() <= max_length
                            && score.is_finite()
                            && score > 0.0
                            && !is_low_complexity(&clean_seq)
                        {
                            peptides.push((clean_seq, score));
                        }
                    } else {
                        invalid += 1;
                    }
                } else {
                    invalid += 1;
                }
            } else {
                invalid += 1;
            }
        }
        if invalid > 0 {
            eprintln!(
                "DIANovo: excluded {invalid} unresolved/invalid prediction rows in {}",
                path.display()
            );
        }
    } else if is_alphanova {
        // Explicit `if let Ok` skips malformed rows; `.flatten()` would drop errors silently.
        #[allow(clippy::manual_flatten)]
        for result in rdr.deserialize::<AlphaNovoRow>() {
            if let Ok(row) = result {
                if let (Some(seq), Some(score)) = (row.sequence, row.score) {
                    let clean_seq: String =
                        seq.chars().filter(|c| c.is_ascii_uppercase()).collect();
                    if clean_seq.len() >= min_length
                        && clean_seq.len() <= max_length
                        && score.is_finite()
                        && score > 0.0
                        && !is_low_complexity(&clean_seq)
                    {
                        peptides.push((clean_seq, score));
                    }
                }
            }
        }
    } else if is_casanovo {
        // Generic format: "sequence" + "score" columns
        #[derive(Deserialize)]
        struct GenericRow {
            sequence: Option<String>,
            score: Option<f64>,
        }
        // Explicit `if let Ok` skips malformed rows; `.flatten()` would drop errors silently.
        #[allow(clippy::manual_flatten)]
        for result in rdr.deserialize::<GenericRow>() {
            if let Ok(row) = result {
                if let (Some(seq), Some(score)) = (row.sequence, row.score) {
                    let clean_seq: String =
                        seq.chars().filter(|c| c.is_ascii_uppercase()).collect();
                    if clean_seq.len() >= min_length
                        && clean_seq.len() <= max_length
                        && score.is_finite()
                        && score > 0.0
                        && !is_low_complexity(&clean_seq)
                    {
                        peptides.push((clean_seq, score));
                    }
                }
            }
        }
    } else {
        eprintln!("WARNING: Unrecognized CSV format in {:?}. Expected columns: pred_seq/pred_prob (DIANovo), peptide_prediction_detokenized_unmodified/score (AlphaNovo), or sequence/score (generic).", path);
    }

    peptides
}

/// Parse a FASTA file into (header_start, header_end, seq_start, seq_end) byte positions
fn parse_fasta_regions(data: &[u8]) -> Vec<(usize, usize, usize, usize)> {
    let mut regions = Vec::with_capacity(1_000_000);
    let mut i = 0;
    let len = data.len();

    while i < len {
        if data[i] == b'>' {
            let header_start = i + 1; // skip '>'

            // Find end of header
            while i < len && data[i] != b'\n' {
                i += 1;
            }
            let header_end = i;
            i += 1; // skip newline

            let seq_start = i;

            // Find end of sequence
            while i < len && data[i] != b'>' {
                i += 1;
            }
            let seq_end = i;

            if seq_end > seq_start {
                regions.push((header_start, header_end, seq_start, seq_end));
            }
        } else {
            i += 1;
        }
    }
    regions
}

/// Extract sequence bytes, removing newlines
fn extract_sequence(data: &[u8], start: usize, end: usize) -> Vec<u8> {
    let mut seq = Vec::with_capacity(end - start);
    for &byte in &data[start..end] {
        if byte != b'\n' && byte != b'\r' {
            seq.push(byte);
        }
    }
    seq
}

/// Normalize I to L for mass spec equivalence
fn normalize_il_seq(seq: &str) -> String {
    seq.replace("I", "L")
}

/// Reject low-complexity peptides: poly-amino acid repeats, dipeptide repeats,
/// and sequences dominated by a single residue.
fn is_low_complexity(seq: &str) -> bool {
    let n = seq.len();
    if n == 0 {
        return true;
    }
    let bytes = seq.as_bytes();

    // Single amino acid dominance: >60% one residue
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

    // Too few unique AAs for longer peptides (e.g. 12-mer with only 2 AAs)
    let unique_aas = counts.iter().filter(|&&c| c > 0).count();
    if n >= 12 && unique_aas <= 2 {
        return true;
    }

    // Dipeptide repeat (ALALAL, AGAGAG, etc.)
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

/// Peptide collection with scores and sample provenance.
struct PeptideCollection {
    /// Unique peptide sequences (for Aho-Corasick)
    sequences: HashSet<String>,
    /// Best de novo score per peptide
    best_scores: HashMap<String, f64>,
    /// Number of samples each peptide was observed in
    sample_counts: HashMap<String, usize>,
    sample_ids: HashMap<String, HashSet<usize>>,
    /// Total CSV files (samples) processed
    total_samples: usize,
}

impl PeptideCollection {
    fn new() -> Self {
        PeptideCollection {
            sequences: HashSet::new(),
            best_scores: HashMap::new(),
            sample_counts: HashMap::new(),
            sample_ids: HashMap::new(),
            total_samples: 0,
        }
    }

    fn add(&mut self, seq: String, score: f64, is_new_sample: bool) {
        let entry = self.best_scores.entry(seq.clone()).or_insert(0.0);
        if score > *entry {
            *entry = score;
        }
        let _ = is_new_sample; // Identity-based membership also handles duplicate CSV rows.
        let ids = self.sample_ids.entry(seq.clone()).or_default();
        ids.insert(self.total_samples);
        self.sample_counts.insert(seq.clone(), ids.len());
        self.sequences.insert(seq);
    }
}

/// Recursively collect all CSV files from a directory
fn collect_csvs(dir: &std::path::Path) -> Vec<PathBuf> {
    let mut out = Vec::new();
    fn walk(dir: &std::path::Path, out: &mut Vec<PathBuf>) {
        if let Ok(entries) = std::fs::read_dir(dir) {
            for entry in entries.filter_map(|e| e.ok()) {
                let p = entry.path();
                if p.is_dir() {
                    walk(&p, out);
                } else if p.extension().is_some_and(|ext| ext == "csv") {
                    out.push(p);
                }
            }
        }
    }
    walk(dir, &mut out);
    out.sort();
    out
}

/// Load peptides from all AlphaNovo CSVs with TOP X% filtering per sample
/// Pool retained prediction queries after ranking rows within each CSV.
///
/// The top-fraction score cutoff includes ties. Sequence deduplication and optional
/// I/L normalization follow eligibility filtering and per-file ranking. Reference
/// chunking calls this same procedure with the complete prediction roster.
fn load_peptides_top_percent(
    dir: &std::path::Path,
    min_length: usize,
    max_length: usize,
    top_fraction: f64,
    do_normalize_il: bool,
) -> PeptideCollection {
    let mut collection = PeptideCollection::new();

    let csv_files = collect_csvs(dir);

    println!(
        "Loading peptides from {} CSV files (top {:.0}% per sample)...",
        csv_files.len(),
        top_fraction * 100.0
    );

    let pb = ProgressBar::new(csv_files.len() as u64);
    pb.set_style(ProgressStyle::default_bar()
        .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} files | {msg}")
        .unwrap());

    let mut total_retained = 0usize;
    let mut total_rows = 0usize;

    for path in &csv_files {
        let file_peptides = parse_prediction_csv(path, min_length, max_length, do_normalize_il);

        total_rows += file_peptides.len();

        // Track which peptides are new for this sample (for sample_counts)
        let mut seen_in_sample: HashSet<String> = HashSet::new();

        if !file_peptides.is_empty() {
            let mut scores: Vec<f64> = file_peptides.iter().map(|(_, s)| *s).collect();
            scores.sort_by(|a, b| b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal));

            let cutoff_idx = ((scores.len() as f64 * top_fraction).ceil() as usize)
                .min(scores.len())
                .max(1)
                - 1;
            let threshold = scores[cutoff_idx];

            for (seq, score) in file_peptides {
                if score >= threshold && !is_low_complexity(&seq) {
                    let final_seq = if do_normalize_il {
                        normalize_il_seq(&seq)
                    } else {
                        seq
                    };
                    let is_new = seen_in_sample.insert(final_seq.clone());
                    collection.add(final_seq, score, is_new);
                    total_retained += 1;
                }
            }
        }

        collection.total_samples += 1;
        pb.inc(1);
        pb.set_message(format!("{} unique", collection.sequences.len()));
    }

    pb.finish_with_message(format!(
        "{} unique peptides from {} retained",
        collection.sequences.len(),
        total_retained
    ));

    println!("  Total rows processed: {}", total_rows);
    println!(
        "  Retained (top {:.0}%): {}",
        top_fraction * 100.0,
        total_retained
    );
    println!("  Unique peptides: {}", collection.sequences.len());
    println!("  Samples processed: {}", collection.total_samples);

    collection
}

/// Load peptides from all AlphaNovo CSVs in directory (flat threshold)
fn load_peptides(
    dir: &std::path::Path,
    min_length: usize,
    max_length: usize,
    min_score: f64,
    do_normalize_il: bool,
) -> PeptideCollection {
    let mut collection = PeptideCollection::new();

    let csv_files = collect_csvs(dir);

    println!(
        "Loading peptides from {} CSV files (score >= {})...",
        csv_files.len(),
        min_score
    );

    let pb = ProgressBar::new(csv_files.len() as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} files")
            .unwrap(),
    );

    for path in csv_files {
        let file_peptides = parse_prediction_csv(&path, min_length, max_length, do_normalize_il);
        let mut seen_in_sample: HashSet<String> = HashSet::new();
        for (seq, score) in file_peptides {
            if score >= min_score {
                let final_seq = if do_normalize_il {
                    normalize_il_seq(&seq)
                } else {
                    seq
                };
                let is_new = seen_in_sample.insert(final_seq.clone());
                collection.add(final_seq, score, is_new);
            }
        }
        collection.total_samples += 1;
        pb.inc(1);
    }
    pb.finish_with_message("Peptides loaded");

    collection
}

/// Load peptides from a single CSV file (flat threshold, auto-detects format)
fn load_peptides_single_csv(
    path: &PathBuf,
    min_length: usize,
    max_length: usize,
    min_score: f64,
    do_normalize_il: bool,
) -> PeptideCollection {
    let mut collection = PeptideCollection::new();

    println!(
        "Loading peptides from single CSV: {:?} (score >= {})...",
        path, min_score
    );

    let file_peptides = parse_prediction_csv(path, min_length, max_length, do_normalize_il);
    for (seq, score) in file_peptides {
        if score >= min_score {
            let final_seq = if do_normalize_il {
                normalize_il_seq(&seq)
            } else {
                seq
            };
            collection.add(final_seq, score, true);
        }
    }
    collection.total_samples = 1;

    println!("  Loaded {} unique peptides", collection.sequences.len());
    collection
}

/// Load peptides from a single CSV with TOP X% filtering (auto-detects format)
fn load_peptides_top_percent_single_csv(
    path: &PathBuf,
    min_length: usize,
    max_length: usize,
    top_fraction: f64,
    do_normalize_il: bool,
) -> PeptideCollection {
    let mut collection = PeptideCollection::new();

    println!(
        "Loading peptides from single CSV: {:?} (top {:.0}%)...",
        path,
        top_fraction * 100.0
    );

    let file_peptides = parse_prediction_csv(path, min_length, max_length, do_normalize_il);

    println!("  Total rows: {}", file_peptides.len());

    if !file_peptides.is_empty() {
        let mut scores: Vec<f64> = file_peptides.iter().map(|(_, s)| *s).collect();
        scores.sort_by(|a, b| b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal));

        let cutoff_idx = ((scores.len() as f64 * top_fraction).ceil() as usize)
            .min(scores.len())
            .max(1)
            - 1;
        let threshold = scores[cutoff_idx];

        let mut retained = 0usize;
        for (seq, score) in file_peptides {
            if score >= threshold && !is_low_complexity(&seq) {
                let final_seq = if do_normalize_il {
                    normalize_il_seq(&seq)
                } else {
                    seq
                };
                collection.add(final_seq, score, true);
                retained += 1;
            }
        }

        println!("  Score threshold: {:.4}", threshold);
        println!(
            "  Retained (top {:.0}%): {}",
            top_fraction * 100.0,
            retained
        );
    }
    collection.total_samples = 1;

    println!("  Unique peptides: {}", collection.sequences.len());
    collection
}

/// Load peptides from a simple text file (one per line).
/// Matches the CSV path: applies the max-length upper bound and I/L normalization.
fn load_peptides_txt(
    path: &PathBuf,
    min_length: usize,
    max_length: usize,
    do_normalize_il: bool,
) -> PeptideCollection {
    let mut collection = PeptideCollection::new();

    if let Ok(file) = File::open(path) {
        let reader = BufReader::new(file);
        // NOTE: kept an explicit `if let Ok(seq) = line` rather than `.flatten()`
        // so read errors are not silently swallowed differently than before
        // (behavior is at least as strict as the original).
        #[allow(clippy::manual_flatten)]
        for line in reader.lines() {
            if let Ok(seq) = line {
                let clean: String = seq
                    .trim()
                    .chars()
                    .filter(|c| c.is_ascii_uppercase())
                    .collect();
                if clean.len() >= min_length && clean.len() <= max_length {
                    let final_seq = if do_normalize_il {
                        normalize_il_seq(&clean)
                    } else {
                        clean
                    };
                    collection.add(final_seq, 1.0, true);
                }
            }
        }
    }
    collection.total_samples = 1;

    collection
}

/// Match retained peptide queries against a prepared protein FASTA.
///
/// Write complete proteins containing exact query substrings, together with the
/// requested evidence/hit tables. I/L equivalence affects matching, while original
/// FASTA sequence and header identities are retained. No razor or clustering is
/// performed by this executable.
fn main() -> std::io::Result<()> {
    let args = Args::parse();
    let start_time = Instant::now();
    let invalid = |message: &str| std::io::Error::new(std::io::ErrorKind::InvalidInput, message);
    if args.min_length == 0 || args.max_length < args.min_length {
        return Err(invalid(
            "Lengths must satisfy 1 <= min-length <= max-length",
        ));
    }
    if args.chunk_size == 0 || args.threads == Some(0) {
        return Err(invalid("chunk-size and threads must be positive"));
    }
    if !args.min_score.is_finite()
        || args
            .top_percent
            .is_some_and(|v| !v.is_finite() || v <= 0.0 || v > 1.0)
    {
        return Err(invalid(
            "Scores must be finite; top-percent is a fraction in (0, 1], e.g. 0.30",
        ));
    }
    if !args.database.is_file() || !args.peptides_dir.exists() {
        return Err(invalid("Database and prediction inputs must exist"));
    }
    for path in [args.include.as_ref(), args.entrapment.as_ref()]
        .into_iter()
        .flatten()
    {
        if !path.is_file() {
            return Err(invalid(&format!(
                "Requested FASTA does not exist: {}",
                path.display()
            )));
        }
    }
    if args.entrapment.is_some() && args.entrapment_method.is_some() {
        return Err(invalid("Choose entrapment OR entrapment-method, not both"));
    }
    if args.entrapment_method.as_deref() == Some("noble") {
        return Err(invalid("Randomly splitting the reference does not establish known-absent entrapments; use an independently justified entrapment FASTA"));
    }
    if args.entrapment_method.is_some() && args.entrapment_n == 0 {
        return Err(invalid("entrapment-n must be positive"));
    }
    let stats_path = args.output.with_extension("stats.json");
    let mut destinations = HashSet::new();
    for path in [
        Some(&args.output),
        Some(&stats_path),
        args.output_hits.as_ref(),
        args.output_evidence.as_ref(),
    ]
    .into_iter()
    .flatten()
    {
        if path.exists() || path.is_symlink() {
            return Err(invalid(&format!(
                "Refusing to overwrite existing output: {}",
                path.display()
            )));
        }
        let parent = path
            .parent()
            .filter(|p| !p.as_os_str().is_empty())
            .unwrap_or(std::path::Path::new("."));
        let resolved = parent.canonicalize()?.join(
            path.file_name()
                .ok_or_else(|| invalid("Output must name a file"))?,
        );
        if !destinations.insert(resolved) {
            return Err(invalid("Output files must have distinct paths"));
        }
    }

    // C11: seed 0 makes xorshift64 a fixed point (all-zero stream -> identity
    // permutation), which silently defeats the shuffle/split entrapment. Reject it.
    if args.entrapment_method.is_some() && args.entrapment_seed == 0 {
        eprintln!("ERROR: --entrapment-seed 0 produces a degenerate all-zero xorshift64 stream");
        eprintln!(
            "       (identity permutation), silently defeating entrapment. Use a nonzero seed."
        );
        std::process::exit(1);
    }

    // Configure thread pool
    let num_threads = args.threads.unwrap_or_else(|| {
        std::thread::available_parallelism()
            .map(|p| p.get())
            .unwrap_or(8)
    });

    rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build_global()
        .ok();

    println!("=======================================================");
    println!("FASTA EXTRACTOR - High Performance Parallel Search");
    println!("=======================================================");
    println!("Threads: {}", num_threads);
    println!("Database: {:?}", args.database);
    println!("Output: {:?}", args.output);
    println!("Source tag: {}", args.source_tag);
    println!("Peptide length: {}-{} aa", args.min_length, args.max_length);
    if let Some(top_pct) = args.top_percent {
        println!("Score filter: top {:.0}% per sample", top_pct * 100.0);
    } else {
        println!("Score filter: >= {}", args.min_score);
    }
    println!("I/L normalization: {}", args.normalize_il);
    println!();

    // Load peptides
    println!("Step 1: Loading peptides...");
    let collection: PeptideCollection = if args.peptides_dir.is_dir() {
        if let Some(top_pct) = args.top_percent {
            load_peptides_top_percent(
                &args.peptides_dir,
                args.min_length,
                args.max_length,
                top_pct,
                args.normalize_il,
            )
        } else {
            load_peptides(
                &args.peptides_dir,
                args.min_length,
                args.max_length,
                args.min_score,
                args.normalize_il,
            )
        }
    } else if args
        .peptides_dir
        .extension()
        .is_some_and(|ext| ext == "csv")
    {
        if let Some(top_pct) = args.top_percent {
            load_peptides_top_percent_single_csv(
                &args.peptides_dir,
                args.min_length,
                args.max_length,
                top_pct,
                args.normalize_il,
            )
        } else {
            load_peptides_single_csv(
                &args.peptides_dir,
                args.min_length,
                args.max_length,
                args.min_score,
                args.normalize_il,
            )
        }
    } else {
        load_peptides_txt(
            &args.peptides_dir,
            args.min_length,
            args.max_length,
            args.normalize_il,
        )
    };

    let peptides = &collection.sequences;
    println!("  Loaded {} unique peptides", peptides.len());

    if peptides.is_empty() {
        eprintln!("ERROR: No peptides loaded from {:?}!", args.peptides_dir);
        std::process::exit(1);
    }

    // Sanity: check peptides look like amino acid sequences
    let non_aa_count = peptides
        .iter()
        .filter(|p| !p.chars().all(|c| "ACDEFGHIKLMNPQRSTVWY".contains(c)))
        .count();
    if non_aa_count > peptides.len() / 10 {
        eprintln!(
            "WARNING: {}% of peptides contain non-standard amino acids. Check input format.",
            non_aa_count * 100 / peptides.len()
        );
    }

    // Sanity: typical peptide length distribution
    let avg_len: f64 = peptides.iter().map(|p| p.len() as f64).sum::<f64>() / peptides.len() as f64;
    if !(7.0..=30.0).contains(&avg_len) {
        eprintln!(
            "WARNING: Average peptide length {:.1} is unusual (expected 8-25). Check input format.",
            avg_len
        );
    }

    // Build Aho-Corasick automaton
    println!("\nStep 2: Building Aho-Corasick automaton...");
    let ac_start = Instant::now();
    let peptide_vec: Vec<&str> = peptides.iter().map(|s| s.as_str()).collect();
    let ac = AhoCorasick::new(&peptide_vec).expect("Failed to build automaton");
    println!("  Automaton built in {:.2?}", ac_start.elapsed());

    // Memory-map the database
    println!("\nStep 3: Memory-mapping database...");
    let file = File::open(&args.database)?;
    let file_size = file.metadata()?.len();
    println!("  Database size: {:.2} GB", file_size as f64 / 1e9);

    // SAFETY: the input FASTA database is opened read-only and is not mutated or
    // truncated by any other process for the duration of this run.
    let mmap = unsafe { MmapOptions::new().map(&file)? };
    let data = &mmap[..];

    // Parse FASTA regions
    println!("\nStep 4: Parsing FASTA structure...");
    let parse_start = Instant::now();
    let regions = parse_fasta_regions(data);
    let total_proteins = regions.len();
    println!(
        "  Found {} proteins in {:.2?}",
        total_proteins,
        parse_start.elapsed()
    );

    // Parallel search
    println!(
        "\nStep 5: Parallel search ({} threads, chunk size {})...",
        num_threads, args.chunk_size
    );
    let search_start = Instant::now();

    let matched_indices: DashMap<usize, Vec<(String, usize)>> = DashMap::new();
    let proteins_searched = AtomicUsize::new(0);
    let matches_found = AtomicUsize::new(0);

    let pb = ProgressBar::new(total_proteins as u64);
    pb.set_style(ProgressStyle::default_bar()
        .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} proteins ({per_sec}) | Matches: {msg}")
        .unwrap());

    // Process in parallel chunks
    let do_normalize = args.normalize_il;
    let chunk_size = args.chunk_size;
    regions
        .par_chunks(chunk_size)
        .enumerate()
        .for_each(|(chunk_idx, chunk)| {
            let chunk_start_idx = chunk_idx * chunk_size;
            proteins_searched.fetch_add(chunk.len(), Ordering::Relaxed);

            for (local_idx, &(_header_start, _header_end, seq_start, seq_end)) in
                chunk.iter().enumerate()
            {
                let idx = chunk_start_idx + local_idx;
                let seq = extract_sequence(data, seq_start, seq_end);

                // Normalize I/L in protein sequence if enabled (same as peptides)
                let search_seq: Vec<u8> = if do_normalize {
                    seq.iter()
                        .map(|&b| if b == b'I' { b'L' } else { b })
                        .collect()
                } else {
                    seq
                };

                // Find all peptide matches
                let mut hits: Vec<(String, usize)> = Vec::new();
                for mat in ac.find_overlapping_iter(&search_seq) {
                    let peptide = peptide_vec[mat.pattern().as_usize()].to_string();
                    hits.push((peptide, mat.start()));
                }

                if !hits.is_empty() {
                    matches_found.fetch_add(1, Ordering::Relaxed);
                    matched_indices.insert(idx, hits);
                }
            }

            pb.inc(chunk.len() as u64);
            pb.set_message(format!("{}", matches_found.load(Ordering::Relaxed)));
        });

    pb.finish_with_message(format!(
        "{} proteins matched",
        matches_found.load(Ordering::Relaxed)
    ));
    println!("  Search completed in {:.2?}", search_start.elapsed());

    let total_matched = matched_indices.len();
    println!(
        "  Matched proteins: {} ({:.2}%)",
        total_matched,
        100.0 * total_matched as f64 / total_proteins as f64
    );

    // Sanity: warn on extremely low match rates
    let match_rate = total_matched as f64 / total_proteins as f64;
    if match_rate < 0.001 && total_proteins > 100 {
        eprintln!(
            "WARNING: Match rate {:.4}% is extremely low. Possible input mismatch.",
            match_rate * 100.0
        );
        eprintln!("  Peptides loaded: {}", peptides.len());
        eprintln!("  Proteins searched: {}", total_proteins);
        let sample: Vec<&String> = peptides.iter().take(5).collect();
        eprintln!("  Sample peptides: {:?}", sample);
    }

    // Step 5b: Entrapment search (if requested)
    // Supports: --entrapment <fasta> (user/phylogenetic)
    //           --entrapment-method shuffled (Varela: shuffle DB sequences)
    //           --entrapment-method noble (Noble 2009: split DB 50/50)
    let mut entrapment_matched: usize = 0;
    let mut entrapment_total: usize = 0;
    let mut entrapment_hits_detail: Vec<(String, usize)> = Vec::new();
    let mut entrapment_method_name = String::new();
    let mut noble_excluded: HashSet<usize> = HashSet::new(); // indices excluded by noble split

    let do_entrapment = args.entrapment.is_some() || args.entrapment_method.is_some();

    if do_entrapment {
        println!("\nStep 5b: Entrapment search...");
        let entrap_start = Instant::now();

        // Build entrapment protein list: Vec<(id, sequence_bytes)>
        let mut entrap_proteins: Vec<(String, Vec<u8>)>;

        if let Some(ref method) = args.entrapment_method {
            entrapment_method_name = method.clone();

            match method.as_str() {
                "shuffled" => {
                    // Varela method: shuffle database sequences preserving AA composition
                    println!(
                        "  Method: shuffled (Varela) — {} proteins, seed {}",
                        args.entrapment_n, args.entrapment_seed
                    );
                    let n = args.entrapment_n.min(total_proteins);
                    let mut rng_state = args.entrapment_seed;

                    // Simple xorshift64 PRNG for reproducibility without external crate
                    let mut next_rand = move || -> u64 {
                        rng_state ^= rng_state << 13;
                        rng_state ^= rng_state >> 7;
                        rng_state ^= rng_state << 17;
                        rng_state
                    };

                    // Sample n proteins from the database, shuffle their sequences
                    let mut sampled: Vec<usize> = (0..total_proteins).collect();
                    // Fisher-Yates partial shuffle to get n random indices
                    for i in 0..n {
                        let j = i + (next_rand() as usize % (total_proteins - i));
                        sampled.swap(i, j);
                    }
                    sampled.truncate(n);

                    let mut shuffled_proteins = Vec::with_capacity(n);
                    for (eidx, &pidx) in sampled.iter().enumerate() {
                        let (_, _, seq_start, seq_end) = regions[pidx];
                        let mut seq = extract_sequence(data, seq_start, seq_end);

                        // Fisher-Yates shuffle of the sequence
                        let slen = seq.len();
                        if slen > 1 {
                            for i in (1..slen).rev() {
                                let j = next_rand() as usize % (i + 1);
                                seq.swap(i, j);
                            }
                        }

                        shuffled_proteins.push((format!("ENTRAP_SHUFFLED_{:05}", eidx), seq));
                    }
                    entrap_proteins = shuffled_proteins;
                }
                "noble" => {
                    // Noble 2009: split database 50/50
                    println!(
                        "  Method: noble paired — split {} proteins 50/50, seed {}",
                        total_proteins, args.entrapment_seed
                    );

                    let mut indices: Vec<usize> = (0..total_proteins).collect();
                    let mut rng_state = args.entrapment_seed;
                    let mut next_rand = move || -> u64 {
                        rng_state ^= rng_state << 13;
                        rng_state ^= rng_state >> 7;
                        rng_state ^= rng_state << 17;
                        rng_state
                    };

                    // Fisher-Yates shuffle
                    for i in (1..indices.len()).rev() {
                        let j = next_rand() as usize % (i + 1);
                        indices.swap(i, j);
                    }

                    let mid = indices.len() / 2;
                    // First half = target (stays in search), second half = entrapment
                    let entrap_indices = &indices[mid..];

                    // Mark entrapment indices for exclusion from output
                    for &idx in entrap_indices {
                        noble_excluded.insert(idx);
                    }

                    let mut noble_proteins = Vec::with_capacity(entrap_indices.len());
                    for &pidx in entrap_indices {
                        let (header_start, header_end, seq_start, seq_end) = regions[pidx];
                        let header =
                            std::str::from_utf8(&data[header_start..header_end]).unwrap_or("?");
                        let protein_id = header.split_whitespace().next().unwrap_or(header);
                        let seq = extract_sequence(data, seq_start, seq_end);
                        noble_proteins.push((format!("ENTRAP_NOBLE_{}", protein_id), seq));
                    }

                    println!(
                        "  Target: {} proteins, Entrapment: {} proteins",
                        mid,
                        entrap_indices.len()
                    );
                    entrap_proteins = noble_proteins;
                }
                _ => {
                    eprintln!("ERROR: Unknown entrapment method: {}", method);
                    entrap_proteins = Vec::new();
                }
            }
        } else if let Some(ref entrap_path) = args.entrapment {
            // User-supplied FASTA (phylogenetic / Kleiner / custom)
            entrapment_method_name = "phylogenetic".to_string();
            if !entrap_path.exists() {
                eprintln!("WARNING: Entrapment FASTA not found: {:?}", entrap_path);
                entrap_proteins = Vec::new();
            } else {
                println!("  Method: phylogenetic (user FASTA: {:?})", entrap_path);
                let entrap_file = File::open(entrap_path)?;
                // SAFETY: the entrapment FASTA is opened read-only and is not mutated
                // or truncated by any other process for the duration of this run.
                let entrap_mmap = unsafe { MmapOptions::new().map(&entrap_file)? };
                let entrap_data = &entrap_mmap[..];
                let entrap_regions = parse_fasta_regions(entrap_data);

                let mut loaded = Vec::with_capacity(entrap_regions.len());
                for &(hs, he, ss, se) in &entrap_regions {
                    let header = std::str::from_utf8(&entrap_data[hs..he]).unwrap_or("?");
                    let protein_id = header.split_whitespace().next().unwrap_or(header);
                    let seq = extract_sequence(entrap_data, ss, se);
                    loaded.push((protein_id.to_string(), seq));
                }
                entrap_proteins = loaded;
            }
        } else {
            entrap_proteins = Vec::new();
        }

        entrap_proteins.retain(|(_, sequence)| !sequence.is_empty());
        entrapment_total = entrap_proteins.len();
        if entrapment_total == 0 {
            return Err(invalid("Entrapment FASTA contains no nonempty sequences"));
        }
        println!("  Entrapment proteins: {}", entrapment_total);

        // Search ALL entrapment proteins with the SAME automaton
        for (pid, seq) in &entrap_proteins {
            let search_seq: Vec<u8> = if args.normalize_il {
                seq.iter()
                    .map(|&b| if b == b'I' { b'L' } else { b })
                    .collect()
            } else {
                seq.clone()
            };

            let mut n_hits = 0usize;
            for _mat in ac.find_overlapping_iter(&search_seq) {
                n_hits += 1;
            }

            if n_hits > 0 {
                entrapment_matched += 1;
                entrapment_hits_detail.push((pid.clone(), n_hits));
            }
        }

        let entrap_rate = if entrapment_total > 0 {
            entrapment_matched as f64 / entrapment_total as f64
        } else {
            0.0
        };

        println!(
            "  Entrapment matched: {} / {} ({:.4}%)",
            entrapment_matched,
            entrapment_total,
            entrap_rate * 100.0
        );
        println!("  Entrapment search: {:.2?}", entrap_start.elapsed());

        if entrapment_matched > 0 {
            println!(
                "  WARNING: {} entrapment proteins had peptide hits!",
                entrapment_matched
            );
            for (pid, n) in entrapment_hits_detail.iter().take(20) {
                println!("    {} ({} peptides)", pid, n);
            }
            if entrapment_hits_detail.len() > 20 {
                println!("    ... and {} more", entrapment_hits_detail.len() - 20);
            }
        } else {
            println!("  No entrapment matches observed; this does not certify identification FDR.");
        }
    }

    // Write output FASTA
    println!("\nStep 6: Writing output FASTA...");
    let write_start = Instant::now();

    let output_file = File::create(&args.output)?;
    let mut writer = BufWriter::with_capacity(16 * 1024 * 1024, output_file);

    // Collect and sort matched indices for deterministic output
    // For noble paired: exclude entrapment-half proteins from output
    let mut sorted_indices: Vec<usize> = matched_indices
        .iter()
        .map(|r| *r.key())
        .filter(|idx| !noble_excluded.contains(idx))
        .collect();
    sorted_indices.sort();

    let mut peptide_counts: Vec<usize> = Vec::new();

    for &idx in &sorted_indices {
        let (header_start, header_end, seq_start, seq_end) = regions[idx];
        let seq = extract_sequence(data, seq_start, seq_end);

        let hits = matched_indices.get(&idx).unwrap();
        let num_peptides = hits.len();
        peptide_counts.push(num_peptides);

        // C2: write the original header BYTES verbatim (no lossy UTF-8 round-trip
        // that would collapse a non-UTF-8 header to the literal ">unknown").
        writer.write_all(b">")?;
        writer.write_all(&data[header_start..header_end])?;
        writer.write_all(b"\n")?;

        // Write sequence in 60-char lines
        for chunk in seq.chunks(60) {
            writer.write_all(chunk)?;
            writer.write_all(b"\n")?;
        }
    }

    // Step 6b: Append included (core proteome) proteins not already matched
    let mut n_included = 0usize;
    let mut n_included_already = 0usize;

    if let Some(ref include_path) = args.include {
        if include_path.exists() {
            println!("  Appending core proteome from {:?}...", include_path);

            // Collect IDs already in the output
            let mut existing_ids: HashSet<String> = HashSet::new();
            for &idx in &sorted_indices {
                let (header_start, header_end, _, _) = regions[idx];
                let header = std::str::from_utf8(&data[header_start..header_end]).unwrap_or("?");
                let protein_id = header.split_whitespace().next().unwrap_or(header);
                existing_ids.insert(protein_id.to_string());
            }

            // Read include FASTA and append new proteins
            let include_file = File::open(include_path)?;
            // SAFETY: the include FASTA is opened read-only and is not mutated or
            // truncated by any other process for the duration of this run.
            let include_mmap = unsafe { MmapOptions::new().map(&include_file)? };
            let include_data = &include_mmap[..];
            let include_regions = parse_fasta_regions(include_data);

            for &(header_start, header_end, seq_start, seq_end) in &include_regions {
                let header =
                    std::str::from_utf8(&include_data[header_start..header_end]).unwrap_or("?");
                let protein_id = header.split_whitespace().next().unwrap_or(header);

                if existing_ids.contains(protein_id) {
                    n_included_already += 1;
                    continue;
                }

                // Write this protein (not already in output).
                // C2: emit the header BYTES verbatim; the lossy `header` above is used
                // only to derive protein_id for dedup logic, not for output.
                let seq = extract_sequence(include_data, seq_start, seq_end);
                writer.write_all(b">")?;
                writer.write_all(&include_data[header_start..header_end])?;
                writer.write_all(b"\n")?;
                for chunk in seq.chunks(60) {
                    writer.write_all(chunk)?;
                    writer.write_all(b"\n")?;
                }
                existing_ids.insert(protein_id.to_string());
                n_included += 1;
            }

            println!(
                "  Core proteome: {} new proteins added, {} already present",
                n_included, n_included_already
            );
        } else {
            eprintln!("WARNING: Include FASTA not found: {:?}", include_path);
        }
    }

    writer.flush()?;
    let total_output = total_matched + n_included;
    println!(
        "  Written {} proteins total ({} matched + {} core) in {:.2?}",
        total_output,
        total_matched,
        n_included,
        write_start.elapsed()
    );

    // Optionally write hits TSV
    if let Some(hits_path) = args.output_hits {
        println!("\nStep 7: Writing peptide hits TSV...");
        let hits_file = File::create(&hits_path)?;
        let mut hits_writer = BufWriter::new(hits_file);

        writeln!(hits_writer, "peptide\tprotein_id\tposition")?;

        for &idx in &sorted_indices {
            let (header_start, header_end, _, _) = regions[idx];
            let header = std::str::from_utf8(&data[header_start..header_end]).unwrap_or("unknown");
            let protein_id = header.split_whitespace().next().unwrap_or(header);

            let hits = matched_indices.get(&idx).unwrap();
            // D7: hit order depends on the AhoCorasick pattern order, which is derived
            // from a HashSet and thus nondeterministic across runs. Sort by a total key
            // (position, peptide) so the hits TSV is byte-reproducible. protein_id is
            // constant within this loop, so it need not be part of the sort key.
            let mut rows: Vec<&(String, usize)> = hits.iter().collect();
            rows.sort_by(|a, b| a.1.cmp(&b.1).then_with(|| a.0.cmp(&b.0)));
            for (peptide, pos) in rows {
                writeln!(hits_writer, "{}\t{}\t{}", peptide, protein_id, pos)?;
            }
        }
        hits_writer.flush()?;
    }

    // Optionally write evidence manifest TSV
    if let Some(ref evidence_path) = args.output_evidence {
        println!("\nStep 7b: Writing evidence manifest...");
        let ev_file = File::create(evidence_path)?;
        let mut ev_writer = BufWriter::new(ev_file);

        writeln!(ev_writer, "protein_id\tn_peptides\tn_unique_peptides\tn_samples\tmean_score\tmax_score\tsupport_score")?;

        let total_samples = collection.total_samples.max(1) as f64;

        for &idx in &sorted_indices {
            let (header_start, header_end, _, _) = regions[idx];
            let header = std::str::from_utf8(&data[header_start..header_end]).unwrap_or("unknown");
            let protein_id = header.split_whitespace().next().unwrap_or(header);

            let hits = matched_indices.get(&idx).unwrap();
            let n_peptide_hits = hits.len();

            // Unique peptides and their scores
            let unique_peps: HashSet<&str> = hits.iter().map(|(p, _)| p.as_str()).collect();
            let n_unique = unique_peps.len();

            let mut score_sum = 0.0f64;
            let mut max_score = 0.0f64;
            let mut sample_ids: HashSet<usize> = HashSet::new();

            let mut ordered_peps: Vec<_> = unique_peps.iter().copied().collect();
            ordered_peps.sort_unstable();
            for pep in &ordered_peps {
                let s = collection.best_scores.get(*pep).copied().unwrap_or(0.5);
                score_sum += s;
                if s > max_score {
                    max_score = s;
                }
                if let Some(ids) = collection.sample_ids.get(*pep) {
                    sample_ids.extend(ids.iter().copied());
                }
            }

            let sample_union = sample_ids.len();
            let mean_score = if n_unique > 0 {
                score_sum / n_unique as f64
            } else {
                0.0
            };

            // Support score: composite of peptide count, sample breadth, and de novo confidence
            // log2(n_peptides + 1) * sqrt(n_samples / total_samples) * mean_score
            let support = (((n_unique as f64) + 1.0).log2())
                * ((sample_union as f64 / total_samples).sqrt())
                * mean_score;
            let support_score = support.min(1.0);

            writeln!(
                ev_writer,
                "{}\t{}\t{}\t{}\t{:.4}\t{:.4}\t{:.4}",
                protein_id,
                n_peptide_hits,
                n_unique,
                sample_union,
                mean_score,
                max_score,
                support_score
            )?;
        }

        ev_writer.flush()?;
        println!(
            "  Written evidence for {} proteins to {:?}",
            sorted_indices.len(),
            evidence_path
        );
    }

    // Summary statistics
    let total_time = start_time.elapsed();
    let avg_peptides: f64 = if !peptide_counts.is_empty() {
        peptide_counts.iter().sum::<usize>() as f64 / peptide_counts.len() as f64
    } else {
        0.0
    };

    println!("\n=======================================================");
    println!("EXTRACTION COMPLETE");
    println!("=======================================================");
    println!("Total time: {:.2?}", total_time);
    println!("Proteins searched: {}", total_proteins);
    println!(
        "Proteins matched: {} ({:.2}%)",
        total_matched,
        100.0 * total_matched as f64 / total_proteins as f64
    );
    println!("Avg peptides per protein: {:.1}", avg_peptides);
    println!("Output: {:?}", args.output);

    // Write stats JSON
    let stats_path = args.output.with_extension("stats.json");
    let mut stats = serde_json::json!({
        "database": args.database.to_string_lossy(),
        "source_tag": args.source_tag,
        "peptide_length_range": [args.min_length, args.max_length],
        "score_filter": if let Some(tp) = args.top_percent {
            format!("top_{:.0}%_per_sample", tp * 100.0)
        } else {
            format!("min_score_{}", args.min_score)
        },
        "normalize_il": args.normalize_il,
        "total_peptides": peptides.len(),
        "proteins_searched": total_proteins,
        "proteins_matched": total_matched,
        "core_proteome_added": n_included,
        "total_output_proteins": total_matched + n_included,
        "match_rate": total_matched as f64 / total_proteins as f64,
        "avg_peptides_per_protein": avg_peptides,
        "total_seconds": total_time.as_secs_f64()
    });

    // Extraction-stage match diagnostics are not identification-stage FDP estimates.
    if do_entrapment {
        let total_hits = total_matched + entrapment_matched;
        let fraction = if total_hits > 0 {
            Some(entrapment_matched as f64 / total_hits as f64)
        } else {
            None
        };
        stats["entrapment"] = serde_json::json!({
            "method": entrapment_method_name,
            "total_proteins": entrapment_total,
            "matched_proteins": entrapment_matched,
            "target_matched": total_matched,
            "fraction_of_matches": fraction,
            "certification": "NOT_APPLICABLE",
            "interpretation": "Extraction-stage negative-control match fraction; not an identification FDR estimate or certificate.",
            "hits": entrapment_hits_detail.iter().map(|(pid, n)| {
                serde_json::json!({"protein_id": pid, "peptide_hits": n})
            }).collect::<Vec<_>>()
        });
    }

    let stats_file = File::create(&stats_path)?;
    serde_json::to_writer_pretty(stats_file, &stats)?;
    println!("Stats: {:?}", stats_path);

    Ok(())
}

#[cfg(test)]
mod release_review_tests {
    use super::*;
    #[test]
    fn sample_breadth_is_union_not_max_and_duplicate_rows_do_not_inflate_it() {
        let mut c = PeptideCollection::new();
        c.add("PEPTIDEAK".into(), 0.9, true);
        c.add("PEPTIDEAK".into(), 0.8, true);
        assert_eq!(c.sample_counts["PEPTIDEAK"], 1);
        c.total_samples = 1;
        c.add("OTHERPEPK".into(), 0.7, true);
        let union: HashSet<usize> = ["PEPTIDEAK", "OTHERPEPK"]
            .iter()
            .flat_map(|p| c.sample_ids[*p].iter().copied())
            .collect();
        assert_eq!(union.len(), 2);
    }
}
