//! Mass k-mer Anchoring: Rapid Validation of De Novo Peptide Sequences Against Protein Databases
//!
//! Rust implementation of Lebedev, Treit, Mann (2026).
//!
//! Given a FASTA database and AlphaNovo prediction CSVs, finds peptide-protein matches
//! using mass-based 5-mer windows instead of exact sequence matching. Tolerates adjacent
//! AA swaps, I/L substitutions, and other mass-preserving de novo errors.
//!
//! Every candidate is verified against a precursor mass at +/-20 ppm, trying +/-2 length
//! offsets (N <-> GG type substitutions), as in the reference implementation. The mass is
//! the measured precursor mass when the CSV carries `prec_mz` and `prec_charge`; otherwise
//! it is the theoretical monoisotopic mass of the de novo sequence itself, which the
//! reference also falls back to. The `mass_source` output column records which was used.
//!
//! Usage:
//!   mass_kmer_anchor -d <database.fasta> -p <predictions_dir> -o <output.tsv> [-t threads]

// Legacy/experimental reference binary — not on the canonical pipeline; clippy relaxed.
#![allow(clippy::all)]
#![allow(dead_code, unused)]

use clap::Parser;
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::PathBuf;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Instant;

/// Monoisotopic amino acid residue masses (Da)
fn aa_mass(aa: u8) -> f64 {
    match aa {
        b'G' => 57.02146,
        b'A' => 71.03711,
        b'V' => 99.06841,
        b'L' | b'I' => 113.08406, // Isobaric
        b'P' => 97.05276,
        b'F' => 147.06841,
        b'W' => 186.07931,
        b'M' => 131.04049,
        b'S' => 87.03203,
        b'T' => 101.04768,
        b'C' => 103.00919,
        b'Y' => 163.06333,
        b'H' => 137.05891,
        b'D' => 115.02694,
        b'E' => 129.04259,
        b'N' => 114.04293,
        b'Q' => 128.05858,
        b'K' => 128.09496,
        b'R' => 156.10111,
        _ => 0.0,
    }
}

const K: usize = 5;
const BIN_WIDTH: f64 = 0.001;
const WATER_MASS: f64 = 18.010565;
const PROTON_MASS: f64 = 1.007276;
const MASS_TOLERANCE_PPM: f64 = 20.0;

#[inline]
fn mass_to_bin(mass: f64) -> i64 {
    (mass / BIN_WIDTH).round() as i64
}

/// A protein in the database
struct Protein {
    id: String,
    sequence: Vec<u8>,
    cumulative_mass: Vec<f64>,
}

impl Protein {
    fn new(id: String, sequence: Vec<u8>) -> Self {
        let mut cum = Vec::with_capacity(sequence.len() + 1);
        cum.push(0.0);
        let mut running = 0.0;
        for &aa in &sequence {
            running += aa_mass(aa);
            cum.push(running);
        }
        Protein {
            id,
            sequence,
            cumulative_mass: cum,
        }
    }

    fn kmer_mass(&self, pos: usize) -> f64 {
        self.cumulative_mass[pos + K] - self.cumulative_mass[pos]
    }

    fn peptide_mass(&self, start: usize, end: usize) -> f64 {
        self.cumulative_mass[end] - self.cumulative_mass[start] + WATER_MASS
    }
}

/// Neutral monoisotopic mass of a peptide sequence (residues + water).
fn theoretical_mass(seq: &[u8]) -> f64 {
    seq.iter().map(|&aa| aa_mass(aa)).sum::<f64>() + WATER_MASS
}

/// A de novo prediction to validate
struct Prediction {
    sequence: Vec<u8>,
    /// Mass every candidate is verified against (never 0).
    precursor_mass: f64,
    /// "measured" (from prec_mz/prec_charge) or "theoretical" (from the sequence).
    mass_source: &'static str,
    score: f64,
}

impl Prediction {
    /// `measured` is the instrument precursor mass if the input carried one.
    fn new(sequence: Vec<u8>, measured: Option<f64>, score: f64) -> Self {
        let (precursor_mass, mass_source) = match measured {
            Some(m) if m > 0.0 => (m, "measured"),
            _ => (theoretical_mass(&sequence), "theoretical"),
        };
        Prediction {
            sequence,
            precursor_mass,
            mass_source,
            score,
        }
    }
}

/// A match result
struct AnchorHit {
    query_seq: String,
    db_seq: String,
    protein_id: String,
    position: usize,
    votes: usize,
    n_kmers: usize,
    ppm_error: f64,
    length_offset: i32,
    score: f64,
    mass_source: &'static str,
}

/// Mass k-mer index: bin → Vec<(protein_idx, position)>
struct MassKmerIndex {
    proteins: Vec<Protein>,
    bins: HashMap<i64, Vec<(u32, u32)>>, // (protein_idx as u32, position as u32)
}

impl MassKmerIndex {
    fn build(proteins: Vec<Protein>) -> Self {
        let t0 = Instant::now();
        let mut bins: HashMap<i64, Vec<(u32, u32)>> = HashMap::new();
        let mut total_kmers: usize = 0;

        for (pidx, prot) in proteins.iter().enumerate() {
            if prot.sequence.len() < K {
                continue;
            }
            for pos in 0..=(prot.sequence.len() - K) {
                let mass = prot.kmer_mass(pos);
                let bin_id = mass_to_bin(mass);
                bins.entry(bin_id)
                    .or_default()
                    .push((pidx as u32, pos as u32));
                total_kmers += 1;
            }

            if (pidx + 1) % 500_000 == 0 {
                eprintln!(
                    "  Indexed {}/{} proteins ({} k-mers)...",
                    pidx + 1,
                    proteins.len(),
                    total_kmers
                );
            }
        }

        let elapsed = t0.elapsed();
        eprintln!(
            "Index built: {} proteins, {} k-mers, {} bins, {:.1}s",
            proteins.len(),
            total_kmers,
            bins.len(),
            elapsed.as_secs_f64()
        );

        MassKmerIndex { proteins, bins }
    }

    fn query(&self, pred: &Prediction, min_votes: usize) -> Vec<AnchorHit> {
        let seq = &pred.sequence;
        if seq.len() < K {
            return vec![];
        }
        let n_kmers = seq.len() - K + 1;

        // Compute query k-mer masses
        let mut query_cum = Vec::with_capacity(seq.len() + 1);
        query_cum.push(0.0);
        let mut running = 0.0;
        for &aa in seq {
            running += aa_mass(aa);
            query_cum.push(running);
        }

        // Vote: (protein_idx, inferred_start) → count
        let mut votes: HashMap<(u32, i32), usize> = HashMap::new();

        for kmer_idx in 0..n_kmers {
            let mass = query_cum[kmer_idx + K] - query_cum[kmer_idx];
            let bin_id = mass_to_bin(mass);

            for b in [bin_id - 1, bin_id, bin_id + 1] {
                if let Some(entries) = self.bins.get(&b) {
                    for &(pidx, pos) in entries {
                        let inferred_start = pos as i32 - kmer_idx as i32;
                        if inferred_start >= 0 {
                            *votes.entry((pidx, inferred_start)).or_insert(0) += 1;
                        }
                    }
                }
            }
        }

        // Filter and verify
        let mut hits = Vec::new();

        for (&(pidx, start), &vote_count) in &votes {
            if vote_count < min_votes {
                continue;
            }
            let prot = &self.proteins[pidx as usize];
            let start = start as usize;

            // Try ±2 length offsets
            for offset in -2i32..=2 {
                let pep_len = seq.len() as i32 + offset;
                if pep_len < K as i32 {
                    continue;
                }
                let end = start + pep_len as usize;
                if end > prot.sequence.len() {
                    continue;
                }

                let candidate_mass = prot.peptide_mass(start, end);

                // The mass check is unconditional (reference behaviour): a candidate
                // window is a hit only if its mass agrees with the precursor mass,
                // measured or theoretical, within MASS_TOLERANCE_PPM.
                let ppm_error =
                    ((candidate_mass - pred.precursor_mass).abs() / pred.precursor_mass) * 1e6;
                if ppm_error <= MASS_TOLERANCE_PPM {
                    let db_seq = String::from_utf8_lossy(&prot.sequence[start..end]).to_string();
                    hits.push(AnchorHit {
                        query_seq: String::from_utf8_lossy(seq).to_string(),
                        db_seq,
                        protein_id: prot.id.clone(),
                        position: start,
                        votes: vote_count,
                        n_kmers,
                        ppm_error,
                        length_offset: offset,
                        score: pred.score,
                        mass_source: pred.mass_source,
                    });
                }
            }
        }

        // Determinism: the `votes` HashMap iterates in a nondeterministic order, so
        // the order is a total key: votes (desc), mass agreement at the printed
        // precision (0.01 ppm, so two identical windows in different proteins do
        // not rank by prefix-sum rounding noise), then protein_id, position and
        // length_offset. The selected rescue is reproducible bit-for-bit and can
        // be reconstructed from the written TSV.
        hits.sort_by(|a, b| {
            let ppm_key = |h: &AnchorHit| (h.ppm_error * 100.0).round() as i64;
            b.votes
                .cmp(&a.votes)
                .then_with(|| ppm_key(a).cmp(&ppm_key(b)))
                .then_with(|| a.protein_id.cmp(&b.protein_id))
                .then_with(|| a.position.cmp(&b.position))
                .then_with(|| a.length_offset.cmp(&b.length_offset))
        });
        hits
    }
}

fn load_fasta(path: &str) -> Vec<Protein> {
    let t0 = Instant::now();
    let file = File::open(path).unwrap_or_else(|e| panic!("Cannot open {}: {}", path, e));
    let reader = BufReader::with_capacity(8 * 1024 * 1024, file);

    let mut proteins = Vec::new();
    let mut current_id = String::new();
    let mut current_seq = Vec::new();

    for line in reader.lines() {
        let line = line.expect("Failed to read line");
        if line.starts_with('>') {
            if !current_id.is_empty() && !current_seq.is_empty() {
                proteins.push(Protein::new(
                    std::mem::take(&mut current_id),
                    std::mem::take(&mut current_seq),
                ));
            }
            current_id = line[1..]
                .split_whitespace()
                .next()
                .unwrap_or("")
                .to_string();
            current_seq.clear();
        } else {
            current_seq.extend(
                line.trim()
                    .as_bytes()
                    .iter()
                    .map(|b| b.to_ascii_uppercase()),
            );
        }
    }
    if !current_id.is_empty() && !current_seq.is_empty() {
        proteins.push(Protein::new(current_id, current_seq));
    }

    eprintln!(
        "Loaded {} proteins from {} in {:.1}s",
        proteins.len(),
        path,
        t0.elapsed().as_secs_f64()
    );
    proteins
}

fn load_predictions(csv_dir: &str) -> Vec<Prediction> {
    let t0 = Instant::now();
    let mut predictions = Vec::new();
    let mut best_scores: HashMap<Vec<u8>, usize> = HashMap::new(); // seq → index in predictions

    // Try *predictions.csv first, then fall back to *.csv
    let pattern = format!("{}/**/*predictions.csv", csv_dir);
    let mut paths: Vec<_> = glob::glob(&pattern)
        .expect("Invalid glob pattern")
        .filter_map(|p| p.ok())
        .collect();

    if paths.is_empty() {
        let pattern2 = format!("{}/**/*.csv", csv_dir);
        paths = glob::glob(&pattern2)
            .expect("Invalid glob pattern")
            .filter_map(|p| p.ok())
            .collect();
    }

    if paths.is_empty() {
        // Try direct file
        let path = PathBuf::from(csv_dir);
        if path.is_file() {
            return load_predictions_from_file(&path);
        }
        eprintln!("No prediction CSVs found in {}", csv_dir);
        return predictions;
    }

    for path in &paths {
        let file_preds = load_predictions_from_file(path);
        for pred in file_preds {
            if let Some(&existing_idx) = best_scores.get(&pred.sequence) {
                if pred.score > predictions[existing_idx].score {
                    predictions[existing_idx] = pred;
                }
            } else {
                let idx = predictions.len();
                best_scores.insert(pred.sequence.clone(), idx);
                predictions.push(pred);
            }
        }
    }

    eprintln!(
        "Loaded {} unique predictions from {} files in {:.1}s",
        predictions.len(),
        paths.len(),
        t0.elapsed().as_secs_f64()
    );
    predictions
}

fn load_predictions_from_file(path: &PathBuf) -> Vec<Prediction> {
    let mut predictions = Vec::new();
    let file = match File::open(path) {
        Ok(f) => f,
        Err(e) => {
            eprintln!("Warning: cannot open {}: {}", path.display(), e);
            return predictions;
        }
    };
    let mut reader = csv::Reader::from_reader(BufReader::new(file));

    // Find column indices from headers
    let headers = match reader.headers() {
        Ok(h) => h.clone(),
        Err(_) => return predictions,
    };

    let seq_idx = headers
        .iter()
        .position(|c| c == "peptide_prediction_detokenized_unmodified")
        .or_else(|| headers.iter().position(|c| c == "pred_seq"))
        .or_else(|| headers.iter().position(|c| c == "sequence"));
    let score_idx = headers
        .iter()
        .position(|c| c == "score")
        .or_else(|| headers.iter().position(|c| c == "pred_prob"));
    let mz_idx = headers.iter().position(|c| c == "prec_mz");
    let charge_idx = headers.iter().position(|c| c == "prec_charge");

    let seq_idx = match seq_idx {
        Some(i) => i,
        None => {
            eprintln!("Warning: no sequence column in {}", path.display());
            return predictions;
        }
    };

    for result in reader.records() {
        let record = match result {
            Ok(r) => r,
            Err(_) => continue,
        };

        let seq_str = record.get(seq_idx).unwrap_or("").trim().to_uppercase();
        let seq: Vec<u8> = seq_str.bytes().collect();

        if seq.len() < 8 || seq.len() > 50 {
            continue;
        }
        if !seq.iter().all(|&b| aa_mass(b) > 0.0) {
            continue;
        }
        // Low-complexity filter: skip if any single AA is >50% of sequence
        {
            let mut counts = [0u32; 26];
            for &b in &seq {
                if b >= b'A' && b <= b'Z' {
                    counts[(b - b'A') as usize] += 1;
                }
            }
            let max_count = *counts.iter().max().unwrap_or(&0) as usize;
            if max_count * 2 > seq.len() {
                continue;
            }
        }

        let score: f64 = score_idx
            .and_then(|i| record.get(i))
            .and_then(|s| s.parse().ok())
            .unwrap_or(0.0);

        let mz: f64 = mz_idx
            .and_then(|i| record.get(i))
            .and_then(|s| s.parse().ok())
            .unwrap_or(0.0);

        let charge: f64 = charge_idx
            .and_then(|i| record.get(i))
            .and_then(|s| s.parse().ok())
            .unwrap_or(0.0);

        let measured = if mz > 0.0 && charge > 0.0 {
            Some(mz * charge - charge * PROTON_MASS)
        } else {
            None
        };

        predictions.push(Prediction::new(seq, measured, score));
    }

    predictions
}

/// Exact substring matching with I/L normalization
fn exact_match(peptides: &[Vec<u8>], proteins: &[Protein]) -> Vec<bool> {
    let t0 = Instant::now();
    let n = peptides.len();
    let matched = AtomicUsize::new(0);

    // Normalize peptides: I→L
    let norm_peptides: Vec<Vec<u8>> = peptides
        .iter()
        .map(|p| {
            p.iter()
                .map(|&b| if b == b'I' { b'L' } else { b })
                .collect()
        })
        .collect();

    // Build Aho-Corasick automaton
    let ac = aho_corasick::AhoCorasick::builder()
        .build(&norm_peptides)
        .expect("Failed to build Aho-Corasick");

    let results: Vec<std::sync::atomic::AtomicBool> = (0..n)
        .map(|_| std::sync::atomic::AtomicBool::new(false))
        .collect();

    proteins.par_iter().for_each(|prot| {
        let norm_seq: Vec<u8> = prot
            .sequence
            .iter()
            .map(|&b| if b == b'I' { b'L' } else { b })
            .collect();

        for mat in ac.find_overlapping_iter(&norm_seq) {
            let idx = mat.pattern().as_usize();
            results[idx].store(true, Ordering::Relaxed);
            matched.fetch_add(1, Ordering::Relaxed);
        }
    });

    let results: Vec<bool> = results.iter().map(|b| b.load(Ordering::Relaxed)).collect();
    let m = results.iter().filter(|&&b| b).count();
    eprintln!(
        "Exact match: {}/{} peptides matched in {:.1}s",
        m,
        n,
        t0.elapsed().as_secs_f64()
    );
    results
}

const VERSION: &str = concat!(env!("CARGO_PKG_VERSION"), "+", env!("GIT_REV"));

#[derive(Parser)]
#[command(name = "mass_kmer_anchor")]
#[command(version = VERSION)]
#[command(about = "Mass k-mer anchoring: validate de novo peptides against protein databases")]
struct Args {
    /// Path to protein FASTA database
    #[arg(short = 'd', long)]
    database: String,

    /// Path to predictions CSV file or directory
    #[arg(short = 'p', long)]
    predictions: String,

    /// Output TSV file
    #[arg(short = 'o', long)]
    output: String,

    /// Minimum votes for a match
    #[arg(long, default_value = "3")]
    min_votes: usize,

    /// Number of threads
    #[arg(short = 't', long, default_value = "0")]
    threads: usize,

    /// Skip exact matching (run mass k-mer on all peptides)
    #[arg(long)]
    skip_exact: bool,

    /// Output stats JSON
    #[arg(long)]
    output_stats: Option<String>,
}

fn main() {
    let args = Args::parse();

    if args.threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global()
            .ok();
    }

    let total_t0 = Instant::now();

    // Load database
    let proteins = load_fasta(&args.database);

    // Load predictions
    let predictions = load_predictions(&args.predictions);
    if predictions.is_empty() {
        eprintln!("No predictions loaded, exiting.");
        return;
    }

    let peptide_seqs: Vec<Vec<u8>> = predictions.iter().map(|p| p.sequence.clone()).collect();

    // Exact matching
    let failed_indices: Vec<usize>;
    let n_exact;

    if args.skip_exact {
        failed_indices = (0..predictions.len()).collect();
        n_exact = 0;
        eprintln!(
            "Skipping exact matching, querying all {} peptides",
            predictions.len()
        );
    } else {
        let matched = exact_match(&peptide_seqs, &proteins);
        n_exact = matched.iter().filter(|&&b| b).count();
        failed_indices = matched
            .iter()
            .enumerate()
            .filter(|&(_, b)| !b)
            .map(|(i, _)| i)
            .collect();
        eprintln!(
            "Exact: {} matched, {} failed → will query with mass k-mer",
            n_exact,
            failed_indices.len()
        );
    }

    // Build mass k-mer index
    let index = MassKmerIndex::build(proteins);

    // Query failed predictions
    eprintln!(
        "Querying {} peptides with mass k-mer anchoring...",
        failed_indices.len()
    );
    let query_t0 = Instant::now();
    let queried = AtomicUsize::new(0);
    let rescued_count = AtomicUsize::new(0);

    let all_hits: Vec<Option<AnchorHit>> = failed_indices
        .par_iter()
        .map(|&idx| {
            let pred = &predictions[idx];
            let hits = index.query(pred, args.min_votes);
            let q = queried.fetch_add(1, Ordering::Relaxed) + 1;
            if q % 10000 == 0 {
                eprintln!(
                    "  Queried {}/{} ({} rescued)",
                    q,
                    failed_indices.len(),
                    rescued_count.load(Ordering::Relaxed)
                );
            }
            if let Some(best) = hits.into_iter().next() {
                rescued_count.fetch_add(1, Ordering::Relaxed);
                Some(best)
            } else {
                None
            }
        })
        .collect();

    let query_time = query_t0.elapsed();
    let hits: Vec<&AnchorHit> = all_hits.iter().filter_map(|h| h.as_ref()).collect();

    eprintln!(
        "Mass k-mer: rescued {} / {} failed peptides in {:.1}s ({:.0} pep/s)",
        hits.len(),
        failed_indices.len(),
        query_time.as_secs_f64(),
        failed_indices.len() as f64 / query_time.as_secs_f64()
    );

    // Write output
    let out_file = File::create(&args.output).expect("Cannot create output file");
    let mut writer = BufWriter::new(out_file);
    writeln!(
        writer,
        "query_sequence\tdb_sequence\tprotein_id\tposition\tvotes\tn_kmers\tppm_error\tlength_offset\tscore\tmass_source"
    )
    .unwrap();

    for hit in &hits {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{:.2}\t{}\t{:.4}\t{}",
            hit.query_seq,
            hit.db_seq,
            hit.protein_id,
            hit.position,
            hit.votes,
            hit.n_kmers,
            hit.ppm_error,
            hit.length_offset,
            hit.score,
            hit.mass_source,
        )
        .unwrap();
    }

    let total_time = total_t0.elapsed();
    eprintln!("\nSUMMARY:");
    eprintln!("  Total peptides:     {}", predictions.len());
    eprintln!("  Exact matches:      {}", n_exact);
    eprintln!("  Failed:             {}", failed_indices.len());
    eprintln!("  Rescued (k-mer):    {}", hits.len());
    let unique_proteins: std::collections::HashSet<&str> =
        hits.iter().map(|h| h.protein_id.as_str()).collect();
    eprintln!("  Rescued proteins:   {}", unique_proteins.len());
    eprintln!("  Output:             {}", args.output);
    eprintln!("  Total time:         {:.1}s", total_time.as_secs_f64());

    // Stats JSON
    if let Some(stats_path) = &args.output_stats {
        let stats = serde_json::json!({
            "total_peptides": predictions.len(),
            "exact_matches": n_exact,
            "failed": failed_indices.len(),
            "rescued_peptides": hits.len(),
            "rescued_proteins": unique_proteins.len(),
            "query_rate": failed_indices.len() as f64 / query_time.as_secs_f64(),
            "total_time_s": total_time.as_secs_f64(),
        });
        std::fs::write(stats_path, serde_json::to_string_pretty(&stats).unwrap())
            .expect("Cannot write stats");
    }
}

// ============================================================================
// TESTS
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_aa_mass_known_values() {
        assert!((aa_mass(b'G') - 57.02146).abs() < 1e-5);
        assert!((aa_mass(b'A') - 71.03711).abs() < 1e-5);
        assert!((aa_mass(b'L') - 113.08406).abs() < 1e-5);
        // I and L are isobaric
        assert_eq!(aa_mass(b'I'), aa_mass(b'L'));
    }

    #[test]
    fn test_il_isobaric() {
        // Core feature: I/L substitution preserves mass exactly
        assert_eq!(aa_mass(b'I'), aa_mass(b'L'));
    }

    #[test]
    fn test_mass_to_bin_roundtrip() {
        let mass = 500.12345;
        let bin = mass_to_bin(mass);
        let recovered = bin as f64 * BIN_WIDTH;
        assert!((mass - recovered).abs() < BIN_WIDTH);
    }

    #[test]
    fn test_protein_cumulative_mass() {
        let seq = b"GACL".to_vec();
        let prot = Protein::new("test".into(), seq);
        // cumulative[0] = 0
        // cumulative[1] = G = 57.02
        // cumulative[2] = G+A = 57.02 + 71.04 = 128.06
        assert_eq!(prot.cumulative_mass[0], 0.0);
        assert!((prot.cumulative_mass[1] - 57.02146).abs() < 1e-5);
        assert!((prot.cumulative_mass[2] - 128.05857).abs() < 1e-4);
    }

    #[test]
    fn test_protein_kmer_mass() {
        // 5-mer: GACLV
        let seq = b"GACLVPF".to_vec();
        let prot = Protein::new("test".into(), seq);
        let kmer_mass = prot.kmer_mass(0); // G+A+C+L+V
        let expected = 57.02146 + 71.03711 + 103.00919 + 113.08406 + 99.06841;
        assert!((kmer_mass - expected).abs() < 1e-4);
    }

    #[test]
    fn test_index_build_small() {
        let proteins = vec![
            Protein::new("P1".into(), b"MKWVTFISLLFLFSSAYSRGVFRRDAHK".to_vec()),
            Protein::new("P2".into(), b"ACDEFGHIKLMNPQRSTVWY".to_vec()),
        ];
        let index = MassKmerIndex::build(proteins);
        assert_eq!(index.proteins.len(), 2);
        assert!(!index.bins.is_empty());
    }

    #[test]
    fn test_query_exact_match() {
        // A peptide that exactly matches a protein segment should get max votes
        let proteins = vec![Protein::new(
            "P1".into(),
            b"MKWVTFISLLFLFSSAYSRGVFRRDAHK".to_vec(),
        )];
        let index = MassKmerIndex::build(proteins);

        // Query with exact substring (10 residues from position 5)
        let sub = b"ISLLFLFSSA".to_vec();
        let pred = Prediction::new(sub, None, 0.95);
        let hits = index.query(&pred, 3);
        assert!(!hits.is_empty(), "Exact substring should produce hits");
        assert!(hits[0].votes >= 3, "Should have at least 3 votes");
    }

    #[test]
    fn test_query_il_swap_rescued() {
        // A peptide with I→L swap should still match (same mass)
        let proteins = vec![Protein::new(
            "P1".into(),
            b"MKWVTFISLLFLFSSAYSRGVFRRDAHK".to_vec(),
        )];
        let index = MassKmerIndex::build(proteins);

        // Original substring: ISLLFLFSSA (positions 5-14)
        // Swap I→L at position 0: LSLLFLFSSA
        let swapped = b"LSLLFLFSSA".to_vec();
        let pred = Prediction::new(swapped, None, 0.9);
        let hits = index.query(&pred, 3);
        assert!(
            !hits.is_empty(),
            "I/L swap should still produce hits (isobaric)"
        );
    }

    #[test]
    fn test_query_no_match_random() {
        let proteins = vec![Protein::new(
            "P1".into(),
            b"MKWVTFISLLFLFSSAYSRGVFRRDAHK".to_vec(),
        )];
        let index = MassKmerIndex::build(proteins);

        // Random peptide with no relation to the protein
        let random = b"GGGGGGGGGG".to_vec();
        let pred = Prediction::new(random, None, 0.5);
        let hits = index.query(&pred, 5);
        assert!(hits.is_empty(), "Random poly-G should not match");
    }

    #[test]
    fn test_query_short_peptide_skipped() {
        let proteins = vec![Protein::new(
            "P1".into(),
            b"MKWVTFISLLFLFSSAYSRGVFRRDAHK".to_vec(),
        )];
        let index = MassKmerIndex::build(proteins);

        // Peptide shorter than K=5
        let short = b"MKWV".to_vec();
        let pred = Prediction::new(short, None, 0.9);
        let hits = index.query(&pred, 1);
        assert!(
            hits.is_empty(),
            "Peptide shorter than K should return no hits"
        );
    }

    fn small_index() -> MassKmerIndex {
        // ...SRGG... so that a de novo N (isobaric with GG) can be tested.
        MassKmerIndex::build(vec![Protein::new(
            "P1".into(),
            b"MKWVTFISLLFLFSSAYSRGGVFRRDAHKQE".to_vec(),
        )])
    }

    #[test]
    fn theoretical_mass_is_the_fallback_and_is_recorded() {
        let pred = Prediction::new(b"ISLLFLFSSA".to_vec(), None, 0.9);
        assert_eq!(pred.mass_source, "theoretical");
        assert!((pred.precursor_mass - theoretical_mass(b"ISLLFLFSSA")).abs() < 1e-9);
        let measured = Prediction::new(b"ISLLFLFSSA".to_vec(), Some(1234.5), 0.9);
        assert_eq!(measured.mass_source, "measured");
        assert_eq!(measured.precursor_mass, 1234.5);
    }

    #[test]
    fn exact_window_has_zero_ppm_and_zero_offset() {
        let index = small_index();
        let hits = index.query(&Prediction::new(b"ISLLFLFSSAYSR".to_vec(), None, 0.9), 3);
        assert!(!hits.is_empty());
        assert_eq!(hits[0].length_offset, 0);
        assert_eq!(hits[0].db_seq, "ISLLFLFSSAYSR");
        assert!(hits[0].ppm_error < 0.01, "ppm={}", hits[0].ppm_error);
        assert_eq!(hits[0].mass_source, "theoretical");
    }

    #[test]
    fn n_for_gg_is_accepted_with_length_offset() {
        // Query carries N where the protein has GG (isobaric, one residue shorter).
        // The substitution sits in the second half so the clean windows before it
        // carry the vote and anchor the true start; windows after it would vote a
        // start shifted by one, exactly as in the reference implementation.
        let index = small_index();
        let hits = index.query(
            &Prediction::new(b"ISLLFLFSSAYSRNVFR".to_vec(), None, 0.9),
            3,
        );
        assert!(!hits.is_empty(), "N<->GG must be rescued");
        assert_eq!(hits[0].length_offset, 1);
        assert_eq!(hits[0].db_seq, "ISLLFLFSSAYSRGGVFR");
        assert!(hits[0].ppm_error <= MASS_TOLERANCE_PPM);
    }

    #[test]
    fn mass_changing_substitution_is_rejected_without_measured_mass() {
        // Last residue V -> A: six clean 5-mers still vote, but the mass is 28 Da off.
        let index = small_index();
        let hits = index.query(&Prediction::new(b"ISLLFLFSSAYSRGA".to_vec(), None, 0.9), 3);
        assert!(
            hits.is_empty(),
            "a candidate failing the ppm check must not be reported; got {:?}",
            hits.iter()
                .map(|h| (&h.db_seq, h.ppm_error))
                .collect::<Vec<_>>()
        );
    }
}
