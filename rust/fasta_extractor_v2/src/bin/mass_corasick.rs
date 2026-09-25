//! mass_corasick: Single-pass mass-tolerant peptide-protein matching
//!
//! Like Aho-Corasick for masses. Indexes peptide k-mer masses, then sweeps
//! through the protein database in one pass. No sharding needed.
//!
//! Architecture:
//!   1. Load all peptide predictions → compute 5-mer mass fingerprints
//!   2. Build inverted index: mass_bin → [(peptide_idx, kmer_position)]
//!   3. Stream proteins one by one:
//!      - Compute protein k-mer masses
//!      - Look up matching peptide k-mers
//!      - Vote per (peptide, inferred_offset)
//!      - Emit rescues where votes ≥ min_votes
//!   4. Output rescued peptide-protein pairs
//!
//! Memory: O(peptides × k-mers_per_peptide) for the index — typically 5-50 GB
//! Time: O(proteins × avg_length) for the sweep — single pass, no sharding

// Legacy/experimental reference binary — not on the canonical pipeline; clippy relaxed.
#![allow(clippy::all)]
#![allow(dead_code, unused)]

use clap::Parser;
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;
use std::time::Instant;

/// Monoisotopic amino acid residue masses (Da)
#[inline]
fn aa_mass(aa: u8) -> f64 {
    match aa {
        b'G' => 57.02146,
        b'A' => 71.03711,
        b'V' => 99.06841,
        b'L' | b'I' => 113.08406,
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
const BIN_WIDTH: f64 = 0.01; // Wider bins than mass_kmer_anchor (0.01 vs 0.001)
const WATER_MASS: f64 = 18.010565;
const PROTON_MASS: f64 = 1.007276;

#[inline]
fn mass_to_bin(mass: f64) -> i64 {
    (mass / BIN_WIDTH).round() as i64
}

/// A peptide query with precomputed k-mer masses
struct PeptideQuery {
    sequence: Vec<u8>,
    kmer_masses: Vec<i64>, // binned masses for each k-mer position
    precursor_mass: f64,
    score: f64,
}

/// Inverted index: mass_bin → Vec<(peptide_idx, kmer_position)>
struct PeptideIndex {
    queries: Vec<PeptideQuery>,
    bins: HashMap<i64, Vec<(u32, u16)>>, // (peptide_idx, kmer_pos)
}

impl PeptideIndex {
    fn build(queries: Vec<PeptideQuery>) -> Self {
        let t0 = Instant::now();
        let mut bins: HashMap<i64, Vec<(u32, u16)>> = HashMap::new();
        let mut total_kmers: usize = 0;

        for (qidx, query) in queries.iter().enumerate() {
            for (kpos, &bin_id) in query.kmer_masses.iter().enumerate() {
                bins.entry(bin_id)
                    .or_default()
                    .push((qidx as u32, kpos as u16));
                total_kmers += 1;
            }

            if (qidx + 1) % 1_000_000 == 0 {
                eprintln!(
                    "  Indexed {}/{} peptides ({} k-mers)...",
                    qidx + 1,
                    queries.len(),
                    total_kmers
                );
            }
        }

        eprintln!(
            "Peptide index built: {} peptides, {} k-mers, {} bins, {:.1}s",
            queries.len(),
            total_kmers,
            bins.len(),
            t0.elapsed().as_secs_f64()
        );

        // Stats on bin sizes
        let max_bin = bins.values().map(|v| v.len()).max().unwrap_or(0);
        let avg_bin = total_kmers as f64 / bins.len().max(1) as f64;
        eprintln!(
            "  Bin stats: avg {:.0} entries, max {} entries",
            avg_bin, max_bin
        );

        PeptideIndex { queries, bins }
    }

    /// Sweep one protein against the peptide index
    /// Returns rescued (peptide_idx, votes, position_in_protein) tuples
    fn sweep_protein(
        &self,
        prot_seq: &[u8],
        prot_id: &str,
        min_votes: usize,
    ) -> Vec<(usize, usize, usize, String)> {
        if prot_seq.len() < K {
            return vec![];
        }

        // Compute protein k-mer masses
        let mut cum = Vec::with_capacity(prot_seq.len() + 1);
        cum.push(0.0f64);
        let mut running = 0.0;
        for &aa in prot_seq {
            running += aa_mass(aa);
            cum.push(running);
        }

        // Vote: (peptide_idx, inferred_prot_start) → count
        let mut votes: HashMap<(u32, i32), u16> = HashMap::new();

        for prot_pos in 0..=(prot_seq.len() - K) {
            let prot_kmer_mass = cum[prot_pos + K] - cum[prot_pos];
            let bin_id = mass_to_bin(prot_kmer_mass);

            for b in [bin_id - 1, bin_id, bin_id + 1] {
                if let Some(entries) = self.bins.get(&b) {
                    for &(qidx, kpos) in entries {
                        let inferred_start = prot_pos as i32 - kpos as i32;
                        if inferred_start >= 0 {
                            let v = votes.entry((qidx, inferred_start)).or_insert(0);
                            *v = v.saturating_add(1);
                        }
                    }
                }
            }
        }

        // Collect rescues
        let mut results = Vec::new();
        for (&(qidx, start), &vote_count) in &votes {
            if (vote_count as usize) < min_votes {
                continue;
            }
            let query = &self.queries[qidx as usize];
            let start = start as usize;

            // Try ±2 length offsets, pick best
            for offset in -2i32..=2 {
                let pep_len = query.sequence.len() as i32 + offset;
                if pep_len < K as i32 {
                    continue;
                }
                let end = start + pep_len as usize;
                if end > prot_seq.len() {
                    continue;
                }

                let db_seq = String::from_utf8_lossy(&prot_seq[start..end]).to_string();
                results.push((qidx as usize, vote_count as usize, start, db_seq));
                break; // Take first valid offset
            }
        }

        results
    }
}

fn load_predictions(csv_dir: &str) -> Vec<PeptideQuery> {
    let t0 = Instant::now();
    let mut raw_preds: Vec<(Vec<u8>, f64, f64)> = Vec::new();
    let mut best_scores: HashMap<Vec<u8>, usize> = HashMap::new();

    let pattern = format!("{}/**/*predictions.csv", csv_dir);
    let mut paths: Vec<_> = glob::glob(&pattern)
        .expect("Invalid glob")
        .filter_map(|p| p.ok())
        .collect();

    if paths.is_empty() {
        let pattern2 = format!("{}/**/*.csv", csv_dir);
        paths = glob::glob(&pattern2)
            .expect("Invalid glob")
            .filter_map(|p| p.ok())
            .collect();
    }

    if paths.is_empty() {
        let path = std::path::PathBuf::from(csv_dir);
        if path.is_file() {
            paths.push(path);
        }
    }

    for path in &paths {
        let file = match File::open(path) {
            Ok(f) => f,
            Err(_) => continue,
        };
        let mut reader = csv::Reader::from_reader(BufReader::new(file));
        let headers = match reader.headers() {
            Ok(h) => h.clone(),
            Err(_) => continue,
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
            None => continue,
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
            // Low-complexity filter
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
            let precursor_mass = if mz > 0.0 && charge > 0.0 {
                mz * charge - charge * PROTON_MASS
            } else {
                0.0
            };

            if let Some(&existing_idx) = best_scores.get(&seq) {
                if score > raw_preds[existing_idx].1 {
                    raw_preds[existing_idx] = (seq.clone(), score, precursor_mass);
                }
            } else {
                let idx = raw_preds.len();
                best_scores.insert(seq.clone(), idx);
                raw_preds.push((seq, score, precursor_mass));
            }
        }
    }

    // Convert to PeptideQuery with precomputed k-mer masses
    let queries: Vec<PeptideQuery> = raw_preds
        .into_iter()
        .filter(|(seq, _, _)| seq.len() >= K)
        .map(|(seq, score, precursor_mass)| {
            let mut cum = Vec::with_capacity(seq.len() + 1);
            cum.push(0.0);
            let mut running = 0.0;
            for &aa in &seq {
                running += aa_mass(aa);
                cum.push(running);
            }
            let kmer_masses: Vec<i64> = (0..=(seq.len() - K))
                .map(|i| mass_to_bin(cum[i + K] - cum[i]))
                .collect();

            PeptideQuery {
                sequence: seq,
                kmer_masses,
                precursor_mass,
                score,
            }
        })
        .collect();

    eprintln!(
        "Loaded {} unique peptides from {} files in {:.1}s",
        queries.len(),
        paths.len(),
        t0.elapsed().as_secs_f64()
    );
    queries
}

#[derive(Parser)]
#[command(name = "mass_corasick")]
#[command(about = "Single-pass mass-tolerant peptide matching — like Aho-Corasick for masses")]
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
    #[arg(long, default_value = "5")]
    min_votes: usize,

    /// Number of threads
    #[arg(short = 't', long, default_value = "0")]
    threads: usize,

    /// Chunk size for parallel protein processing
    #[arg(long, default_value = "10000")]
    chunk_size: usize,

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

    // Step 1: Load and index peptides
    eprintln!("=== mass_corasick: single-pass mass-tolerant matching ===");
    eprintln!("");
    let queries = load_predictions(&args.predictions);
    if queries.is_empty() {
        eprintln!("No predictions loaded, exiting.");
        return;
    }
    let index = PeptideIndex::build(queries);

    // Step 2: Stream proteins and sweep
    eprintln!("");
    eprintln!("Sweeping database: {}", args.database);

    let file = File::open(&args.database)
        .unwrap_or_else(|e| panic!("Cannot open {}: {}", args.database, e));
    let reader = BufReader::with_capacity(8 * 1024 * 1024, file);

    // Load all proteins into chunks for parallel processing
    let load_t0 = Instant::now();
    let mut proteins: Vec<(String, Vec<u8>)> = Vec::new();
    let mut current_id = String::new();
    let mut current_seq = Vec::new();

    for line in reader.lines() {
        let line = line.expect("Failed to read line");
        if line.starts_with('>') {
            if !current_id.is_empty() && !current_seq.is_empty() {
                proteins.push((
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

        if proteins.len() % 5_000_000 == 0 && proteins.len() > 0 && current_seq.is_empty() {
            eprintln!("  Loaded {} proteins...", proteins.len());
        }
    }
    if !current_id.is_empty() && !current_seq.is_empty() {
        proteins.push((current_id, current_seq));
    }
    eprintln!(
        "Loaded {} proteins in {:.1}s",
        proteins.len(),
        load_t0.elapsed().as_secs_f64()
    );

    // Parallel sweep
    let sweep_t0 = Instant::now();
    let processed = AtomicUsize::new(0);
    let rescued_total = AtomicUsize::new(0);

    // Best hit per peptide: peptide_idx → (votes, protein_id, position, db_seq)
    let best_hits: Vec<Mutex<Option<(usize, String, usize, String)>>> =
        (0..index.queries.len()).map(|_| Mutex::new(None)).collect();

    proteins.par_chunks(args.chunk_size).for_each(|chunk| {
        for (prot_id, prot_seq) in chunk {
            let hits = index.sweep_protein(prot_seq, prot_id, args.min_votes);

            for (qidx, votes, pos, db_seq) in hits {
                let mut best = best_hits[qidx].lock().unwrap();
                if best.is_none() || votes > best.as_ref().unwrap().0 {
                    *best = Some((votes, prot_id.clone(), pos, db_seq));
                }
            }
        }

        let p = processed.fetch_add(chunk.len(), Ordering::Relaxed) + chunk.len();
        if p % 500_000 < args.chunk_size {
            let r = best_hits
                .iter()
                .filter(|h| h.lock().unwrap().is_some())
                .count();
            rescued_total.store(r, Ordering::Relaxed);
            eprintln!(
                "  Swept {}/{} proteins ({} peptides rescued)",
                p,
                proteins.len(),
                r
            );
        }
    });

    let sweep_time = sweep_t0.elapsed();
    let final_rescued: Vec<_> = best_hits
        .iter()
        .enumerate()
        .filter_map(|(qidx, h)| {
            let guard = h.lock().unwrap();
            guard.as_ref().map(|(votes, prot_id, pos, db_seq)| {
                (qidx, *votes, prot_id.clone(), *pos, db_seq.clone())
            })
        })
        .collect();

    eprintln!(
        "\nSweep complete: {} proteins in {:.1}s ({:.0} prot/s)",
        proteins.len(),
        sweep_time.as_secs_f64(),
        proteins.len() as f64 / sweep_time.as_secs_f64()
    );
    eprintln!("Rescued: {} peptides", final_rescued.len());

    let unique_proteins: std::collections::HashSet<&str> = final_rescued
        .iter()
        .map(|(_, _, pid, _, _)| pid.as_str())
        .collect();
    eprintln!("Unique rescued proteins: {}", unique_proteins.len());

    // Write output
    let out_file = File::create(&args.output).expect("Cannot create output file");
    let mut writer = BufWriter::new(out_file);
    writeln!(
        writer,
        "query_sequence\tdb_sequence\tprotein_id\tposition\tvotes\tn_kmers\tscore"
    )
    .unwrap();

    for (qidx, votes, prot_id, pos, db_seq) in &final_rescued {
        let query = &index.queries[*qidx];
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{:.4}",
            String::from_utf8_lossy(&query.sequence),
            db_seq,
            prot_id,
            pos,
            votes,
            query.kmer_masses.len(),
            query.score,
        )
        .unwrap();
    }

    let total_time = total_t0.elapsed();
    eprintln!("\nSUMMARY:");
    eprintln!("  Total peptides:     {}", index.queries.len());
    eprintln!("  Rescued:            {}", final_rescued.len());
    eprintln!("  Rescued proteins:   {}", unique_proteins.len());
    eprintln!("  Database proteins:  {}", proteins.len());
    eprintln!(
        "  Sweep rate:         {:.0} proteins/s",
        proteins.len() as f64 / sweep_time.as_secs_f64()
    );
    eprintln!("  Total time:         {:.1}s", total_time.as_secs_f64());

    // Stats JSON
    if let Some(stats_path) = &args.output_stats {
        let stats = serde_json::json!({
            "total_peptides": index.queries.len(),
            "rescued_peptides": final_rescued.len(),
            "rescued_proteins": unique_proteins.len(),
            "database_proteins": proteins.len(),
            "sweep_rate": proteins.len() as f64 / sweep_time.as_secs_f64(),
            "total_time_s": total_time.as_secs_f64(),
        });
        std::fs::write(stats_path, serde_json::to_string_pretty(&stats).unwrap())
            .expect("Cannot write stats");
    }
}
