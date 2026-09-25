//! mass_corasick_v2: High-performance single-pass mass-tolerant peptide matching
//!
//! Optimized over v1:
//!   - Flat sorted arrays instead of HashMap bins → cache-friendly, SIMD-friendly
//!   - Two-level filtering: precursor mass range → k-mer mass bins
//!   - Tighter bin width (0.002 Da) for lower collision rate
//!   - Streaming protein processing — constant memory per protein
//!   - Chunk-parallel with per-thread vote buffers (no Mutex contention)

// Legacy/experimental reference binary — not on the canonical pipeline; clippy relaxed.
#![allow(clippy::all)]
#![allow(dead_code, unused)]

use clap::Parser;
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Instant;

#[inline(always)]
fn aa_mass(aa: u8) -> f64 {
    // Lookup table for speed — avoid match branching in hot loop
    const TABLE: [f64; 26] = [
        71.03711,  // A
        0.0,       // B
        103.00919, // C
        115.02694, // D
        129.04259, // E
        147.06841, // F
        57.02146,  // G
        137.05891, // H
        113.08406, // I
        0.0,       // J
        128.09496, // K
        113.08406, // L
        131.04049, // M
        114.04293, // N
        0.0,       // O
        97.05276,  // P
        128.05858, // Q
        156.10111, // R
        87.03203,  // S
        101.04768, // T
        0.0,       // U
        99.06841,  // V
        186.07931, // W
        0.0,       // X
        163.06333, // Y
        0.0,       // Z
    ];
    if aa >= b'A' && aa <= b'Z' {
        TABLE[(aa - b'A') as usize]
    } else {
        0.0
    }
}

const K: usize = 5;
const BIN_WIDTH: f64 = 0.001; // Matched to mass_kmer_anchor for identical results
const WATER_MASS: f64 = 18.010565;
const PROTON_MASS: f64 = 1.007276;

// Mass range for 5-mers: 5×Gly(57) = 285 to 5×Trp(186) = 932
const MIN_KMER_MASS: f64 = 280.0;
const MAX_KMER_MASS: f64 = 940.0;
const N_BINS: usize = ((MAX_KMER_MASS - MIN_KMER_MASS) / BIN_WIDTH) as usize + 100; // +margin for safety

#[inline(always)]
fn mass_to_bin(mass: f64) -> usize {
    let raw = ((mass - MIN_KMER_MASS) / BIN_WIDTH) as i64;
    if raw < 0 {
        0
    } else if raw >= N_BINS as i64 {
        N_BINS - 1
    } else {
        raw as usize
    }
}

/// Peptide with precomputed data
struct PeptideQuery {
    sequence: Vec<u8>,
    kmer_bins: Vec<usize>,
    precursor_mass: f64,
    score: f64,
}

/// Flat-array peptide index for cache-friendly lookups
struct PeptideIndex {
    queries: Vec<PeptideQuery>,
    /// For each bin: a sorted slice of (peptide_idx, kmer_position) pairs
    /// Stored as one contiguous Vec with offset/length arrays
    entries: Vec<(u32, u16)>,
    bin_offsets: Vec<u32>,
    bin_lengths: Vec<u32>,
    total_kmers: usize,
}

impl PeptideIndex {
    fn build(queries: Vec<PeptideQuery>) -> Self {
        let t0 = Instant::now();

        // Count entries per bin
        let mut bin_counts = vec![0u32; N_BINS];
        let mut total_kmers = 0usize;

        for (qidx, query) in queries.iter().enumerate() {
            for (kpos, &bin) in query.kmer_bins.iter().enumerate() {
                bin_counts[bin] += 1;
                total_kmers += 1;
                // Also count adjacent bins for ±1 tolerance
                if bin > 0 {
                    bin_counts[bin - 1] += 1;
                    total_kmers += 1;
                }
                if bin + 1 < N_BINS {
                    bin_counts[bin + 1] += 1;
                    total_kmers += 1;
                }
            }

            if (qidx + 1) % 1_000_000 == 0 {
                eprintln!(
                    "  Counting bins: {}/{} peptides...",
                    qidx + 1,
                    queries.len()
                );
            }
        }

        // Compute offsets
        let mut bin_offsets = vec![0u32; N_BINS];
        let mut running = 0u32;
        for i in 0..N_BINS {
            bin_offsets[i] = running;
            running += bin_counts[i];
        }

        // Fill entries
        let mut entries = vec![(0u32, 0u16); total_kmers];
        let mut write_pos = vec![0u32; N_BINS]; // current write position per bin

        for (qidx, query) in queries.iter().enumerate() {
            for (kpos, &bin) in query.kmer_bins.iter().enumerate() {
                let entry = (qidx as u32, kpos as u16);

                // Write to bin and adjacent bins
                for b in [bin.wrapping_sub(1), bin, bin + 1] {
                    if b < N_BINS {
                        let offset = bin_offsets[b] as usize + write_pos[b] as usize;
                        entries[offset] = entry;
                        write_pos[b] += 1;
                    }
                }
            }

            if (qidx + 1) % 1_000_000 == 0 {
                eprintln!(
                    "  Filling index: {}/{} peptides...",
                    qidx + 1,
                    queries.len()
                );
            }
        }

        // Convert write_pos to bin_lengths
        let bin_lengths = write_pos;

        let elapsed = t0.elapsed();
        let non_empty = bin_counts.iter().filter(|&&c| c > 0).count();
        let max_bin = bin_counts.iter().max().unwrap_or(&0);
        let avg_bin = if non_empty > 0 {
            total_kmers as f64 / non_empty as f64
        } else {
            0.0
        };

        eprintln!(
            "Peptide index built: {} peptides, {} k-mer entries, {:.1}s",
            queries.len(),
            total_kmers,
            elapsed.as_secs_f64()
        );
        eprintln!(
            "  {} non-empty bins, avg {:.0} entries, max {} entries",
            non_empty, avg_bin, max_bin
        );
        eprintln!("  Index memory: {:.1} GB", (total_kmers * 6) as f64 / 1e9);

        PeptideIndex {
            queries,
            entries,
            bin_offsets,
            bin_lengths,
            total_kmers,
        }
    }

    /// Get entries for a mass bin — returns a slice into the flat array
    #[inline(always)]
    fn get_bin(&self, bin: usize) -> &[(u32, u16)] {
        if bin >= N_BINS {
            return &[];
        }
        let offset = self.bin_offsets[bin] as usize;
        let length = self.bin_lengths[bin] as usize;
        &self.entries[offset..offset + length]
    }

    /// Sweep one protein against the peptide index
    fn sweep_protein(&self, prot_seq: &[u8], min_votes: usize) -> Vec<(u32, u16, u16)> {
        // Returns: Vec<(peptide_idx, votes, position_in_protein)>
        if prot_seq.len() < K {
            return vec![];
        }

        // Compute cumulative mass for this protein
        let mut cum = Vec::with_capacity(prot_seq.len() + 1);
        cum.push(0.0f64);
        let mut running = 0.0;
        for &aa in prot_seq {
            running += aa_mass(aa);
            cum.push(running);
        }

        // Vote buffer: (peptide_idx, inferred_start) → vote_count
        // Use a flat HashMap — cleared per protein
        let mut votes: HashMap<u64, u16> = HashMap::with_capacity(256);

        for prot_pos in 0..=(prot_seq.len() - K) {
            let mass = cum[prot_pos + K] - cum[prot_pos];
            let bin = mass_to_bin(mass);

            // Tolerance already baked into index (±1 bin during build)
            let entries = self.get_bin(bin);

            for &(qidx, kpos) in entries {
                let inferred_start = prot_pos as i32 - kpos as i32;
                if inferred_start < 0 {
                    continue;
                }
                // Pack (qidx, inferred_start) into u64 key
                let key = (qidx as u64) << 32 | (inferred_start as u64);
                let v = votes.entry(key).or_insert(0);
                *v = v.saturating_add(1);
            }
        }

        // Collect hits above threshold
        let mut results = Vec::new();
        for (&key, &vote_count) in &votes {
            if (vote_count as usize) < min_votes {
                continue;
            }
            let qidx = (key >> 32) as u32;
            let start = (key & 0xFFFFFFFF) as u16;
            results.push((qidx, vote_count, start));
        }

        results
    }
}

fn load_predictions(csv_dir: &str) -> Vec<PeptideQuery> {
    let t0 = Instant::now();
    let mut raw: Vec<(Vec<u8>, f64, f64)> = Vec::new();
    let mut best: HashMap<Vec<u8>, usize> = HashMap::new();

    let pattern = format!("{}/**/*predictions.csv", csv_dir);
    let mut paths: Vec<_> = glob::glob(&pattern)
        .expect("Invalid glob")
        .filter_map(|p| p.ok())
        .collect();
    if paths.is_empty() {
        let p2 = format!("{}/**/*.csv", csv_dir);
        paths = glob::glob(&p2)
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
                if *counts.iter().max().unwrap_or(&0) as usize * 2 > seq.len() {
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

            if let Some(&eidx) = best.get(&seq) {
                if score > raw[eidx].1 {
                    raw[eidx] = (seq.clone(), score, precursor_mass);
                }
            } else {
                let idx = raw.len();
                best.insert(seq.clone(), idx);
                raw.push((seq, score, precursor_mass));
            }
        }
    }

    let queries: Vec<PeptideQuery> = raw
        .into_iter()
        .filter(|(seq, _, _)| seq.len() >= K)
        .map(|(seq, score, precursor_mass)| {
            let mut cum = Vec::with_capacity(seq.len() + 1);
            cum.push(0.0);
            let mut r = 0.0;
            for &aa in &seq {
                r += aa_mass(aa);
                cum.push(r);
            }
            let kmer_bins: Vec<usize> = (0..=(seq.len() - K))
                .map(|i| mass_to_bin(cum[i + K] - cum[i]))
                .collect();
            PeptideQuery {
                sequence: seq,
                kmer_bins,
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
#[command(name = "mass_corasick_v2")]
#[command(
    about = "High-performance single-pass mass-tolerant peptide matching (v2: flat index, no sharding)"
)]
struct Args {
    #[arg(short = 'd', long)]
    database: String,
    #[arg(short = 'p', long)]
    predictions: String,
    #[arg(short = 'o', long)]
    output: String,
    #[arg(long, default_value = "5")]
    min_votes: usize,
    #[arg(short = 't', long, default_value = "0")]
    threads: usize,
    #[arg(long, default_value = "5000")]
    chunk_size: usize,
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

    eprintln!("=== mass_corasick v2: flat-index single-pass matching ===");
    eprintln!("");

    // Load and index peptides
    let queries = load_predictions(&args.predictions);
    if queries.is_empty() {
        eprintln!("No predictions loaded, exiting.");
        return;
    }
    let index = PeptideIndex::build(queries);

    // Load proteins
    eprintln!("");
    eprintln!("Loading database: {}", args.database);
    let load_t0 = Instant::now();
    let file = File::open(&args.database)
        .unwrap_or_else(|e| panic!("Cannot open {}: {}", args.database, e));
    let reader = BufReader::with_capacity(16 * 1024 * 1024, file);

    let mut proteins: Vec<(String, Vec<u8>)> = Vec::new();
    let mut current_id = String::new();
    let mut current_seq: Vec<u8> = Vec::new();

    for line in reader.lines() {
        let line = line.expect("read error");
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
        if proteins.len() % 10_000_000 == 0 && !proteins.is_empty() && current_seq.is_empty() {
            eprintln!("  Loaded {}M proteins...", proteins.len() / 1_000_000);
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
    eprintln!("");
    eprintln!(
        "Sweeping {} proteins with {} peptide queries...",
        proteins.len(),
        index.queries.len()
    );
    let sweep_t0 = Instant::now();
    let processed = AtomicUsize::new(0);

    // Per-peptide best hit: use atomic-friendly approach
    // Each thread collects local results, then merge
    let all_results: Vec<Vec<(u32, u16, String, u16)>> = proteins
        .par_chunks(args.chunk_size)
        .map(|chunk| {
            let mut local_results = Vec::new();
            for (prot_id, prot_seq) in chunk {
                let hits = index.sweep_protein(prot_seq, args.min_votes);
                for (qidx, votes, start) in hits {
                    // Extract db subsequence
                    let query = &index.queries[qidx as usize];
                    let start_u = start as usize;
                    let end = (start_u + query.sequence.len()).min(prot_seq.len());
                    if end <= start_u {
                        continue;
                    }
                    local_results.push((qidx, votes, prot_id.clone(), start));
                }
            }

            let p = processed.fetch_add(chunk.len(), Ordering::Relaxed) + chunk.len();
            if p % 1_000_000 < args.chunk_size {
                let elapsed = sweep_t0.elapsed().as_secs_f64();
                let rate = p as f64 / elapsed;
                let remaining = (proteins.len() - p) as f64 / rate;
                eprintln!(
                    "  {}/{}  ({:.0} prot/s, ETA {:.0}m)",
                    p,
                    proteins.len(),
                    rate,
                    remaining / 60.0
                );
            }

            local_results
        })
        .collect();

    let sweep_time = sweep_t0.elapsed();

    // Merge: keep best hit per peptide across all chunks
    eprintln!("Merging results...");
    let mut best_per_peptide: HashMap<u32, (u16, String, u16)> = HashMap::new();
    for chunk_results in &all_results {
        for (qidx, votes, prot_id, start) in chunk_results {
            let entry = best_per_peptide
                .entry(*qidx)
                .or_insert((0, String::new(), 0));
            if *votes > entry.0 {
                *entry = (*votes, prot_id.clone(), *start);
            }
        }
    }

    let rescued_count = best_per_peptide.len();
    let unique_proteins: std::collections::HashSet<&str> = best_per_peptide
        .values()
        .map(|(_, pid, _)| pid.as_str())
        .collect();

    eprintln!(
        "\nSweep complete: {} proteins in {:.1}s ({:.0} prot/s)",
        proteins.len(),
        sweep_time.as_secs_f64(),
        proteins.len() as f64 / sweep_time.as_secs_f64()
    );
    eprintln!(
        "Rescued: {} peptides → {} proteins",
        rescued_count,
        unique_proteins.len()
    );

    // Write output
    let out_file = File::create(&args.output).expect("Cannot create output");
    let mut writer = BufWriter::new(out_file);
    writeln!(
        writer,
        "query_sequence\tprotein_id\tposition\tvotes\tn_kmers\tscore"
    )
    .unwrap();

    let mut sorted: Vec<_> = best_per_peptide.iter().collect();
    sorted.sort_by(|a, b| b.1 .0.cmp(&a.1 .0));

    for (qidx, (votes, prot_id, start)) in &sorted {
        let query = &index.queries[**qidx as usize];
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{:.4}",
            String::from_utf8_lossy(&query.sequence),
            prot_id,
            start,
            votes,
            query.kmer_bins.len(),
            query.score,
        )
        .unwrap();
    }

    let total_time = total_t0.elapsed();
    eprintln!("\n=== SUMMARY ===");
    eprintln!("  Peptides:         {}", index.queries.len());
    eprintln!("  Database:         {} proteins", proteins.len());
    eprintln!("  Rescued:          {} peptides", rescued_count);
    eprintln!("  Proteins found:   {}", unique_proteins.len());
    eprintln!(
        "  Sweep rate:       {:.0} prot/s",
        proteins.len() as f64 / sweep_time.as_secs_f64()
    );
    eprintln!(
        "  Index memory:     {:.1} GB",
        (index.total_kmers * 6) as f64 / 1e9
    );
    eprintln!("  Total time:       {:.1}s", total_time.as_secs_f64());

    if let Some(stats_path) = &args.output_stats {
        let stats = serde_json::json!({
            "total_peptides": index.queries.len(),
            "rescued_peptides": rescued_count,
            "rescued_proteins": unique_proteins.len(),
            "database_proteins": proteins.len(),
            "sweep_rate_prot_s": proteins.len() as f64 / sweep_time.as_secs_f64(),
            "index_memory_gb": (index.total_kmers * 6) as f64 / 1e9,
            "total_time_s": total_time.as_secs_f64(),
        });
        std::fs::write(stats_path, serde_json::to_string_pretty(&stats).unwrap())
            .expect("Cannot write stats");
    }
}
