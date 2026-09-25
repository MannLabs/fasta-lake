//! Bayesian protein parsimony — evidence-driven FASTA reduction.
//!
//! Runs EM to compute per-protein probabilities from de novo peptide evidence,
//! then filters to proteins above a threshold. Designed for multi-million
//! protein evidence lakes where Python is too slow.
//!
//! Usage:
//!   parsimony --fasta evidence_lake.fasta \
//!             --predictions-dir csvs/ \
//!             --output scored_lake.fasta \
//!             --min-probability 0.01 \
//!             --iterations 10
//!
//! Also outputs a TSV of protein_id → probability for analysis.

// Legacy/experimental reference binary — not on the canonical pipeline; clippy relaxed.
#![allow(clippy::all)]
#![allow(dead_code, unused)]

#[path = "../../../common/dianovo.rs"]
mod dianovo;

use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::PathBuf;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Instant;

use clap::Parser;
use dashmap::DashMap;
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;

#[derive(Parser)]
#[command(name = "parsimony", about = "Bayesian protein parsimony with EM")]
struct Args {
    /// FASTA database (evidence lake or per-sample)
    #[arg(short, long)]
    fasta: PathBuf,

    /// Directory containing AlphaNovo/DIANovo prediction CSVs
    #[arg(short, long)]
    predictions_dir: PathBuf,

    /// Output FASTA (filtered to proteins above threshold)
    #[arg(short, long)]
    output: PathBuf,

    /// Output scores TSV (protein_id \t probability)
    #[arg(long)]
    output_scores: Option<PathBuf>,

    /// Output stats JSON
    #[arg(long)]
    output_stats: Option<PathBuf>,

    /// Minimum probability threshold (default: 0.01)
    #[arg(long, default_value_t = 0.01)]
    min_probability: f64,

    /// Number of EM iterations (default: 10)
    #[arg(long, default_value_t = 10)]
    iterations: usize,

    /// Top percent of peptides to keep per sample (default: 0.30)
    #[arg(long, default_value_t = 0.30)]
    top_percent: f64,

    /// Minimum peptide length (default: 8)
    #[arg(long, default_value_t = 8)]
    min_length: usize,

    /// Maximum peptide length (default: 50)
    #[arg(long, default_value_t = 50)]
    max_length: usize,

    /// Number of threads (default: all available)
    #[arg(short, long, default_value_t = 0)]
    threads: usize,

    /// Strategy: bayesian, info_score, razor (default: bayesian)
    #[arg(long, default_value = "bayesian")]
    strategy: String,

    /// Normalize I/L (default: true)
    #[arg(long, default_value_t = true)]
    normalize_il: bool,
}

/// A protein: id, sequence, header line
struct Protein {
    id: String,
    sequence: Vec<u8>,
    header: String,
}

/// Parsed peptide with score
struct Peptide {
    sequence: String,
    score: f64,
    sample_idx: usize,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();
    let t_start = Instant::now();

    if args.threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global()
            .ok();
    }
    let n_threads = rayon::current_num_threads();

    println!("=======================================================");
    println!("PARSIMONY — Bayesian Protein Evidence Scoring");
    println!("=======================================================");
    println!("Strategy: {}", args.strategy);
    println!("Threads: {}", n_threads);
    println!("Min probability: {}", args.min_probability);
    println!("EM iterations: {}", args.iterations);
    println!();

    // ---------------------------------------------------------------
    // Step 1: Load FASTA
    // ---------------------------------------------------------------
    println!("Step 1: Loading FASTA...");
    let t1 = Instant::now();
    let proteins = load_fasta(&args.fasta)?;
    let n_proteins = proteins.len();
    println!(
        "  {} proteins in {:.2}s",
        n_proteins,
        t1.elapsed().as_secs_f64()
    );

    // Build protein index: id → index
    let prot_idx: HashMap<&str, usize> = proteins
        .iter()
        .enumerate()
        .map(|(i, p)| (p.id.as_str(), i))
        .collect();

    // Protein lengths (sqrt for normalization)
    let sqrt_lengths: Vec<f64> = proteins
        .iter()
        .map(|p| (p.sequence.len().max(50) as f64).sqrt())
        .collect();

    // ---------------------------------------------------------------
    // Step 2: Load predictions
    // ---------------------------------------------------------------
    println!("\nStep 2: Loading predictions...");
    let t2 = Instant::now();
    let csv_files = collect_csvs(&args.predictions_dir);
    println!("  {} CSV files", csv_files.len());

    // Parse all CSVs in parallel, collect (peptide, score, sample_idx)
    let all_peptides: Vec<Vec<Peptide>> = csv_files
        .par_iter()
        .enumerate()
        .map(|(sample_idx, path)| {
            parse_prediction_csv(
                path,
                args.min_length,
                args.max_length,
                args.normalize_il,
                args.top_percent,
                sample_idx,
            )
        })
        .collect();

    let n_samples = csv_files.len();
    let total_peptides: usize = all_peptides.iter().map(|v| v.len()).sum();
    println!(
        "  {} peptides from {} samples in {:.2}s",
        total_peptides,
        n_samples,
        t2.elapsed().as_secs_f64()
    );

    // ---------------------------------------------------------------
    // Step 3: Build peptide → protein map (Aho-Corasick)
    // ---------------------------------------------------------------
    println!("\nStep 3: Building peptide-protein map...");
    let t3 = Instant::now();

    // Collect unique peptides with best score
    let mut best_scores: HashMap<String, f64> = HashMap::new();
    let mut pep_samples: HashMap<String, HashSet<usize>> = HashMap::new();
    for sample_peps in &all_peptides {
        for pep in sample_peps {
            let entry = best_scores.entry(pep.sequence.clone()).or_insert(0.0);
            if pep.score > *entry {
                *entry = pep.score;
            }
            pep_samples
                .entry(pep.sequence.clone())
                .or_default()
                .insert(pep.sample_idx);
        }
    }
    let unique_peps: Vec<String> = best_scores.keys().cloned().collect();
    println!("  {} unique peptides", unique_peps.len());

    // Build Aho-Corasick automaton
    let ac = aho_corasick::AhoCorasick::new(&unique_peps).expect("Failed to build automaton");

    // Search all proteins in parallel
    let pb = ProgressBar::new(n_proteins as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("[{elapsed_precise}] {bar:50} {pos}/{len} proteins ({per_sec})")
            .unwrap(),
    );

    // pep2prot: peptide_index → Vec<protein_index>
    let pep2prot: DashMap<usize, Vec<usize>> = DashMap::new();
    // prot2pep: protein_index → Vec<peptide_index>
    let prot2pep: DashMap<usize, Vec<usize>> = DashMap::new();

    proteins
        .par_iter()
        .enumerate()
        .for_each(|(prot_i, protein)| {
            let matches: HashSet<usize> = ac
                .find_overlapping_iter(&protein.sequence)
                .map(|m| m.pattern().as_usize())
                .collect();

            if !matches.is_empty() {
                let mut peps_for_prot = Vec::with_capacity(matches.len());
                for pep_i in matches {
                    pep2prot.entry(pep_i).or_default().push(prot_i);
                    peps_for_prot.push(pep_i);
                }
                prot2pep.insert(prot_i, peps_for_prot);
            }
            pb.inc(1);
        });
    pb.finish();

    let n_mapped_peps = pep2prot.len();
    let n_mapped_prots = prot2pep.len();
    println!(
        "  {} peptides mapped to {} proteins in {:.2}s",
        n_mapped_peps,
        n_mapped_prots,
        t3.elapsed().as_secs_f64()
    );

    // Convert DashMaps to Vecs for fast iteration
    let pep2prot_vec: Vec<(usize, Vec<usize>)> = pep2prot.into_iter().collect();
    let prot2pep_vec: HashMap<usize, Vec<usize>> = prot2pep.into_iter().collect();

    // Count unique peptides
    let n_unique = pep2prot_vec
        .iter()
        .filter(|(_, prots)| prots.len() == 1)
        .count();
    let n_shared = pep2prot_vec.len() - n_unique;
    println!(
        "  {} unique ({:.1}%), {} shared ({:.1}%)",
        n_unique,
        n_unique as f64 / pep2prot_vec.len() as f64 * 100.0,
        n_shared,
        n_shared as f64 / pep2prot_vec.len() as f64 * 100.0
    );

    // Pre-compute sample counts per protein
    let prot_sample_count: Vec<f64> = (0..n_proteins)
        .into_par_iter()
        .map(|prot_i| {
            let peps = match prot2pep_vec.get(&prot_i) {
                Some(p) => p,
                None => return 0.0,
            };
            let mut samples: HashSet<usize> = HashSet::new();
            for &pep_i in peps {
                if let Some(s) = pep_samples.get(&unique_peps[pep_i]) {
                    samples.extend(s);
                }
            }
            (samples.len() as f64).sqrt()
        })
        .collect();

    // ---------------------------------------------------------------
    // Step 4: Run strategy
    // ---------------------------------------------------------------
    println!(
        "\nStep 4: Running {} ({} iterations)...",
        args.strategy, args.iterations
    );
    let t4 = Instant::now();

    let probabilities = match args.strategy.as_str() {
        "bayesian" => run_bayesian_em(
            n_proteins,
            &pep2prot_vec,
            &prot2pep_vec,
            &unique_peps,
            &best_scores,
            &sqrt_lengths,
            &prot_sample_count,
            args.iterations,
        ),
        "info_score" => run_info_score(
            n_proteins,
            &pep2prot_vec,
            &prot2pep_vec,
            &unique_peps,
            &best_scores,
            &sqrt_lengths,
        ),
        "razor" => run_razor(
            n_proteins,
            &pep2prot_vec,
            &prot2pep_vec,
            &unique_peps,
            &best_scores,
        ),
        other => {
            eprintln!(
                "Unknown strategy: {}. Use bayesian, info_score, or razor.",
                other
            );
            std::process::exit(1);
        }
    };

    println!("  Done in {:.2}s", t4.elapsed().as_secs_f64());

    // Scale to [0, 1]
    let max_prob = probabilities
        .iter()
        .cloned()
        .fold(f64::NEG_INFINITY, f64::max);
    let scaled: Vec<f64> = if max_prob > 0.0 {
        probabilities.iter().map(|p| p / max_prob).collect()
    } else {
        probabilities.clone()
    };

    // ---------------------------------------------------------------
    // Step 5: Filter and write output
    // ---------------------------------------------------------------
    println!("\nStep 5: Filtering and writing output...");

    // Threshold analysis
    println!("  Threshold analysis:");
    for &thresh in &[0.0001, 0.001, 0.005, 0.01, 0.02, 0.05, 0.10, 0.25, 0.50] {
        let n = scaled.iter().filter(|&&p| p >= thresh).count();
        let marker = if (thresh - args.min_probability).abs() < 1e-6 {
            " <-- selected"
        } else {
            ""
        };
        println!("    P >= {:.4}: {:>10}{}", thresh, n, marker);
    }

    let selected: Vec<usize> = (0..n_proteins)
        .filter(|&i| scaled[i] >= args.min_probability)
        .collect();
    println!(
        "\n  Selected: {} / {} proteins ({:.1}x reduction)",
        selected.len(),
        n_proteins,
        n_proteins as f64 / selected.len().max(1) as f64
    );

    // Write output FASTA
    let mut writer = BufWriter::new(File::create(&args.output)?);
    for &i in &selected {
        writeln!(writer, "{}", proteins[i].header)?;
        for chunk in proteins[i].sequence.chunks(60) {
            writer.write_all(chunk)?;
            writeln!(writer)?;
        }
    }
    drop(writer);
    println!("  Wrote {}", args.output.display());

    // Write scores TSV
    if let Some(ref scores_path) = args.output_scores {
        let mut w = BufWriter::new(File::create(scores_path)?);
        writeln!(
            w,
            "protein_id\tprobability\tn_peptides\tn_unique\tn_samples"
        )?;
        let mut scored: Vec<(usize, f64)> = (0..n_proteins)
            .filter(|&i| scaled[i] > 0.0)
            .map(|i| (i, scaled[i]))
            .collect();
        scored.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
        for (i, prob) in scored {
            let peps = prot2pep_vec.get(&i).map(|v| v.len()).unwrap_or(0);
            let n_uniq = prot2pep_vec
                .get(&i)
                .map(|v| {
                    v.iter()
                        .filter(|&&pi| {
                            pep2prot_vec
                                .iter()
                                .find(|(k, _)| *k == pi)
                                .map(|(_, prots)| prots.len() == 1)
                                .unwrap_or(false)
                        })
                        .count()
                })
                .unwrap_or(0);
            let sc = prot_sample_count[i].powi(2) as usize; // undo sqrt
            writeln!(
                w,
                "{}\t{:.8}\t{}\t{}\t{}",
                proteins[i].id, prob, peps, n_uniq, sc
            )?;
        }
        println!("  Wrote {}", scores_path.display());
    }

    // Write stats JSON
    if let Some(ref stats_path) = args.output_stats {
        let stats = serde_json::json!({
            "strategy": args.strategy,
            "fasta_proteins": n_proteins,
            "mapped_proteins": n_mapped_prots,
            "mapped_peptides": n_mapped_peps,
            "unique_peptides": n_unique,
            "shared_peptides": n_shared,
            "n_samples": n_samples,
            "total_predictions": total_peptides,
            "min_probability": args.min_probability,
            "iterations": args.iterations,
            "selected_proteins": selected.len(),
            "reduction_factor": n_proteins as f64 / selected.len().max(1) as f64,
            "runtime_seconds": t_start.elapsed().as_secs_f64(),
        });
        let w = BufWriter::new(File::create(stats_path)?);
        serde_json::to_writer_pretty(w, &stats)?;
        println!("  Wrote {}", stats_path.display());
    }

    println!("\nTotal runtime: {:.2}s", t_start.elapsed().as_secs_f64());
    Ok(())
}

// ---------------------------------------------------------------
// Bayesian EM
// ---------------------------------------------------------------
fn run_bayesian_em(
    n_proteins: usize,
    pep2prot: &[(usize, Vec<usize>)],
    prot2pep: &HashMap<usize, Vec<usize>>,
    unique_peps: &[String],
    best_scores: &HashMap<String, f64>,
    sqrt_lengths: &[f64],
    sample_counts: &[f64],
    n_iterations: usize,
) -> Vec<f64> {
    // Pre-build pep_idx → n_proteins lookup for O(1) access
    let pep_n_prots: HashMap<usize, usize> = pep2prot
        .iter()
        .map(|(pi, prots)| (*pi, prots.len()))
        .collect();

    // Initialize: proportional to score-weighted evidence / sqrt(length) * sqrt(n_samples)
    let mut prob: Vec<f64> = (0..n_proteins)
        .into_par_iter()
        .map(|i| {
            let peps = match prot2pep.get(&i) {
                Some(p) => p,
                None => return 0.0,
            };
            let raw: f64 = peps
                .iter()
                .map(|&pi| {
                    let n_match = *pep_n_prots.get(&pi).unwrap_or(&1) as f64;
                    best_scores.get(&unique_peps[pi]).unwrap_or(&1.0) / n_match
                })
                .sum();
            raw * sample_counts[i] / sqrt_lengths[i]
        })
        .collect();

    normalize(&mut prob);

    for iter in 0..n_iterations {
        let t_iter = Instant::now();

        // E-step: compute new evidence per protein (parallel over peptides)
        let new_evidence: Vec<f64> = {
            let evidence = DashMap::with_capacity(n_proteins);

            pep2prot.par_iter().for_each(|(pep_i, prot_ids)| {
                let score = *best_scores.get(&unique_peps[*pep_i]).unwrap_or(&1.0);

                if prot_ids.len() == 1 {
                    // Unique peptide: full evidence, boosted by sample count
                    let pid = prot_ids[0];
                    let boost = sample_counts[pid].max(1.0);
                    *evidence.entry(pid).or_insert(0.0) += score * boost;
                } else {
                    // Shared peptide: distribute proportional to current probabilities
                    let prob_sum: f64 = prot_ids.iter().map(|&p| prob[p]).sum();
                    if prob_sum > 0.0 {
                        for &pid in prot_ids {
                            let weight = prob[pid] / prob_sum;
                            *evidence.entry(pid).or_insert(0.0) += score * weight;
                        }
                    } else {
                        let share = score / prot_ids.len() as f64;
                        for &pid in prot_ids {
                            *evidence.entry(pid).or_insert(0.0) += share;
                        }
                    }
                }
            });

            // Convert to Vec
            let mut ev = vec![0.0f64; n_proteins];
            for entry in evidence.iter() {
                ev[*entry.key()] = *entry.value();
            }
            ev
        };

        // M-step: update probabilities with length normalization + sample boost
        let old_prob = prob.clone();
        prob = (0..n_proteins)
            .into_par_iter()
            .map(|i| new_evidence[i] * sample_counts[i] / sqrt_lengths[i])
            .collect();
        normalize(&mut prob);

        // Convergence check: max delta AFTER normalization
        let max_delta = AtomicUsize::new(0);
        prob.par_iter().enumerate().for_each(|(i, &new_p)| {
            let delta_int = ((new_p - old_prob[i]).abs() * 1e15) as usize;
            max_delta.fetch_max(delta_int, Ordering::Relaxed);
        });
        let delta = max_delta.load(Ordering::Relaxed) as f64 / 1e15;
        println!(
            "  Iteration {}: max_delta={:.2e} ({:.2}s)",
            iter + 1,
            delta,
            t_iter.elapsed().as_secs_f64()
        );

        if delta < 1e-10 {
            println!("  Converged at iteration {}", iter + 1);
            break;
        }
    }

    prob
}

// ---------------------------------------------------------------
// Info Score (non-iterative)
// ---------------------------------------------------------------
fn run_info_score(
    n_proteins: usize,
    pep2prot: &[(usize, Vec<usize>)],
    prot2pep: &HashMap<usize, Vec<usize>>,
    unique_peps: &[String],
    best_scores: &HashMap<String, f64>,
    sqrt_lengths: &[f64],
) -> Vec<f64> {
    // Build pep_idx → n_proteins lookup
    let pep_n_prots: HashMap<usize, usize> = pep2prot
        .iter()
        .map(|(pi, prots)| (*pi, prots.len()))
        .collect();

    (0..n_proteins)
        .into_par_iter()
        .map(|i| {
            let peps = match prot2pep.get(&i) {
                Some(p) => p,
                None => return 0.0,
            };
            let info: f64 = peps
                .iter()
                .map(|&pi| {
                    let n = *pep_n_prots.get(&pi).unwrap_or(&1) as f64;
                    best_scores.get(&unique_peps[pi]).unwrap_or(&1.0) / n
                })
                .sum();
            info / sqrt_lengths[i]
        })
        .collect()
}

// ---------------------------------------------------------------
// Razor (assign shared to best protein)
// ---------------------------------------------------------------
fn run_razor(
    n_proteins: usize,
    pep2prot: &[(usize, Vec<usize>)],
    prot2pep: &HashMap<usize, Vec<usize>>,
    _unique_peps: &[String],
    _best_scores: &HashMap<String, f64>,
) -> Vec<f64> {
    // Count unique peptides per protein
    let unique_count: HashMap<usize, usize> = {
        let mut counts = HashMap::new();
        for (_, prot_ids) in pep2prot {
            if prot_ids.len() == 1 {
                *counts.entry(prot_ids[0]).or_insert(0usize) += 1;
            }
        }
        counts
    };

    // Assign each peptide to the protein with most unique peptides
    let assigned: DashMap<usize, bool> = DashMap::new();

    pep2prot.par_iter().for_each(|(_, prot_ids)| {
        if prot_ids.len() == 1 {
            assigned.insert(prot_ids[0], true);
        } else {
            // Pick protein with most unique peptides (tie-break: most total peptides, then ID)
            let best = prot_ids
                .iter()
                .max_by_key(|&&p| {
                    let u = *unique_count.get(&p).unwrap_or(&0);
                    let t = prot2pep.get(&p).map(|v| v.len()).unwrap_or(0);
                    (u, t)
                })
                .unwrap();
            assigned.insert(*best, true);
        }
    });

    // Binary output: 1.0 if assigned, 0.0 if not
    (0..n_proteins)
        .into_par_iter()
        .map(|i| if assigned.contains_key(&i) { 1.0 } else { 0.0 })
        .collect()
}

// ---------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------
fn normalize(v: &mut [f64]) {
    let sum: f64 = v.iter().sum();
    if sum > 0.0 {
        v.iter_mut().for_each(|x| *x /= sum);
    }
}

fn load_fasta(path: &std::path::Path) -> Result<Vec<Protein>, Box<dyn std::error::Error>> {
    let file = File::open(path)?;
    let reader = BufReader::with_capacity(64 * 1024 * 1024, file);
    let mut proteins = Vec::new();
    let mut current_header = String::new();
    let mut current_id = String::new();
    let mut current_seq = Vec::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('>') {
            if !current_id.is_empty() {
                proteins.push(Protein {
                    id: std::mem::take(&mut current_id),
                    sequence: std::mem::take(&mut current_seq),
                    header: std::mem::take(&mut current_header),
                });
            }
            current_header = line.clone();
            current_id = line[1..]
                .split_whitespace()
                .next()
                .unwrap_or("")
                .to_string();
        } else {
            current_seq.extend(line.trim().as_bytes());
        }
    }
    if !current_id.is_empty() {
        proteins.push(Protein {
            id: current_id,
            sequence: current_seq,
            header: current_header,
        });
    }
    Ok(proteins)
}

fn collect_csvs(dir: &std::path::Path) -> Vec<PathBuf> {
    let mut out = Vec::new();
    fn walk(dir: &std::path::Path, out: &mut Vec<PathBuf>) {
        if let Ok(entries) = std::fs::read_dir(dir) {
            for entry in entries.flatten() {
                let path = entry.path();
                if path.is_dir() {
                    walk(&path, out);
                } else if path.extension().map(|e| e == "csv").unwrap_or(false) {
                    out.push(path);
                }
            }
        }
    }
    walk(dir, &mut out);
    out.sort();
    out
}

fn parse_prediction_csv(
    path: &std::path::Path,
    min_length: usize,
    max_length: usize,
    normalize_il: bool,
    top_percent: f64,
    sample_idx: usize,
) -> Vec<Peptide> {
    let mut rdr = match csv::Reader::from_path(path) {
        Ok(r) => r,
        Err(_) => return Vec::new(),
    };

    let headers: Vec<String> = match rdr.headers() {
        Ok(h) => h.iter().map(|s| s.to_string()).collect(),
        Err(_) => return Vec::new(),
    };

    // Detect format
    let is_dianovo = headers.iter().any(|h| h == "pred_seq");
    let is_alphanovo = headers
        .iter()
        .any(|h| h == "peptide_prediction_detokenized_unmodified");

    let (seq_col, score_col) = if is_dianovo {
        ("pred_seq", "pred_prob")
    } else if is_alphanovo {
        ("peptide_prediction_detokenized_unmodified", "score")
    } else if headers.contains(&"sequence".to_string()) {
        ("sequence", "score")
    } else {
        return Vec::new();
    };

    let seq_idx = headers.iter().position(|h| h == seq_col);
    let score_idx = headers.iter().position(|h| h == score_col);

    if seq_idx.is_none() || score_idx.is_none() {
        return Vec::new();
    }
    let seq_idx = seq_idx.unwrap();
    let score_idx = score_idx.unwrap();

    let mut all_peps: Vec<(String, f64)> = Vec::new();

    for result in rdr.records() {
        let record = match result {
            Ok(r) => r,
            Err(_) => continue,
        };
        let seq_raw = record.get(seq_idx).unwrap_or("").to_string();
        let score_str = record.get(score_idx).unwrap_or("0");

        let (mut seq, score) = if is_dianovo {
            match dianovo::parse_prediction(&seq_raw, score_str) {
                Some(value) => value,
                None => continue,
            }
        } else {
            (
                seq_raw
                    .chars()
                    .filter(|c| c.is_ascii_uppercase())
                    .collect::<String>(),
                score_str.parse::<f64>().unwrap_or(0.0),
            )
        };

        if normalize_il {
            seq = seq.replace("I", "L");
        }

        if seq.len() >= min_length && seq.len() <= max_length && score > 0.0 {
            all_peps.push((seq, score));
        }
    }

    // Top percent filter
    if top_percent < 1.0 && !all_peps.is_empty() {
        let mut scores: Vec<f64> = all_peps.iter().map(|(_, s)| *s).collect();
        scores.sort_by(|a, b| b.partial_cmp(a).unwrap());
        let cutoff_idx = ((scores.len() as f64 * top_percent) as usize)
            .max(1)
            .min(scores.len() - 1);
        let threshold = scores[cutoff_idx];
        all_peps.retain(|(_, s)| *s >= threshold);
    }

    all_peps
        .into_iter()
        .map(|(seq, score)| Peptide {
            sequence: seq,
            score,
            sample_idx,
        })
        .collect()
}
