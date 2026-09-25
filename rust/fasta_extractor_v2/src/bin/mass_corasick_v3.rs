//! mass_corasick_v3: Per-sample parallel mass k-mer rescue
//!
//! Three modes:
//!   prep:   Aho-Corasick exact match → matched peptide set + sample manifest
//!   rescue: Per-sample k-mer rescue (SLURM array, one task per sample)
//!   merge:  Union all rescued proteins → enhanced evidence lake
//!
//! With 306 samples × 273M lake: 45 minutes total (30 min wall clock for rescue array)

// Legacy/experimental reference binary — not on the canonical pipeline; clippy relaxed.
#![allow(clippy::all)]
#![allow(dead_code, unused)]

use clap::{Parser, Subcommand};
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write as IoWrite};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Instant;

// ============================================================
// Shared constants and utilities
// ============================================================

#[inline(always)]
fn aa_mass(aa: u8) -> f64 {
    const TABLE: [f64; 26] = [
        71.03711, 0.0, 103.00919, 115.02694, 129.04259, 147.06841, 57.02146, 137.05891, 113.08406,
        0.0, 128.09496, 113.08406, 131.04049, 114.04293, 0.0, 97.05276, 128.05858, 156.10111,
        87.03203, 101.04768, 0.0, 99.06841, 186.07931, 0.0, 163.06333, 0.0,
    ];
    if aa >= b'A' && aa <= b'Z' {
        TABLE[(aa - b'A') as usize]
    } else {
        0.0
    }
}

const K: usize = 5;
const BIN_WIDTH: f64 = 0.001;
const MIN_KMER_MASS: f64 = 280.0;
const MAX_KMER_MASS: f64 = 940.0;
const N_BINS: usize = ((MAX_KMER_MASS - MIN_KMER_MASS) / BIN_WIDTH) as usize + 100;

#[inline(always)]
fn mass_to_bin(mass: f64) -> usize {
    let raw = ((mass - MIN_KMER_MASS) / BIN_WIDTH) as i64;
    raw.max(0).min(N_BINS as i64 - 1) as usize
}

fn is_valid_peptide(seq: &[u8]) -> bool {
    if seq.len() < 8 || seq.len() > 50 {
        return false;
    }
    if !seq.iter().all(|&b| aa_mass(b) > 0.0) {
        return false;
    }
    // Low-complexity filter
    let mut counts = [0u32; 26];
    for &b in seq {
        if b >= b'A' && b <= b'Z' {
            counts[(b - b'A') as usize] += 1;
        }
    }
    let max_count = *counts.iter().max().unwrap_or(&0) as usize;
    max_count * 2 <= seq.len()
}

fn compute_kmer_bins(seq: &[u8]) -> Vec<usize> {
    if seq.len() < K {
        return vec![];
    }
    let mut cum = Vec::with_capacity(seq.len() + 1);
    cum.push(0.0f64);
    let mut r = 0.0;
    for &aa in seq {
        r += aa_mass(aa);
        cum.push(r);
    }
    (0..=(seq.len() - K))
        .map(|i| mass_to_bin(cum[i + K] - cum[i]))
        .collect()
}

// ============================================================
// Peptide index (small, per-sample)
// ============================================================

struct PeptideQuery {
    sequence: Vec<u8>,
    kmer_bins: Vec<usize>,
    score: f64,
}

struct PeptideIndex {
    queries: Vec<PeptideQuery>,
    bin_start: Vec<u32>,
    bin_count: Vec<u32>,
    entries: Vec<u64>, // packed: (peptide_idx:u32 << 16) | kmer_pos:u16
    total_entries: usize,
}

impl PeptideIndex {
    fn build(queries: Vec<PeptideQuery>) -> Self {
        let mut counts = vec![0u32; N_BINS];
        for query in &queries {
            for &bin in &query.kmer_bins {
                if bin > 0 {
                    counts[bin - 1] += 1;
                }
                counts[bin] += 1;
                if bin + 1 < N_BINS {
                    counts[bin + 1] += 1;
                }
            }
        }

        let total: usize = counts.iter().map(|&c| c as usize).sum();
        let mut bin_start = vec![0u32; N_BINS];
        let mut running = 0u32;
        for i in 0..N_BINS {
            bin_start[i] = running;
            running += counts[i];
        }

        let mut entries = vec![0u64; total];
        let mut write_pos = vec![0u32; N_BINS];

        for (qidx, query) in queries.iter().enumerate() {
            for (kpos, &bin) in query.kmer_bins.iter().enumerate() {
                let packed = ((qidx as u64) << 16) | (kpos as u64);
                for target in [bin.wrapping_sub(1), bin, bin + 1] {
                    if target < N_BINS {
                        let offset = bin_start[target] as usize + write_pos[target] as usize;
                        entries[offset] = packed;
                        write_pos[target] += 1;
                    }
                }
            }
        }

        let non_empty = counts.iter().filter(|&&c| c > 0).count();
        let max_bin = *counts.iter().max().unwrap_or(&0);
        eprintln!(
            "  Index: {} peptides, {} entries, {} non-empty bins (max {}), {:.1} MB",
            queries.len(),
            total,
            non_empty,
            max_bin,
            (total * 8) as f64 / 1e6
        );

        PeptideIndex {
            queries,
            bin_start,
            bin_count: write_pos,
            entries,
            total_entries: total,
        }
    }

    #[inline(always)]
    fn get_bin(&self, bin: usize) -> &[u64] {
        if bin >= N_BINS {
            return &[];
        }
        let s = self.bin_start[bin] as usize;
        let n = self.bin_count[bin] as usize;
        &self.entries[s..s + n]
    }

    fn sweep_protein(&self, prot_seq: &[u8], min_votes: usize) -> Vec<(u32, u16)> {
        if prot_seq.len() < K {
            return vec![];
        }

        let mut cum = Vec::with_capacity(prot_seq.len() + 1);
        cum.push(0.0f64);
        let mut r = 0.0;
        for &aa in prot_seq {
            r += aa_mass(aa);
            cum.push(r);
        }

        let mut votes: HashMap<u64, u16> = HashMap::with_capacity(64);

        for prot_pos in 0..=(prot_seq.len() - K) {
            let mass = cum[prot_pos + K] - cum[prot_pos];
            let bin = mass_to_bin(mass);
            let entries = self.get_bin(bin);

            for &packed in entries {
                let qidx = (packed >> 16) as u32;
                let kpos = (packed & 0xFFFF) as i32;
                let inferred_start = prot_pos as i32 - kpos;
                if inferred_start < 0 {
                    continue;
                }
                let key = ((qidx as u64) << 32) | (inferred_start as u64);
                let v = votes.entry(key).or_insert(0);
                *v = v.saturating_add(1);
            }
        }

        let mut results = Vec::new();
        for (&key, &vote_count) in &votes {
            if (vote_count as usize) >= min_votes {
                let qidx = (key >> 32) as u32;
                results.push((qidx, vote_count));
            }
        }
        results
    }
}

// ============================================================
// CSV loading
// ============================================================

fn load_csv_peptides(path: &Path) -> Vec<(Vec<u8>, f64)> {
    let mut result = Vec::new();
    let file = match File::open(path) {
        Ok(f) => f,
        Err(_) => return result,
    };
    let mut reader = csv::Reader::from_reader(BufReader::new(file));
    let headers = match reader.headers() {
        Ok(h) => h.clone(),
        Err(_) => return result,
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
    let seq_idx = match seq_idx {
        Some(i) => i,
        None => return result,
    };

    for rec in reader.records() {
        let rec = match rec {
            Ok(r) => r,
            Err(_) => continue,
        };
        let seq: Vec<u8> = rec
            .get(seq_idx)
            .unwrap_or("")
            .trim()
            .to_uppercase()
            .bytes()
            .collect();
        if !is_valid_peptide(&seq) {
            continue;
        }
        let score: f64 = score_idx
            .and_then(|i| rec.get(i))
            .and_then(|s| s.parse().ok())
            .unwrap_or(0.0);
        result.push((seq, score));
    }
    result
}

fn find_csv_files(dir: &str) -> Vec<PathBuf> {
    let pattern = format!("{}/**/*predictions.csv", dir);
    let mut paths: Vec<_> = glob::glob(&pattern)
        .expect("glob")
        .filter_map(|p| p.ok())
        .collect();
    if paths.is_empty() {
        let p2 = format!("{}/**/*.csv", dir);
        paths = glob::glob(&p2)
            .expect("glob")
            .filter_map(|p| p.ok())
            .collect();
    }
    if paths.is_empty() {
        let p = PathBuf::from(dir);
        if p.is_file() {
            paths.push(p);
        }
    }
    paths.sort();
    paths
}

fn load_proteins(path: &str) -> Vec<(String, Vec<u8>)> {
    let t0 = Instant::now();
    let file = File::open(path).unwrap_or_else(|e| panic!("{}: {}", path, e));
    let reader = BufReader::with_capacity(16 * 1024 * 1024, file);
    let mut proteins = Vec::new();
    let mut cid = String::new();
    let mut cseq: Vec<u8> = Vec::new();

    for line in reader.lines() {
        let line = line.expect("read");
        if line.starts_with('>') {
            if !cid.is_empty() && !cseq.is_empty() {
                proteins.push((std::mem::take(&mut cid), std::mem::take(&mut cseq)));
            }
            cid = line[1..]
                .split_whitespace()
                .next()
                .unwrap_or("")
                .to_string();
            cseq.clear();
        } else {
            cseq.extend(
                line.trim()
                    .as_bytes()
                    .iter()
                    .map(|b| b.to_ascii_uppercase()),
            );
        }
        if proteins.len() % 10_000_000 == 0 && !proteins.is_empty() && cseq.is_empty() {
            eprintln!("  {}M proteins...", proteins.len() / 1_000_000);
        }
    }
    if !cid.is_empty() && !cseq.is_empty() {
        proteins.push((cid, cseq));
    }
    eprintln!(
        "Loaded {} proteins in {:.1}s",
        proteins.len(),
        t0.elapsed().as_secs_f64()
    );
    proteins
}

// ============================================================
// Mode: prep
// ============================================================

fn run_prep(database: &str, predictions: &str, output: &str, top_percent: f64, threads: usize) {
    let t0 = Instant::now();
    let outdir = Path::new(output);
    std::fs::create_dir_all(outdir).ok();

    eprintln!("=== mass_corasick v3: PREP ===\n");

    // Find all CSVs
    let csv_files = find_csv_files(predictions);
    eprintln!("Found {} CSV files", csv_files.len());

    // Write sample manifest
    let manifest_path = outdir.join("sample_manifest.tsv");
    {
        let mut mf = BufWriter::new(File::create(&manifest_path).expect("manifest"));
        writeln!(mf, "sample_id\tcsv_path\tpeptide_count").unwrap();
        for (i, path) in csv_files.iter().enumerate() {
            let peps = load_csv_peptides(path);
            writeln!(mf, "{}\t{}\t{}", i, path.display(), peps.len()).unwrap();
            if (i + 1) % 50 == 0 {
                eprintln!("  Scanned {}/{} CSVs...", i + 1, csv_files.len());
            }
        }
    }
    eprintln!(
        "Manifest: {} samples → {}",
        csv_files.len(),
        manifest_path.display()
    );

    // Pool all peptides, deduplicate, apply score filter
    eprintln!(
        "\nPooling peptides (top {:.0}% per sample)...",
        top_percent * 100.0
    );
    let mut all_peptides: HashMap<Vec<u8>, f64> = HashMap::new();

    for path in &csv_files {
        let mut peps = load_csv_peptides(path);
        // Top percent by score
        peps.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
        let keep = (peps.len() as f64 * top_percent).ceil() as usize;
        peps.truncate(keep);

        for (seq, score) in peps {
            let entry = all_peptides.entry(seq).or_insert(0.0);
            if score > *entry {
                *entry = score;
            }
        }
    }
    eprintln!("Total unique peptides: {}", all_peptides.len());

    // Run Aho-Corasick exact match
    eprintln!("\nAho-Corasick exact match against {}...", database);
    let proteins = load_proteins(database);

    let peptide_list: Vec<Vec<u8>> = all_peptides.keys().cloned().collect();

    // Normalize I→L for matching
    let norm_peptides: Vec<Vec<u8>> = peptide_list
        .iter()
        .map(|p| {
            p.iter()
                .map(|&b| if b == b'I' { b'L' } else { b })
                .collect()
        })
        .collect();

    let ac = aho_corasick::AhoCorasick::builder()
        .build(&norm_peptides)
        .expect("AC build");

    let matched_flags: Vec<std::sync::atomic::AtomicBool> = (0..peptide_list.len())
        .map(|_| std::sync::atomic::AtomicBool::new(false))
        .collect();

    let match_t0 = Instant::now();
    let evidence_ids: dashmap::DashSet<String> = dashmap::DashSet::new();

    proteins.par_iter().for_each(|(pid, pseq)| {
        let norm_seq: Vec<u8> = pseq
            .iter()
            .map(|&b| if b == b'I' { b'L' } else { b })
            .collect();
        let mut found = false;
        for mat in ac.find_overlapping_iter(&norm_seq) {
            matched_flags[mat.pattern().as_usize()].store(true, Ordering::Relaxed);
            found = true;
        }
        if found {
            evidence_ids.insert(pid.clone());
        }
    });

    let matched_set: HashSet<Vec<u8>> = peptide_list
        .iter()
        .enumerate()
        .filter(|(i, _)| matched_flags[*i].load(Ordering::Relaxed))
        .map(|(_, seq)| seq.clone())
        .collect();

    eprintln!(
        "Exact match: {} peptides matched, {} proteins, {:.1}s",
        matched_set.len(),
        evidence_ids.len(),
        match_t0.elapsed().as_secs_f64()
    );

    // Write exact evidence FASTA
    let evidence_path = outdir.join("exact_evidence.fasta");
    {
        let mut ef = BufWriter::new(File::create(&evidence_path).expect("evidence fasta"));
        for (pid, pseq) in &proteins {
            if evidence_ids.contains(pid) {
                writeln!(ef, ">{}", pid).unwrap();
                for chunk in pseq.chunks(60) {
                    ef.write_all(chunk).unwrap();
                    ef.write_all(b"\n").unwrap();
                }
            }
        }
    }
    eprintln!(
        "Evidence: {} → {}",
        evidence_ids.len(),
        evidence_path.display()
    );

    // Write matched peptides (binary: sorted list of sequences)
    let matched_path = outdir.join("matched_peptides.tsv");
    {
        let mut mf = BufWriter::new(File::create(&matched_path).expect("matched"));
        for seq in &matched_set {
            mf.write_all(seq).unwrap();
            mf.write_all(b"\n").unwrap();
        }
    }
    eprintln!(
        "Matched peptides: {} → {}",
        matched_set.len(),
        matched_path.display()
    );

    // Stats
    let stats = serde_json::json!({
        "total_peptides": all_peptides.len(),
        "matched_peptides": matched_set.len(),
        "evidence_proteins": evidence_ids.len(),
        "database_proteins": proteins.len(),
        "csv_files": csv_files.len(),
        "top_percent": top_percent,
        "time_s": t0.elapsed().as_secs_f64(),
    });
    let stats_path = outdir.join("prep_stats.json");
    std::fs::write(&stats_path, serde_json::to_string_pretty(&stats).unwrap()).ok();
    eprintln!("\nPrep complete in {:.1}s", t0.elapsed().as_secs_f64());
}

// ============================================================
// Mode: rescue
// ============================================================

fn run_rescue(database: &str, output: &str, sample_id: usize, min_votes: usize, threads: usize) {
    let t0 = Instant::now();
    let outdir = Path::new(output);

    eprintln!("=== mass_corasick v3: RESCUE sample {} ===\n", sample_id);

    // Load manifest
    let manifest_path = outdir.join("sample_manifest.tsv");
    let manifest = BufReader::new(File::open(&manifest_path).expect("manifest"));
    let mut csv_path = String::new();
    for line in manifest.lines().skip(1) {
        let line = line.expect("manifest line");
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() >= 2 {
            let sid: usize = parts[0].parse().unwrap_or(usize::MAX);
            if sid == sample_id {
                csv_path = parts[1].to_string();
                break;
            }
        }
    }
    if csv_path.is_empty() {
        eprintln!("Sample {} not found in manifest", sample_id);
        return;
    }
    eprintln!("Sample {}: {}", sample_id, csv_path);

    // Load matched peptides
    let matched_path = outdir.join("matched_peptides.tsv");
    let mut matched: HashSet<Vec<u8>> = HashSet::new();
    if let Ok(f) = File::open(&matched_path) {
        for line in BufReader::new(f).lines() {
            if let Ok(seq) = line {
                matched.insert(seq.trim().as_bytes().to_vec());
            }
        }
    }
    eprintln!("Loaded {} matched peptides to exclude", matched.len());

    // Load this sample's peptides (top 30% by score, matching prep filter)
    let mut peps = load_csv_peptides(Path::new(&csv_path));
    peps.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
    let keep = (peps.len() as f64 * 0.30).ceil() as usize;
    peps.truncate(keep);
    eprintln!(
        "Loaded {} peptides (top 30% of {})",
        peps.len(),
        peps.len() * 10 / 3
    );

    let mut queries: Vec<PeptideQuery> = Vec::new();
    let mut seen: HashSet<Vec<u8>> = HashSet::new();

    for (seq, score) in peps {
        if matched.contains(&seq) {
            continue;
        }
        if seen.contains(&seq) {
            continue;
        }
        seen.insert(seq.clone());
        let kmer_bins = compute_kmer_bins(&seq);
        if kmer_bins.is_empty() {
            continue;
        }
        queries.push(PeptideQuery {
            sequence: seq,
            kmer_bins,
            score,
        });
    }
    eprintln!("Failed peptides for rescue: {}", queries.len());

    if queries.is_empty() {
        eprintln!("No peptides to rescue");
        // Write empty output
        let out_path = outdir.join(format!("rescued_sample_{:04}.tsv", sample_id));
        let mut w = BufWriter::new(File::create(&out_path).expect("output"));
        writeln!(w, "query_sequence\tprotein_id\tvotes\tn_kmers\tscore").unwrap();
        return;
    }

    // Build peptide index
    let index = PeptideIndex::build(queries);

    // Load proteins
    let proteins = load_proteins(database);

    // Sweep
    eprintln!("\nSweeping {} proteins...", proteins.len());
    let sweep_t0 = Instant::now();
    let processed = AtomicUsize::new(0);

    let chunk_size = 5000;
    let all_results: Vec<Vec<(u32, u16, String)>> = proteins
        .par_chunks(chunk_size)
        .map(|chunk| {
            let mut local = Vec::new();
            for (pid, pseq) in chunk {
                let hits = index.sweep_protein(pseq, min_votes);
                for (qidx, votes) in hits {
                    local.push((qidx, votes, pid.clone()));
                }
            }
            let p = processed.fetch_add(chunk.len(), Ordering::Relaxed) + chunk.len();
            if p % 1_000_000 < chunk_size {
                let elapsed = sweep_t0.elapsed().as_secs_f64();
                let rate = p as f64 / elapsed;
                let eta = (proteins.len() - p) as f64 / rate;
                eprintln!(
                    "  {}/{}  {:.0} prot/s  ETA {:.0}m  rescued: {}",
                    p,
                    proteins.len(),
                    rate,
                    eta / 60.0,
                    local.len()
                );
            }
            local
        })
        .collect();

    let sweep_time = sweep_t0.elapsed();

    // Merge best per peptide
    let mut best: HashMap<u32, (u16, String)> = HashMap::new();
    for chunk in &all_results {
        for (qidx, votes, pid) in chunk {
            let e = best.entry(*qidx).or_insert((0, String::new()));
            if *votes > e.0 {
                *e = (*votes, pid.clone());
            }
        }
    }

    eprintln!(
        "\nSweep: {} proteins in {:.1}s ({:.0} prot/s)",
        proteins.len(),
        sweep_time.as_secs_f64(),
        proteins.len() as f64 / sweep_time.as_secs_f64()
    );
    eprintln!("Rescued: {} peptides", best.len());

    // Write output
    let out_path = outdir.join(format!("rescued_sample_{:04}.tsv", sample_id));
    let mut w = BufWriter::new(File::create(&out_path).expect("output"));
    writeln!(w, "query_sequence\tprotein_id\tvotes\tn_kmers\tscore").unwrap();

    let mut sorted: Vec<_> = best.iter().collect();
    sorted.sort_by(|a, b| b.1 .0.cmp(&a.1 .0));
    for (qidx, (votes, pid)) in &sorted {
        let q = &index.queries[**qidx as usize];
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{:.4}",
            String::from_utf8_lossy(&q.sequence),
            pid,
            votes,
            q.kmer_bins.len(),
            q.score
        )
        .unwrap();
    }

    let stats = serde_json::json!({
        "sample_id": sample_id,
        "csv_path": csv_path,
        "total_peptides_in_csv": seen.len() + matched.len(),
        "failed_peptides": index.queries.len(),
        "rescued_peptides": best.len(),
        "rescued_proteins": best.values().map(|(_, p)| p.as_str()).collect::<HashSet<_>>().len(),
        "sweep_rate": proteins.len() as f64 / sweep_time.as_secs_f64(),
        "total_time_s": t0.elapsed().as_secs_f64(),
    });
    let stats_path = outdir.join(format!("rescued_sample_{:04}_stats.json", sample_id));
    std::fs::write(&stats_path, serde_json::to_string_pretty(&stats).unwrap()).ok();

    eprintln!("Output: {}", out_path.display());
    eprintln!("Total: {:.1}s", t0.elapsed().as_secs_f64());
}

// ============================================================
// Mode: merge
// ============================================================

fn run_merge(database: &str, output: &str, _threads: usize) {
    let t0 = Instant::now();
    let outdir = Path::new(output);

    eprintln!("=== mass_corasick v3: MERGE ===\n");

    // Collect all rescued TSVs
    let pattern = format!("{}/rescued_sample_*.tsv", outdir.display());
    let paths: Vec<_> = glob::glob(&pattern)
        .expect("glob")
        .filter_map(|p| p.ok())
        .collect();
    eprintln!("Found {} sample rescue files", paths.len());

    // Union: best hit per peptide across all samples
    let mut best_per_peptide: HashMap<String, (u16, String, f64)> = HashMap::new();
    let mut per_sample_stats: Vec<(usize, usize)> = Vec::new();

    for path in &paths {
        let file = File::open(path).expect("open rescued tsv");
        let mut reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .from_reader(BufReader::new(file));
        let mut count = 0;
        for rec in reader.records() {
            let rec = match rec {
                Ok(r) => r,
                Err(_) => continue,
            };
            let seq = rec.get(0).unwrap_or("").to_string();
            let pid = rec.get(1).unwrap_or("").to_string();
            let votes: u16 = rec.get(2).and_then(|s| s.parse().ok()).unwrap_or(0);
            let score: f64 = rec.get(4).and_then(|s| s.parse().ok()).unwrap_or(0.0);

            let e = best_per_peptide
                .entry(seq)
                .or_insert((0, String::new(), 0.0));
            if votes > e.0 {
                *e = (votes, pid, score);
            }
            count += 1;
        }
        // Extract sample id from filename
        let sample_num = path
            .file_stem()
            .and_then(|s| s.to_str())
            .and_then(|s| s.strip_prefix("rescued_sample_"))
            .and_then(|s| s.parse::<usize>().ok())
            .unwrap_or(0);
        per_sample_stats.push((sample_num, count));
    }

    let rescued_proteins: HashSet<&str> = best_per_peptide
        .values()
        .map(|(_, pid, _)| pid.as_str())
        .collect();

    eprintln!(
        "Merged: {} rescued peptides → {} unique proteins",
        best_per_peptide.len(),
        rescued_proteins.len()
    );

    // Write merged rescued TSV
    let merged_path = outdir.join("kmer_rescued_merged.tsv");
    {
        let mut w = BufWriter::new(File::create(&merged_path).expect("merged tsv"));
        writeln!(w, "query_sequence\tprotein_id\tvotes\tscore").unwrap();
        let mut sorted: Vec<_> = best_per_peptide.iter().collect();
        sorted.sort_by(|a, b| (b.1).0.cmp(&(a.1).0));
        for (seq, (votes, pid, score)) in &sorted {
            writeln!(w, "{}\t{}\t{}\t{:.4}", seq, pid, votes, score).unwrap();
        }
    }

    // Load exact evidence protein IDs
    let exact_path = outdir.join("exact_evidence.fasta");
    let mut exact_ids: HashSet<String> = HashSet::new();
    if let Ok(f) = File::open(&exact_path) {
        for line in BufReader::new(f).lines() {
            if let Ok(l) = line {
                if l.starts_with('>') {
                    exact_ids.insert(l[1..].split_whitespace().next().unwrap_or("").to_string());
                }
            }
        }
    }

    let new_proteins: HashSet<&str> = rescued_proteins
        .difference(&exact_ids.iter().map(|s| s.as_str()).collect())
        .cloned()
        .collect();

    eprintln!(
        "Already in evidence: {}",
        rescued_proteins.len() - new_proteins.len()
    );
    eprintln!("New proteins to add: {}", new_proteins.len());

    // Build enhanced evidence: copy exact + extract new from lake
    let enhanced_path = outdir.join("enhanced_evidence.fasta");
    std::fs::copy(&exact_path, &enhanced_path).expect("copy exact evidence");

    let mut added = 0;
    let lake_file = File::open(database).expect("lake");
    let lake_reader = BufReader::with_capacity(16 * 1024 * 1024, lake_file);
    let mut enhanced = std::fs::OpenOptions::new()
        .append(true)
        .open(&enhanced_path)
        .expect("append");
    let mut writing = false;

    for line in lake_reader.lines() {
        let line = line.expect("read");
        if line.starts_with('>') {
            let pid = line[1..].split_whitespace().next().unwrap_or("");
            writing = new_proteins.contains(pid);
            if writing {
                added += 1;
            }
        }
        if writing {
            enhanced.write_all(line.as_bytes()).unwrap();
            enhanced.write_all(b"\n").unwrap();
        }
    }

    eprintln!("Added {} new proteins to evidence lake", added);

    // Summary
    let summary = serde_json::json!({
        "exact_evidence_proteins": exact_ids.len(),
        "rescued_peptides": best_per_peptide.len(),
        "rescued_proteins": rescued_proteins.len(),
        "new_proteins_from_kmer": added,
        "enhanced_evidence_proteins": exact_ids.len() + added,
        "improvement_pct": if exact_ids.len() > 0 { 100.0 * added as f64 / exact_ids.len() as f64 } else { 0.0 },
        "samples_processed": paths.len(),
        "total_time_s": t0.elapsed().as_secs_f64(),
    });
    let summary_path = outdir.join("summary.json");
    std::fs::write(
        &summary_path,
        serde_json::to_string_pretty(&summary).unwrap(),
    )
    .ok();

    eprintln!("\n{}", serde_json::to_string_pretty(&summary).unwrap());
    eprintln!("\nMerge complete in {:.1}s", t0.elapsed().as_secs_f64());
}

// ============================================================
// CLI
// ============================================================

#[derive(Parser)]
#[command(name = "mass_corasick_v3")]
#[command(about = "Per-sample parallel mass k-mer rescue (prep → rescue array → merge)")]
struct Args {
    #[command(subcommand)]
    mode: Mode,
}

#[derive(Subcommand)]
enum Mode {
    /// Step 1: Exact match + build manifest + matched peptide set
    Prep {
        #[arg(short = 'd', long)]
        database: String,
        #[arg(short = 'p', long)]
        predictions: String,
        #[arg(short = 'o', long)]
        output: String,
        #[arg(long, default_value = "0.30")]
        top_percent: f64,
        #[arg(short = 't', long, default_value = "0")]
        threads: usize,
    },
    /// Step 2: Per-sample k-mer rescue (run as SLURM array)
    Rescue {
        #[arg(short = 'd', long)]
        database: String,
        #[arg(short = 'o', long)]
        output: String,
        #[arg(long)]
        sample_id: usize,
        #[arg(long, default_value = "5")]
        min_votes: usize,
        #[arg(short = 't', long, default_value = "0")]
        threads: usize,
    },
    /// Step 3: Merge all sample rescues → enhanced evidence lake
    Merge {
        #[arg(short = 'd', long)]
        database: String,
        #[arg(short = 'o', long)]
        output: String,
        #[arg(short = 't', long, default_value = "0")]
        threads: usize,
    },
}

fn main() {
    let args = Args::parse();

    match args.mode {
        Mode::Prep {
            database,
            predictions,
            output,
            top_percent,
            threads,
        } => {
            if threads > 0 {
                rayon::ThreadPoolBuilder::new()
                    .num_threads(threads)
                    .build_global()
                    .ok();
            }
            run_prep(&database, &predictions, &output, top_percent, threads);
        }
        Mode::Rescue {
            database,
            output,
            sample_id,
            min_votes,
            threads,
        } => {
            if threads > 0 {
                rayon::ThreadPoolBuilder::new()
                    .num_threads(threads)
                    .build_global()
                    .ok();
            }
            run_rescue(&database, &output, sample_id, min_votes, threads);
        }
        Mode::Merge {
            database,
            output,
            threads,
        } => {
            run_merge(&database, &output, threads);
        }
    }
}
