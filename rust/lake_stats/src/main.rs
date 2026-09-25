//! lake_stats — Sequence-space saturation (D/K) of tryptic peptides in a FASTA.
//!
//! Mirrors the Python reference (`sequence_space_saturation.py`) but in
//! parallel Rust for ~50-100× speedup.
//!
//! Tryptic cleavage:  cut after K/R unless followed by P.
//! Filter:            drop peptides containing B,J,O,U,X,Z,*,- (ambiguous/non-standard).
//! Output:            TSV with database, length, K=20^N, D=unique peptides, D/K%.

use ahash::AHashSet;
use clap::Parser;
use rayon::prelude::*;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::PathBuf;
use std::sync::Mutex;
use std::time::Instant;

const CHUNK_SIZE: usize = 50_000;

const VERSION: &str = concat!(env!("CARGO_PKG_VERSION"), "+", env!("GIT_REV"));

#[derive(Parser)]
#[command(version = VERSION, about = "Compute D/K saturation for tryptic peptides per length")]
struct Args {
    #[arg(short, long)]
    fasta: PathBuf,
    #[arg(short, long)]
    output: PathBuf,
    /// Lengths to count, comma-separated
    #[arg(short, long, default_value = "6,7,8,9,10,11,12")]
    lengths: String,
    /// Database name for output
    #[arg(short, long, default_value = "lake")]
    name: String,
    /// Threads (0 = use all)
    #[arg(short, long, default_value_t = 0)]
    threads: usize,
}

#[inline]
fn is_clean(pep: &[u8]) -> bool {
    for &b in pep {
        match b {
            b'A' | b'C' | b'D' | b'E' | b'F' | b'G' | b'H' | b'I' | b'K' | b'L' | b'M' | b'N'
            | b'P' | b'Q' | b'R' | b'S' | b'T' | b'V' | b'W' | b'Y' => {}
            _ => return false,
        }
    }
    true
}

/// Yield tryptic peptides from a protein sequence as byte slices.
/// Cleave after K/R unless followed by P. Last peptide is always emitted.
fn trypsin_into(seq: &[u8], lengths: &[usize], local: &mut [AHashSet<Vec<u8>>]) {
    let n = seq.len();
    if n == 0 {
        return;
    }
    let mut last = 0usize;
    let emit = |pep: &[u8], local: &mut [AHashSet<Vec<u8>>]| {
        if !is_clean(pep) {
            return;
        }
        let len = pep.len();
        for (idx, &target_len) in lengths.iter().enumerate() {
            if len == target_len {
                local[idx].insert(pep.to_vec());
                return;
            }
        }
    };
    for i in 0..n.saturating_sub(1) {
        let c = seq[i];
        if (c == b'K' || c == b'R') && seq[i + 1] != b'P' {
            emit(&seq[last..=i], local);
            last = i + 1;
        }
    }
    if last < n {
        emit(&seq[last..], local);
    }
}

/// Stream FASTA, yielding owned protein byte vecs (whitespace/digits stripped,
/// letters uppercased). `*` (stop codon) and `-` (gap) are PRESERVED so that a
/// tryptic peptide spanning them reaches `is_clean` and is dropped — never
/// silently concatenated across the boundary (see spec in module docstring).
struct FastaReader<R: BufRead> {
    reader: R,
    next_header: Option<Vec<u8>>,
}

impl<R: BufRead> FastaReader<R> {
    fn new(reader: R) -> Self {
        Self {
            reader,
            next_header: None,
        }
    }
}

impl<R: BufRead> Iterator for FastaReader<R> {
    type Item = Vec<u8>;
    fn next(&mut self) -> Option<Self::Item> {
        let mut seq: Vec<u8> = Vec::with_capacity(512);
        let mut line = Vec::with_capacity(256);
        // If we left a header from previous call, we're already inside a record
        let mut have_record = self.next_header.take().is_some();
        loop {
            line.clear();
            match self.reader.read_until(b'\n', &mut line) {
                Ok(0) => {
                    if have_record && !seq.is_empty() {
                        return Some(seq);
                    }
                    return None;
                }
                Ok(_) => {}
                Err(_) => return None,
            }
            // Strip trailing \n and \r
            while matches!(line.last(), Some(b'\n') | Some(b'\r')) {
                line.pop();
            }
            if line.is_empty() {
                continue;
            }
            if line[0] == b'>' {
                if have_record {
                    self.next_header = Some(line.clone());
                    return Some(seq);
                }
                have_record = true;
            } else {
                // Append uppercased letters; preserve `*`/`-` so peptides
                // spanning a stop codon or gap are dropped (not concatenated)
                // by `is_clean`. Whitespace, digits and other punctuation are
                // stripped as before.
                for &b in &line {
                    if b.is_ascii_alphabetic() {
                        seq.push(b.to_ascii_uppercase());
                    } else if b == b'*' || b == b'-' {
                        seq.push(b);
                    }
                }
            }
        }
    }
}

fn main() {
    let args = Args::parse();
    let lengths: Vec<usize> = args
        .lengths
        .split(',')
        .map(|s| s.trim().parse::<usize>().expect("lengths must be integers"))
        .collect();

    if args.threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global()
            .unwrap();
    }
    let nthreads = rayon::current_num_threads();
    eprintln!(
        "lake_stats — input: {}, lengths: {:?}, threads: {}",
        args.fasta.display(),
        lengths,
        nthreads
    );

    let file = File::open(&args.fasta).expect("Cannot open FASTA");
    let reader = BufReader::with_capacity(1 << 20, file);
    let fasta = FastaReader::new(reader);

    let global: Vec<Mutex<AHashSet<Vec<u8>>>> = lengths
        .iter()
        .map(|_| Mutex::new(AHashSet::new()))
        .collect();

    let t0 = Instant::now();
    let mut total_proteins: u64 = 0;
    let mut chunk: Vec<Vec<u8>> = Vec::with_capacity(CHUNK_SIZE);

    let process_chunk = |chunk: &[Vec<u8>]| {
        // Each rayon worker builds local sets, then merges into global with one lock.
        // Splitting chunk into sub-chunks to reduce per-protein overhead.
        let sub_size = (chunk.len() / nthreads).max(1);
        let local_per_worker: Vec<Vec<AHashSet<Vec<u8>>>> = chunk
            .par_chunks(sub_size)
            .map(|sub| {
                let mut local: Vec<AHashSet<Vec<u8>>> = lengths
                    .iter()
                    .map(|_| AHashSet::with_capacity(1024))
                    .collect();
                for protein in sub {
                    trypsin_into(protein, &lengths, &mut local);
                }
                local
            })
            .collect();
        // Merge each worker's sets into the global sets
        for local in local_per_worker {
            for (idx, set) in local.into_iter().enumerate() {
                let mut g = global[idx].lock().unwrap();
                g.extend(set);
            }
        }
    };

    for protein in fasta {
        chunk.push(protein);
        total_proteins += 1;
        if chunk.len() >= CHUNK_SIZE {
            process_chunk(&chunk);
            chunk.clear();
            if total_proteins.is_multiple_of(200_000) {
                let sizes: Vec<usize> = global.iter().map(|g| g.lock().unwrap().len()).collect();
                eprintln!(
                    "  {} proteins, unique sizes: {:?}, {:.0}s",
                    total_proteins,
                    sizes,
                    t0.elapsed().as_secs_f64()
                );
            }
        }
    }
    if !chunk.is_empty() {
        process_chunk(&chunk);
    }
    let elapsed = t0.elapsed().as_secs_f64();

    // Write output
    let mut out = File::create(&args.output).expect("cannot create output");
    writeln!(
        out,
        "database\tmin_length\tK_theoretical_20pow_N\tD_unique_peptides_observed\tsaturation_percent\tn_proteins\telapsed_seconds"
    )
    .unwrap();
    for (i, &len) in lengths.iter().enumerate() {
        let d = global[i].lock().unwrap().len() as u64;
        let k: u64 = 20u64.checked_pow(len as u32).unwrap_or_else(|| {
            eprintln!(
                "error: min_length {} overflows u64 in K_theoretical = 20^{} (maximum supported length is 14)",
                len, len
            );
            std::process::exit(1);
        });
        let sat = (d as f64 / k as f64) * 100.0;
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{:.6}\t{}\t{:.1}",
            args.name, len, k, d, sat, total_proteins, elapsed
        )
        .unwrap();
        eprintln!("  L={}: K={:.3e}  D={}  D/K={:.6}%", len, k as f64, d, sat);
    }
    eprintln!(
        "Done. {} proteins in {:.1}s ({:.0} prot/s)",
        total_proteins,
        elapsed,
        total_proteins as f64 / elapsed
    );
}
