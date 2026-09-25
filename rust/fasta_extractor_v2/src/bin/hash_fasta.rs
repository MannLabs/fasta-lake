//! Cross-language hash verification utility.
//!
//! Reads a FASTA and emits one line per sequence:
//!   ``<sha256_hex>\t<accession>\t<normalized_length>``
//!
//! The normalization + hash implemented here is byte-for-byte identical to
//! ``lake_builder.rs`` and to Python's ``fasta_lake.hashing.sequence_hash``.
//! ``tests/test_hashing_crosslang.py`` runs this binary on several real
//! FASTA fixtures and asserts every line equals the Python computation.
//!
//! If the three implementations ever diverge, the cross-language test fails
//! loudly — exactly the signal we want, because silent dedup drift between
//! pipeline stages is a category of bug that otherwise goes unnoticed for
//! years.

use clap::Parser;
use sha2::{Digest, Sha256};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[command(about = "Emit SHA256 hash per FASTA record (cross-language test utility)")]
struct Args {
    /// Input FASTA file
    #[arg(short, long)]
    input: PathBuf,

    /// Output TSV (stdout if omitted)
    #[arg(short, long)]
    output: Option<PathBuf>,
}

/// Uppercase A-Z letters only. Drops everything else — whitespace, digits,
/// stop codons, punctuation. MUST stay byte-identical to
/// ``fasta_lake.hashing.normalize_sequence`` and to the in-line body of
/// ``extract_clean_sequence`` in ``lake_builder.rs``.
fn normalize(seq: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(seq.len());
    for &b in seq {
        match b {
            b'A'..=b'Z' => out.push(b),
            b'a'..=b'z' => out.push(b - 32),
            _ => {}
        }
    }
    out
}

fn hash_hex(bytes: &[u8]) -> String {
    let digest: [u8; 32] = Sha256::digest(bytes).into();
    digest.iter().map(|b| format!("{:02x}", b)).collect()
}

fn emit<W: Write>(out: &mut W, header: &str, seq: &[u8]) -> std::io::Result<()> {
    let norm = normalize(seq);
    let hex = hash_hex(&norm);
    let accession = header.split_whitespace().next().unwrap_or("");
    writeln!(out, "{}\t{}\t{}", hex, accession, norm.len())
}

fn main() -> std::io::Result<()> {
    let args = Args::parse();
    let inp = BufReader::new(File::open(&args.input)?);
    let mut out: Box<dyn Write> = match args.output {
        Some(p) => Box::new(BufWriter::new(File::create(p)?)),
        None => Box::new(BufWriter::new(std::io::stdout())),
    };

    let mut header: Option<String> = None;
    let mut seq_bytes: Vec<u8> = Vec::new();

    for line in inp.lines() {
        let line = line?;
        if let Some(rest) = line.strip_prefix('>') {
            if let Some(h) = header.take() {
                emit(&mut out, &h, &seq_bytes)?;
            }
            header = Some(rest.to_string());
            seq_bytes.clear();
        } else {
            seq_bytes.extend(line.trim_end().as_bytes());
        }
    }
    if let Some(h) = header.take() {
        emit(&mut out, &h, &seq_bytes)?;
    }
    Ok(())
}
