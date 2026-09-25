//! FASTALAKE Stage 0: Reference Protein Lake Builder
//!
//! Two-phase memory-efficient SHA256-deduplicated FASTA lake builder.
//! Handles 1B+ sequences using a compact index (~91 bytes/entry).
//!
//! Hashing contract: uppercased A-Z only (drops stop codons, whitespace,
//! digits, punctuation); SHA256 over the normalized bytes. This is shared
//! byte-for-byte with `fasta_lake/hashing.py` — both sides are pinned by
//! tests/test_hashing_crosslang.py. v5.3+ lakes use SHA256; pre-v5.3 lakes
//! used MD5 and must be rebuilt.
//!
//! Phase 1 (Index):
//!   Streams through all source files, computing SHA256 hashes.
//!   Stores only a compact {SHA256 -> (source_idx, n_sources)} index.
//!   Tracks full joint headers only for reference proteins (P1-P3).
//!
//! Phase 2 (Write):
//!   Re-reads source files, writes unique sequences to output FASTA.
//!   Each sequence written once, by the highest-priority source that owns it.
//!
//! Priority system:
//!   Lower number = higher priority = keeps header on collision.
//!   When a smORF matches a reference protein, the reference header is kept
//!   and the smORF ID is recorded in the joint headers TSV.
//!
//! Usage:
//!   lake_builder \
//!     -s "PRODIGAL:1:/path/to/dir:*_proteins.faa" \
//!     -s "POXA48_NR:2:/path/to/file.fasta" \
//!     -s "UNIPROT_ECOLI:3:/path/to/ecoli.fasta" \
//!     -s "GMSC:4:/path/to/GMSC10.100AA.faa" \
//!     -o /path/to/output.fasta \
//!     --output-headers /path/to/headers.tsv \
//!     --output-stats /path/to/stats.txt \
//!     -t 24

use clap::Parser;
use dashmap::DashMap;
use indicatif::{ProgressBar, ProgressStyle};
use memmap2::MmapOptions;
use parking_lot::Mutex;
use rayon::prelude::*;
use sha2::{Digest, Sha256};
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use std::sync::atomic::{AtomicU32, AtomicUsize, Ordering};
use std::time::Instant;

// ============================================================================
// Hashing contract provenance
// ============================================================================
// Frozen sequence-hashing contract, kept byte-for-byte in sync with
// fasta_lake/hashing.py (HASH_VERSION / HASH_ALGORITHM) and pinned by
// tests/test_hashing_crosslang.py. v5.3+ lakes use SHA256; pre-v5.3 lakes used
// MD5 and must be rebuilt. Recorded in the stats sidecar for provenance.
const HASH_VERSION: &str = "5.3";
const HASH_ALGORITHM: &str = "sha256";

// ============================================================================
// CLI
// ============================================================================

#[derive(Parser, Debug)]
#[command(name = "lake_builder")]
#[command(about = "Build deduplicated reference protein lake from multiple FASTA sources")]
struct Args {
    /// Source: TAG:PRIORITY:PATH[:GLOB]
    /// Lower priority = keeps header on collision.
    /// Examples:
    ///   -s "PRODIGAL:1:/data/prodigal:*_proteins.faa"
    ///   -s "GMSC:4:/data/GMSC10.100AA.faa"
    #[arg(short = 's', long = "source")]
    sources: Vec<String>,

    /// Output FASTA file
    #[arg(short = 'o', long)]
    output: PathBuf,

    /// Output joint headers TSV (sha256 -> all source headers)
    #[arg(long)]
    output_headers: Option<PathBuf>,

    /// Output statistics file
    #[arg(long)]
    output_stats: Option<PathBuf>,

    /// Number of threads (default: all available)
    #[arg(short = 't', long)]
    threads: Option<usize>,

    /// Chunk size for parallel processing
    #[arg(long, default_value = "100000")]
    chunk_size: usize,

    /// Initial hash-index capacity; grows with input. Avoids multi-GB toy-run reservation.
    #[arg(long, default_value = "65536")]
    initial_capacity: usize,

    /// Priority threshold for full header tracking.
    /// Sources with priority <= this value get full joint header tracking.
    /// Higher-priority sources only get counted.
    #[arg(long, default_value = "3")]
    header_tracking_priority: u8,

    /// Comma-separated source TAGs whose FASTA headers are REASSIGNED to a
    /// hash-derived stable id `{TAG}_{sha256hex}` instead of being kept
    /// verbatim. Use for metaG / assembly sources (MEGAHIT/Prodigal) whose
    /// native contig ids (e.g. k141_1234_5) are unique only WITHIN one
    /// assembly and therefore collide across samples — producing
    /// cross-sample protein-identity conflation downstream. Catalogue sources
    /// (SwissProt / UHGP / GMGC / GMSC / ...) are left untouched because their
    /// accessions are already globally unique. Because the id is derived from
    /// the sequence sha256, the SAME sequence receives the SAME id in every
    /// sample and every study (cross-study reproducibility). The original
    /// header is preserved in the joint-headers sidecar for provenance.
    #[arg(long)]
    reassign_header_tags: Option<String>,

    /// Reassign the accession of EVERY sequence, from EVERY source, to
    /// `{TAG}_{sha256hex}` using this ONE shared tag.
    ///
    /// WHY THIS EXISTS. `--reassign-header-tags` assigns a *per-source* tag, so
    /// reassigning all sources yields `GMGC_<hash>`, `GMSC_<hash>`, `UHGG_<hash>`.
    /// Those still encode the source, and they sort GMGC < GMSC < UHGG — which
    /// re-creates the very bias we are removing, because razor's default `alpha`
    /// tie-break picks the lexicographically smallest accession. Measured on
    /// healthy stool: alpha placed 57.96% of the database in GMGC where three
    /// mechanistically unrelated rules all placed ~40%.
    ///
    /// With a single shared tag the accession carries NO source signal, so no
    /// catalogue can be preferred by any downstream ordering. Measured on PREDICT,
    /// which already has hash-only accessions, the largest alpha-vs-hash-acc
    /// difference is 0.016 pp — i.e. with no prefix there is nothing left to bias.
    ///
    /// Provenance is NOT lost: every original header, from every source that owns
    /// the sequence, is written to the joint-headers sidecar. Use
    /// `--header-tracking-priority` high enough to cover all sources, or provenance
    /// for the uncovered ones is silently dropped.
    ///
    /// Takes precedence over `--reassign-header-tags` when both are given.
    #[arg(long)]
    uniform_tag: Option<String>,

    /// Track joint headers for EVERY source, discarding no provenance.
    ///
    /// Equivalent to setting `--header-tracking-priority` above the highest source
    /// number, but says what it means. That flag is a TRACKING THRESHOLD, not a
    /// ranking: sources numbered above it are counted but their headers are thrown
    /// away. The healthy-stool build ran the default 3 with GMGC at 4, so GMGC's
    /// headers were discarded -- 187,204,001 of 233,497,702 sequences tracked --
    /// which destroyed exactly the provenance needed to report where a sequence
    /// came from. Prefer this flag; it cannot be got wrong by miscounting sources.
    #[arg(long, default_value_t = false)]
    track_all_sources: bool,
}

// ============================================================================
// Types
// ============================================================================

type Sha256Hash = [u8; 32];

struct SourceSpec {
    tag: String,
    priority: u8,
    files: Vec<PathBuf>,
    description: String,
}

/// Compact index entry: 5 bytes of payload + DashMap overhead.
/// No sequence stored — sequences are re-read in Phase 2.
struct IndexEntry {
    source_idx: u8,
    n_sources: AtomicU32,
}

#[derive(Default)]
struct SourceStats {
    tag: String,
    priority: u8,
    n_files: usize,
    total_proteins: usize,
    new_unique: usize,
    duplicates: usize,
    elapsed_secs: f64,
}

// ============================================================================
// Helpers
// ============================================================================

fn fmt_num(n: usize) -> String {
    let s = n.to_string();
    let mut result = String::new();
    for (i, c) in s.chars().rev().enumerate() {
        if i > 0 && i % 3 == 0 {
            result.push(',');
        }
        result.push(c);
    }
    result.chars().rev().collect()
}

fn sha256_hex(hash: &Sha256Hash) -> String {
    hash.iter().map(|b| format!("{:02x}", b)).collect()
}

/// Deterministically order the collected joint headers for one sequence hash.
///
/// Sort by source priority (lower number = higher priority, kept first), then
/// lexicographically by the `"tag\theader"` string as a total secondary
/// tie-break. This is order-independent: the same multiset of headers always
/// produces the same ordering regardless of the parallel push order that
/// produced it, so `kept_header` (== `first()`) and `all_source_headers` are
/// reproducible bit-for-bit across runs and thread counts.
///
/// Does NOT change the SET of headers — only their order.
fn sort_joint_headers(mut headers: Vec<(u8, String)>) -> Vec<(u8, String)> {
    headers.sort_by(|a, b| a.0.cmp(&b.0).then_with(|| a.1.cmp(&b.1)));
    headers
}

/// Deterministic order that does NOT consult source priority: sorts on the
/// "TAG\theader" string alone.
///
/// WHY. `sort_joint_headers` sorts by (priority, string), so the FIRST element --
/// written to the `kept_header` column -- is whichever source was given the lowest
/// number by hand. That makes the provenance file itself carry the ranking, even
/// when the FASTA accession is source-neutral. Two cohorts numbered the same
/// catalogues differently (CAPSCAN ranked GMGC 2nd, healthy stool ranked it 4th),
/// so the same shared sequence was designated GMGC in one and GMSC in the other.
/// Sorting on the string alone removes the ranking from the sidecar too.
fn sort_joint_headers_source_blind(mut headers: Vec<(u8, String)>) -> Vec<(u8, String)> {
    headers.sort_by(|a, b| a.1.cmp(&b.1));
    headers
}

fn parse_source_spec(spec: &str) -> Result<SourceSpec, String> {
    let parts: Vec<&str> = spec.splitn(4, ':').collect();
    if parts.len() < 3 {
        return Err(format!(
            "Invalid source spec '{}': need TAG:PRIORITY:PATH[:GLOB]",
            spec
        ));
    }

    let tag = parts[0].to_string();
    let priority: u8 = parts[1]
        .parse()
        .map_err(|_| format!("Invalid priority '{}' in '{}'", parts[1], spec))?;
    let path = std::path::Path::new(parts[2]);
    let glob_pattern = parts.get(3).copied();

    let files = if let Some(pattern) = glob_pattern {
        if !path.is_dir() {
            return Err(format!(
                "'{}' is not a directory (glob requires directory)",
                parts[2]
            ));
        }
        let full_pattern = format!("{}/{}", parts[2], pattern);
        let mut matched: Vec<PathBuf> = glob::glob(&full_pattern)
            .map_err(|e| format!("Glob error: {}", e))?
            .filter_map(|r| r.ok())
            .collect();
        matched.sort();
        if matched.is_empty() {
            return Err(format!("No files matched '{}'", full_pattern));
        }
        matched
    } else if path.is_file() {
        vec![path.to_path_buf()]
    } else {
        return Err(format!("'{}' not found or not a file", parts[2]));
    };

    let description = if let Some(pattern) = glob_pattern {
        format!("{} ({} files, {})", tag, files.len(), pattern)
    } else {
        format!(
            "{} ({})",
            tag,
            path.file_name().unwrap_or_default().to_string_lossy()
        )
    };

    Ok(SourceSpec {
        tag,
        priority,
        files,
        description,
    })
}

/// Parse FASTA regions from memory-mapped data.
fn parse_fasta_regions(data: &[u8]) -> Vec<(usize, usize, usize, usize)> {
    let est = data.len() / 150;
    let mut regions = Vec::with_capacity(est.min(200_000_000));
    let mut i = 0;
    let len = data.len();

    while i < len {
        if data[i] == b'>' {
            let header_start = i + 1;
            while i < len && data[i] != b'\n' {
                i += 1;
            }
            let header_end = i;
            if i < len {
                i += 1;
            }
            let seq_start = i;
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

/// Extract and clean protein sequence (uppercase, strip *, whitespace).
fn extract_clean_sequence(data: &[u8], start: usize, end: usize) -> Vec<u8> {
    let mut seq = Vec::with_capacity(end - start);
    for &b in &data[start..end] {
        match b {
            b'A'..=b'Z' => seq.push(b),
            b'a'..=b'z' => seq.push(b - 32),
            _ => {}
        }
    }
    seq
}

fn progress_bar(total: u64, label: &str) -> ProgressBar {
    let pb = ProgressBar::new(total);
    pb.set_style(
        ProgressStyle::default_bar()
            .template(&format!(
                "  {} [{{elapsed_precise}}] [{{bar:50.cyan/blue}}] {{pos}}/{{len}} ({{per_sec}}) {{msg}}",
                label
            ))
            .unwrap(),
    );
    pb
}

// ============================================================================
// Main
// ============================================================================

/// Build a sequence-keyed reservoir from the declared FASTA sources.
///
/// Clean and hash complete sequences, merge exact duplicates and retain original
/// source headers in sidecars. I and L remain distinct in the protein hash.
/// This preparatory command does not perform peptide lookup or similarity clustering.
fn main() -> std::io::Result<()> {
    let args = Args::parse();
    let t0 = Instant::now();
    let invalid = |m: &str| std::io::Error::new(std::io::ErrorKind::InvalidInput, m);
    if args.sources.is_empty() || args.chunk_size == 0 || args.threads == Some(0) {
        return Err(invalid(
            "Require at least one source, positive chunk-size and positive threads",
        ));
    }
    if args.uniform_tag.as_deref().is_some_and(|s| {
        s.is_empty()
            || !s
                .bytes()
                .all(|c| c.is_ascii_alphanumeric() || b"._-".contains(&c))
    }) {
        return Err(invalid(
            "uniform-tag must contain only ASCII letters, digits, dot, underscore or hyphen",
        ));
    }
    let mut destinations = std::collections::HashSet::new();
    for p in [
        Some(&args.output),
        args.output_headers.as_ref(),
        args.output_stats.as_ref(),
    ]
    .into_iter()
    .flatten()
    {
        if p.exists() || p.is_symlink() {
            return Err(invalid("Refusing to overwrite an existing output"));
        }
        let parent = p
            .parent()
            .filter(|p| !p.as_os_str().is_empty())
            .unwrap_or(std::path::Path::new("."));
        fs::create_dir_all(parent)?;
        let dest = parent.canonicalize()?.join(
            p.file_name()
                .ok_or_else(|| invalid("Output must name a file"))?,
        );
        if !destinations.insert(dest) {
            return Err(invalid("Output paths must be distinct"));
        }
    }

    let num_threads = args.threads.unwrap_or_else(|| {
        std::thread::available_parallelism()
            .map(|p| p.get())
            .unwrap_or(8)
    });
    rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build_global()
        .ok();

    println!("======================================================================");
    println!("FASTALAKE STAGE 0: Reference Protein Lake Builder (Rust)");
    println!("======================================================================");
    println!("Threads: {}", num_threads);
    match &args.reassign_header_tags {
        Some(t) => println!(
            "Headers: catalogue sources kept verbatim; reassigned to hash ids for tags [{}]",
            t
        ),
        None => println!("Headers: original (kept verbatim for all sources)"),
    }
    println!(
        "Header tracking for priority <= {}",
        args.header_tracking_priority
    );
    println!();

    // --- Parse sources, sort by priority ---
    let mut sources: Vec<SourceSpec> = Vec::new();
    for spec in &args.sources {
        match parse_source_spec(spec) {
            Ok(s) => {
                println!("  [P{}] {}", s.priority, s.description);
                sources.push(s);
            }
            Err(e) => {
                eprintln!("ERROR: {}", e);
                std::process::exit(1);
            }
        }
    }
    sources.sort_by_key(|s| s.priority);

    // Source indices are stored as u8 (IndexEntry.source_idx). With more than
    // 256 sources the `src_idx as u8` casts below would wrap and silently
    // conflate distinct sources, corrupting ownership/dedup. Fail loudly.
    if sources.len() > 256 {
        eprintln!(
            "ERROR: {} sources exceeds the 256-source limit (source_idx is u8)",
            sources.len()
        );
        std::process::exit(1);
    }
    println!();

    // --- Data structures ---
    // Compact index: SHA256 -> (source_idx, n_sources). No sequence stored.
    let index: DashMap<Sha256Hash, IndexEntry> = DashMap::with_capacity(args.initial_capacity);

    // Full headers only for reference proteins (priority <= threshold).
    // Value = (source_priority, "tag\theader"); priority is stored so the
    // joint-header list can be ordered deterministically (see sort_joint_headers)
    // instead of in the race-dependent parallel push order.
    let ref_headers: DashMap<Sha256Hash, Vec<(u8, String)>> = DashMap::new();

    // --track-all-sources wins: track every source, discard no provenance.
    let header_prio_threshold = if args.track_all_sources {
        u8::MAX
    } else {
        args.header_tracking_priority
    };
    if args.track_all_sources {
        println!("Tracking joint headers for ALL sources (no provenance discarded).");
    }
    let mut all_stats: Vec<SourceStats> = Vec::new();
    let mut grand_total: usize = 0;

    // ====================================================================
    // PHASE 1: INDEX
    // ====================================================================
    println!("========== PHASE 1: Indexing all sources ==========");
    println!();

    for (src_idx, source) in sources.iter().enumerate() {
        let t_src = Instant::now();
        let cnt_new = AtomicUsize::new(0);
        let cnt_dup = AtomicUsize::new(0);
        let cnt_total = AtomicUsize::new(0);

        println!("--- [P{}] {} ---", source.priority, source.description);

        for file_path in &source.files {
            let fsize = fs::metadata(file_path)?.len();
            let fname = file_path.file_name().unwrap_or_default().to_string_lossy();

            if source.files.len() > 1 {
                print!("\r  {}: ", fname);
                std::io::stdout().flush()?;
            } else {
                println!("  File: {} ({:.2} GB)", fname, fsize as f64 / 1e9);
            }

            let file = File::open(file_path)?;
            // SAFETY: source FASTA files must not be truncated or modified for
            // the duration of the build; concurrent external mutation of a
            // mmap'd file is UB (typically SIGBUS on read).
            let mmap = unsafe { MmapOptions::new().map(&file)? };
            let data = &mmap[..];
            let regions = parse_fasta_regions(data);
            let n = regions.len();

            if source.files.len() > 1 {
                print!("{} proteins", fmt_num(n));
                std::io::stdout().flush()?;
            } else {
                println!("  Proteins: {}", fmt_num(n));
            }

            let pb = if n > 1_000_000 {
                Some(progress_bar(n as u64, "IDX"))
            } else {
                None
            };

            let src_i = src_idx as u8;
            let prio = source.priority;
            let tag = &source.tag;

            regions.par_chunks(args.chunk_size).for_each(|chunk| {
                for &(hs, he, ss, se) in chunk {
                    let seq = extract_clean_sequence(data, ss, se);
                    if seq.is_empty() {
                        continue;
                    }

                    let digest: Sha256Hash = Sha256::digest(&seq).into();
                    cnt_total.fetch_add(1, Ordering::Relaxed);

                    // Try insert into index
                    let is_new = {
                        let mut new = false;
                        index
                            .entry(digest)
                            .and_modify(|e| {
                                e.n_sources.fetch_add(1, Ordering::Relaxed);
                            })
                            .or_insert_with(|| {
                                new = true;
                                IndexEntry {
                                    source_idx: src_i,
                                    n_sources: AtomicU32::new(1),
                                }
                            });
                        new
                    };

                    if is_new {
                        cnt_new.fetch_add(1, Ordering::Relaxed);
                    } else {
                        cnt_dup.fetch_add(1, Ordering::Relaxed);
                    }

                    // Track headers for reference proteins + any smORF overlaps
                    let dominated_by_ref = if !is_new {
                        index
                            .get(&digest)
                            .map(|e| e.source_idx < src_i)
                            .unwrap_or(false)
                    } else {
                        false
                    };

                    if prio <= header_prio_threshold || dominated_by_ref {
                        let header = std::str::from_utf8(&data[hs..he])
                            .unwrap_or("unknown")
                            .to_string();
                        // Store source tag and original header separately, keyed
                        // with the source priority so the list can be ordered
                        // deterministically at write time (parallel push order
                        // here is race-dependent and must not leak into output).
                        // "TAG header": a space, not a tab, so the sidecar row keeps
                        // exactly the four columns its header line declares.
                        let entry_str = format!("{} {}", tag, header);
                        ref_headers
                            .entry(digest)
                            .and_modify(|v| v.push((prio, entry_str.clone())))
                            .or_insert_with(|| vec![(prio, entry_str)]);
                    }
                }

                if let Some(ref pb) = pb {
                    pb.inc(chunk.len() as u64);
                    pb.set_message(format!("new: {}", fmt_num(cnt_new.load(Ordering::Relaxed))));
                }
            });

            if let Some(pb) = pb {
                pb.finish_and_clear();
            }
        }

        if source.files.len() > 1 {
            println!();
        }

        let total = cnt_total.load(Ordering::Relaxed);
        let new = cnt_new.load(Ordering::Relaxed);
        let dup = cnt_dup.load(Ordering::Relaxed);
        let elapsed = t_src.elapsed().as_secs_f64();
        grand_total += total;

        println!(
            "  Total: {}  New: {}  Dup: {}  Time: {:.1}s  Index: {}",
            fmt_num(total),
            fmt_num(new),
            fmt_num(dup),
            elapsed,
            fmt_num(index.len())
        );
        println!();

        all_stats.push(SourceStats {
            tag: source.tag.clone(),
            priority: source.priority,
            n_files: source.files.len(),
            total_proteins: total,
            new_unique: new,
            duplicates: dup,
            elapsed_secs: elapsed,
        });
    }

    let lake_size = index.len();
    let ref_tracked = ref_headers.len();

    println!("======================================================================");
    println!("PHASE 1 COMPLETE");
    println!("  Total input:     {}", fmt_num(grand_total));
    println!("  Unique (index):  {}", fmt_num(lake_size));
    println!(
        "  Ref-tracked:     {} (full joint headers)",
        fmt_num(ref_tracked)
    );
    if grand_total > 0 {
        println!(
            "  Reduction:       {:.1}%",
            (1.0 - lake_size as f64 / grand_total as f64) * 100.0
        );
    }
    println!(
        "  Index memory:    ~{:.1} GB estimated",
        lake_size as f64 * 72.0 / 1e9
    );
    println!("======================================================================");
    println!();

    // ====================================================================
    // PHASE 2: WRITE FASTA
    // ====================================================================
    println!("========== PHASE 2: Writing output FASTA ==========");
    println!();

    if let Some(parent) = args.output.parent() {
        fs::create_dir_all(parent)?;
    }

    let output_file = File::create(&args.output)?;
    let writer = Mutex::new(BufWriter::with_capacity(64 * 1024 * 1024, output_file));
    let written_count = AtomicUsize::new(0);

    // Track written SHA256s to prevent duplicate writes within multi-file sources
    let written_hashes: DashMap<Sha256Hash, ()> = DashMap::with_capacity(lake_size);
    // Only provenance-tracked entries need retained header metadata.
    let mut emitted_headers: std::collections::HashMap<String, String> =
        std::collections::HashMap::new();

    // Source tags whose headers are reassigned to hash-derived stable ids.
    let reassign_tags: std::collections::HashSet<String> = args
        .reassign_header_tags
        .as_deref()
        .map(|s| {
            s.split(',')
                .map(|t| t.trim().to_string())
                .filter(|t| !t.is_empty())
                .collect()
        })
        .unwrap_or_default();
    if !reassign_tags.is_empty() {
        println!(
            "Reassigning hash-derived ids for source tags: {:?}",
            reassign_tags
        );
    }

    // Source-neutral accessions: ONE shared tag for every source (see --uniform-tag).
    // Fail loudly on an empty or whitespace tag rather than silently emitting ">_<hash>".
    let uniform_tag: Option<String> = match args.uniform_tag.as_deref() {
        None => None,
        Some(t) if t.trim().is_empty() => {
            eprintln!("ERROR: --uniform-tag was given an empty tag.");
            std::process::exit(2);
        }
        Some(t) => Some(t.trim().to_string()),
    };
    if let Some(ref ut) = uniform_tag {
        println!(
            "SOURCE-NEUTRAL ACCESSIONS: every sequence from every source -> {}_<sha256>",
            ut
        );
        if !reassign_tags.is_empty() {
            println!(
                "  NOTE: --uniform-tag overrides --reassign-header-tags ({:?}); \
                 a per-source tag would re-encode the source.",
                reassign_tags
            );
        }
        if (args.header_tracking_priority as usize) < sources.len() {
            println!(
                "  WARNING: --header-tracking-priority {} does not cover all {} sources. \
                 Accessions will carry no source, and provenance for the uncovered sources \
                 is NOT tracked -- their contribution becomes unrecoverable.",
                args.header_tracking_priority,
                sources.len()
            );
        }
    }

    for (src_idx, source) in sources.iter().enumerate() {
        let src_i = src_idx as u8;
        let tag = &source.tag;
        let _tag = tag; // keep for stats/logging
        let reassign = reassign_tags.contains(tag);

        // Count how many entries belong to this source
        let expected: usize = index
            .iter()
            .filter(|e| e.value().source_idx == src_i)
            .count();

        if expected == 0 {
            println!(
                "  [P{}] {}: 0 sequences to write, skipping re-read",
                source.priority, tag
            );
            continue;
        }

        println!(
            "  [P{}] {}: writing {} sequences",
            source.priority,
            tag,
            fmt_num(expected)
        );

        for file_path in &source.files {
            let file = File::open(file_path)?;
            // SAFETY: source FASTA files must not be truncated or modified for
            // the duration of the build; concurrent external mutation of a
            // mmap'd file is UB (typically SIGBUS on read).
            let mmap = unsafe { MmapOptions::new().map(&file)? };
            let data = &mmap[..];
            let regions = parse_fasta_regions(data);
            let n = regions.len();

            let pb = if n > 1_000_000 {
                Some(progress_bar(n as u64, "WRT"))
            } else {
                None
            };

            // Sequential write to maintain deterministic output per-file
            for (batch_i, chunk) in regions.chunks(args.chunk_size).enumerate() {
                // Parallel: compute SHA256 + extract data for every region this
                // source OWNS. Deliberately does NOT dedup here — deduplication
                // (which region's header survives for a duplicate sequence) is
                // done sequentially below so it cannot depend on thread/completion
                // order. par_iter().filter_map().collect() preserves region order.
                let owned_entries: Vec<_> = chunk
                    .par_iter()
                    .filter_map(|&(hs, he, ss, se)| {
                        let seq = extract_clean_sequence(data, ss, se);
                        if seq.is_empty() {
                            return None;
                        }

                        let digest: Sha256Hash = Sha256::digest(&seq).into();

                        // Keep only sequences this source owns.
                        if let Some(entry) = index.get(&digest) {
                            if entry.source_idx == src_i {
                                let header =
                                    std::str::from_utf8(&data[hs..he]).unwrap_or("unknown");
                                let n_sources = entry.n_sources.load(Ordering::Relaxed);
                                return Some((digest, header.to_string(), seq, n_sources));
                            }
                        }
                        None
                    })
                    .collect();

                // Sequential write in region order (no contention). The dedup
                // check-and-insert runs HERE, so when several regions in this
                // source carry the same sequence the LOWEST region index wins
                // deterministically (chunks are processed in order, and within a
                // chunk owned_entries is in region order). Only the surviving
                // HEADER / output position is affected — the set of unique
                // sequences written is identical to before.
                if !owned_entries.is_empty() {
                    let mut w = writer.lock();
                    let mut n_written_here = 0usize;
                    for (_digest, header, seq, _n_sources) in &owned_entries {
                        if written_hashes.insert(*_digest, ()).is_some() {
                            continue; // already written by an earlier (lower-index) region
                        }
                        // Catalogue sources: keep original (globally-unique) header.
                        // metaG/assembly sources (in --reassign-header-tags): emit a
                        // hash-derived stable id so identical sequences share one id
                        // across all samples/studies (the original is kept in the
                        // joint-headers sidecar for provenance).
                        let emitted = if let Some(ref ut) = uniform_tag {
                            format!("{}_{}", ut, sha256_hex(_digest))
                        } else if reassign {
                            format!("{}_{}", tag, sha256_hex(_digest))
                        } else {
                            header.clone()
                        };
                        writeln!(w, ">{}", emitted)?;
                        if args.output_headers.is_some() && ref_headers.contains_key(_digest) {
                            emitted_headers.insert(sha256_hex(_digest), emitted);
                        }
                        for line in seq.chunks(60) {
                            w.write_all(line).unwrap();
                            w.write_all(b"\n").unwrap();
                        }
                        n_written_here += 1;
                    }
                    written_count.fetch_add(n_written_here, Ordering::Relaxed);
                }

                if let Some(ref pb) = pb {
                    pb.set_position(((batch_i + 1) * args.chunk_size).min(n) as u64);
                    pb.set_message(format!(
                        "written: {}",
                        fmt_num(written_count.load(Ordering::Relaxed))
                    ));
                }
            }

            if let Some(pb) = pb {
                pb.finish_and_clear();
            }
        }

        // Remove written entries from index to avoid duplicate writes
        // (handles sequences appearing in multiple files of the same source)
        index.retain(|_, v| v.source_idx != src_i);
    }

    writer.lock().flush()?;
    let total_written = written_count.load(Ordering::Relaxed);

    println!();
    println!("  Written: {} proteins", fmt_num(total_written));
    println!();

    // ====================================================================
    // WRITE HEADERS TSV
    // ====================================================================
    if let Some(ref hdr_path) = args.output_headers {
        println!("Writing joint headers: {:?}", hdr_path);
        let t_h = Instant::now();
        let hf = File::create(hdr_path)?;
        let mut hw = BufWriter::with_capacity(32 * 1024 * 1024, hf);

        // all_source_headers: "TAG header" entries joined by "|||", one per source
        // that carried the sequence, ordered by (priority, string).
        writeln!(
            hw,
            "sha256_hash\tn_sources\tkept_header\tall_source_headers"
        )?;
        let uniform_tag_for_sidecar: Option<String> = args
            .uniform_tag
            .as_deref()
            .map(|t| t.trim().to_string())
            .filter(|t| !t.is_empty());
        // Collect and sort rows by sha256 hex (unique + total) so TSV row order
        // is independent of DashMap iteration order.
        let mut rows: Vec<(String, Vec<(u8, String)>)> = ref_headers
            .iter()
            .map(|entry| (sha256_hex(entry.key()), entry.value().clone()))
            .collect();
        rows.sort_by(|a, b| a.0.cmp(&b.0));
        for (hex, headers) in &rows {
            // Order the headers deterministically (priority, then string) so that
            // kept_header, all_source_headers, and first() are reproducible and
            // mutually consistent across runs / thread counts.
            // With a uniform tag there is no "winning" source, so the sidecar must
            // not designate one. Order the list source-blind and set kept_header to
            // the accession actually written to the FASTA. Every source header stays
            // in all_source_headers, so provenance is a SET rather than a winner.
            let neutral = uniform_tag_for_sidecar.is_some();
            let ordered = if neutral {
                sort_joint_headers_source_blind(headers.clone())
            } else {
                sort_joint_headers(headers.clone())
            };
            let kept = emitted_headers.get(hex).ok_or_else(|| {
                std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!("No emitted FASTA header for tracked sequence {}", hex),
                )
            })?;
            let all = ordered
                .iter()
                .map(|(_, s)| s.as_str())
                .collect::<Vec<_>>()
                .join("|||");
            writeln!(hw, "{}\t{}\t{}\t{}", hex, ordered.len(), kept, all)?;
        }
        hw.flush()?;
        println!(
            "  {} entries in {:.1}s",
            fmt_num(ref_headers.len()),
            t_h.elapsed().as_secs_f64()
        );
    }

    // ====================================================================
    // WRITE STATISTICS
    // ====================================================================
    if let Some(ref stats_path) = args.output_stats {
        println!("Writing statistics: {:?}", stats_path);
        let sf = File::create(stats_path)?;
        let mut sw = BufWriter::new(sf);
        let total_time = t0.elapsed().as_secs_f64();

        writeln!(sw, "FASTALAKE STAGE 0: Reference Protein Lake Statistics")?;
        writeln!(sw, "{}", "=".repeat(90))?;
        writeln!(sw, "Builder: lake_builder (Rust, two-phase)")?;
        // Provenance: record the frozen hashing contract so downstream consumers
        // can tell a SHA256 v5.3 lake from a pre-v5.3 MD5 lake without re-reading
        // the FASTA. Kept byte-for-byte in sync with fasta_lake/hashing.py
        // (HASH_VERSION / HASH_ALGORITHM), pinned by tests/test_hashing_crosslang.py.
        writeln!(sw, "Hash algorithm: {}", HASH_ALGORITHM)?;
        writeln!(sw, "Hash version: {}", HASH_VERSION)?;
        writeln!(sw, "Threads: {}", num_threads)?;
        writeln!(sw, "Total time: {:.0}s", total_time)?;
        writeln!(sw)?;
        writeln!(sw, "INPUT SUMMARY")?;
        writeln!(sw, "  Total input proteins:    {}", fmt_num(grand_total))?;
        writeln!(sw, "  Unique sequences (lake): {}", fmt_num(lake_size))?;
        writeln!(sw, "  Written to FASTA:        {}", fmt_num(total_written))?;
        if grand_total > 0 {
            writeln!(
                sw,
                "  Reduction:               {:.1}%",
                (1.0 - lake_size as f64 / grand_total as f64) * 100.0
            )?;
        }
        writeln!(sw)?;

        writeln!(sw, "SOURCES BY PRIORITY")?;
        writeln!(
            sw,
            "{:<25} {:>5} {:>6} {:>14} {:>14} {:>14} {:>8}",
            "Source", "Prio", "Files", "Input", "Unique", "Dupes", "Time"
        )?;
        writeln!(sw, "{}", "-".repeat(90))?;
        for s in &all_stats {
            writeln!(
                sw,
                "{:<25} {:>5} {:>6} {:>14} {:>14} {:>14} {:>7.1}s",
                s.tag,
                s.priority,
                s.n_files,
                fmt_num(s.total_proteins),
                fmt_num(s.new_unique),
                fmt_num(s.duplicates),
                s.elapsed_secs
            )?;
        }
        writeln!(sw, "{}", "-".repeat(90))?;
        writeln!(
            sw,
            "{:<25} {:>5} {:>6} {:>14} {:>14}",
            "TOTAL",
            "",
            "",
            fmt_num(grand_total),
            fmt_num(lake_size)
        )?;
        writeln!(sw)?;

        writeln!(sw, "JOINT HEADER TRACKING")?;
        writeln!(sw, "  Reference proteins tracked: {}", fmt_num(ref_tracked))?;
        let multi = ref_headers.iter().filter(|r| r.value().len() > 1).count();
        writeln!(
            sw,
            "  With >1 source header:      {} ({:.1}%)",
            fmt_num(multi),
            if ref_tracked > 0 {
                multi as f64 / ref_tracked as f64 * 100.0
            } else {
                0.0
            }
        )?;

        sw.flush()?;
    }

    // --- Summary ---
    let total_time = t0.elapsed();
    println!();
    println!("======================================================================");
    println!("REFERENCE LAKE COMPLETE");
    println!("  Unique proteins: {}", fmt_num(lake_size));
    println!("  Written:         {}", fmt_num(total_written));
    println!("  Output:          {:?}", args.output);
    println!(
        "  Time:            {:.1}s ({:.1} min)",
        total_time.as_secs_f64(),
        total_time.as_secs_f64() / 60.0
    );
    println!("======================================================================");

    Ok(())
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    /// Duplicate sequence under different accessions/priorities: the kept
    /// header must be the priority winner (lowest priority number), with a
    /// lexicographic tie-break, and the result must be identical no matter what
    /// order the parallel producers pushed the headers in.
    #[test]
    fn joint_header_tiebreak_is_deterministic_and_priority_first() {
        // Same sequence hash seen from 4 sources, two of which share priority 2.
        let variants = vec![
            vec![
                (2u8, "UHGP\tMGYG_9".to_string()),
                (1u8, "SP\tsp|P12345|".to_string()),
                (2u8, "GMGC\tGMGC_1".to_string()),
                (3u8, "TREMBL\ttr|Q0000|".to_string()),
            ],
            // Same multiset, scrambled (simulates a different parallel push order).
            vec![
                (3u8, "TREMBL\ttr|Q0000|".to_string()),
                (2u8, "GMGC\tGMGC_1".to_string()),
                (2u8, "UHGP\tMGYG_9".to_string()),
                (1u8, "SP\tsp|P12345|".to_string()),
            ],
        ];

        let expected = vec![
            (1u8, "SP\tsp|P12345|".to_string()), // highest priority wins
            (2u8, "GMGC\tGMGC_1".to_string()),   // prio 2, lex-smaller string
            (2u8, "UHGP\tMGYG_9".to_string()),   // prio 2, lex-larger string
            (3u8, "TREMBL\ttr|Q0000|".to_string()),
        ];

        for v in variants {
            let ordered = sort_joint_headers(v);
            assert_eq!(ordered, expected, "ordering must be order-independent");
            // kept_header = first element's header field.
            let kept = ordered
                .first()
                .map(|(_, s)| s.split_once('\t').map(|(_, h)| h).unwrap_or(s.as_str()))
                .unwrap();
            assert_eq!(kept, "sp|P12345|");
        }
    }

    /// The tie-break must never change the SET of headers, only their order.
    #[test]
    fn joint_header_sort_preserves_the_set() {
        let input = vec![
            (5u8, "B\tb".to_string()),
            (1u8, "A\ta".to_string()),
            (5u8, "A\ta2".to_string()),
        ];
        let mut before = input.clone();
        before.sort();
        let mut after = sort_joint_headers(input);
        after.sort();
        assert_eq!(before, after);
    }

    /// The normalization contract (uppercase A-Z only, drop everything else)
    /// must remain byte-stable — guards the shared hashing contract.
    #[test]
    fn extract_clean_sequence_normalizes_deterministically() {
        let data = b"mk*L 12\nRq";
        let out = extract_clean_sequence(data, 0, data.len());
        assert_eq!(out, b"MKLRQ");
    }
}
