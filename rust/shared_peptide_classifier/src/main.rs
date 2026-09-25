//! shared_peptide_classifier — FastaLake v5.4
//!
//! Projects k-mer-rescued peptides from cluster representatives onto member
//! proteins and classifies each UNIQUE or SHARED (>=80% of members) — the
//! conservation call `parsimony_engine --weights` consumes.
//!
//! Given peptides anchored to cluster representatives by `mass_kmer_anchor`,
//! lift each peptide back to the unclustered member proteins that actually
//! contain it. Tags each (peptide, rep) pair as UNIQUE or SHARED based on
//! the fraction of cluster members containing the peptide.
//!
//! Inputs
//!   --rescued        TSV with header. Columns: peptide, cluster_rep, ... (extra cols ignored)
//!   --cluster-tsv    2-col TSV: representative<TAB>member (one row per member; reps may repeat)
//!   --lake           unclustered FASTA; headers are `>id [optional desc]` — only first whitespace-delimited token used as id
//!
//! Outputs
//!   --output         TSV: peptide, cluster_rep, n_cluster_members, n_with_hit, members_with_hit (csv), is_shared
//!
//! Algorithm
//!   1. Parse rescued.tsv → set of cluster_reps we care about.
//!   2. Stream cluster.tsv once → for needed reps, build Vec<member_id>.
//!      Also build set of all member ids we'll need to read from the lake.
//!   3. Stream lake.fasta once → for each header whose id ∈ needed_members,
//!      store (id → seq_bytes) in memory.
//!   4. For each (peptide, rep): scan each member's seq for the peptide via
//!      memmem::find. Record members containing it.
//!   5. is_shared = members_with_hit / RESOLVABLE_members_in_cluster ≥ threshold,
//!      where resolvable = cluster members whose sequence was found in the lake.
//!      Members absent from the lake cannot be assessed and are excluded from the
//!      denominator (see C9 note at the classification loop). The output still
//!      reports the true total member count in n_cluster_members.

use ahash::{AHashMap, AHashSet};
use clap::Parser;
use memchr::memmem;
use rayon::prelude::*;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

/// Build-time version string: <Cargo version>+<git short hash>[-dirty]
const VERSION: &str = concat!(env!("CARGO_PKG_VERSION"), "+", env!("GIT_REV"));

#[derive(Parser)]
#[command(version = VERSION, about = "Lift k-mer rescued peptides from cluster reps to unclustered members")]
struct Args {
    #[arg(long)]
    rescued: PathBuf,
    #[arg(long, value_name = "PATH")]
    cluster_tsv: PathBuf,
    #[arg(long)]
    lake: PathBuf,
    #[arg(long)]
    output: PathBuf,
    /// Members-with-hit fraction at or above this is tagged SHARED [default: 0.8]
    /// 0.8 means: a peptide hitting ≥80% of cluster members is treated as shared
    /// across the cluster rather than uniquely mapping to one representative. Tuned
    /// on the Masters healthy-stool benchmark — drives v5.4 pooled-evidence weights
    /// (see Supp Note 3 / shared_peptide_classifier README).
    #[arg(long, default_value_t = 0.8)]
    shared_threshold: f64,
    /// Column name in rescued.tsv holding the peptide [default: peptide]
    #[arg(long, default_value = "peptide")]
    peptide_col: String,
    /// Column name in rescued.tsv holding the cluster representative [default: cluster_rep]
    #[arg(long, default_value = "cluster_rep")]
    rep_col: String,
    /// Threads (0 = all) [default: 0]
    #[arg(short, long, default_value_t = 0)]
    threads: usize,
}

/// Read rescued.tsv → Vec<(peptide, rep)>. Picks columns by header name.
fn load_rescued(
    path: &Path,
    peptide_col: &str,
    rep_col: &str,
) -> std::io::Result<(Vec<(String, String)>, usize)> {
    let f = File::open(path)?;
    let r = BufReader::new(f);
    let mut iter = r.lines();
    let header = iter
        .next()
        .ok_or_else(|| std::io::Error::other("empty rescued file"))??;
    let cols: Vec<&str> = header.split('\t').collect();
    let pep_idx = cols
        .iter()
        .position(|c| *c == peptide_col)
        .ok_or_else(|| std::io::Error::other(format!("missing column: {}", peptide_col)))?;
    let rep_idx = cols
        .iter()
        .position(|c| *c == rep_col)
        .ok_or_else(|| std::io::Error::other(format!("missing column: {}", rep_col)))?;
    let mut out = Vec::new();
    let mut skipped = 0usize;
    for line in iter {
        let line = line?;
        if line.is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() <= pep_idx.max(rep_idx) {
            skipped += 1;
            continue;
        }
        out.push((parts[pep_idx].to_string(), parts[rep_idx].to_string()));
    }
    Ok((out, skipped))
}

/// Read 2-col cluster.tsv (rep<TAB>member) for the subset of reps we need.
/// Returns rep → Vec<member>.
fn load_cluster_subset(
    path: &Path,
    needed_reps: &AHashSet<String>,
) -> std::io::Result<(AHashMap<String, Vec<String>>, usize)> {
    let f = File::open(path)?;
    let r = BufReader::with_capacity(1 << 20, f);
    let mut clusters: AHashMap<String, Vec<String>> = AHashMap::new();
    let mut skipped = 0usize;
    for line in r.lines() {
        let line = line?;
        if line.is_empty() {
            continue;
        }
        let mut parts = line.split('\t');
        let rep = match parts.next() {
            Some(s) => s,
            None => {
                skipped += 1;
                continue;
            }
        };
        let member = match parts.next() {
            Some(s) => s,
            None => {
                skipped += 1;
                continue;
            }
        };
        if needed_reps.contains(rep) {
            clusters
                .entry(rep.to_string())
                .or_default()
                .push(member.to_string());
        }
    }
    Ok((clusters, skipped))
}

/// Stream unclustered lake.fasta once, keep only sequences whose id ∈ needed_members.
/// Returns id → uppercase byte sequence.
fn load_lake_subset(
    path: &Path,
    needed_members: &AHashSet<String>,
) -> std::io::Result<AHashMap<String, Vec<u8>>> {
    let f = File::open(path)?;
    let r = BufReader::with_capacity(1 << 20, f);
    let mut out: AHashMap<String, Vec<u8>> = AHashMap::with_capacity(needed_members.len());
    let mut current_id: Option<String> = None;
    let mut current_seq: Vec<u8> = Vec::with_capacity(512);
    let flush = |id: Option<String>, seq: &mut Vec<u8>, out: &mut AHashMap<String, Vec<u8>>| {
        if let Some(id) = id {
            if !seq.is_empty() {
                out.insert(id, std::mem::take(seq));
            }
        } else {
            seq.clear();
        }
    };
    for line in r.lines() {
        let line = line?;
        if line.is_empty() {
            continue;
        }
        if line.as_bytes()[0] == b'>' {
            flush(current_id.take(), &mut current_seq, &mut out);
            // Pick first whitespace-delimited token after '>'
            let header_body = &line[1..];
            let id_end = header_body
                .find(|c: char| c.is_whitespace())
                .unwrap_or(header_body.len());
            let id = &header_body[..id_end];
            if needed_members.contains(id) {
                current_id = Some(id.to_string());
            } else {
                current_id = None;
            }
        } else if current_id.is_some() {
            for &b in line.as_bytes() {
                if b.is_ascii_alphabetic() {
                    current_seq.push(b.to_ascii_uppercase());
                }
            }
        }
    }
    flush(current_id, &mut current_seq, &mut out);
    Ok(out)
}

fn main() {
    let args = Args::parse();

    if args.threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global()
            .unwrap();
    }

    let t0 = Instant::now();

    eprintln!("[1/4] Reading rescued.tsv: {}", args.rescued.display());
    let (rescued, rescued_skipped) = load_rescued(&args.rescued, &args.peptide_col, &args.rep_col)
        .expect("failed reading rescued");
    eprintln!("       {} (peptide, rep) pairs", rescued.len());

    let needed_reps: AHashSet<String> = rescued.iter().map(|(_, r)| r.clone()).collect();
    eprintln!(
        "       {} distinct cluster representatives",
        needed_reps.len()
    );

    eprintln!(
        "[2/4] Streaming cluster.tsv: {}",
        args.cluster_tsv.display()
    );
    let (clusters, cluster_skipped) =
        load_cluster_subset(&args.cluster_tsv, &needed_reps).expect("failed reading cluster_tsv");
    let n_members: usize = clusters.values().map(|v| v.len()).sum();
    eprintln!(
        "       {} clusters resolved, {} total members ({:.1} avg)",
        clusters.len(),
        n_members,
        n_members as f64 / clusters.len().max(1) as f64
    );

    let needed_members: AHashSet<String> =
        clusters.values().flat_map(|v| v.iter().cloned()).collect();
    eprintln!(
        "       {} distinct member proteins to fetch",
        needed_members.len()
    );

    eprintln!(
        "[3/4] Streaming lake.fasta: {} ({:.1}s elapsed)",
        args.lake.display(),
        t0.elapsed().as_secs_f64()
    );
    let member_seqs = load_lake_subset(&args.lake, &needed_members).expect("failed reading lake");
    eprintln!(
        "       {} member sequences loaded into memory ({:.1}s)",
        member_seqs.len(),
        t0.elapsed().as_secs_f64()
    );

    if member_seqs.len() < needed_members.len() {
        let missing = needed_members.len() - member_seqs.len();
        eprintln!(
            "       WARNING: {} members not found in lake (will be skipped)",
            missing
        );
    }

    eprintln!("[4/4] Lifting peptides ({} pairs)...", rescued.len());

    // D5 (determinism): classify every (peptide, rep) pair in parallel into an
    // owned result row, then sort by the total key (peptide, cluster_rep) and
    // write sequentially. Previously rows were streamed straight from
    // par_iter().for_each() through a shared Mutex<BufWriter>, so the on-disk row
    // order depended on rayon thread scheduling — non-deterministic and
    // thread-count-variant. Sorting on a total key makes output byte-identical
    // regardless of thread count. (The per-row members_with_hit CSV is built from
    // cluster.tsv Vec order, which is already deterministic — left unchanged.)
    let mut rows: Vec<(String, String, String)> = rescued
        .par_iter()
        .filter_map(|(peptide, rep)| {
            let members = clusters.get(rep)?;
            let pep_bytes = peptide.as_bytes();
            let finder = memmem::Finder::new(pep_bytes);
            let mut hits: Vec<&str> = Vec::new();
            // C9 (correctness): the SHARED decision divides by RESOLVABLE members
            // only — members whose sequence we actually loaded from the lake.
            // Members absent from the lake cannot be assessed (no sequence to
            // search), so they can never contribute a hit; counting them in the
            // denominator would deflate the hit fraction and bias the call toward
            // UNIQUE. We therefore use the resolvable-member count as the
            // denominator while still reporting the true total member count in the
            // n_cluster_members output column.
            let mut n_resolvable = 0usize;
            for m in members {
                if let Some(seq) = member_seqs.get(m) {
                    n_resolvable += 1;
                    if finder.find(seq).is_some() {
                        hits.push(m.as_str());
                    }
                }
            }
            let n_total = members.len();
            let n_hits = hits.len();
            // Zero resolvable members => conservation is unassessable; default to
            // UNIQUE and avoid a divide-by-zero.
            let is_shared =
                n_resolvable > 0 && (n_hits as f64 / n_resolvable as f64) >= args.shared_threshold;
            let members_csv = hits.join(",");
            let row = format!(
                "{}\t{}\t{}\t{}\t{}\t{}\n",
                peptide, rep, n_total, n_hits, members_csv, is_shared
            );
            Some((peptide.clone(), rep.clone(), row))
        })
        .collect();

    rows.sort_by(|a, b| a.0.cmp(&b.0).then_with(|| a.1.cmp(&b.1)));

    let out_file = File::create(&args.output).expect("cannot create output");
    let mut w = BufWriter::new(out_file);
    writeln!(
        w,
        "peptide\tcluster_rep\tn_cluster_members\tn_with_hit\tmembers_with_hit\tis_shared"
    )
    .unwrap();
    for (_, _, row) in &rows {
        w.write_all(row.as_bytes()).unwrap();
    }
    w.flush().unwrap();
    drop(w);

    let total_skipped = rescued_skipped + cluster_skipped;
    if total_skipped > 0 {
        eprintln!(
            "WARN: skipped {} malformed rows ({} in rescued.tsv, {} in cluster.tsv)",
            total_skipped, rescued_skipped, cluster_skipped
        );
    }

    eprintln!(
        "Done. Output: {} ({:.1}s total)",
        args.output.display(),
        t0.elapsed().as_secs_f64()
    );
}
