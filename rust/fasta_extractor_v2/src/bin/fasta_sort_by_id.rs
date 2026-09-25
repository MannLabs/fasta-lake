//! Sort a FASTA deterministically by accession, so record order is a function of
//! sequence content and nothing else.
//!
//! WHY THIS EXISTS
//!
//! `lake_builder` writes records in the order it reads sources, and sources are read
//! in the order they are listed. So the lake's record order encodes which catalogue
//! was named first. Permuting the `-s` specs produces the same sequences with the
//! same accessions in a DIFFERENT order (control gate G9a).
//!
//! That order is not cosmetic. It becomes the clustering input, and MMseqs2 is run
//! with `--shuffle 0` so it respects input order; in greedy set-cover, tied candidate
//! representatives are resolved by the lower database key. Whichever source was
//! listed first therefore wins representatives systematically — and the
//! representative's accession is what `db_breakdown` reports for the whole cluster.
//! Measured consequence of changing this order: ~18% depth (defect D23).
//!
//! Sorting by accession removes it. With `lake_builder --uniform-tag TAG`, every
//! accession is `TAG_<sha256hex>`, so sorting by accession IS sorting by sequence
//! content: a total order that no source ordering can influence.
//!
//! Deliberately a separate tool rather than a `lake_builder` flag: the lake is tens
//! of gigabytes and sorting it in place would need the whole record index in memory,
//! whereas the file that actually feeds clustering (the evidence lake) is far
//! smaller. Sorting is applied where the bias acts.
//!
//! Guarantees:
//!   * output is a permutation of the input — no record added, dropped or altered;
//!   * output order depends only on accessions, so it is stable across runs, thread
//!     counts, and any permutation of the upstream source list;
//!   * duplicate accessions are an ERROR, not a silent merge.

use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::PathBuf;

use clap::Parser;

#[derive(Parser, Debug)]
#[command(about = "Sort a FASTA by accession so record order is content-determined")]
struct Args {
    /// Input FASTA
    #[arg(short = 'i', long)]
    input: PathBuf,

    /// Output FASTA (sorted by accession, ascending byte order)
    #[arg(short = 'o', long)]
    output: PathBuf,

    /// Permit duplicate accessions instead of failing. Off by default: two records
    /// sharing an accession means the identifier is not unique, and silently keeping
    /// both would let a downstream join match the wrong one.
    #[arg(long, default_value_t = false)]
    allow_duplicate_ids: bool,
}

fn main() -> std::io::Result<()> {
    let args = Args::parse();

    let f = File::open(&args.input)?;
    let mut rdr = BufReader::with_capacity(32 * 1024 * 1024, f);

    // (accession, full record including its '>' header line and sequence lines)
    let mut records: Vec<(String, String)> = Vec::new();
    let mut cur_id: Option<String> = None;
    let mut cur = String::new();
    let mut line = String::new();

    let flush = |records: &mut Vec<(String, String)>, id: &mut Option<String>, buf: &mut String| {
        if let Some(i) = id.take() {
            records.push((i, std::mem::take(buf)));
        }
    };

    loop {
        line.clear();
        if rdr.read_line(&mut line)? == 0 {
            break;
        }
        if let Some(header) = line.strip_prefix('>') {
            flush(&mut records, &mut cur_id, &mut cur);
            // Accession is the first whitespace-delimited token after '>'.
            let acc = header.split_whitespace().next().unwrap_or("").to_string();
            cur_id = Some(acc);
            cur.push_str(&line);
        } else {
            if cur_id.is_none() {
                eprintln!("ERROR: sequence data before the first '>' header line.");
                std::process::exit(2);
            }
            cur.push_str(&line);
        }
    }
    flush(&mut records, &mut cur_id, &mut cur);

    if records.is_empty() {
        eprintln!("ERROR: no FASTA records found in {}", args.input.display());
        std::process::exit(2);
    }

    if !args.allow_duplicate_ids {
        let mut seen: HashSet<&str> = HashSet::with_capacity(records.len());
        let mut dupes = 0usize;
        let mut first_dupe = String::new();
        for (id, _) in &records {
            if !seen.insert(id.as_str()) {
                dupes += 1;
                if first_dupe.is_empty() {
                    first_dupe = id.clone();
                }
            }
        }
        if dupes > 0 {
            eprintln!(
                "ERROR: {} duplicate accession(s), first '{}'. Refusing to sort: a \
                 non-unique identifier cannot give a well-defined order, and a \
                 downstream join could match the wrong record. Pass \
                 --allow-duplicate-ids only if that is genuinely acceptable.",
                dupes, first_dupe
            );
            std::process::exit(3);
        }
    }

    let n_in = records.len();
    // Sort on the accession alone. With TAG_<sha256hex> accessions this is a sort on
    // sequence content, so no upstream source ordering can influence it.
    records.sort_by(|a, b| a.0.cmp(&b.0));

    let out = File::create(&args.output)?;
    let mut w = BufWriter::with_capacity(32 * 1024 * 1024, out);
    for (_, rec) in &records {
        w.write_all(rec.as_bytes())?;
    }
    w.flush()?;

    eprintln!(
        "sorted {} records by accession -> {}",
        n_in,
        args.output.display()
    );
    Ok(())
}
