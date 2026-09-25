use clap::{Parser, ValueEnum};
use std::fs;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::PathBuf;

// ============================================================================
// CLI
// ============================================================================

const VERSION: &str = concat!(env!("CARGO_PKG_VERSION"), "+", env!("GIT_REV"));

#[derive(Parser)]
#[command(name = "fasta_export")]
#[command(version = VERSION)]
#[command(about = "Engine-aware FASTA export for FastaLake")]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(clap::Subcommand)]
enum Command {
    /// Prepare FASTA for a specific search engine
    Export {
        #[arg(short, long)]
        input: PathBuf,
        #[arg(short, long)]
        output: PathBuf,
        #[arg(short, long, value_enum)]
        engine: Engine,
    },
    /// Compute entrapment FDP from SAGE results
    Entrapment {
        #[arg(
            short,
            long,
            help = "Directory containing sample subdirs with results.sage.tsv"
        )]
        input: PathBuf,
        #[arg(short, long, help = "Entrapment prefix in protein IDs")]
        prefix: String,
        #[arg(long, help = "Number of target proteins in database")]
        target_db_size: usize,
        #[arg(long, help = "Number of entrapment proteins in database")]
        entrap_db_size: usize,
        #[arg(short, long, default_value = "0.01")]
        q_threshold: f64,
    },
}

#[derive(Clone, Debug, ValueEnum)]
enum Engine {
    Sage,
    Diann,
    Msfragger,
}

// ============================================================================
// Header parsing
// ============================================================================

/// Parsed FASTA header with extracted metadata.
#[derive(Debug, Clone)]
struct ParsedHeader {
    accession: String,
    gene: Option<String>,
    description: String,
    // Parsed from the OS= tag and asserted in tests; retained for completeness even
    // though the current exporter does not emit it. Kept intentionally, not dead.
    #[allow(dead_code)]
    organism: Option<String>,
    source_db: HeaderFormat,
    raw: String,
}

#[derive(Debug, Clone, PartialEq)]
enum HeaderFormat {
    SwissProt,
    TrEMBL,
    Mgyg,
    Gmgc,
    Prodigal,
    Ensembl,
    GnomAD,
    Unknown,
}

/// Parse a FASTA header line (without leading '>').
fn parse_header(raw: &str) -> ParsedHeader {
    let parts: Vec<&str> = raw.splitn(2, ' ').collect();
    let accession = parts[0].to_string();
    let description = if parts.len() > 1 {
        parts[1].to_string()
    } else {
        String::new()
    };

    // Extract GN= gene name if present
    let gn_gene = extract_gn(&description);

    // Detect format and extract metadata
    if accession.starts_with("sp|") {
        let gene = gn_gene.or_else(|| extract_entry_gene(&accession));
        let organism = extract_os(&description);
        ParsedHeader {
            accession,
            gene,
            description,
            organism,
            source_db: HeaderFormat::SwissProt,
            raw: raw.to_string(),
        }
    } else if accession.starts_with("tr|") {
        let gene = gn_gene.or_else(|| extract_entry_gene(&accession));
        let organism = extract_os(&description);
        ParsedHeader {
            accession,
            gene,
            description,
            organism,
            source_db: HeaderFormat::TrEMBL,
            raw: raw.to_string(),
        }
    } else if accession.starts_with("var-") || accession.starts_with("var_") {
        // gnomAD variant: var-nfe|P02768|ALB_Glu177Asp|...
        let gene = gn_gene.or_else(|| extract_gnomad_gene(&description));
        ParsedHeader {
            accession,
            gene,
            description,
            organism: None,
            source_db: HeaderFormat::GnomAD,
            raw: raw.to_string(),
        }
    } else if accession.starts_with("MGYG") {
        // UHGP: MGYG000001373_00501
        let gene = gn_gene.or_else(|| Some(accession.clone()));
        ParsedHeader {
            accession,
            gene,
            description,
            organism: None,
            source_db: HeaderFormat::Mgyg,
            raw: raw.to_string(),
        }
    } else if accession.starts_with("GMGC") {
        // GMGC: GMGC10.000_200_695.UNKNOWN
        let gene = gn_gene.or_else(|| extract_gmgc_gene(&accession));
        ParsedHeader {
            accession,
            gene,
            description,
            organism: None,
            source_db: HeaderFormat::Gmgc,
            raw: raw.to_string(),
        }
    } else if accession.starts_with("ENSP") || accession.starts_with("ENST") {
        let gene = gn_gene;
        ParsedHeader {
            accession,
            gene,
            description,
            organism: None,
            source_db: HeaderFormat::Ensembl,
            raw: raw.to_string(),
        }
    } else if accession.contains("_") && description.contains("# ") {
        // Prodigal: k141_130475_1 # 1 # 1404 # 1 # ID=...
        let gene = gn_gene.or_else(|| Some(accession.clone()));
        ParsedHeader {
            accession,
            gene,
            description,
            organism: None,
            source_db: HeaderFormat::Prodigal,
            raw: raw.to_string(),
        }
    } else {
        let gene = gn_gene.or_else(|| Some(accession.clone()));
        ParsedHeader {
            accession,
            gene,
            description,
            organism: None,
            source_db: HeaderFormat::Unknown,
            raw: raw.to_string(),
        }
    }
}

fn extract_gn(desc: &str) -> Option<String> {
    if let Some(pos) = desc.find("GN=") {
        let rest = &desc[pos + 3..];
        let end = rest.find(|c: char| c.is_whitespace()).unwrap_or(rest.len());
        let gene = &rest[..end];
        if !gene.is_empty() {
            Some(gene.to_string())
        } else {
            None
        }
    } else {
        None
    }
}

fn extract_os(desc: &str) -> Option<String> {
    if let Some(pos) = desc.find("OS=") {
        let rest = &desc[pos + 3..];
        let end = rest
            .find(" OX=")
            .or_else(|| rest.find(" GN="))
            .unwrap_or(rest.len());
        Some(rest[..end].trim().to_string())
    } else {
        None
    }
}

fn extract_entry_gene(accession: &str) -> Option<String> {
    // sp|P02768|ALBU_HUMAN → gene from entry name (ALBU)
    let parts: Vec<&str> = accession.split('|').collect();
    if parts.len() >= 3 {
        let entry = parts[2];
        if let Some(pos) = entry.find('_') {
            return Some(entry[..pos].to_string());
        }
    }
    None
}

fn extract_gnomad_gene(desc: &str) -> Option<String> {
    // Try GN= from description first, then parent sp| accession
    extract_gn(desc).or_else(|| {
        if let Some(sp_pos) = desc.find("sp|") {
            let rest = &desc[sp_pos..];
            let parts: Vec<&str> = rest.splitn(2, ' ').collect();
            extract_entry_gene(parts[0])
        } else {
            None
        }
    })
}

fn extract_gmgc_gene(accession: &str) -> Option<String> {
    // GMGC10.000_200_695.UNKNOWN → gene after last dot
    if let Some(pos) = accession.rfind('.') {
        let gene = &accession[pos + 1..];
        if gene != "UNKNOWN" && !gene.is_empty() {
            return Some(gene.to_string());
        }
    }
    Some(accession.to_string())
}

// ============================================================================
// Header healing (engine-specific formatting)
// ============================================================================

/// Heal a header for SAGE: add GN= if missing.
fn heal_for_sage(parsed: &ParsedHeader) -> String {
    if parsed.raw.contains("GN=") {
        return parsed.raw.clone();
    }
    let gene = parsed.gene.as_deref().unwrap_or(&parsed.accession);
    format!("{} GN={}", parsed.raw, gene)
}

/// Heal a header for DIA-NN: gene must be FIRST WORD of description
/// for non-sp/tr entries. sp|/tr| are handled correctly by DIA-NN.
fn heal_for_diann(parsed: &ParsedHeader) -> String {
    match parsed.source_db {
        HeaderFormat::SwissProt | HeaderFormat::TrEMBL => {
            // DIA-NN handles sp|/tr| natively — just ensure GN= present
            if parsed.raw.contains("GN=") {
                parsed.raw.clone()
            } else {
                let gene = parsed.gene.as_deref().unwrap_or(&parsed.accession);
                format!("{} GN={}", parsed.raw, gene)
            }
        }
        _ => {
            // All non-sp/tr: gene must be FIRST WORD of description
            let gene = parsed.gene.as_deref().unwrap_or(&parsed.accession);
            if parsed.description.is_empty() {
                format!("{} {} GN={}", parsed.accession, gene, gene)
            } else {
                format!(
                    "{} {} {} GN={}",
                    parsed.accession, gene, parsed.description, gene
                )
            }
        }
    }
}

/// Heal a header for MSFragger: add GN= if missing (decoys handled separately).
fn heal_for_msfragger(parsed: &ParsedHeader) -> String {
    // MSFragger is tolerant like SAGE
    heal_for_sage(parsed)
}

/// Heal a header for any supported engine.
fn heal_header(raw: &str, engine: &Engine) -> String {
    let parsed = parse_header(raw);
    match engine {
        Engine::Sage => heal_for_sage(&parsed),
        Engine::Diann => heal_for_diann(&parsed),
        Engine::Msfragger => heal_for_msfragger(&parsed),
    }
}

// ============================================================================
// Decoy generation
// ============================================================================

/// Generate reversed decoy sequences with a prefix.
/// Returns (n_targets, n_decoys).
fn generate_decoys(
    input: &std::path::Path,
    output: &std::path::Path,
    prefix: &str,
) -> std::io::Result<(usize, usize)> {
    let mut entries: Vec<(String, String)> = Vec::new(); // (header, sequence)
    let reader = BufReader::new(fs::File::open(input)?);
    let mut current_header = String::new();
    let mut current_seq = String::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('>') {
            if !current_header.is_empty() {
                entries.push((current_header.clone(), current_seq.clone()));
                current_seq.clear();
            }
            current_header = line;
        } else {
            current_seq.push_str(line.trim());
        }
    }
    if !current_header.is_empty() {
        entries.push((current_header, current_seq));
    }

    let n_targets = entries.len();
    let file = fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(output)?;
    let mut w = BufWriter::new(file);

    // Write targets
    for (header, seq) in &entries {
        writeln!(w, "{}", header)?;
        writeln!(w, "{}", seq)?;
    }

    // Write reversed decoys
    for (header, seq) in &entries {
        let acc_rest = &header[1..]; // strip >
        let parts: Vec<&str> = acc_rest.splitn(2, ' ').collect();
        let acc = parts[0];
        let desc = if parts.len() > 1 { parts[1] } else { "" };
        if desc.is_empty() {
            writeln!(w, ">{}{}", prefix, acc)?;
        } else {
            writeln!(w, ">{}{} {}", prefix, acc, desc)?;
        }
        let reversed: String = seq.chars().rev().collect();
        writeln!(w, "{}", reversed)?;
    }

    w.flush()?;
    Ok((n_targets, n_targets))
}

// ============================================================================
// FASTA export: heal + optionally add decoys
// ============================================================================

fn export_fasta(
    input: &std::path::Path,
    output: &std::path::Path,
    engine: &Engine,
) -> std::io::Result<(usize, usize)> {
    match engine {
        Engine::Msfragger => {
            // A private directory owns its intermediate; never touch a predictable
            // sibling path that may belong to the caller.
            let parent = output
                .parent()
                .filter(|p| !p.as_os_str().is_empty())
                .unwrap_or_else(|| std::path::Path::new("."));
            let temporary = tempfile::Builder::new()
                .prefix(".fasta-export-")
                .tempdir_in(parent)?;
            let tmp = temporary.path().join("healed.fasta");
            heal_fasta(input, &tmp, engine)?;
            let (n_t, n_d) = generate_decoys(&tmp, output, "rev_")?;
            Ok((n_t, n_d))
        }
        _ => {
            // SAGE, DIA-NN: heal headers only, no decoys
            let n = heal_fasta(input, output, engine)?;
            Ok((n, 0))
        }
    }
}

/// Heal all headers in a FASTA file. Returns count of proteins written.
fn heal_fasta(
    input: &std::path::Path,
    output: &std::path::Path,
    engine: &Engine,
) -> std::io::Result<usize> {
    let reader = BufReader::new(fs::File::open(input)?);
    let file = fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(output)?;
    let mut w = BufWriter::new(file);
    let mut count = 0;

    for line in reader.lines() {
        let line = line?;
        if let Some(raw) = line.strip_prefix('>') {
            let healed = heal_header(raw, engine);
            writeln!(w, ">{}", healed)?;
            count += 1;
        } else {
            writeln!(w, "{}", line)?;
        }
    }

    w.flush()?;
    Ok(count)
}

// ============================================================================
// Entrapment FDP
// ============================================================================

/// Compute False Discovery Proportion using the combined method
/// (Eq 1, Wen, Keich & Noble 2025 Nature Methods).
///
/// r = entrap_db_size / target_db_size (Noble convention)
/// FDP = N_E × (1 + 1/r) / (N_T + N_E)
fn compute_fdp_combined(
    n_entrapment: usize,
    n_target: usize,
    target_db_size: usize,
    entrap_db_size: usize,
) -> f64 {
    let n_total = n_entrapment + n_target;
    if n_total == 0 || target_db_size == 0 {
        return 0.0;
    }
    let r = entrap_db_size as f64 / target_db_size as f64;
    if r <= 0.0 {
        return 0.0;
    }
    (n_entrapment as f64 * (1.0 + 1.0 / r)) / n_total as f64
}

/// Lower bound FDP (insufficient alone — only proves failure, not success).
fn compute_fdp_lower(n_entrapment: usize, n_target: usize) -> f64 {
    let n_total = n_entrapment + n_target;
    if n_total == 0 {
        0.0
    } else {
        n_entrapment as f64 / n_total as f64
    }
}

/// Parse a SAGE results.sage.tsv and count entrapment hits at a given q threshold.
fn count_entrapment_hits(
    results_path: &std::path::Path,
    entrap_prefix: &str,
    q_threshold: f64,
) -> std::io::Result<(usize, usize)> {
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(results_path)?;

    let headers = reader.headers()?.clone();
    let proteins_idx = headers.iter().position(|h| h == "proteins");
    let spectrum_q_idx = headers.iter().position(|h| h == "spectrum_q");
    let label_idx = headers.iter().position(|h| h == "label");

    let (proteins_idx, spectrum_q_idx) = match (proteins_idx, spectrum_q_idx) {
        (Some(p), Some(q)) => (p, q),
        _ => return Ok((0, 0)),
    };

    let mut n_target = 0usize;
    let mut n_entrap = 0usize;

    for record in reader.records() {
        let record = record?;

        // Skip decoys
        if let Some(li) = label_idx {
            if let Some(label) = record.get(li) {
                if label == "-1" {
                    continue;
                }
            }
        }

        let q: f64 = record
            .get(spectrum_q_idx)
            .and_then(|v| v.parse().ok())
            .unwrap_or(1.0);
        if q > q_threshold {
            continue;
        }

        let proteins = record.get(proteins_idx).unwrap_or("");
        if proteins
            .split(';')
            .any(|p| p.trim().starts_with(entrap_prefix))
        {
            n_entrap += 1;
        } else {
            n_target += 1;
        }
    }

    Ok((n_target, n_entrap))
}

// ============================================================================
// Main
// ============================================================================

fn main() {
    let cli = Cli::parse();

    match cli.command {
        Command::Export {
            input,
            output,
            engine,
        } => {
            eprintln!("FastaLake fasta_export v1.0.0");
            eprintln!("Engine: {:?}", engine);
            eprintln!("Input:  {}", input.display());
            eprintln!("Output: {}", output.display());

            match export_fasta(&input, &output, &engine) {
                Ok((n_targets, n_decoys)) => {
                    eprintln!("Targets: {}", n_targets);
                    if n_decoys > 0 {
                        eprintln!("Decoys:  {}", n_decoys);
                    }
                    eprintln!("Total:   {}", n_targets + n_decoys);
                    eprintln!("Done.");
                }
                Err(e) => {
                    eprintln!("ERROR: {}", e);
                    std::process::exit(1);
                }
            }
        }
        Command::Entrapment {
            input,
            prefix,
            target_db_size,
            entrap_db_size,
            q_threshold,
        } => {
            eprintln!("FastaLake entrapment FDP v1.0.0");
            eprintln!("Prefix: {}", prefix);
            eprintln!(
                "DB sizes: {} target, {} entrapment",
                target_db_size, entrap_db_size
            );
            eprintln!("q threshold: {}", q_threshold);

            // Guard: zero DB sizes would silently certify FDP=0.0 (mis-set CLI arg).
            if target_db_size == 0 {
                eprintln!("ERROR: --target-db-size must be > 0 (got 0)");
                std::process::exit(2);
            }
            if entrap_db_size == 0 {
                eprintln!("ERROR: --entrap-db-size must be > 0 (got 0)");
                std::process::exit(2);
            }

            // Find all results.sage.tsv files
            let mut results: Vec<(String, usize, usize, f64, f64)> = Vec::new();
            // read_dir order is filesystem-dependent; the per-sample table and the
            // summary mean (a float sum) must not change with it.
            let mut entries: Vec<_> = fs::read_dir(&input)
                .expect("Cannot read input directory")
                .filter_map(|e| e.ok())
                .filter(|e| e.path().is_dir())
                .collect();
            entries.sort_by_key(|e| e.file_name());

            for entry in &entries {
                let sage_tsv = entry.path().join("results.sage.tsv");
                if !sage_tsv.exists() {
                    continue;
                }

                let sample = entry.file_name().to_string_lossy().to_string();
                match count_entrapment_hits(&sage_tsv, &prefix, q_threshold) {
                    Ok((n_t, n_e)) => {
                        let fdp_lower = compute_fdp_lower(n_e, n_t);
                        let fdp_combined =
                            compute_fdp_combined(n_e, n_t, target_db_size, entrap_db_size);
                        results.push((sample, n_t, n_e, fdp_lower, fdp_combined));
                    }
                    Err(e) => eprintln!("  WARN: {} — {}", sample, e),
                }
            }

            // Print results
            println!("sample\tn_target\tn_entrapment\tfdp_lower\tfdp_combined");
            for (s, n_t, n_e, fl, fc) in &results {
                println!("{}\t{}\t{}\t{:.6}\t{:.6}", s, n_t, n_e, fl, fc);
            }

            // Summary
            if !results.is_empty() {
                let fdps: Vec<f64> = results.iter().map(|r| r.4).collect();
                let mean = fdps.iter().sum::<f64>() / fdps.len() as f64;
                let mut sorted = fdps.clone();
                sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
                let median = sorted[sorted.len() / 2];
                let max = sorted.last().copied().unwrap_or(0.0);
                eprintln!("\nSummary ({} samples):", results.len());
                eprintln!(
                    "  FDP combined: mean={:.4}, median={:.4}, max={:.4}",
                    mean, median, max
                );
            }
        }
    }
}

// ============================================================================
// TESTS
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    // --- Header parsing ---

    #[test]
    fn test_parse_swissprot() {
        let h =
            parse_header("sp|P02768|ALBU_HUMAN Albumin OS=Homo sapiens OX=9606 GN=ALB PE=1 SV=2");
        assert_eq!(h.source_db, HeaderFormat::SwissProt);
        assert_eq!(h.accession, "sp|P02768|ALBU_HUMAN");
        assert_eq!(h.gene, Some("ALB".into()));
        assert_eq!(h.organism, Some("Homo sapiens".into()));
    }

    #[test]
    fn test_parse_swissprot_no_gn() {
        let h = parse_header(
            "sp|P0DOX5|IGG1_HUMAN Immunoglobulin gamma-1 OS=Homo sapiens OX=9606 PE=1 SV=2",
        );
        assert_eq!(h.source_db, HeaderFormat::SwissProt);
        assert_eq!(h.gene, Some("IGG1".into())); // from entry name
    }

    #[test]
    fn test_parse_trembl() {
        let h = parse_header("tr|A0A218KGR2|A0A218KGR2_HUMAN Amyloid OS=Homo sapiens GN=APP");
        assert_eq!(h.source_db, HeaderFormat::TrEMBL);
        assert_eq!(h.gene, Some("APP".into()));
    }

    #[test]
    fn test_parse_mgyg() {
        let h = parse_header("MGYG000001373_00501 Uncharacterized protein");
        assert_eq!(h.source_db, HeaderFormat::Mgyg);
        assert_eq!(h.gene, Some("MGYG000001373_00501".into()));
    }

    #[test]
    fn test_parse_gmgc() {
        let h = parse_header("GMGC10.000_200_695.DNAK hypothetical protein");
        assert_eq!(h.source_db, HeaderFormat::Gmgc);
        assert_eq!(h.gene, Some("DNAK".into()));
    }

    #[test]
    fn test_parse_gmgc_unknown() {
        let h = parse_header("GMGC10.000_200_695.UNKNOWN hypothetical protein");
        assert_eq!(h.source_db, HeaderFormat::Gmgc);
        // UNKNOWN → falls back to full accession
        assert_eq!(h.gene, Some("GMGC10.000_200_695.UNKNOWN".into()));
    }

    #[test]
    fn test_parse_gnomad() {
        let h = parse_header(
            "var-nfe|P02768|ALB_Glu177Asp|4-73409403-A-T sp|P02768|ALBU_HUMAN Albumin GN=ALB",
        );
        assert_eq!(h.source_db, HeaderFormat::GnomAD);
        assert_eq!(h.gene, Some("ALB".into()));
    }

    // --- Header healing: SAGE ---

    #[test]
    fn test_heal_sage_adds_gn() {
        let healed = heal_header("MGYG000001373_00501 Uncharacterized protein", &Engine::Sage);
        assert!(healed.contains("GN=MGYG000001373_00501"));
    }

    #[test]
    fn test_heal_sage_preserves_existing_gn() {
        let raw = "sp|P02768|ALBU_HUMAN Albumin GN=ALB";
        let healed = heal_header(raw, &Engine::Sage);
        assert_eq!(healed, raw); // unchanged
    }

    // --- Header healing: DIA-NN ---

    #[test]
    fn test_heal_diann_swissprot_unchanged() {
        let raw = "sp|P02768|ALBU_HUMAN Albumin OS=Homo sapiens GN=ALB";
        let healed = heal_header(raw, &Engine::Diann);
        assert_eq!(healed, raw); // sp| entries unchanged
    }

    #[test]
    fn test_heal_diann_mgyg_gene_first_word() {
        let healed = heal_header(
            "MGYG000001373_00501 Uncharacterized protein",
            &Engine::Diann,
        );
        let parts: Vec<&str> = healed.splitn(3, ' ').collect();
        assert_eq!(parts[0], "MGYG000001373_00501"); // accession
        assert_eq!(parts[1], "MGYG000001373_00501"); // gene as first word
        assert!(healed.contains("GN=MGYG000001373_00501"));
    }

    #[test]
    fn test_heal_diann_gmgc_gene_first_word() {
        let healed = heal_header(
            "GMGC10.050_359_667.DNAK hypothetical protein",
            &Engine::Diann,
        );
        let parts: Vec<&str> = healed.splitn(3, ' ').collect();
        assert_eq!(parts[1], "DNAK"); // gene as first word of description
    }

    // --- Header healing: MSFragger ---

    #[test]
    fn test_heal_msfragger_adds_gn() {
        let healed = heal_header(
            "MGYG000001373_00501 Uncharacterized protein",
            &Engine::Msfragger,
        );
        assert!(healed.contains("GN=MGYG000001373_00501"));
    }

    // --- Decoy generation ---

    #[test]
    fn test_generate_decoys() {
        let dir = tempfile::tempdir().unwrap();
        let input = dir.path().join("input.fasta");
        let output = dir.path().join("output.fasta");

        fs::write(
            &input,
            ">sp|P02768|ALBU_HUMAN Albumin\nMKWVTFISL\n>MGYG_001 test\nACDEFG\n",
        )
        .unwrap();

        let (n_t, n_d) = generate_decoys(&input, &output, "rev_").unwrap();
        assert_eq!(n_t, 2);
        assert_eq!(n_d, 2);

        let content = fs::read_to_string(&output).unwrap();
        let headers: Vec<&str> = content.lines().filter(|l| l.starts_with('>')).collect();
        assert_eq!(headers.len(), 4);
        assert!(headers[2].starts_with(">rev_sp|P02768|ALBU_HUMAN"));
        assert!(headers[3].starts_with(">rev_MGYG_001"));

        // Check sequence reversal
        let lines: Vec<&str> = content.lines().collect();
        assert_eq!(lines[1], "MKWVTFISL"); // target
        assert_eq!(lines[5], "LSIFTVWKM"); // reversed
        assert_eq!(lines[3], "ACDEFG"); // target
        assert_eq!(lines[7], "GFEDCA"); // reversed
    }

    // --- FDP computation ---

    #[test]
    fn test_fdp_combined_r_equals_1() {
        // r=1 (equal DBs): combined = 2 × lower bound
        let fdp = compute_fdp_combined(10, 990, 1000, 1000);
        let expected = 10.0 * 2.0 / 1000.0; // 0.02
        assert!((fdp - expected).abs() < 1e-10);
    }

    #[test]
    fn test_fdp_combined_r_less_1() {
        // r=0.1: entrapment much smaller than target
        // combined = 10 * (1 + 10) / 1000 = 0.11
        let fdp = compute_fdp_combined(10, 990, 50000, 5000);
        let r = 5000.0 / 50000.0; // 0.1
        let expected = 10.0 * (1.0 + 1.0 / r) / 1000.0; // 0.11
        assert!((fdp - expected).abs() < 1e-10);
    }

    #[test]
    fn test_fdp_combined_large_r() {
        // r=10: entrapment much larger than target
        // combined = 10 * (1 + 0.1) / 1000 = 0.011
        let fdp = compute_fdp_combined(10, 990, 1000, 10000);
        let r = 10000.0 / 1000.0; // 10
        let expected = 10.0 * (1.0 + 1.0 / r) / 1000.0;
        assert!((fdp - expected).abs() < 1e-10);
    }

    #[test]
    fn test_fdp_zero_entrapment() {
        let fdp = compute_fdp_combined(0, 1000, 50000, 5000);
        assert_eq!(fdp, 0.0);
    }

    #[test]
    fn test_fdp_zero_total() {
        let fdp = compute_fdp_combined(0, 0, 50000, 5000);
        assert_eq!(fdp, 0.0);
    }

    #[test]
    fn test_fdp_lower_bound() {
        let fdp = compute_fdp_lower(10, 990);
        assert!((fdp - 0.01).abs() < 1e-10);
    }
}
