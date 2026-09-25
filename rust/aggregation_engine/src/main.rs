use clap::{Parser, ValueEnum};
use csv::ReaderBuilder;
use indexmap::IndexMap;
use rayon::prelude::*;
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::fs;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

// ============================================================================
// CLI
// ============================================================================

/// Default functional levels = all eggNOG-emapper columns we commonly want.
/// Any missing column silently produces a trivial matrix (see load_functional_lookup).
const DEFAULT_FUNCTIONAL_LEVELS: &str =
    "OG,KEGG_ko,KEGG_Pathway,KEGG_Module,KEGG_Reaction,COG_category,CAZy,Pfam,GOs,BRITE,EC";

const VERSION: &str = concat!(env!("CARGO_PKG_VERSION"), "+", env!("GIT_REV"));

#[derive(Parser)]
#[command(name = "aggregation_engine")]
#[command(version = VERSION)]
#[command(about = "Cross-sample aggregation for metaproteomics (FastaLake)")]
struct Cli {
    #[arg(short, long)]
    input: PathBuf,

    #[arg(short, long)]
    output: PathBuf,

    #[arg(
        short,
        long,
        help = "Dataset/project identifier (e.g., MicrobPredict, CAPSCAN)"
    )]
    dataset: String,

    #[arg(
        short,
        long,
        default_value = "FastaLake",
        help = "Method name for output files (e.g., FastaLake, MetaG, MetaLab)"
    )]
    method: String,

    #[arg(long, value_enum)]
    data_type: DataType,

    #[arg(short, long, default_value_t = 1)]
    threads: usize,

    #[arg(
        long,
        default_value_t = false,
        help = "Unsupported legacy option: q-values cannot be treated as p-values for global BH"
    )]
    bh_correction: bool,

    #[arg(
        long,
        help = "Path to functional lookup TSV with columns: accession + functional columns. First line must be header."
    )]
    functional_lookup: Option<PathBuf>,

    /// Host accessions for a hash-keyed build.
    ///
    /// Stage 0 under `--uniform-tag` renames every protein to `<COHORT>_<sha256>`,
    /// which carries no source, so `classify_group` cannot tell host from
    /// microbial and every group reads Unknown. This supplies the missing fact.
    ///
    /// The file lists ONLY the host accessions, because they are a rounding error
    /// of the lake: CapScan has 20,342 host sequences against 264,925,317 total.
    /// Absence from the list therefore means microbial. That is a contract backed
    /// by a complete pass over the stage-0 sidecar, not the silent default the
    /// classifier used to make from an unreadable accession.
    #[arg(
        long,
        help = "TSV of host accessions (accession, classification) for hash-keyed builds. Absence means microbial."
    )]
    host_accessions: Option<PathBuf>,

    #[arg(
        long,
        help = "Comma-separated functional level column names. Default: all eggNOG columns we support. Pass '' to disable."
    )]
    functional_levels: Option<String>,

    #[arg(
        long,
        default_value_t = false,
        help = "Disable functional aggregation entirely."
    )]
    no_functional: bool,

    /// Pre-compute every functional roll-up matrix.
    ///
    /// OFF BY DEFAULT, and that default is the point. A roll-up is the core matrix
    /// grouped by one functional annotation, so it is fully determined by two files that
    /// are always written: the core matrix and `<dataset>_group_to_<level>_map.tsv`.
    /// Pre-computing them writes 11 levels x 2 modes x 2 granularities x 4 thresholds =
    /// 176 extra files.
    ///
    /// Measured on PREDICT, 2026-08-03: those 176 files are 155.1 GB of a 178.3 GB
    /// output (87%), and writing them IS the runtime. The write phase is serial and
    /// formats billions of f64s into decimal text -- one peptide matrix is 4,281,281
    /// rows x 306 columns = 1.31 billion cells, and 92 of the 198 files are at peptide
    /// level. Sustained throughput was 19 MB/s on one core at 100%, giving 2.6 h.
    /// Dropping the roll-ups leaves 22 files and ~17 GB, so a full cohort aggregates in
    /// roughly 20 minutes.
    ///
    /// Nothing in the figure layer reads a roll-up (verified 2026-08-03: the only
    /// reference to `_by_` in the extraction code is the clause that EXCLUDES them,
    /// added after a plain glob silently fed a figure a 2,440-row BRITE roll-up in place
    /// of its 14,629-group matrix). Four near-identical files matching one pattern is
    /// how that bug happened, so writing fewer of them removes the bug class as well as
    /// the bytes.
    #[arg(
        long,
        default_value_t = false,
        help = "Pre-compute all functional roll-up matrices (11 levels x 2 modes x 2 \
                  granularities x 4 thresholds = 176 extra files, ~87% of output size and \
                  ~8x the runtime). Off by default: every roll-up is derivable from the \
                  core matrix plus the group_to_<level> map, both of which are always \
                  written."
    )]
    functional_rollups: bool,

    #[arg(long, value_enum, default_value_t = GroupKeyMode::MemberSet,
          help = "Protein-group key mode. 'member-set' (default): key by the full \
                  sorted, deduped set of member accessions (exact legacy behavior). \
                  'representative': collapse co-occurring near-identical orthologs by \
                  keying on a single canonical anchor (the lexicographically-smallest \
                  member accession), improving cross-sample comparability.")]
    group_key_mode: GroupKeyMode,

    #[arg(
        long,
        help = "Comma-separated list of protein/peptide FDR q-value thresholds at \
                  which to emit matrices (e.g. '0.01,0.05'). The literal token \
                  'unfiltered' (or 'none') requests the unfiltered matrix. \
                  Default (when omitted): '0.01,0.05,unfiltered', which reproduces \
                  the historical output exactly."
    )]
    q_thresholds: Option<String>,
}

/// Parse the `--q-thresholds` value into the (Option<f64>, label) pairs the
/// matrix writers consume. `None` (flag omitted) reproduces the historical
/// default set: q<=0.01, q<=0.05, and unfiltered — in that exact order, so
/// default output is byte-identical.
///
/// Accepted tokens: a numeric value in (0, 1] (e.g. `0.01`) emits a q<=value
/// matrix labelled `q<value>`; the literal `unfiltered` or `none` (case
/// insensitive) emits the unfiltered matrix.
fn parse_q_thresholds(spec: Option<&str>) -> Result<Vec<(Option<f64>, String)>, String> {
    let spec = match spec {
        None => {
            return Ok(vec![
                (Some(0.01), "q0.01".to_string()),
                (Some(0.05), "q0.05".to_string()),
                (None, "unfiltered".to_string()),
            ])
        }
        Some(s) => s,
    };
    let mut out: Vec<(Option<f64>, String)> = Vec::new();
    let mut seen: HashSet<String> = HashSet::new();
    for raw in spec.split(',') {
        let tok = raw.trim();
        if tok.is_empty() {
            continue;
        }
        let (thresh, label): (Option<f64>, String) =
            if tok.eq_ignore_ascii_case("unfiltered") || tok.eq_ignore_ascii_case("none") {
                (None, "unfiltered".to_string())
            } else {
                let v: f64 = tok.parse().map_err(|_| {
                    format!("invalid q-threshold token '{tok}' (expected a number or 'unfiltered')")
                })?;
                if !(v > 0.0 && v <= 1.0) {
                    return Err(format!("q-threshold {v} out of range (0, 1]"));
                }
                (Some(v), format!("q{v}"))
            };
        // Preserve order, drop duplicates by label.
        if seen.insert(label.clone()) {
            out.push((thresh, label));
        }
    }
    if out.is_empty() {
        return Err("--q-thresholds resolved to an empty set".to_string());
    }
    Ok(out)
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
enum GroupKeyMode {
    /// Key a protein group by the full sorted, deduped set of member accessions.
    /// This is the original (legacy) behavior and the default; it must not change
    /// any existing output.
    MemberSet,
    /// Key a protein group by a single CANONICAL anchor that is stable across
    /// samples: the lexicographically-smallest stripped member accession. This
    /// collapses co-occurring near-identical orthologs (which differ only in
    /// which members happen to be observed in a given sample) into one stable key,
    /// optimizing cross-sample comparability.
    Representative,
}

#[derive(Clone, Debug, ValueEnum)]
enum DataType {
    Lfq,
    Psm,
}

#[derive(Clone, Copy, Debug)]
enum AggregationMode {
    /// Assign the group's full intensity to its plurality-winner value.
    /// Conserves total intensity but does not apportion annotation ambiguity.
    Plurality,
    /// Split the group's intensity equally across the N unique values its
    /// members collectively point to (mass-balance correct).
    SplitIntensity,
}

impl AggregationMode {
    fn short(self) -> &'static str {
        match self {
            AggregationMode::Plurality => "plurality",
            AggregationMode::SplitIntensity => "splitintensity",
        }
    }
}

// ============================================================================
// Data structures
// ============================================================================

/// One peptide observation in one sample
struct PeptideObs {
    peptide: String,
    protein_group: String, // normalized: sorted, deduped, semicolon-joined
    intensity: f64,
    q_value: f64,    // LFQ: q_value; PSM: unused (use q_values below)
    spectrum_q: f64, // PSM only
    peptide_q: f64,  // PSM only
    protein_q: f64,  // PSM only
}

struct SampleData {
    sample_name: String,
    observations: Vec<PeptideObs>,
    n_peptides: usize,
    n_groups: usize,
}

// ============================================================================
// Source tag stripping & normalization
// ============================================================================

fn strip_source_tags(protein_id: &str) -> &str {
    let parts: Vec<&str> = protein_id.splitn(3, '|').collect();
    if parts.len() >= 3 {
        let first = parts[0];
        let second = parts[1];
        let is_tag = |s: &str| {
            !s.is_empty()
                && s.chars()
                    .all(|c| c.is_ascii_uppercase() || c.is_ascii_digit() || c == '_')
        };
        if is_tag(first) && is_tag(second) {
            return parts[2];
        }
    }
    protein_id
}

/// Normalize a raw `proteins` field into the canonical protein-group key.
///
/// `mode` selects how the key is formed from the deduped, source-tag-stripped
/// member accessions:
/// - `MemberSet` (default): the full set, sorted and semicolon-joined. This is
///   the original behavior; output is byte-for-byte identical to the legacy code.
/// - `Representative`: a single canonical anchor — the lexicographically-smallest
///   member accession. Because the smallest member is stable regardless of which
///   other near-identical orthologs happen to be observed in a given sample,
///   co-occurring orthologs collapse to one comparable key across samples.
///
/// `BTreeSet` iterates in ascending (lexicographic) order, so for `Representative`
/// the first element is exactly the lexicographically-smallest accession.
fn normalize_protein_group(proteins_field: &str, mode: GroupKeyMode) -> String {
    let mut members: BTreeSet<&str> = BTreeSet::new();
    for p in proteins_field.split(';') {
        let stripped = strip_source_tags(p.trim());
        if !stripped.is_empty() {
            members.insert(stripped);
        }
    }
    match mode {
        GroupKeyMode::MemberSet => members.into_iter().collect::<Vec<_>>().join(";"),
        GroupKeyMode::Representative => {
            // BTreeSet's first element is the lexicographically-smallest accession.
            members.into_iter().next().unwrap_or("").to_string()
        }
    }
}

/// Can this accession be read for provenance at all?
///
/// Every host test below keys on a prefix — `sp|`, `tr|`, `GENCODE|` and so on.
/// Stage 0 run with `--uniform-tag` renames every protein to `<COHORT>_<sha256>`,
/// which carries no prefix and no source. Such a member is not evidence of a
/// microbial origin; it is an absence of evidence, and the `else` branch of a
/// host test silently converts one into the other.
///
/// Found 2026-08-02, with v27 aggregation still queued: v26 CAPSCAN labelled all
/// 2,041,914 groups "Microbial", against v25's Human 20,596 / Mixed 4,336 /
/// Microbial 2,060,551 on the same cohort. The host fraction became zero and
/// nothing complained, because "Microbial" is a plausible-looking answer.
fn is_classifiable_accession(member: &str) -> bool {
    if let Some((_, rest)) = member.rsplit_once('_') {
        if rest.len() == 64 && rest.bytes().all(|b| b.is_ascii_hexdigit()) {
            return false;
        }
    }
    true
}

fn classify_group(group_key: &str) -> &'static str {
    let mut has_human = false;
    let mut has_micro = false;
    let mut any_readable = false;
    for member in group_key.split(';') {
        if !is_classifiable_accession(member) {
            continue;
        }
        any_readable = true;
        if is_host_member(member) {
            has_human = true;
        } else {
            has_micro = true;
        }
        if has_human && has_micro {
            return "Mixed";
        }
    }
    // Report the absence rather than guessing past it. Host/microbial for a
    // hash-keyed build must be derived from the stage-0 joint_headers sidecar,
    // which is the only place the source survives.
    if !any_readable {
        return "Unknown";
    }
    if has_human {
        "Human"
    } else {
        "Microbial"
    }
}

/// Host/human membership by SOURCE, not a blanket `sp|`/`tr|` prefix (microbial
/// UniProt entries are also `tr|`). Two host provenances:
///   1. UniProt human canonical/isoforms — `sp|`/`tr|` accession containing `_HUMAN`.
///   2. Tissue-isoform union (header `SOURCE|...`), sources GENCODE / RefSeq / CHESS /
///      UniProt-varsplic — human by construction. These prefixes are unique to the
///      isoform union; no metagenomic catalogue (MGYG/GMGC/GMSC/WP_) collides with them.
fn is_host_member(member: &str) -> bool {
    if (member.starts_with("sp|") || member.starts_with("tr|")) && member.contains("_HUMAN") {
        return true;
    }
    member.starts_with("GENCODE|")
        || member.starts_with("REFSEQ|")
        || member.starts_with("CHESS|")
        || member.starts_with("UNIPROT|")
}

/// An ARBITRARY, DETERMINISTIC DISPLAY LABEL for a protein group.
///
/// ⚠ THIS VALUE ENCODES NO ORIGIN INFORMATION. ⚠
///
/// It is the lexicographically-smallest member accession and nothing more. It is
/// NOT a "best" member, NOT a preferred database entry, and it MUST NEVER be used
/// for host/microbial attribution, for kingdom assignment, for source-catalogue
/// attribution, or for taxonomic assignment. A group's kingdom is the
/// `classification` column (`classify_group`), which is Human / Mixed / Microbial
/// and is the ONLY sanctioned basis for host-vs-microbial attribution. Taxonomy
/// and source catalogue must be derived over ALL members, never from this label.
///
/// History (2026-07-28, interviews I5/I7/I8): this function previously ranked
/// members Human SwissProt > SwissProt > Human TrEMBL > TrEMBL > other and took
/// the first. That resolved every ambiguous group towards human — the one
/// direction that inflates apparent host content — and it implemented a Methods
/// rule that has since been deleted. The ranking is gone. The replacement
/// tie-break is deliberately meaningless: because `sp|`/`tr|` are lowercase they
/// sort AFTER catalogue prefixes, so this label leans away from UniProt exactly
/// as arbitrarily as the old one leaned towards it. That is acceptable ONLY
/// because nothing is permitted to read meaning into it.
///
/// Note this is a different thing from `GroupKeyMode::Representative`, which
/// selects how the group KEY itself is formed. This function only produces a
/// label written into the `representative` column of `*_protein_group_metadata.tsv`.
///
/// Pure function of `group_key`, so shipped metadata files can be corrected by a
/// single streaming pass over column 1 — no re-aggregation is required.
fn representative_protein(group_key: &str) -> &str {
    // Lexicographically smallest member. Computed explicitly rather than assuming
    // the key is sorted, so this holds for any GroupKeyMode and any legacy file.
    group_key
        .split(';')
        .filter(|m| !m.is_empty())
        .min()
        .unwrap_or("")
}

// ============================================================================
// Parsing
// ============================================================================

/// Parse a q-value cell, treating BOTH a parse failure AND a NaN result as the
/// conservative default 1.0 (fails FDR). This is critical: `"NaN".parse::<f64>()`
/// SUCCEEDS (returns NaN), so a plain `unwrap_or(1.0)` would let a NaN slip
/// through — and NaN corrupts the BH `partial_cmp` sort (compares Equal) and the
/// `q <= threshold` gate. Valid numeric cells are returned unchanged.
fn parse_q(s: &str) -> f64 {
    let v: f64 = s.parse().unwrap_or(f64::INFINITY);
    if !v.is_finite() || !(0.0..=1.0).contains(&v) {
        f64::INFINITY
    } else {
        v
    }
}

/// Parse an intensity cell; a parse failure OR a NaN result becomes 0.0, which
/// the surrounding `intensity <= 0.0` guard treats as absent (same as a missing
/// cell). Guards against a literal "NaN" intensity. Valid numbers pass through.
fn parse_intensity(s: &str) -> f64 {
    let v: f64 = s.parse().unwrap_or(0.0);
    if !v.is_finite() || v < 0.0 {
        0.0
    } else {
        v
    }
}

fn find_sample_dirs(input: &Path) -> Vec<(String, PathBuf)> {
    let mut samples = Vec::new();
    if let Ok(entries) = fs::read_dir(input) {
        for entry in entries.flatten() {
            let path = entry.path();
            if path.is_dir() {
                let name = entry.file_name().to_string_lossy().to_string();
                samples.push((name, path));
            }
        }
    }
    samples.sort_by(|a, b| a.0.cmp(&b.0));
    samples
}

/// Load per-peptide search-level FDR (peptide_q, protein_q) from results.sage.tsv
/// in the same sample directory.
///
/// CRITICAL: the lfq.tsv `q_value` column is an LFQ-internal quantity (peak-trace
/// requant confidence) that is inflated — and for some samples pathologically so
/// (best value 0.11-0.20, zeroing the whole sample) — and must NOT be used for
/// identification FDR. The canonical per-sample FDR is the search-level q from
/// results.sage.tsv, matching the manuscript ("protein-level FDR" / "peptide-q").
/// Returns peptide -> (min peptide_q, min protein_q) over target PSMs (label==1).
fn load_search_q(sample_dir: &Path) -> std::io::Result<HashMap<String, (f64, f64)>> {
    let mut map: HashMap<String, (f64, f64)> = HashMap::new();
    let file = sample_dir.join("results.sage.tsv");
    let mut rdr = match ReaderBuilder::new()
        .delimiter(b'\t')
        .flexible(false)
        .from_path(&file)
    {
        Ok(r) => r,
        Err(e) => return Err(e.into()),
    };
    let headers = match rdr.headers() {
        Ok(h) => h.clone(),
        Err(e) => return Err(e.into()),
    };
    let pep_i = headers.iter().position(|h| h == "peptide");
    let lab_i = headers.iter().position(|h| h == "label");
    let pq_i = headers.iter().position(|h| h == "peptide_q");
    let rq_i = headers.iter().position(|h| h == "protein_q");
    let (pep_i, pq_i, rq_i) = match (pep_i, pq_i, rq_i) {
        (Some(a), Some(b), Some(c)) => (a, b, c),
        _ => {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "search table is missing q-value columns",
            ))
        }
    };
    for result in rdr.records() {
        let record = match result {
            Ok(r) => r,
            Err(e) => return Err(e.into()),
        };
        if let Some(li) = lab_i {
            if record.get(li).unwrap_or("0") != "1" {
                continue;
            }
        }
        let pep = record.get(pep_i).unwrap_or("").to_string();
        if pep.is_empty() {
            continue;
        }
        let pq = parse_q(record.get(pq_i).unwrap_or(""));
        let rq = parse_q(record.get(rq_i).unwrap_or(""));
        let e = map.entry(pep).or_insert((f64::INFINITY, f64::INFINITY));
        if pq < e.0 {
            e.0 = pq;
        }
        if rq < e.1 {
            e.1 = rq;
        }
    }
    Ok(map)
}

fn validate_search_schema(file: &Path) -> std::io::Result<()> {
    let mut reader = ReaderBuilder::new().delimiter(b'\t').from_path(file)?;
    let headers = reader.headers()?;
    for required in ["peptide", "proteins", "label", "peptide_q", "protein_q"] {
        if !headers.iter().any(|name| name == required) {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!("{} is missing required column {required}", file.display()),
            ));
        }
    }
    Ok(())
}

fn parse_lfq_sample(
    sample_name: &str,
    sample_dir: &Path,
    group_key_mode: GroupKeyMode,
) -> Option<SampleData> {
    let file = sample_dir.join("lfq.tsv");
    if !file.exists() {
        return None;
    }

    // A missing or malformed search table is an input failure, not a zero sample.
    if let Err(error) = validate_search_schema(&sample_dir.join("results.sage.tsv")) {
        eprintln!("ERROR {sample_name}: {error}");
        return None;
    }
    // Join real search-level FDR from results.sage.tsv (see load_search_q).
    let search_q = match load_search_q(sample_dir) {
        Ok(map) => map,
        Err(error) => {
            eprintln!("ERROR {sample_name}: invalid search table: {error}");
            return None;
        }
    };

    let mut rdr = ReaderBuilder::new()
        .delimiter(b'\t')
        .flexible(false)
        .from_path(&file)
        .ok()?;

    let headers = rdr.headers().ok()?.clone();

    let peptide_idx = headers.iter().position(|h| h == "peptide")?;
    let proteins_idx = headers.iter().position(|h| h == "proteins")?;
    let q_value_idx = headers.iter().position(|h| h == "q_value")?;

    // THE INTENSITY COLUMN IS NAMED AFTER THE RAW FILE, AND IT IS NOT ALWAYS mzML.
    //
    // This was `headers.iter().position(|h| h.ends_with(".mzML"))?`. SAGE names the
    // quantification column after whatever raw input it was given, so a Bruker run
    // yields `<sample>.d` and the `?` returned None -- dropping the sample with no
    // reason recorded beyond "produced no usable observations".
    //
    // Measured 2026-08-04 on ChemProt: all 102 timsTOF samples were dropped and the
    // engine reported "No valid data!", after a search that had in fact quantified
    // 112,327 LFQ features. Nothing in the message pointed at the file extension.
    //
    // So: identify the intensity column by ELIMINATION rather than by extension. Every
    // other column SAGE writes to lfq.tsv is a fixed, known name; whatever remains is
    // the per-file intensity column, whichever vendor format produced it. That cannot
    // be broken by a format this code has never heard of, which is exactly how `.d`
    // broke it.
    const LFQ_FIXED_COLUMNS: [&str; 6] = [
        "peptide",
        "charge",
        "proteins",
        "q_value",
        "score",
        "spectral_angle",
    ];
    let intensity_columns: Vec<usize> = headers
        .iter()
        .enumerate()
        .filter(|(_, h)| !LFQ_FIXED_COLUMNS.contains(h) && !h.is_empty())
        .map(|(i, _)| i)
        .collect();
    let intensity_idx = match intensity_columns.as_slice() {
        [i] => *i,
        _ => {
            // FAIL LOUD. The silent `?` here is what made 102 samples disappear
            // without a diagnosis. If the column cannot be found, say what was seen.
            eprintln!(
                "!! {sample_name}: lfq.tsv must have exactly one intensity column; found an ambiguous schema.\n\
                 !!   headers: {}\n\
                 !!   expected one column beyond {LFQ_FIXED_COLUMNS:?}, named after the raw file.",
                headers.iter().collect::<Vec<_>>().join(", ")
            );
            return None;
        }
    };

    let mut observations = Vec::new();
    let mut peptides: HashSet<String> = HashSet::new();
    let mut groups: HashSet<String> = HashSet::new();

    for result in rdr.records() {
        let record = match result {
            Ok(r) => r,
            Err(error) => {
                eprintln!("ERROR {sample_name}: invalid table row: {error}");
                return None;
            }
        };

        let peptide = record.get(peptide_idx).unwrap_or("").to_string();
        let proteins = record.get(proteins_idx).unwrap_or("");
        let q_value = parse_q(record.get(q_value_idx).unwrap_or(""));
        let intensity = parse_intensity(record.get(intensity_idx).unwrap_or("0.0"));

        if intensity <= 0.0 || peptide.is_empty() {
            continue;
        }

        let protein_group = normalize_protein_group(proteins, group_key_mode);
        if protein_group.is_empty() {
            continue;
        }

        peptides.insert(peptide.clone());
        groups.insert(protein_group.clone());

        // Identification FDR comes from the search (results.sage.tsv), NOT from
        // lfq.tsv's q_value. Missing or invalid search q-values fail every numeric
        // threshold, including 1.0. Explicit unfiltered outputs retain those features.
        let (pep_q, prot_q) = search_q
            .get(&peptide)
            .copied()
            .unwrap_or((f64::INFINITY, f64::INFINITY));

        observations.push(PeptideObs {
            peptide,
            protein_group,
            intensity,
            q_value, // lfq.tsv internal q — retained for reference, NOT used for FDR
            spectrum_q: 0.0,
            peptide_q: pep_q,
            protein_q: prot_q,
        });
    }

    if observations.is_empty() {
        return None;
    }

    Some(SampleData {
        sample_name: sample_name.to_string(),
        n_peptides: peptides.len(),
        n_groups: groups.len(),
        observations,
    })
}

fn parse_psm_sample(
    sample_name: &str,
    sample_dir: &Path,
    group_key_mode: GroupKeyMode,
) -> Option<SampleData> {
    let file = sample_dir.join("results.sage.tsv");
    if !file.exists() {
        return None;
    }

    let mut rdr = ReaderBuilder::new()
        .delimiter(b'\t')
        .flexible(false)
        .from_path(&file)
        .ok()?;

    let headers = rdr.headers().ok()?.clone();

    let peptide_idx = headers.iter().position(|h| h == "peptide")?;
    let proteins_idx = headers.iter().position(|h| h == "proteins")?;
    let label_idx = headers.iter().position(|h| h == "label")?;
    let spectrum_q_idx = headers.iter().position(|h| h == "spectrum_q")?;
    let peptide_q_idx = headers.iter().position(|h| h == "peptide_q")?;
    let protein_q_idx = headers.iter().position(|h| h == "protein_q")?;
    let intensity_idx = headers.iter().position(|h| h == "ms2_intensity")?;

    let mut observations = Vec::new();
    let mut peptides: HashSet<String> = HashSet::new();
    let mut groups: HashSet<String> = HashSet::new();

    for result in rdr.records() {
        let record = match result {
            Ok(r) => r,
            Err(error) => {
                eprintln!("ERROR {sample_name}: invalid table row: {error}");
                return None;
            }
        };

        // Filter targets only
        if record.get(label_idx).unwrap_or("0") != "1" {
            continue;
        }

        let peptide = record.get(peptide_idx).unwrap_or("").to_string();
        let proteins = record.get(proteins_idx).unwrap_or("");
        let spectrum_q = parse_q(record.get(spectrum_q_idx).unwrap_or(""));
        let peptide_q = parse_q(record.get(peptide_q_idx).unwrap_or(""));
        let protein_q = parse_q(record.get(protein_q_idx).unwrap_or(""));
        let intensity = parse_intensity(record.get(intensity_idx).unwrap_or("0.0"));

        if intensity <= 0.0 || peptide.is_empty() {
            continue;
        }

        let protein_group = normalize_protein_group(proteins, group_key_mode);
        if protein_group.is_empty() {
            continue;
        }

        peptides.insert(peptide.clone());
        groups.insert(protein_group.clone());

        observations.push(PeptideObs {
            peptide,
            protein_group,
            intensity,
            q_value: protein_q, // default q for unfiltered
            spectrum_q,
            peptide_q,
            protein_q,
        });
    }

    if observations.is_empty() {
        return None;
    }

    Some(SampleData {
        sample_name: sample_name.to_string(),
        n_peptides: peptides.len(),
        n_groups: groups.len(),
        observations,
    })
}

// ============================================================================
// Matrix building — peptide level
// ============================================================================

/// Key for peptide-level matrix: (protein_group, peptide)
type PeptideKey = (String, String);

fn build_peptide_matrix(
    sample_data: &[SampleData],
    sample_names: &[String],
    q_threshold: Option<f64>,
    q_field: &str, // "q_value", "spectrum_q", "peptide_q", "protein_q"
) -> (Vec<PeptideKey>, HashMap<PeptideKey, Vec<Option<f64>>>) {
    let n_samples = sample_names.len();
    let sample_idx: HashMap<&str, usize> = sample_names
        .iter()
        .enumerate()
        .map(|(i, s)| (s.as_str(), i))
        .collect();

    let mut matrix: IndexMap<PeptideKey, Vec<Option<f64>>> = IndexMap::new();

    for sd in sample_data {
        let sidx = match sample_idx.get(sd.sample_name.as_str()) {
            Some(&i) => i,
            None => continue,
        };

        // Group observations by (protein_group, peptide) and sum intensity
        let mut local: HashMap<PeptideKey, (f64, f64)> = HashMap::new(); // (sum_int, best_q)

        for obs in &sd.observations {
            let q = match q_field {
                "spectrum_q" => obs.spectrum_q,
                "peptide_q" => obs.peptide_q,
                "protein_q" => obs.protein_q,
                _ => obs.q_value,
            };

            // Filter observations before summation. A high-confidence observation
            // cannot admit intensity from another observation failing this cutoff.
            if q_threshold.is_some_and(|threshold| q > threshold) {
                continue;
            }
            let key = (obs.protein_group.clone(), obs.peptide.clone());
            let entry = local.entry(key).or_insert((0.0, 1.0));
            entry.0 += obs.intensity;
            if q < entry.1 {
                entry.1 = q;
            }
        }

        for (key, (sum_int, best_q)) in local {
            let passes = match q_threshold {
                Some(thresh) => best_q <= thresh,
                None => true,
            };
            if passes {
                let row = matrix.entry(key).or_insert_with(|| vec![None; n_samples]);
                row[sidx] = Some(sum_int);
            }
        }
    }

    // Sort keys for deterministic, thread-invariant row order. Insertion order
    // into `matrix` is seeded by the per-sample `local` HashMap's randomized
    // iteration order, so without this sort two runs on identical input emit the
    // same rows in a different order. (String,String) is a total order.
    let mut keys: Vec<PeptideKey> = matrix.keys().cloned().collect();
    keys.sort();
    let map: HashMap<PeptideKey, Vec<Option<f64>>> = matrix.into_iter().collect();
    (keys, map)
}

/// Protein-group-level matrix: sum peptide intensities per group per sample
fn build_protein_group_matrix(
    sample_data: &[SampleData],
    sample_names: &[String],
    q_threshold: Option<f64>,
    q_field: &str,
) -> (Vec<String>, HashMap<String, Vec<Option<f64>>>) {
    let n_samples = sample_names.len();
    let sample_idx: HashMap<&str, usize> = sample_names
        .iter()
        .enumerate()
        .map(|(i, s)| (s.as_str(), i))
        .collect();

    let mut matrix: IndexMap<String, Vec<Option<f64>>> = IndexMap::new();

    for sd in sample_data {
        let sidx = match sample_idx.get(sd.sample_name.as_str()) {
            Some(&i) => i,
            None => continue,
        };

        // Sum intensity per protein group, tracking best q
        let mut local: HashMap<&str, (f64, f64)> = HashMap::new();

        for obs in &sd.observations {
            let q = match q_field {
                "spectrum_q" => obs.spectrum_q,
                "peptide_q" => obs.peptide_q,
                "protein_q" => obs.protein_q,
                _ => obs.q_value,
            };

            if q_threshold.is_some_and(|threshold| q > threshold) {
                continue;
            }
            let entry = local
                .entry(obs.protein_group.as_str())
                .or_insert((0.0, 1.0));
            entry.0 += obs.intensity;
            if q < entry.1 {
                entry.1 = q;
            }
        }

        for (group, (sum_int, best_q)) in local {
            let passes = match q_threshold {
                Some(thresh) => best_q <= thresh,
                None => true,
            };
            if passes {
                let row = matrix
                    .entry(group.to_string())
                    .or_insert_with(|| vec![None; n_samples]);
                // Sum if multiple peptides pass (accumulate)
                match row[sidx] {
                    Some(v) => row[sidx] = Some(v + sum_int),
                    None => row[sidx] = Some(sum_int),
                }
            }
        }
    }

    // Sort keys for deterministic, thread-invariant row order. Insertion order
    // into `matrix` is seeded by the per-sample `local` HashMap's randomized
    // iteration order, so without this sort two runs on identical input emit the
    // same rows in a different order. String is a total order.
    let mut keys: Vec<String> = matrix.keys().cloned().collect();
    keys.sort();
    let map: HashMap<String, Vec<Option<f64>>> = matrix.into_iter().collect();
    (keys, map)
}

// ============================================================================
// Benjamini-Hochberg global FDR correction
// ============================================================================

/// Apply Benjamini-Hochberg FDR correction to q-values globally across
/// all samples. This is the standard approach from the Python new_aggregation
/// reference implementation:
///   adjusted_q = q_value * m / rank(q_value), capped at 1.0
///
/// Returns a mapping from (sample_index, observation_index) -> adjusted q-value.
fn apply_bh_correction(sample_data: &[SampleData], q_field: &str) -> HashMap<(usize, usize), f64> {
    // Collect all (q_value, sample_idx, obs_idx)
    let mut entries: Vec<(f64, usize, usize)> = Vec::new();
    for (si, sd) in sample_data.iter().enumerate() {
        for (oi, obs) in sd.observations.iter().enumerate() {
            let q = match q_field {
                "spectrum_q" => obs.spectrum_q,
                "peptide_q" => obs.peptide_q,
                "protein_q" => obs.protein_q,
                _ => obs.q_value,
            };
            entries.push((q, si, oi));
        }
    }

    // Sort by q-value ascending
    entries.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

    let m = entries.len() as f64;
    let mut adjusted: HashMap<(usize, usize), f64> = HashMap::with_capacity(entries.len());

    // Forward pass: compute raw adjusted values
    let mut raw_adjusted: Vec<f64> = Vec::with_capacity(entries.len());
    for (rank_minus_1, (q, _, _)) in entries.iter().enumerate() {
        let rank = (rank_minus_1 + 1) as f64;
        let adj = (q * m / rank).min(1.0);
        raw_adjusted.push(adj);
    }

    // Backward pass: ensure monotonicity (BH step-up)
    for i in (0..raw_adjusted.len().saturating_sub(1)).rev() {
        if raw_adjusted[i] > raw_adjusted[i + 1] {
            raw_adjusted[i] = raw_adjusted[i + 1];
        }
    }

    // Store results
    for (i, (_, si, oi)) in entries.iter().enumerate() {
        adjusted.insert((*si, *oi), raw_adjusted[i]);
    }

    eprintln!(
        "  BH correction ({q_field}): {m:.0} observations, median adjusted q = {:.6}",
        {
            let mut vals: Vec<f64> = adjusted.values().copied().collect();
            vals.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            if vals.is_empty() {
                0.0
            } else {
                vals[vals.len() / 2]
            }
        }
    );

    adjusted
}

// ============================================================================
// Writing
// ============================================================================

/// Shared context for the matrix writers. Bundles the parameters that are
/// identical across every write call (output location, naming components, and
/// the sample column order) so the per-call signatures carry only what varies
/// (keys, matrix, and for functional writers the level + aggregation mode).
/// Pure plumbing — no behavior change.
struct WriteCtx<'a> {
    output_dir: &'a Path,
    method: &'a str,
    dataset: &'a str,
    q_label: &'a str,
    threshold_str: &'a str,
    sample_names: &'a [String],
}

fn write_peptide_matrix(
    ctx: &WriteCtx,
    keys: &[PeptideKey],
    matrix: &HashMap<PeptideKey, Vec<Option<f64>>>,
) -> std::io::Result<usize> {
    let (output_dir, method, dataset, q_label, threshold_str, sample_names) = (
        ctx.output_dir,
        ctx.method,
        ctx.dataset,
        ctx.q_label,
        ctx.threshold_str,
        ctx.sample_names,
    );
    let n_rows = keys.len();
    let n_samples = sample_names.len();

    let filename = format!(
        "{method}_{dataset}_{q_label}_peptidelevel_{n_rows}rows_{n_samples}samples_{threshold_str}.tsv",
    );

    let path = output_dir.join(&filename);
    let file = fs::File::create(&path)?;
    let mut w = BufWriter::new(file);

    // Header: protein_groups, peptide, then sample columns
    write!(w, "protein_groups\tpeptide")?;
    for s in sample_names {
        write!(w, "\t{}", s)?;
    }
    writeln!(w)?;

    for key in keys {
        write!(w, "{}\t{}", key.0, key.1)?;
        if let Some(row) = matrix.get(key) {
            for val in row {
                match val {
                    Some(v) => write!(w, "\t{:.4}", v)?,
                    None => write!(w, "\t")?,
                }
            }
        }
        writeln!(w)?;
    }

    let file_size = fs::metadata(&path)?.len();
    eprintln!(
        "  Written: {} ({} rows x {} samples, {:.1} MB)",
        filename,
        n_rows,
        n_samples,
        file_size as f64 / 1_048_576.0
    );
    Ok(n_rows)
}

fn write_protein_group_matrix(
    ctx: &WriteCtx,
    keys: &[String],
    matrix: &HashMap<String, Vec<Option<f64>>>,
) -> std::io::Result<usize> {
    let (output_dir, method, dataset, q_label, threshold_str, sample_names) = (
        ctx.output_dir,
        ctx.method,
        ctx.dataset,
        ctx.q_label,
        ctx.threshold_str,
        ctx.sample_names,
    );
    let n_groups = keys.len();
    let n_samples = sample_names.len();

    let filename = format!(
        "{method}_{dataset}_{q_label}_proteingroup_{n_groups}grp_{n_samples}samples_{threshold_str}.tsv",
    );

    let path = output_dir.join(&filename);
    let file = fs::File::create(&path)?;
    let mut w = BufWriter::new(file);

    write!(w, "protein_groups")?;
    for s in sample_names {
        write!(w, "\t{}", s)?;
    }
    writeln!(w)?;

    for key in keys {
        write!(w, "{}", key)?;
        if let Some(row) = matrix.get(key) {
            for val in row {
                match val {
                    Some(v) => write!(w, "\t{:.4}", v)?,
                    None => write!(w, "\t")?,
                }
            }
        }
        writeln!(w)?;
    }

    let file_size = fs::metadata(&path)?.len();
    eprintln!(
        "  Written: {} ({} groups x {} samples, {:.1} MB)",
        filename,
        n_groups,
        n_samples,
        file_size as f64 / 1_048_576.0
    );
    Ok(n_groups)
}

/// Load the host-accession list for a hash-keyed build.
///
/// Only host accessions are listed; absence means microbial. See the
/// `--host-accessions` doc comment for why that is a contract rather than a guess.
fn load_host_accessions(path: &Path) -> std::io::Result<HashMap<String, String>> {
    let mut rdr = ReaderBuilder::new()
        .delimiter(b'\t')
        .flexible(true)
        .from_path(path)?;
    let headers = rdr.headers()?.clone();
    let acc_idx = headers
        .iter()
        .position(|h| h == "accession")
        .ok_or_else(|| {
            std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "host accessions: missing 'accession' column",
            )
        })?;
    let cls_idx = headers.iter().position(|h| h == "classification");
    let mut map = HashMap::new();
    for rec in rdr.records().flatten() {
        if let Some(a) = rec.get(acc_idx) {
            if a.is_empty() {
                continue;
            }
            let c = cls_idx
                .and_then(|i| rec.get(i))
                .filter(|s| !s.is_empty())
                .unwrap_or("Human");
            map.insert(a.to_string(), c.to_string());
        }
    }
    eprintln!("  Host accessions loaded: {} entries", map.len());
    Ok(map)
}

/// Classify using the host list. A member present in the list is host; a member
/// absent is microbial, which is sound only because the list is complete.
fn classify_group_with_hosts(group_key: &str, hosts: &HashMap<String, String>) -> &'static str {
    let mut has_host = false;
    let mut has_micro = false;
    for member in group_key.split(';') {
        match hosts.get(member).map(String::as_str) {
            Some("Mixed") => {
                // A single sequence produced by both a host and a non-host source.
                return "Mixed";
            }
            Some(_) => has_host = true,
            None => has_micro = true,
        }
        if has_host && has_micro {
            return "Mixed";
        }
    }
    if has_host {
        "Human"
    } else {
        "Microbial"
    }
}

fn write_group_metadata(
    output_dir: &Path,
    dataset: &str,
    all_groups: &HashSet<String>,
    hosts: Option<&HashMap<String, String>>,
) -> std::io::Result<()> {
    let path = output_dir.join(format!("{}_protein_group_metadata.tsv", dataset));
    let file = fs::File::create(&path)?;
    let mut w = BufWriter::new(file);

    writeln!(
        w,
        "protein_group\trepresentative\tn_members\tclassification"
    )?;

    let mut sorted: Vec<&String> = all_groups.iter().collect();
    sorted.sort();

    let mut n_human = 0u64;
    let mut n_micro = 0u64;
    let mut n_mixed = 0u64;
    let mut n_unknown = 0u64;

    for g in &sorted {
        let rep = representative_protein(g);
        let n_members = g.split(';').count();
        let class = match hosts {
            Some(h) => classify_group_with_hosts(g, h),
            None => classify_group(g),
        };
        match class {
            "Human" => n_human += 1,
            "Microbial" => n_micro += 1,
            "Mixed" => n_mixed += 1,
            "Unknown" => n_unknown += 1,
            _ => {}
        }
        writeln!(w, "{}\t{}\t{}\t{}", g, rep, n_members, class)?;
    }

    eprintln!(
        "  Group metadata: {} groups (Human={}, Microbial={}, Mixed={}, Unknown={})",
        sorted.len(),
        n_human,
        n_micro,
        n_mixed,
        n_unknown
    );
    if n_unknown > 0 {
        let pct = 100.0 * n_unknown as f64 / sorted.len().max(1) as f64;
        eprintln!(
            "  WARNING: {n_unknown} groups ({pct:.1}%) are NOT CLASSIFIABLE from the \
             matrix. Their members are uniform hash accessions carrying no source, so \
             host vs microbial cannot be determined here. Pass --host-accessions with \
             the list derived from the stage-0 joint_headers sidecar."
        );
    }
    // A host list that matches nothing is a broken join, not a microbe-only cohort.
    // Say so, because a silent zero here is exactly the failure this flag exists to
    // prevent.
    if hosts.is_some() && n_human == 0 && n_mixed == 0 {
        eprintln!(
            "  WARNING: --host-accessions was supplied but matched NO group member. \
             Either the accession namespaces disagree or this cohort's lake genuinely \
             contains no host sequences (SIHUMI by design). Check before trusting a \
             host fraction of zero."
        );
    }
    Ok(())
}

// ============================================================================
// Functional lookup + aggregation
// ============================================================================

/// accession -> { column_name -> Vec<value> (decomposed from semicolon source) }
/// Only levels requested at CLI parse time are stored, to keep memory bounded.
type FunctionalTable = HashMap<String, HashMap<String, Vec<String>>>;

/// Load functional lookup TSV. Keeps only requested columns.
///
/// The TSV is expected to have `accession` as one column (exact name) and
/// any of the requested functional-level columns. Values may contain
/// semicolon- OR comma-separated sub-values which we split further.
/// Sentinel values like `-`, `--`, `` (empty), `NA`, `None` become empty.
/// Missing columns are allowed: the loader logs a warning and those levels
/// will produce trivial (all-Unannotated) matrices downstream.
fn load_functional_lookup(path: &Path, levels: &[String]) -> std::io::Result<FunctionalTable> {
    let t0 = Instant::now();
    let mut rdr = ReaderBuilder::new()
        .delimiter(b'\t')
        .flexible(true)
        .from_path(path)?;

    let headers = rdr.headers()?.clone();
    let acc_idx = headers
        .iter()
        .position(|h| h == "accession")
        .ok_or_else(|| {
            std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "functional lookup: missing 'accession' column",
            )
        })?;

    // Map requested level -> column index (None if missing; we warn once)
    let mut level_idx: Vec<(String, Option<usize>)> = Vec::new();
    for lvl in levels {
        let idx = headers.iter().position(|h| h == lvl.as_str());
        if idx.is_none() {
            eprintln!("  WARNING: functional lookup missing column '{lvl}' — will produce empty matrix for that level");
        }
        level_idx.push((lvl.clone(), idx));
    }

    let mut table: FunctionalTable = HashMap::with_capacity(1_000_000);
    let mut n_rows = 0usize;

    for result in rdr.records() {
        let record = match result {
            Ok(r) => r,
            Err(_) => continue,
        };
        let acc = match record.get(acc_idx) {
            Some(s) if !s.is_empty() => s.to_string(),
            _ => continue,
        };

        let mut per_level: HashMap<String, Vec<String>> = HashMap::new();
        for (lvl, idx_opt) in &level_idx {
            if let Some(idx) = idx_opt {
                let raw = record.get(*idx).unwrap_or("");
                let vals = split_functional_value(raw);
                if !vals.is_empty() {
                    per_level.insert(lvl.clone(), vals);
                }
            }
        }
        table.insert(acc, per_level);
        n_rows += 1;

        if n_rows.is_multiple_of(500_000) {
            eprintln!("  loaded {} rows...", n_rows);
        }
    }

    eprintln!(
        "  Functional lookup loaded: {} rows, {} levels, {:.1}s",
        n_rows,
        level_idx.len(),
        t0.elapsed().as_secs_f64()
    );
    Ok(table)
}

/// Split a raw functional-annotation cell into discrete values.
/// Handles semicolons AND commas (eggNOG uses commas for KEGG_Pathway),
/// drops sentinel null values.
fn split_functional_value(raw: &str) -> Vec<String> {
    let trimmed = raw.trim();
    if trimmed.is_empty() {
        return Vec::new();
    }
    // Common sentinels for missing annotations in eggNOG / KEGG exports
    match trimmed {
        "-" | "--" | "NA" | "None" | "null" | "NULL" => return Vec::new(),
        _ => {}
    }
    // Split by ; or , — eggNOG uses both depending on field
    let mut out = Vec::new();
    for part in trimmed.split([';', ',']) {
        let p = part.trim();
        if p.is_empty() {
            continue;
        }
        match p {
            "-" | "--" | "NA" | "None" | "null" | "NULL" => continue,
            _ => out.push(p.to_string()),
        }
    }
    out
}

/// For a given protein_group (semicolon-joined members) at a given level,
/// collect the per-value votes.
///
/// Returns a sorted Vec<(value, count)> where count = number of member proteins
/// in the group that carry that value. A member contributes one vote per
/// DISTINCT value it carries (prevents one semicolon-heavy member from dominating).
///
/// An empty return vector means NO member has an annotation for this level
/// (i.e. the group is unannotated at this level).
fn group_value_votes(
    group_key: &str,
    level: &str,
    table: &FunctionalTable,
) -> Vec<(String, usize)> {
    let mut counts: BTreeMap<String, usize> = BTreeMap::new();

    for member in group_key.split(';') {
        if let Some(per_level) = table.get(member) {
            if let Some(vals) = per_level.get(level) {
                if !vals.is_empty() {
                    let mut seen_here: BTreeSet<&str> = BTreeSet::new();
                    for v in vals {
                        seen_here.insert(v.as_str());
                    }
                    for v in seen_here {
                        *counts.entry(v.to_string()).or_insert(0) += 1;
                    }
                }
            }
        }
    }

    counts.into_iter().collect()
}

/// Canonical mapping for one protein group at one functional level.
/// Resolved once per aggregate so that rows stay consistent across FDR thresholds.
///
/// - `winner` is the plurality-vote winner (alphabetical tie-break), or None if
///   the group has no annotation at this level (we assign an Unannotated_N id
///   in that case at a higher level).
/// - `values` is the full distinct list of functional values contributed by
///   members — used by SplitIntensity mode (len = N_unique for splitting).
struct GroupCanonical {
    winner: Option<String>,
    values: Vec<String>,
}

fn canonical_for_group(group_key: &str, level: &str, table: &FunctionalTable) -> GroupCanonical {
    let votes = group_value_votes(group_key, level, table);
    if votes.is_empty() {
        return GroupCanonical {
            winner: None,
            values: Vec::new(),
        };
    }

    // BTreeMap iteration was alphabetical; pick highest count, tie-break alphabetical first
    let mut best: Option<(&str, usize)> = None;
    for (val, cnt) in &votes {
        match best {
            None => best = Some((val.as_str(), *cnt)),
            Some((_, bc)) if *cnt > bc => best = Some((val.as_str(), *cnt)),
            _ => {}
        }
    }
    let winner = best.map(|(v, _)| v.to_string());
    let values: Vec<String> = votes.into_iter().map(|(v, _)| v).collect();
    GroupCanonical { winner, values }
}

/// Per-group resolved metadata for a single functional level.
/// Built ONCE per aggregate; reused for every FDR threshold and for both
/// aggregation modes. The `row_id` is the stable row identifier in output
/// matrices — either the canonical functional value (annotated groups) or
/// `Unannotated_{NNNNNN}` (unannotated groups, stable numbering by
/// alphabetical group_key ordering).
struct LevelResolution {
    /// group_key -> per-level resolution
    per_group: HashMap<String, GroupResolution>,
}

struct GroupResolution {
    /// For plurality mode: the row id this group's full intensity maps to.
    /// For annotated groups = winner value; for unannotated = "Unannotated_{N}".
    plurality_row: String,
    /// For split mode: the set of (row_id, weight) this group contributes to.
    /// - annotated: one entry per distinct value, weight = 1/N
    /// - unannotated: single entry `("Unannotated_{N}", 1.0)`
    split_rows: Vec<(String, f64)>,
}

/// Resolve every group for one functional level in a single pass.
/// Assigns stable `Unannotated_{NNNNNN}` ids to unannotated groups in the
/// deterministic order of `group_keys` (caller typically passes sorted keys).
fn resolve_level(level: &str, table: &FunctionalTable, group_keys: &[String]) -> LevelResolution {
    let mut per_group: HashMap<String, GroupResolution> = HashMap::with_capacity(group_keys.len());
    let mut unannot_counter: u64 = 0;

    for gk in group_keys {
        let canon = canonical_for_group(gk, level, table);
        let resolution = match canon.winner {
            Some(winner) => {
                let n = canon.values.len().max(1) as f64;
                let weight = 1.0 / n;
                let split_rows: Vec<(String, f64)> =
                    canon.values.into_iter().map(|v| (v, weight)).collect();
                GroupResolution {
                    plurality_row: winner,
                    split_rows,
                }
            }
            None => {
                unannot_counter += 1;
                let id = format!("Unannotated_{:06}", unannot_counter);
                GroupResolution {
                    plurality_row: id.clone(),
                    split_rows: vec![(id, 1.0)],
                }
            }
        };
        per_group.insert(gk.clone(), resolution);
    }

    LevelResolution { per_group }
}

/// Build a protein-group-level functional matrix given a pre-computed
/// per-group resolution and a q-filtered protein-group intensity matrix.
///
/// Returns (sorted keys, matrix).
fn build_functional_group_matrix(
    resolution: &LevelResolution,
    group_keys: &[String],
    group_matrix: &HashMap<String, Vec<Option<f64>>>,
    n_samples: usize,
    mode: AggregationMode,
) -> (Vec<String>, HashMap<String, Vec<Option<f64>>>) {
    let mut func: IndexMap<String, Vec<Option<f64>>> = IndexMap::new();

    for gk in group_keys {
        let res = match resolution.per_group.get(gk) {
            Some(r) => r,
            None => continue, // shouldn't happen if resolution was built from all groups
        };
        let g_row = match group_matrix.get(gk) {
            Some(r) => r,
            None => continue, // group not in this q-filtered matrix
        };

        match mode {
            AggregationMode::Plurality => {
                let row = func
                    .entry(res.plurality_row.clone())
                    .or_insert_with(|| vec![None; n_samples]);
                for (i, val) in g_row.iter().enumerate() {
                    if let Some(v) = val {
                        match row[i] {
                            Some(existing) => row[i] = Some(existing + v),
                            None => row[i] = Some(*v),
                        }
                    }
                }
            }
            AggregationMode::SplitIntensity => {
                for (rid, weight) in &res.split_rows {
                    let row = func
                        .entry(rid.clone())
                        .or_insert_with(|| vec![None; n_samples]);
                    for (i, val) in g_row.iter().enumerate() {
                        if let Some(v) = val {
                            let contrib = v * weight;
                            match row[i] {
                                Some(existing) => row[i] = Some(existing + contrib),
                                None => row[i] = Some(contrib),
                            }
                        }
                    }
                }
            }
        }
    }

    let mut keys: Vec<String> = func.keys().cloned().collect();
    keys.sort();
    let map: HashMap<String, Vec<Option<f64>>> = func.into_iter().collect();
    (keys, map)
}

/// Dual-index row key for the peptide-level functional matrix:
/// (functional_value, peptide).
type FuncPeptideKey = (String, String);

/// Peptide-level functional matrix: dual-index (functional_value, peptide).
/// Rows propagate the SAME per-group resolution used for the protein-group
/// functional matrix, so the two views are consistent.
///
/// For plurality mode: one row per (winner_value, peptide).
/// For split mode: one row per (value_i, peptide), intensities scaled by 1/N.
fn build_functional_peptide_matrix(
    resolution: &LevelResolution,
    peptide_keys: &[PeptideKey],
    peptide_matrix: &HashMap<PeptideKey, Vec<Option<f64>>>,
    n_samples: usize,
    mode: AggregationMode,
) -> (
    Vec<FuncPeptideKey>,
    HashMap<FuncPeptideKey, Vec<Option<f64>>>,
) {
    let mut func: IndexMap<(String, String), Vec<Option<f64>>> = IndexMap::new();

    for pk in peptide_keys {
        let (gk, pep) = pk;
        let res = match resolution.per_group.get(gk) {
            Some(r) => r,
            None => continue,
        };
        let p_row = match peptide_matrix.get(pk) {
            Some(r) => r,
            None => continue,
        };

        match mode {
            AggregationMode::Plurality => {
                let key = (res.plurality_row.clone(), pep.clone());
                let row = func.entry(key).or_insert_with(|| vec![None; n_samples]);
                for (i, val) in p_row.iter().enumerate() {
                    if let Some(v) = val {
                        match row[i] {
                            Some(existing) => row[i] = Some(existing + v),
                            None => row[i] = Some(*v),
                        }
                    }
                }
            }
            AggregationMode::SplitIntensity => {
                for (rid, weight) in &res.split_rows {
                    let key = (rid.clone(), pep.clone());
                    let row = func.entry(key).or_insert_with(|| vec![None; n_samples]);
                    for (i, val) in p_row.iter().enumerate() {
                        if let Some(v) = val {
                            let contrib = v * weight;
                            match row[i] {
                                Some(existing) => row[i] = Some(existing + contrib),
                                None => row[i] = Some(contrib),
                            }
                        }
                    }
                }
            }
        }
    }

    let mut keys: Vec<(String, String)> = func.keys().cloned().collect();
    keys.sort();
    let map: HashMap<(String, String), Vec<Option<f64>>> = func.into_iter().collect();
    (keys, map)
}

// ============================================================================
// Functional writers
// ============================================================================

fn sanitize_level(level: &str) -> String {
    level
        .chars()
        .map(|c| {
            if c.is_alphanumeric() || c == '_' {
                c
            } else {
                '_'
            }
        })
        .collect()
}

fn write_functional_group_matrix(
    ctx: &WriteCtx,
    level: &str,
    mode: AggregationMode,
    keys: &[String],
    matrix: &HashMap<String, Vec<Option<f64>>>,
) -> std::io::Result<usize> {
    let (output_dir, method, dataset, q_label, threshold_str, sample_names) = (
        ctx.output_dir,
        ctx.method,
        ctx.dataset,
        ctx.q_label,
        ctx.threshold_str,
        ctx.sample_names,
    );
    let n_rows = keys.len();
    let n_samples = sample_names.len();
    let level_safe = sanitize_level(level);
    let mode_s = mode.short();

    let filename = format!(
        "{method}_{dataset}_{q_label}_by_{level_safe}_{mode_s}_proteingroup_{n_rows}rows_{n_samples}samples_{threshold_str}.tsv"
    );
    let path = output_dir.join(&filename);
    let file = fs::File::create(&path)?;
    let mut w = BufWriter::new(file);

    write!(w, "{}", level_safe)?;
    for s in sample_names {
        write!(w, "\t{}", s)?;
    }
    writeln!(w)?;

    for key in keys {
        write!(w, "{}", key)?;
        if let Some(row) = matrix.get(key) {
            for val in row {
                match val {
                    Some(v) => write!(w, "\t{:.4}", v)?,
                    None => write!(w, "\t")?,
                }
            }
        }
        writeln!(w)?;
    }

    let file_size = fs::metadata(&path)?.len();
    eprintln!(
        "  Written: {} ({} rows x {} samples, {:.1} MB)",
        filename,
        n_rows,
        n_samples,
        file_size as f64 / 1_048_576.0
    );
    Ok(n_rows)
}

fn write_functional_peptide_matrix(
    ctx: &WriteCtx,
    level: &str,
    mode: AggregationMode,
    keys: &[(String, String)],
    matrix: &HashMap<(String, String), Vec<Option<f64>>>,
) -> std::io::Result<usize> {
    let (output_dir, method, dataset, q_label, threshold_str, sample_names) = (
        ctx.output_dir,
        ctx.method,
        ctx.dataset,
        ctx.q_label,
        ctx.threshold_str,
        ctx.sample_names,
    );
    let n_rows = keys.len();
    let n_samples = sample_names.len();
    let level_safe = sanitize_level(level);
    let mode_s = mode.short();

    let filename = format!(
        "{method}_{dataset}_{q_label}_by_{level_safe}_{mode_s}_peptidelevel_{n_rows}rows_{n_samples}samples_{threshold_str}.tsv"
    );
    let path = output_dir.join(&filename);
    let file = fs::File::create(&path)?;
    let mut w = BufWriter::new(file);

    write!(w, "{}\tpeptide", level_safe)?;
    for s in sample_names {
        write!(w, "\t{}", s)?;
    }
    writeln!(w)?;

    for key in keys {
        write!(w, "{}\t{}", key.0, key.1)?;
        if let Some(row) = matrix.get(key) {
            for val in row {
                match val {
                    Some(v) => write!(w, "\t{:.4}", v)?,
                    None => write!(w, "\t")?,
                }
            }
        }
        writeln!(w)?;
    }

    let file_size = fs::metadata(&path)?.len();
    eprintln!(
        "  Written: {} ({} rows x {} samples, {:.1} MB)",
        filename,
        n_rows,
        n_samples,
        file_size as f64 / 1_048_576.0
    );
    Ok(n_rows)
}

/// Emit one group -> level resolution map per (dataset, level). This table
/// is identical across FDR thresholds (resolution is precomputed) and is
/// useful for downstream QC and joining.
fn write_group_to_level_map(
    output_dir: &Path,
    dataset: &str,
    level: &str,
    resolution: &LevelResolution,
) -> std::io::Result<()> {
    let level_safe = sanitize_level(level);
    let path = output_dir.join(format!("{dataset}_group_to_{level_safe}_map.tsv"));
    let file = fs::File::create(&path)?;
    let mut w = BufWriter::new(file);

    writeln!(
        w,
        "protein_group\tplurality_value\tn_unique_values\tsplit_values"
    )?;

    let mut entries: Vec<(&String, &GroupResolution)> = resolution.per_group.iter().collect();
    entries.sort_by(|a, b| a.0.cmp(b.0));

    for (gk, res) in entries {
        let split_join: Vec<String> = res
            .split_rows
            .iter()
            .map(|(rid, w)| format!("{}:{:.6}", rid, w))
            .collect();
        writeln!(
            w,
            "{}\t{}\t{}\t{}",
            gk,
            res.plurality_row,
            res.split_rows.len(),
            split_join.join(";")
        )?;
    }
    Ok(())
}

// ============================================================================
// Main
// ============================================================================

fn resolve_functional_levels(cli: &Cli) -> Vec<String> {
    if cli.no_functional {
        return Vec::new();
    }
    let spec = cli
        .functional_levels
        .clone()
        .unwrap_or_else(|| DEFAULT_FUNCTIONAL_LEVELS.to_string());
    spec.split(',')
        .map(|x| x.trim().to_string())
        .filter(|x| !x.is_empty())
        .collect()
}

fn main() -> std::io::Result<()> {
    let cli = Cli::parse();
    let start = Instant::now();
    if cli.bh_correction {
        return Err(std::io::Error::new(std::io::ErrorKind::InvalidInput,
            "BH requires valid p-values; search q-values cannot be used to establish global FDR. Omit --bh-correction."));
    }
    for value in [&cli.dataset, &cli.method] {
        if value.is_empty()
            || value == "."
            || value == ".."
            || !value
                .chars()
                .all(|c| c.is_ascii_alphanumeric() || "._-".contains(c))
        {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                "dataset and method must contain only letters, digits, '.', '_' or '-'",
            ));
        }
    }
    if cli.threads == 0 {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            "threads must be positive",
        ));
    }
    if cli.output.exists() && fs::read_dir(&cli.output)?.next().is_some() {
        return Err(std::io::Error::new(
            std::io::ErrorKind::AlreadyExists,
            "Output directory is nonempty; choose a fresh directory. Existing output is preserved.",
        ));
    }
    // Validate threshold syntax before creating output files.
    parse_q_thresholds(cli.q_thresholds.as_deref())
        .map_err(|error| std::io::Error::new(std::io::ErrorKind::InvalidInput, error))?;

    rayon::ThreadPoolBuilder::new()
        .num_threads(cli.threads)
        .build_global()
        .ok();

    fs::create_dir_all(&cli.output).expect("Failed to create output directory");

    eprintln!("================================================================================");
    eprintln!("AGGREGATION ENGINE v2 — FUNCTIONAL PLURALITY + SPLITINTENSITY, DUAL INDEX");
    eprintln!("Dataset: {} — {:?}", cli.dataset, cli.data_type);
    eprintln!("================================================================================");
    eprintln!("Input:   {}", cli.input.display());
    eprintln!("Output:  {}", cli.output.display());
    eprintln!("Threads: {}", cli.threads);

    let sample_dirs = find_sample_dirs(&cli.input);
    eprintln!("Sample directories found: {}", sample_dirs.len());

    if sample_dirs.is_empty() {
        eprintln!("ERROR: No sample directories found!");
        std::process::exit(1);
    }

    // Parse all samples in parallel
    eprintln!("Parsing samples...");
    eprintln!("Group-key mode: {:?}", cli.group_key_mode);
    let data_type_clone = cli.data_type.clone();
    let group_key_mode = cli.group_key_mode;
    let mut sample_data: Vec<SampleData> = sample_dirs
        .par_iter()
        .filter_map(|(name, dir)| match data_type_clone {
            DataType::Lfq => parse_lfq_sample(name, dir, group_key_mode),
            DataType::Psm => parse_psm_sample(name, dir, group_key_mode),
        })
        .collect();

    eprintln!(
        "Parsed: {}/{} samples",
        sample_data.len(),
        sample_dirs.len()
    );

    // LOUD DROP REPORT (added 2026-07-27). Previously a sample that parsed to None was
    // discarded by the filter_map above and the ONLY trace was this count. SIHUMIx F5
    // vanished that way and nobody noticed until the numbers were audited months later;
    // the "164samples"/"282samples"/"3samples" output labels are the same thing.
    // A dropped sample is now named, so a silently shrinking cohort is visible in the log.
    if sample_data.len() < sample_dirs.len() {
        let kept: std::collections::HashSet<&str> =
            sample_data.iter().map(|s| s.sample_name.as_str()).collect();
        let dropped: Vec<&str> = sample_dirs
            .iter()
            .map(|(n, _)| n.as_str())
            .filter(|n| !kept.contains(n))
            .collect();
        eprintln!(
            "!! WARNING: {} sample(s) produced no usable observations and were DROPPED.",
            dropped.len()
        );
        eprintln!("!! A dropped sample contributes to NO reported value. Verify each is a genuine");
        eprintln!("!! failure and not a pipeline error before using these outputs.");
        for n in &dropped {
            eprintln!("!!   dropped: {}", n);
        }
        return Err(std::io::Error::new(std::io::ErrorKind::InvalidData,
            "Incomplete cohort: inspect the named samples and correct their inputs before aggregation."));
    }

    if sample_data.is_empty() {
        eprintln!("ERROR: No valid data!");
        std::process::exit(1);
    }

    // Stats
    let total_obs: usize = sample_data.iter().map(|s| s.observations.len()).sum();
    let total_peptides: usize = sample_data.iter().map(|s| s.n_peptides).sum();
    let total_groups: usize = sample_data.iter().map(|s| s.n_groups).sum();
    eprintln!("Total observations: {}", total_obs);
    eprintln!("Total peptides (sum): {}", total_peptides);
    eprintln!("Total groups (sum): {}", total_groups);

    // Collect all unique groups — sorted so Unannotated_N ids are stable
    let mut all_groups_set: HashSet<String> = HashSet::new();
    for sd in &sample_data {
        for obs in &sd.observations {
            all_groups_set.insert(obs.protein_group.clone());
        }
    }
    let mut all_groups_sorted: Vec<String> = all_groups_set.iter().cloned().collect();
    all_groups_sorted.sort();
    eprintln!("Unique protein groups: {}", all_groups_sorted.len());

    let sample_names: Vec<String> = sample_data.iter().map(|s| s.sample_name.clone()).collect();
    let n_samples = sample_names.len();

    // Apply global Benjamini-Hochberg FDR correction if requested
    if cli.bh_correction {
        eprintln!("\nApplying global Benjamini-Hochberg correction...");
        match cli.data_type {
            DataType::Lfq => {
                // C5 fix: the LFQ matrices gate on the SEARCH-derived peptide_q and
                // protein_q columns (see the build_peptide_matrix / build_protein_group_matrix
                // calls in the DataType::Lfq write path below, which pass "peptide_q" /
                // "protein_q"), NOT the lfq.tsv-internal `q_value`. Correcting only
                // `q_value` here (the previous behavior) left the gated columns untouched,
                // so --bh-correction was a silent no-op for LFQ. Apply the global BH
                // correction to the SAME columns the matrices gate on, mirroring the PSM
                // path, so the flag actually takes effect. (The BH math in
                // apply_bh_correction is unchanged.)
                for q_field in &["peptide_q", "protein_q"] {
                    let adjusted = apply_bh_correction(&sample_data, q_field);
                    for (si, sd) in sample_data.iter_mut().enumerate() {
                        for (oi, obs) in sd.observations.iter_mut().enumerate() {
                            if let Some(&adj) = adjusted.get(&(si, oi)) {
                                match *q_field {
                                    "peptide_q" => obs.peptide_q = adj,
                                    "protein_q" => obs.protein_q = adj,
                                    _ => {}
                                }
                            }
                        }
                    }
                }
            }
            DataType::Psm => {
                for q_field in &["spectrum_q", "peptide_q", "protein_q"] {
                    let adjusted = apply_bh_correction(&sample_data, q_field);
                    for (si, sd) in sample_data.iter_mut().enumerate() {
                        for (oi, obs) in sd.observations.iter_mut().enumerate() {
                            if let Some(&adj) = adjusted.get(&(si, oi)) {
                                match *q_field {
                                    "spectrum_q" => obs.spectrum_q = adj,
                                    "peptide_q" => obs.peptide_q = adj,
                                    "protein_q" => obs.protein_q = adj,
                                    _ => {}
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // Optionally load functional lookup
    let functional_levels: Vec<String> = resolve_functional_levels(&cli);

    let functional_table: Option<FunctionalTable> = if !functional_levels.is_empty() {
        match &cli.functional_lookup {
            Some(p) => {
                eprintln!("\nLoading functional lookup: {}", p.display());
                match load_functional_lookup(p, &functional_levels) {
                    Ok(t) => Some(t),
                    Err(e) => {
                        eprintln!("ERROR loading functional lookup: {e}");
                        std::process::exit(2);
                    }
                }
            }
            None => {
                eprintln!("INFO: no --functional-lookup provided; skipping functional matrices");
                None
            }
        }
    } else {
        None
    };

    // Resolve each functional level ONCE (stable Unannotated_N ids across FDR thresholds)
    let level_resolutions: HashMap<String, LevelResolution> =
        if let Some(ref ftable) = functional_table {
            eprintln!(
                "\nResolving functional levels for {} groups...",
                all_groups_sorted.len()
            );
            let mut map = HashMap::new();
            for lvl in &functional_levels {
                let t0 = Instant::now();
                let res = resolve_level(lvl, ftable, &all_groups_sorted);
                let n_unannot = res
                    .per_group
                    .values()
                    .filter(|r| r.plurality_row.starts_with("Unannotated_"))
                    .count();
                eprintln!(
                    "  {}: {} groups resolved, {} unannotated ({:.1}%, {:.1}s)",
                    lvl,
                    res.per_group.len(),
                    n_unannot,
                    100.0 * n_unannot as f64 / res.per_group.len().max(1) as f64,
                    t0.elapsed().as_secs_f64()
                );
                map.insert(lvl.clone(), res);
            }
            map
        } else {
            HashMap::new()
        };

    // Write group -> level maps (one per level, threshold-independent)
    if !level_resolutions.is_empty() {
        eprintln!("\nWriting group -> level maps...");
        for (lvl, res) in &level_resolutions {
            write_group_to_level_map(&cli.output, &cli.dataset, lvl, res)?;
        }
    }

    // Build and write matrices
    eprintln!("\nCreating matrices...");

    let q_thresholds: Vec<(Option<f64>, String)> =
        match parse_q_thresholds(cli.q_thresholds.as_deref()) {
            Ok(v) => v,
            Err(e) => {
                eprintln!("ERROR parsing --q-thresholds: {e}");
                std::process::exit(2);
            }
        };
    eprintln!(
        "Q-value thresholds: {}",
        q_thresholds
            .iter()
            .map(|(_, l)| l.as_str())
            .collect::<Vec<_>>()
            .join(", ")
    );

    let data_type_str = match cli.data_type {
        DataType::Lfq => "lfq",
        DataType::Psm => "psm",
    };

    let modes = [AggregationMode::Plurality, AggregationMode::SplitIntensity];

    match cli.data_type {
        DataType::Lfq => {
            for (thresh, label) in &q_thresholds {
                let ctx = WriteCtx {
                    output_dir: &cli.output,
                    method: &cli.method,
                    dataset: &cli.dataset,
                    q_label: "lfq",
                    threshold_str: label,
                    sample_names: &sample_names,
                };
                // Peptide-level — gate on search-level peptide_q (NOT lfq.tsv q_value)
                let (pkeys, pmatrix) =
                    build_peptide_matrix(&sample_data, &sample_names, *thresh, "peptide_q");
                write_peptide_matrix(&ctx, &pkeys, &pmatrix)?;

                // Protein-group-level — gate on search-level protein_q (NOT lfq.tsv q_value)
                let (gkeys, gmatrix) =
                    build_protein_group_matrix(&sample_data, &sample_names, *thresh, "protein_q");
                write_protein_group_matrix(&ctx, &gkeys, &gmatrix)?;

                // Functional roll-ups (per level × mode × granularity), opt-in.
                // The group_to_<level> maps are written unconditionally above, so a
                // consumer can always derive any of these; see --functional-rollups.
                if cli.functional_rollups && !level_resolutions.is_empty() {
                    for lvl in &functional_levels {
                        let resolution = match level_resolutions.get(lvl) {
                            Some(r) => r,
                            None => continue,
                        };
                        for mode in &modes {
                            // Protein-group level
                            let (fkeys, fmatrix) = build_functional_group_matrix(
                                resolution, &gkeys, &gmatrix, n_samples, *mode,
                            );
                            write_functional_group_matrix(&ctx, lvl, *mode, &fkeys, &fmatrix)?;

                            // Peptide level
                            let (pfkeys, pfmatrix) = build_functional_peptide_matrix(
                                resolution, &pkeys, &pmatrix, n_samples, *mode,
                            );
                            write_functional_peptide_matrix(&ctx, lvl, *mode, &pfkeys, &pfmatrix)?;
                        }
                    }
                }
            }
        }
        DataType::Psm => {
            let q_types = ["spectrum_q", "peptide_q", "protein_q"];
            for q_type in &q_types {
                let q_label = q_type.replace("_q", "");
                for (thresh, label) in &q_thresholds {
                    let ctx = WriteCtx {
                        output_dir: &cli.output,
                        method: &cli.method,
                        dataset: &cli.dataset,
                        q_label: &q_label,
                        threshold_str: label,
                        sample_names: &sample_names,
                    };
                    let (pkeys, pmatrix) =
                        build_peptide_matrix(&sample_data, &sample_names, *thresh, q_type);
                    write_peptide_matrix(&ctx, &pkeys, &pmatrix)?;

                    let (gkeys, gmatrix) =
                        build_protein_group_matrix(&sample_data, &sample_names, *thresh, q_type);
                    write_protein_group_matrix(&ctx, &gkeys, &gmatrix)?;

                    // Functional matrices — only for protein_q to keep output size reasonable
                    if cli.functional_rollups
                        && q_type == &"protein_q"
                        && !level_resolutions.is_empty()
                    {
                        for lvl in &functional_levels {
                            let resolution = match level_resolutions.get(lvl) {
                                Some(r) => r,
                                None => continue,
                            };
                            for mode in &modes {
                                let (fkeys, fmatrix) = build_functional_group_matrix(
                                    resolution, &gkeys, &gmatrix, n_samples, *mode,
                                );
                                write_functional_group_matrix(&ctx, lvl, *mode, &fkeys, &fmatrix)?;

                                let (pfkeys, pfmatrix) = build_functional_peptide_matrix(
                                    resolution, &pkeys, &pmatrix, n_samples, *mode,
                                );
                                write_functional_peptide_matrix(
                                    &ctx, lvl, *mode, &pfkeys, &pfmatrix,
                                )?;
                            }
                        }
                    }
                }
            }
        }
    }

    // Write metadata
    eprintln!("\nWriting metadata...");
    let host_map = match &cli.host_accessions {
        Some(p) => {
            eprintln!("Loading host accessions from {}", p.display());
            Some(load_host_accessions(p)?)
        }
        None => None,
    };
    write_group_metadata(
        &cli.output,
        &cli.dataset,
        &all_groups_set,
        host_map.as_ref(),
    )?;

    // Summary
    let summary_path = cli.output.join(format!(
        "{}_{}_{}_{n_samples}samples_summary.tsv",
        cli.method, cli.dataset, data_type_str
    ));
    let mut summary = fs::File::create(&summary_path)?;
    writeln!(summary, "sample\tn_groups\tn_peptides\tn_observations")?;
    for sd in &sample_data {
        writeln!(
            summary,
            "{}\t{}\t{}\t{}",
            sd.sample_name,
            sd.n_groups,
            sd.n_peptides,
            sd.observations.len()
        )?;
    }

    let elapsed = start.elapsed();
    eprintln!("\n================================================================================");
    eprintln!("AGGREGATION v2 COMPLETE");
    eprintln!("================================================================================");
    eprintln!("Samples: {}", n_samples);
    eprintln!("Unique protein groups: {}", all_groups_sorted.len());
    if cli.functional_rollups {
        eprintln!(
            "Functional levels: {} (roll-ups WRITTEN, modes: plurality + splitintensity)",
            functional_levels.len()
        );
    } else {
        eprintln!(
            "Functional levels: {} (group_to_<level> maps written; roll-up matrices NOT \
             written -- derive them from a core matrix + the map, or pass \
             --functional-rollups)",
            functional_levels.len()
        );
    }
    eprintln!("Runtime: {:.1}s", elapsed.as_secs_f64());

    Ok(())
}

// ============================================================================
// TESTS
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    /// Regression guard for the lfq-q FDR bug: parse_lfq_sample must take its
    /// identification FDR (peptide_q / protein_q) from results.sage.tsv, NOT from
    /// lfq.tsv's inflated `q_value`. If anyone reverts the gate, this fails.
    #[test]
    fn test_lfq_fdr_from_search_not_lfq_q() {
        use std::fs;
        let dir = std::env::temp_dir().join(format!("agg_lfqfdr_{}", std::process::id()));
        let _ = fs::create_dir_all(&dir);
        // lfq.tsv: real intensity but INFLATED internal q_value (0.5)
        fs::write(
            dir.join("lfq.tsv"),
            "peptide\tcharge\tproteins\tq_value\tscore\tspectral_angle\tS.mzML\n\
             PEPTIDEK\t2\tsp|P1|A_HUMAN\t0.5\t10\t0.9\t1000.0\n",
        )
        .unwrap();
        // results.sage.tsv: the SAME peptide is a confident target ID (q = 0.001)
        fs::write(
            dir.join("results.sage.tsv"),
            "peptide\tproteins\tlabel\tspectrum_q\tpeptide_q\tprotein_q\n\
             PEPTIDEK\tsp|P1|A_HUMAN\t1\t0.001\t0.001\t0.001\n",
        )
        .unwrap();
        let sd = parse_lfq_sample("S", &dir, GroupKeyMode::MemberSet).expect("parse");
        assert_eq!(sd.observations.len(), 1);
        let obs = &sd.observations[0];
        assert!(
            obs.peptide_q <= 0.01,
            "peptide_q must come from search, got {}",
            obs.peptide_q
        );
        assert!(
            obs.protein_q <= 0.01,
            "protein_q must come from search, got {}",
            obs.protein_q
        );
        // and the matrix must KEEP this peptide at q<=0.01 (it would be dropped on lfq_q=0.5)
        let names = vec!["S".to_string()];
        let (gk, _) =
            build_protein_group_matrix(std::slice::from_ref(&sd), &names, Some(0.01), "protein_q");
        assert_eq!(gk.len(), 1, "group must survive at protein_q<=0.01");
        let _ = &dir; // Keep test artifacts for release review.
    }

    /// Regression guard: the LFQ intensity column is named after the RAW FILE, and
    /// Bruker raw is a `.d` directory, not an `.mzML`.
    ///
    /// This code used to locate that column with `h.ends_with(".mzML")`. On 2026-08-04
    /// that silently dropped all 102 samples of a timsTOF cohort -- the engine printed
    /// "No valid data!" for a search that had quantified 112,327 LFQ features, and
    /// nothing in the output pointed at the file extension. The column is now found by
    /// elimination, so an unfamiliar vendor extension cannot repeat it.
    #[test]
    fn test_lfq_intensity_column_any_vendor_extension() {
        use std::fs;
        for (tag, col) in [
            ("mzml", "S.mzML"),
            ("bruker_d", "20230916_TIMS03_PeTr_SA_T760_100.d"),
            ("thermo_raw", "S.raw"),
            ("lowercase", "S.mzml"),
        ] {
            let dir =
                std::env::temp_dir().join(format!("agg_lfqcol_{}_{}", tag, std::process::id()));
            let _ = fs::create_dir_all(&dir);
            fs::write(
                dir.join("lfq.tsv"),
                format!(
                    "peptide\tcharge\tproteins\tq_value\tscore\tspectral_angle\t{col}\n\
                     PEPTIDEK\t2\tsp|P1|A_HUMAN\t0.01\t10\t0.9\t1234.0\n"
                ),
            )
            .unwrap();
            fs::write(
                dir.join("results.sage.tsv"),
                "peptide\tproteins\tlabel\tspectrum_q\tpeptide_q\tprotein_q\n\
                 PEPTIDEK\tsp|P1|A_HUMAN\t1\t0.001\t0.001\t0.001\n",
            )
            .unwrap();

            let sd = parse_lfq_sample("S", &dir, GroupKeyMode::MemberSet)
                .unwrap_or_else(|| panic!("{tag}: sample was DROPPED for column {col}"));
            assert_eq!(sd.observations.len(), 1, "{tag}: expected one observation");
            assert_eq!(
                sd.observations[0].intensity, 1234.0,
                "{tag}: intensity must be read from column {col}"
            );
            let _ = &dir; // Keep test artifacts for release review.
        }
    }

    /// A lfq.tsv carrying ONLY the fixed columns has no intensity column at all.
    /// That must be reported and the sample skipped -- not treated as zero intensity.
    #[test]
    fn test_lfq_missing_intensity_column_is_rejected() {
        use std::fs;
        let dir = std::env::temp_dir().join(format!("agg_lfqnocol_{}", std::process::id()));
        let _ = fs::create_dir_all(&dir);
        fs::write(
            dir.join("lfq.tsv"),
            "peptide\tcharge\tproteins\tq_value\tscore\tspectral_angle\n\
             PEPTIDEK\t2\tsp|P1|A_HUMAN\t0.01\t10\t0.9\n",
        )
        .unwrap();
        fs::write(
            dir.join("results.sage.tsv"),
            "peptide\tproteins\tlabel\tspectrum_q\tpeptide_q\tprotein_q\n\
             PEPTIDEK\tsp|P1|A_HUMAN\t1\t0.001\t0.001\t0.001\n",
        )
        .unwrap();
        assert!(
            parse_lfq_sample("S", &dir, GroupKeyMode::MemberSet).is_none(),
            "a lfq.tsv with no intensity column must be rejected"
        );
        let _ = &dir; // Keep test artifacts for release review.
    }

    /// Determinism guard (D1): matrix ROW ORDER must be sorted and independent
    /// of observation insertion order, which seeds the per-sample `local`
    /// HashMap's randomized iteration order. Feeding identical observations in
    /// two different orders must yield the identical, sorted key sequence for
    /// both the protein-group and peptide matrices. Cell contents are unchanged.
    #[test]
    fn test_matrix_row_order_is_deterministic_and_sorted() {
        let mk = |pg: &str, pep: &str, intensity: f64| PeptideObs {
            peptide: pep.to_string(),
            protein_group: pg.to_string(),
            intensity,
            q_value: 0.001,
            spectrum_q: 0.001,
            peptide_q: 0.001,
            protein_q: 0.001,
        };
        // Same content, two different insertion orders across two samples.
        let obs_a = vec![
            mk("zebra", "PEPK", 10.0),
            mk("alpha", "AAAK", 20.0),
            mk("mid", "MMMK", 30.0),
        ];
        let obs_b = vec![
            mk("mid", "MMMK", 30.0),
            mk("zebra", "PEPK", 10.0),
            mk("alpha", "AAAK", 20.0),
        ];
        let names = vec!["S1".to_string(), "S2".to_string()];
        let sd_a = SampleData {
            sample_name: "S1".to_string(),
            observations: obs_a,
            n_peptides: 3,
            n_groups: 3,
        };
        let sd_b = SampleData {
            sample_name: "S2".to_string(),
            observations: obs_b,
            n_peptides: 3,
            n_groups: 3,
        };
        let data = vec![sd_a, sd_b];

        // Protein-group matrix: keys sorted and stable.
        let (gk, gmap) = build_protein_group_matrix(&data, &names, Some(0.01), "protein_q");
        assert_eq!(
            gk,
            vec!["alpha".to_string(), "mid".to_string(), "zebra".to_string()]
        );
        let mut gk_sorted = gk.clone();
        gk_sorted.sort();
        assert_eq!(gk, gk_sorted, "protein-group row order must be sorted");
        // Content preserved: alpha seen in both samples at 20.0, none dropped.
        assert_eq!(gmap.len(), 3);
        assert_eq!(gmap["alpha"], vec![Some(20.0), Some(20.0)]);

        // Peptide matrix: (group, peptide) keys sorted and stable.
        let (pk, pmap) = build_peptide_matrix(&data, &names, Some(0.01), "peptide_q");
        assert_eq!(
            pk,
            vec![
                ("alpha".to_string(), "AAAK".to_string()),
                ("mid".to_string(), "MMMK".to_string()),
                ("zebra".to_string(), "PEPK".to_string()),
            ]
        );
        let mut pk_sorted = pk.clone();
        pk_sorted.sort();
        assert_eq!(pk, pk_sorted, "peptide row order must be sorted");
        assert_eq!(pmap.len(), 3);
    }

    /// Default (flag omitted) MUST reproduce the historical threshold set, in
    /// the exact order and with the exact labels, so default output is
    /// byte-identical to before the flag existed.
    #[test]
    fn test_q_thresholds_default_is_historical() {
        let got = parse_q_thresholds(None).expect("default parse");
        assert_eq!(
            got,
            vec![
                (Some(0.01), "q0.01".to_string()),
                (Some(0.05), "q0.05".to_string()),
                (None, "unfiltered".to_string()),
            ]
        );
    }

    /// Passing the historical spec explicitly must yield the identical set.
    #[test]
    fn test_q_thresholds_explicit_matches_default() {
        let got = parse_q_thresholds(Some("0.01,0.05,unfiltered")).expect("explicit parse");
        assert_eq!(got, parse_q_thresholds(None).unwrap());
    }

    #[test]
    fn test_q_thresholds_custom_subset_and_order() {
        let got = parse_q_thresholds(Some("0.01")).expect("single parse");
        assert_eq!(got, vec![(Some(0.01), "q0.01".to_string())]);

        // Order is preserved; 'none' is an alias for unfiltered.
        let got = parse_q_thresholds(Some("0.05, none")).expect("parse");
        assert_eq!(
            got,
            vec![
                (Some(0.05), "q0.05".to_string()),
                (None, "unfiltered".to_string()),
            ]
        );
    }

    #[test]
    fn test_q_thresholds_dedup_and_whitespace() {
        // Duplicate labels collapse; surrounding whitespace and empty tokens are ignored.
        let got = parse_q_thresholds(Some(" 0.01 , 0.01 ,, unfiltered ")).expect("parse");
        assert_eq!(
            got,
            vec![
                (Some(0.01), "q0.01".to_string()),
                (None, "unfiltered".to_string()),
            ]
        );
    }

    #[test]
    fn test_q_thresholds_rejects_bad_input() {
        assert!(parse_q_thresholds(Some("abc")).is_err());
        assert!(parse_q_thresholds(Some("0")).is_err());
        assert!(parse_q_thresholds(Some("1.5")).is_err());
        assert!(parse_q_thresholds(Some("")).is_err());
        assert!(parse_q_thresholds(Some(" , ")).is_err());
    }

    #[test]
    fn test_bh_correction_small() {
        let sample_data = vec![SampleData {
            sample_name: "test".into(),
            observations: vec![
                PeptideObs {
                    peptide: "A".into(),
                    protein_group: "P1".into(),
                    intensity: 100.0,
                    q_value: 0.01,
                    spectrum_q: 0.0,
                    peptide_q: 0.0,
                    protein_q: 0.0,
                },
                PeptideObs {
                    peptide: "B".into(),
                    protein_group: "P2".into(),
                    intensity: 200.0,
                    q_value: 0.04,
                    spectrum_q: 0.0,
                    peptide_q: 0.0,
                    protein_q: 0.0,
                },
                PeptideObs {
                    peptide: "C".into(),
                    protein_group: "P3".into(),
                    intensity: 300.0,
                    q_value: 0.03,
                    spectrum_q: 0.0,
                    peptide_q: 0.0,
                    protein_q: 0.0,
                },
                PeptideObs {
                    peptide: "D".into(),
                    protein_group: "P4".into(),
                    intensity: 400.0,
                    q_value: 0.005,
                    spectrum_q: 0.0,
                    peptide_q: 0.0,
                    protein_q: 0.0,
                },
                PeptideObs {
                    peptide: "E".into(),
                    protein_group: "P5".into(),
                    intensity: 500.0,
                    q_value: 0.10,
                    spectrum_q: 0.0,
                    peptide_q: 0.0,
                    protein_q: 0.0,
                },
            ],
            n_peptides: 5,
            n_groups: 5,
        }];

        let adjusted = apply_bh_correction(&sample_data, "q_value");
        assert_eq!(adjusted.len(), 5);
        for v in adjusted.values() {
            assert!(*v <= 1.0);
            assert!(*v >= 0.0);
        }
    }

    #[test]
    fn test_bh_monotonic() {
        let sample_data = vec![SampleData {
            sample_name: "test".into(),
            observations: vec![
                PeptideObs {
                    peptide: "A".into(),
                    protein_group: "P1".into(),
                    intensity: 1.0,
                    q_value: 0.001,
                    spectrum_q: 0.0,
                    peptide_q: 0.0,
                    protein_q: 0.0,
                },
                PeptideObs {
                    peptide: "B".into(),
                    protein_group: "P2".into(),
                    intensity: 1.0,
                    q_value: 0.50,
                    spectrum_q: 0.0,
                    peptide_q: 0.0,
                    protein_q: 0.0,
                },
                PeptideObs {
                    peptide: "C".into(),
                    protein_group: "P3".into(),
                    intensity: 1.0,
                    q_value: 0.003,
                    spectrum_q: 0.0,
                    peptide_q: 0.0,
                    protein_q: 0.0,
                },
                PeptideObs {
                    peptide: "D".into(),
                    protein_group: "P4".into(),
                    intensity: 1.0,
                    q_value: 0.002,
                    spectrum_q: 0.0,
                    peptide_q: 0.0,
                    protein_q: 0.0,
                },
            ],
            n_peptides: 4,
            n_groups: 4,
        }];

        let adjusted = apply_bh_correction(&sample_data, "q_value");

        let mut pairs: Vec<(f64, f64)> = Vec::new();
        for (oi, obs) in sample_data[0].observations.iter().enumerate() {
            let adj = adjusted[&(0, oi)];
            pairs.push((obs.q_value, adj));
        }
        pairs.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

        for i in 1..pairs.len() {
            assert!(pairs[i].1 >= pairs[i - 1].1);
        }
    }

    /// C5 regression guard: for LFQ, --bh-correction must adjust the SEARCH-derived
    /// peptide_q / protein_q columns — the fields the LFQ matrices actually gate on —
    /// NOT the lfq.tsv-internal q_value. This mirrors the exact loop in main()'s
    /// DataType::Lfq branch. Before the fix, only q_value was overwritten, so the
    /// correction was a silent no-op that never touched gating.
    #[test]
    fn test_bh_correction_lfq_adjusts_gated_columns() {
        // Ten observations with small raw peptide_q/protein_q values. BH scales
        // q up by m/rank, so at least some corrected values must exceed the raw.
        let mk = |i: usize, q: f64| PeptideObs {
            peptide: format!("PEP{i}"),
            protein_group: format!("P{i}"),
            intensity: 100.0,
            q_value: 0.5, // inflated lfq-internal q — must be IGNORED by the fix
            spectrum_q: 0.0,
            peptide_q: q,
            protein_q: q,
        };
        let obs: Vec<PeptideObs> = (0..10).map(|i| mk(i, 0.001 * (i as f64 + 1.0))).collect();
        let raw_pep: Vec<f64> = obs.iter().map(|o| o.peptide_q).collect();
        let raw_prot: Vec<f64> = obs.iter().map(|o| o.protein_q).collect();
        let mut sample_data = vec![SampleData {
            sample_name: "S".into(),
            observations: obs,
            n_peptides: 10,
            n_groups: 10,
        }];

        // Exact replica of main()'s DataType::Lfq BH loop.
        for q_field in &["peptide_q", "protein_q"] {
            let adjusted = apply_bh_correction(&sample_data, q_field);
            for (si, sd) in sample_data.iter_mut().enumerate() {
                for (oi, o) in sd.observations.iter_mut().enumerate() {
                    if let Some(&adj) = adjusted.get(&(si, oi)) {
                        match *q_field {
                            "peptide_q" => o.peptide_q = adj,
                            "protein_q" => o.protein_q = adj,
                            _ => {}
                        }
                    }
                }
            }
        }

        // The gated columns MUST have changed (BH inflates the higher-rank q-values).
        let new_pep: Vec<f64> = sample_data[0]
            .observations
            .iter()
            .map(|o| o.peptide_q)
            .collect();
        let new_prot: Vec<f64> = sample_data[0]
            .observations
            .iter()
            .map(|o| o.protein_q)
            .collect();
        assert_ne!(new_pep, raw_pep, "BH must change peptide_q for LFQ");
        assert_ne!(new_prot, raw_prot, "BH must change protein_q for LFQ");
        // The lfq-internal q_value must be left untouched (not used for gating).
        assert!(
            sample_data[0].observations.iter().all(|o| o.q_value == 0.5),
            "q_value must be ignored by the LFQ BH fix"
        );
    }

    #[test]
    fn test_strip_source_tags() {
        assert_eq!(strip_source_tags("REF|SP|P02768"), "P02768");
        assert_eq!(
            strip_source_tags("sp|P02768|ALBU_HUMAN"),
            "sp|P02768|ALBU_HUMAN"
        );
        assert_eq!(
            strip_source_tags("MGYG000001373_00501"),
            "MGYG000001373_00501"
        );
    }

    #[test]
    fn test_normalize_protein_group() {
        // Default MemberSet mode: full sorted, deduped set — legacy behavior.
        let pg = normalize_protein_group("REF|SP|B;REF|SP|A;REF|SP|A", GroupKeyMode::MemberSet);
        assert_eq!(pg, "A;B");
    }

    // ========================================================================
    // Group-key-mode: member_set (default) vs representative
    // ========================================================================

    #[test]
    fn test_group_key_mode_member_set_is_legacy() {
        // MemberSet must reproduce the exact legacy output: sorted, deduped,
        // semicolon-joined full set. This guards default behavior from drifting.
        assert_eq!(
            normalize_protein_group("MGYG002;MGYG001;MGYG003", GroupKeyMode::MemberSet),
            "MGYG001;MGYG002;MGYG003"
        );
        // dedup + source-tag stripping preserved
        assert_eq!(
            normalize_protein_group(
                "CLUSTER97|CAPSCAN|B;CLUSTER97|CAPSCAN|A;CLUSTER97|CAPSCAN|A",
                GroupKeyMode::MemberSet
            ),
            "A;B"
        );
        // single member is itself
        assert_eq!(
            normalize_protein_group("sp|P02768|ALBU_HUMAN", GroupKeyMode::MemberSet),
            "sp|P02768|ALBU_HUMAN"
        );
    }

    #[test]
    fn test_group_key_mode_representative_picks_lexicographically_smallest() {
        // Representative collapses the group to the lexicographically-smallest
        // stripped accession.
        assert_eq!(
            normalize_protein_group("MGYG002;MGYG001;MGYG003", GroupKeyMode::Representative),
            "MGYG001"
        );
        // Source tags stripped before comparison.
        assert_eq!(
            normalize_protein_group(
                "CLUSTER97|CAPSCAN|B;CLUSTER97|CAPSCAN|A",
                GroupKeyMode::Representative
            ),
            "A"
        );
        // Single member: same as that member.
        assert_eq!(
            normalize_protein_group("sp|P02768|ALBU_HUMAN", GroupKeyMode::Representative),
            "sp|P02768|ALBU_HUMAN"
        );
    }

    #[test]
    fn test_group_key_mode_representative_collapses_cooccurring_orthologs() {
        // Two samples observe different SUBSETS of the same near-identical ortholog
        // cluster. Under MemberSet they form DIFFERENT keys (not comparable);
        // under Representative they collapse to the SAME stable anchor.
        let sample_a = "MGYG001;MGYG002;MGYG003";
        let sample_b = "MGYG001;MGYG002"; // MGYG003 not observed in sample B

        // MemberSet: keys differ → orthologs do NOT line up across samples
        assert_ne!(
            normalize_protein_group(sample_a, GroupKeyMode::MemberSet),
            normalize_protein_group(sample_b, GroupKeyMode::MemberSet),
        );

        // Representative: both collapse to MGYG001 → comparable across samples
        assert_eq!(
            normalize_protein_group(sample_a, GroupKeyMode::Representative),
            normalize_protein_group(sample_b, GroupKeyMode::Representative),
        );
        assert_eq!(
            normalize_protein_group(sample_a, GroupKeyMode::Representative),
            "MGYG001"
        );
    }

    #[test]
    fn test_group_key_mode_empty_field_is_empty_in_both_modes() {
        assert_eq!(normalize_protein_group("", GroupKeyMode::MemberSet), "");
        assert_eq!(
            normalize_protein_group("", GroupKeyMode::Representative),
            ""
        );
    }

    #[test]
    fn test_classify_group() {
        assert_eq!(classify_group("sp|P02768|ALBU_HUMAN"), "Human");
        assert_eq!(classify_group("MGYG000001373_00501"), "Microbial");
        assert_eq!(
            classify_group("sp|P02768|ALBU_HUMAN;MGYG000001373_00501"),
            "Mixed"
        );

        // A uniform hash accession carries no source. It must read Unknown, not
        // Microbial: v26 CAPSCAN labelled all 2,041,914 groups Microbial and the
        // host fraction silently became zero.
        let h1 = format!("CAPSCAN_{}", "a".repeat(64));
        let h2 = format!("CAPSCAN_{}", "b".repeat(64));
        assert_eq!(classify_group(&h1), "Unknown");
        assert_eq!(classify_group(&format!("{h1};{h2}")), "Unknown");
        // A readable member alongside hashes still decides; hashes must not vote.
        assert_eq!(
            classify_group(&format!("{h1};sp|P02768|ALBU_HUMAN")),
            "Human"
        );
        assert_eq!(
            classify_group(&format!("{h1};MGYG000001373_00501")),
            "Microbial"
        );
        assert_eq!(
            classify_group(&format!("{h1};sp|P02768|ALBU_HUMAN;MGYG000001373_00501")),
            "Mixed"
        );

        // The hash detector must not blind the classifier to ordinary accessions
        // that merely contain an underscore.
        assert!(!is_classifiable_accession(&format!(
            "HSTOOL_{}",
            "f".repeat(64)
        )));
        assert!(is_classifiable_accession("MGYG000001373_00501"));
        assert!(is_classifiable_accession("GMGC10.210_551_413.ATPD"));
        assert!(is_classifiable_accession("GENCODE|ENSP00000493376.2"));
        assert!(is_classifiable_accession(&format!("X_{}", "a".repeat(63))));
        assert!(is_classifiable_accession(&format!("X_{}", "a".repeat(65))));
        assert!(is_classifiable_accession(&format!("X_{}z", "a".repeat(63))));
    }

    /// With a host list, a hash-keyed group gets the RIGHT answer, not Unknown.
    /// This is the whole point of --host-accessions: the classification column
    /// stays correct so nothing downstream has to change.
    #[test]
    fn test_classify_group_with_hosts() {
        let h_host = format!("CAPSCAN_{}", "a".repeat(64));
        let h_mic1 = format!("CAPSCAN_{}", "b".repeat(64));
        let h_mic2 = format!("CAPSCAN_{}", "c".repeat(64));
        let h_both = format!("CAPSCAN_{}", "d".repeat(64));
        let mut hosts: HashMap<String, String> = HashMap::new();
        hosts.insert(h_host.clone(), "Human".to_string());
        hosts.insert(h_both.clone(), "Mixed".to_string());

        // listed -> host; absent -> microbial, because the list is complete
        assert_eq!(classify_group_with_hosts(&h_host, &hosts), "Human");
        assert_eq!(classify_group_with_hosts(&h_mic1, &hosts), "Microbial");
        // a group spanning both is Mixed
        assert_eq!(
            classify_group_with_hosts(&format!("{h_host};{h_mic1}"), &hosts),
            "Mixed"
        );
        // several unlisted members stay microbial, they do not accumulate into host
        assert_eq!(
            classify_group_with_hosts(&format!("{h_mic1};{h_mic2}"), &hosts),
            "Microbial"
        );
        // a single sequence produced by host AND non-host sources is Mixed on its own
        assert_eq!(classify_group_with_hosts(&h_both, &hosts), "Mixed");

        // An empty list must yield an all-microbial answer, which is what the
        // "matched nothing" warning is there to catch.
        let empty: HashMap<String, String> = HashMap::new();
        assert_eq!(classify_group_with_hosts(&h_host, &empty), "Microbial");
        // Source-based host detection: tissue-isoform union members are Human, not Microbial.
        assert_eq!(classify_group("GENCODE|ENSP00000493376.2|ENST00000641515.2|ENSG00000186092.7|OTTHUMG00000001094.4|OTTHUMT00000003223.4|OR4F5-201|OR4F5|326"), "Human");
        assert_eq!(classify_group("REFSEQ|NP_001005484.2|OR4F5|305"), "Human");
        assert_eq!(classify_group("CHESS|CHS.1.1|GENE|100"), "Human");
        assert_eq!(classify_group("UNIPROT|P02768-2|ALB|500"), "Human");
        // Isoform + microbial => Mixed (previously the isoform would have read as Microbial).
        assert_eq!(
            classify_group("GENCODE|ENSP00000493376.2|OR4F5|326;MGYG000001373_00501"),
            "Mixed"
        );
        // A microbial tr| (not _HUMAN) stays Microbial.
        assert_eq!(classify_group("tr|A0A1234|A0A1234_ECOLI"), "Microbial");
    }

    /// The `representative` column is an arbitrary deterministic display label and
    /// must carry NO origin information (I5/I7/I8, 2026-07-28). These assertions
    /// exist to make a reintroduction of the human-first ranking fail the suite.
    #[test]
    fn test_representative_protein_is_unbiased_display_label() {
        // Lexicographically smallest member, full stop.
        assert_eq!(representative_protein("MGYG002;MGYG001;MGYG003"), "MGYG001");
        assert_eq!(
            representative_protein("sp|P02768|ALBU_HUMAN"),
            "sp|P02768|ALBU_HUMAN"
        );

        // NEGATIVE CONTROL for the deleted rule: a Mixed group must NOT resolve to
        // its human member. Under the old ranking these returned the sp|/tr| _HUMAN
        // accession; under the display-label rule they return the smallest member.
        assert_eq!(
            representative_protein("sp|P02768|ALBU_HUMAN;MGYG000001373_00501"),
            "MGYG000001373_00501"
        );
        assert_eq!(
            representative_protein("tr|A0A1|A0A1_HUMAN;GMGC10.001_01;sp|P0A7|P0A7_ECOLI"),
            "GMGC10.001_01"
        );
        // Human SwissProt no longer outranks non-human SwissProt.
        assert_eq!(
            representative_protein("sp|P02768|ALBU_HUMAN;sp|P0A799|PGK_ECOLI"),
            "sp|P02768|ALBU_HUMAN"
        );
        assert_eq!(
            representative_protein("sp|Z99999|ZZZZ_HUMAN;sp|P0A799|PGK_ECOLI"),
            "sp|P0A799|PGK_ECOLI"
        );

        // Deterministic under member reordering (order-independence).
        assert_eq!(
            representative_protein("C;A;B"),
            representative_protein("B;C;A")
        );

        // Degenerate inputs behave as before.
        assert_eq!(representative_protein(""), "");
        assert_eq!(representative_protein("A;;B"), "A");
    }

    // ========================================================================
    // Functional-lookup tests
    // ========================================================================

    fn make_table(pairs: &[(&str, &[(&str, &str)])]) -> FunctionalTable {
        let mut t: FunctionalTable = HashMap::new();
        for (acc, lvls) in pairs {
            let mut m: HashMap<String, Vec<String>> = HashMap::new();
            for (lvl, raw) in *lvls {
                let vals = split_functional_value(raw);
                if !vals.is_empty() {
                    m.insert(lvl.to_string(), vals);
                }
            }
            t.insert(acc.to_string(), m);
        }
        t
    }

    #[test]
    fn test_split_functional_value_basic() {
        assert_eq!(split_functional_value(""), Vec::<String>::new());
        assert_eq!(split_functional_value("-"), Vec::<String>::new());
        assert_eq!(split_functional_value("NA"), Vec::<String>::new());
        assert_eq!(
            split_functional_value("COG0107"),
            vec!["COG0107".to_string()]
        );
        assert_eq!(
            split_functional_value("COG0107;COG0200"),
            vec!["COG0107".to_string(), "COG0200".to_string()]
        );
        assert_eq!(
            split_functional_value("ko00010,ko00020,-"),
            vec!["ko00010".to_string(), "ko00020".to_string()]
        );
        assert_eq!(
            split_functional_value(" COG0107 ; COG0200 "),
            vec!["COG0107".to_string(), "COG0200".to_string()]
        );
    }

    #[test]
    fn test_group_value_votes_plurality() {
        let table = make_table(&[
            ("A1", &[("OG", "COG0107")]),
            ("A2", &[("OG", "COG0107")]),
            ("A3", &[("OG", "COG0200")]),
        ]);
        let votes = group_value_votes("A1;A2;A3", "OG", &table);
        let as_map: HashMap<&str, usize> = votes.iter().map(|(v, c)| (v.as_str(), *c)).collect();
        assert_eq!(as_map["COG0107"], 2);
        assert_eq!(as_map["COG0200"], 1);
    }

    #[test]
    fn test_group_value_votes_empty_when_unannotated() {
        let table = make_table(&[("A", &[("OG", "-")]), ("B", &[])]);
        let votes = group_value_votes("A;B;NOT_IN_TABLE", "OG", &table);
        assert!(votes.is_empty());
    }

    #[test]
    fn test_canonical_tie_alphabetical() {
        let table = make_table(&[("X", &[("OG", "COG0200")]), ("Y", &[("OG", "COG0107")])]);
        let canon = canonical_for_group("X;Y", "OG", &table);
        assert_eq!(canon.winner.as_deref(), Some("COG0107"));
        // Values sorted by BTreeMap ordering (alphabetical)
        assert_eq!(
            canon.values,
            vec!["COG0107".to_string(), "COG0200".to_string()]
        );
    }

    #[test]
    fn test_canonical_single_member() {
        let table = make_table(&[("SOLO", &[("OG", "COG9999")])]);
        let canon = canonical_for_group("SOLO", "OG", &table);
        assert_eq!(canon.winner.as_deref(), Some("COG9999"));
        assert_eq!(canon.values, vec!["COG9999".to_string()]);
    }

    #[test]
    fn test_canonical_all_missing() {
        let table = make_table(&[("A", &[("OG", "-")]), ("B", &[])]);
        let canon = canonical_for_group("A;B;NOT_IN_TABLE", "OG", &table);
        assert!(canon.winner.is_none());
        assert!(canon.values.is_empty());
    }

    #[test]
    fn test_canonical_semicolon_source_counts_distinct() {
        let table = make_table(&[
            ("A", &[("OG", "COG0107;COG0200")]),
            ("B", &[("OG", "COG0107")]),
        ]);
        let canon = canonical_for_group("A;B", "OG", &table);
        assert_eq!(canon.winner.as_deref(), Some("COG0107"));
        assert_eq!(canon.values.len(), 2);
    }

    // ========================================================================
    // Resolution: stable Unannotated_N ids
    // ========================================================================

    #[test]
    fn test_resolve_level_unannotated_numbering_stable() {
        let table = make_table(&[
            ("ANNOT1", &[("OG", "COG0001")]),
            ("ANNOT2", &[("OG", "COG0002")]),
            // The rest are deliberately missing
        ]);

        // Feed groups in sorted order — Unannotated_N must be 1,2,3... in
        // deterministic sorted-key order.
        let mut group_keys = vec![
            "ANNOT1".to_string(),
            "UNK_A".to_string(),
            "UNK_B".to_string(),
            "ANNOT2".to_string(),
            "UNK_C;UNK_D".to_string(),
        ];
        group_keys.sort();

        let res = resolve_level("OG", &table, &group_keys);

        // Annotated groups keep their canonical value
        assert_eq!(res.per_group["ANNOT1"].plurality_row, "COG0001");
        assert_eq!(res.per_group["ANNOT2"].plurality_row, "COG0002");

        // Unannotated groups in sorted order: UNK_A, UNK_B, UNK_C;UNK_D
        // → Unannotated_000001, 000002, 000003
        assert_eq!(res.per_group["UNK_A"].plurality_row, "Unannotated_000001");
        assert_eq!(res.per_group["UNK_B"].plurality_row, "Unannotated_000002");
        assert_eq!(
            res.per_group["UNK_C;UNK_D"].plurality_row,
            "Unannotated_000003"
        );
    }

    #[test]
    fn test_resolve_level_split_weights() {
        let table = make_table(&[
            // Group with 2 values → 0.5 weight each
            ("A1", &[("OG", "COG_X")]),
            ("A2", &[("OG", "COG_Y")]),
            // Group with 1 value → 1.0 weight
            ("B1", &[("OG", "COG_Z")]),
            ("B2", &[("OG", "COG_Z")]),
        ]);
        let group_keys = vec!["A1;A2".to_string(), "B1;B2".to_string()];
        let res = resolve_level("OG", &table, &group_keys);

        let a = &res.per_group["A1;A2"];
        assert_eq!(a.split_rows.len(), 2);
        for (_, w) in &a.split_rows {
            assert!((w - 0.5).abs() < 1e-9);
        }

        let b = &res.per_group["B1;B2"];
        assert_eq!(b.split_rows.len(), 1);
        assert!((b.split_rows[0].1 - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_resolve_level_unannotated_split_is_full_weight() {
        let table = make_table(&[]);
        let group_keys = vec!["X".to_string()];
        let res = resolve_level("OG", &table, &group_keys);
        let r = &res.per_group["X"];
        assert_eq!(r.plurality_row, "Unannotated_000001");
        assert_eq!(r.split_rows.len(), 1);
        assert_eq!(r.split_rows[0].0, "Unannotated_000001");
        assert!((r.split_rows[0].1 - 1.0).abs() < 1e-9);
    }

    // ========================================================================
    // Matrix building — plurality
    // ========================================================================

    #[test]
    fn test_build_functional_group_matrix_plurality() {
        let table = make_table(&[
            ("P1", &[("OG", "COG0107")]),
            ("P2", &[("OG", "COG0107")]),
            ("P3", &[("OG", "COG0200")]),
        ]);
        let group_keys = vec!["P1;P2".to_string(), "P3".to_string()];
        let mut gm: HashMap<String, Vec<Option<f64>>> = HashMap::new();
        gm.insert("P1;P2".to_string(), vec![Some(100.0), Some(200.0)]);
        gm.insert("P3".to_string(), vec![Some(50.0), None]);

        let res = resolve_level("OG", &table, &group_keys);
        let (keys, m) =
            build_functional_group_matrix(&res, &group_keys, &gm, 2, AggregationMode::Plurality);

        assert!(keys.contains(&"COG0107".to_string()));
        assert!(keys.contains(&"COG0200".to_string()));

        assert_eq!(m["COG0107"][0], Some(100.0));
        assert_eq!(m["COG0107"][1], Some(200.0));
        assert_eq!(m["COG0200"][0], Some(50.0));
        assert_eq!(m["COG0200"][1], None);
    }

    #[test]
    fn test_build_functional_group_matrix_splitintensity() {
        // P1;P2 has TWO distinct OGs (COG_X via P1, COG_Y via P2)
        // With splitintensity: each OG gets HALF the group's intensity
        let table = make_table(&[("P1", &[("OG", "COG_X")]), ("P2", &[("OG", "COG_Y")])]);
        let group_keys = vec!["P1;P2".to_string()];
        let mut gm: HashMap<String, Vec<Option<f64>>> = HashMap::new();
        gm.insert("P1;P2".to_string(), vec![Some(100.0), Some(200.0)]);

        let res = resolve_level("OG", &table, &group_keys);
        let (keys, m) = build_functional_group_matrix(
            &res,
            &group_keys,
            &gm,
            2,
            AggregationMode::SplitIntensity,
        );

        assert_eq!(keys, vec!["COG_X".to_string(), "COG_Y".to_string()]);
        assert_eq!(m["COG_X"], vec![Some(50.0), Some(100.0)]);
        assert_eq!(m["COG_Y"], vec![Some(50.0), Some(100.0)]);

        // Sum across rows = 2x the input (each value duplicated by 1/N * N) — WRONG
        // Actually: 1/2 * 100 + 1/2 * 100 = 100 → conservation holds
        let s0: f64 = m.values().filter_map(|r| r[0]).sum();
        let s1: f64 = m.values().filter_map(|r| r[1]).sum();
        assert!((s0 - 100.0).abs() < 1e-9);
        assert!((s1 - 200.0).abs() < 1e-9);
    }

    #[test]
    fn test_unannotated_kept_individual_not_lumped() {
        // Two distinct unannotated groups must produce TWO distinct rows.
        let table = make_table(&[("KNOWN", &[("OG", "COG0001")])]);
        let group_keys = vec![
            "KNOWN".to_string(),
            "UNK_ALPHA".to_string(),
            "UNK_BETA".to_string(),
        ];
        let mut gm: HashMap<String, Vec<Option<f64>>> = HashMap::new();
        gm.insert("KNOWN".to_string(), vec![Some(1000.0)]);
        gm.insert("UNK_ALPHA".to_string(), vec![Some(500.0)]);
        gm.insert("UNK_BETA".to_string(), vec![Some(300.0)]);

        let res = resolve_level("OG", &table, &group_keys);
        let (keys, m) =
            build_functional_group_matrix(&res, &group_keys, &gm, 1, AggregationMode::Plurality);

        // Should have THREE rows, NOT two
        assert_eq!(keys.len(), 3);
        assert!(keys.iter().any(|k| k == "COG0001"));
        // Two distinct Unannotated_ rows
        let unannot_rows: Vec<&String> = keys
            .iter()
            .filter(|k| k.starts_with("Unannotated_"))
            .collect();
        assert_eq!(unannot_rows.len(), 2);

        // Intensity correctly placed per row
        assert_eq!(m["COG0001"][0], Some(1000.0));
        let sum_unannot: f64 = unannot_rows.iter().map(|k| m[*k][0].unwrap_or(0.0)).sum();
        assert!((sum_unannot - 800.0).abs() < 1e-9);
    }

    #[test]
    fn test_intensity_conservation_plurality() {
        // Plurality assigns ALL intensity to one row per group → conservation holds exactly.
        let table = make_table(&[
            ("P1", &[("OG", "COG_A")]),
            ("P2", &[("OG", "COG_B")]),
            ("P3", &[("OG", "COG_A")]),
        ]);
        let group_keys = vec![
            "P1".to_string(),
            "P2".to_string(),
            "P3".to_string(),
            "P4_UNK".to_string(),
        ];
        let mut g: HashMap<String, Vec<Option<f64>>> = HashMap::new();
        g.insert("P1".to_string(), vec![Some(10.0), Some(20.0)]);
        g.insert("P2".to_string(), vec![Some(30.0), None]);
        g.insert("P3".to_string(), vec![None, Some(40.0)]);
        g.insert("P4_UNK".to_string(), vec![Some(5.0), Some(5.0)]);

        let group_total: f64 = g.values().flatten().filter_map(|x| x.as_ref()).sum();
        let res = resolve_level("OG", &table, &group_keys);
        let (_k, m) =
            build_functional_group_matrix(&res, &group_keys, &g, 2, AggregationMode::Plurality);
        let func_total: f64 = m.values().flatten().filter_map(|x| x.as_ref()).sum();

        assert!((group_total - func_total).abs() < 1e-9);
    }

    #[test]
    fn test_intensity_conservation_splitintensity() {
        // SplitIntensity divides each group's intensity across its N values → conservation holds.
        let table = make_table(&[
            ("P1", &[("OG", "COG_A")]),
            ("P2", &[("OG", "COG_B")]),
            ("P3", &[("OG", "COG_C")]),
        ]);
        // One merged group with 3 distinct OGs
        let group_keys = vec!["P1;P2;P3".to_string()];
        let mut g: HashMap<String, Vec<Option<f64>>> = HashMap::new();
        g.insert("P1;P2;P3".to_string(), vec![Some(90.0), Some(300.0)]);

        let group_total: f64 = g.values().flatten().filter_map(|x| x.as_ref()).sum();
        let res = resolve_level("OG", &table, &group_keys);
        let (_k, m) = build_functional_group_matrix(
            &res,
            &group_keys,
            &g,
            2,
            AggregationMode::SplitIntensity,
        );
        let func_total: f64 = m.values().flatten().filter_map(|x| x.as_ref()).sum();

        assert!(
            (group_total - func_total).abs() < 1e-9,
            "splitintensity must conserve total: {} vs {}",
            group_total,
            func_total
        );

        // 3 rows, each with 90/3=30 and 300/3=100
        assert_eq!(m.len(), 3);
        for row in m.values() {
            assert_eq!(row[0], Some(30.0));
            assert_eq!(row[1], Some(100.0));
        }
    }

    // ========================================================================
    // Peptide-level functional matrix
    // ========================================================================

    #[test]
    fn test_build_functional_peptide_matrix_plurality() {
        let table = make_table(&[
            ("P1", &[("OG", "COG_A")]),
            ("P2", &[("OG", "COG_A")]),
            ("P3", &[("OG", "COG_B")]),
        ]);
        let group_keys = vec!["P1;P2".to_string(), "P3".to_string()];
        let res = resolve_level("OG", &table, &group_keys);

        // Two peptides under P1;P2, one under P3
        let pkeys: Vec<PeptideKey> = vec![
            ("P1;P2".to_string(), "PEPA".to_string()),
            ("P1;P2".to_string(), "PEPB".to_string()),
            ("P3".to_string(), "PEPC".to_string()),
        ];
        let mut pm: HashMap<PeptideKey, Vec<Option<f64>>> = HashMap::new();
        pm.insert(pkeys[0].clone(), vec![Some(10.0), Some(20.0)]);
        pm.insert(pkeys[1].clone(), vec![Some(5.0), None]);
        pm.insert(pkeys[2].clone(), vec![Some(100.0), Some(200.0)]);

        let (keys, m) =
            build_functional_peptide_matrix(&res, &pkeys, &pm, 2, AggregationMode::Plurality);

        // Expect 3 dual-index rows: (COG_A, PEPA), (COG_A, PEPB), (COG_B, PEPC)
        assert_eq!(keys.len(), 3);
        assert_eq!(
            m[&("COG_A".to_string(), "PEPA".to_string())],
            vec![Some(10.0), Some(20.0)]
        );
        assert_eq!(
            m[&("COG_A".to_string(), "PEPB".to_string())],
            vec![Some(5.0), None]
        );
        assert_eq!(
            m[&("COG_B".to_string(), "PEPC".to_string())],
            vec![Some(100.0), Some(200.0)]
        );
    }

    #[test]
    fn test_build_functional_peptide_matrix_splitintensity() {
        // A peptide under a merged group with 2 OGs → appears under BOTH OGs,
        // each at half intensity.
        let table = make_table(&[("P1", &[("OG", "COG_X")]), ("P2", &[("OG", "COG_Y")])]);
        let group_keys = vec!["P1;P2".to_string()];
        let res = resolve_level("OG", &table, &group_keys);

        let pkeys: Vec<PeptideKey> = vec![("P1;P2".to_string(), "PEPX".to_string())];
        let mut pm: HashMap<PeptideKey, Vec<Option<f64>>> = HashMap::new();
        pm.insert(pkeys[0].clone(), vec![Some(100.0)]);

        let (keys, m) =
            build_functional_peptide_matrix(&res, &pkeys, &pm, 1, AggregationMode::SplitIntensity);

        assert_eq!(keys.len(), 2);
        assert_eq!(
            m[&("COG_X".to_string(), "PEPX".to_string())],
            vec![Some(50.0)]
        );
        assert_eq!(
            m[&("COG_Y".to_string(), "PEPX".to_string())],
            vec![Some(50.0)]
        );
    }

    #[test]
    fn test_peptide_functional_matrix_unannotated_per_group() {
        // Two unannotated groups, each with a peptide — peptides must appear
        // under DISTINCT Unannotated_ rows, not merged.
        let table = make_table(&[]);
        let group_keys = vec!["GRP_ALPHA".to_string(), "GRP_BETA".to_string()];
        let res = resolve_level("OG", &table, &group_keys);

        let pkeys: Vec<PeptideKey> = vec![
            ("GRP_ALPHA".to_string(), "PEP1".to_string()),
            ("GRP_BETA".to_string(), "PEP2".to_string()),
        ];
        let mut pm: HashMap<PeptideKey, Vec<Option<f64>>> = HashMap::new();
        pm.insert(pkeys[0].clone(), vec![Some(10.0)]);
        pm.insert(pkeys[1].clone(), vec![Some(20.0)]);

        let (keys, _m) =
            build_functional_peptide_matrix(&res, &pkeys, &pm, 1, AggregationMode::Plurality);

        // Distinct rows per (unannotated-group, peptide)
        assert_eq!(keys.len(), 2);
        let functional_values: HashSet<&String> = keys.iter().map(|(f, _)| f).collect();
        assert_eq!(functional_values.len(), 2);
    }
}

#[cfg(test)]
mod release_review_tests {
    use super::*;

    fn observation(peptide: &str, intensity: f64, q: f64) -> PeptideObs {
        PeptideObs {
            peptide: peptide.into(),
            protein_group: "P1".into(),
            intensity,
            q_value: q,
            spectrum_q: q,
            peptide_q: q,
            protein_q: q,
        }
    }

    #[test]
    fn failing_observations_do_not_contribute_intensity() {
        let data = vec![SampleData {
            sample_name: "S".into(),
            n_peptides: 2,
            n_groups: 1,
            observations: vec![
                observation("PEPTIDEAA", 10.0, 0.001),
                observation("PEPTIDEBB", 1000.0, 0.9),
            ],
        }];
        let names = vec!["S".into()];
        let (_, groups) = build_protein_group_matrix(&data, &names, Some(0.01), "protein_q");
        assert_eq!(groups["P1"][0], Some(10.0));
        let (_, raw) = build_protein_group_matrix(&data, &names, None, "protein_q");
        assert_eq!(raw["P1"][0], Some(1010.0));
    }

    #[test]
    fn repeated_peptide_does_not_rescue_failing_psm_intensity() {
        let data = vec![SampleData {
            sample_name: "S".into(),
            n_peptides: 1,
            n_groups: 1,
            observations: vec![
                observation("PEPTIDEAA", 10.0, 0.001),
                observation("PEPTIDEAA", 1000.0, 0.9),
            ],
        }];
        let (_, matrix) = build_peptide_matrix(&data, &["S".into()], Some(0.01), "spectrum_q");
        assert_eq!(matrix[&("P1".into(), "PEPTIDEAA".into())][0], Some(10.0));
    }

    #[test]
    fn invalid_numeric_cells_cannot_pass_a_confidence_filter() {
        for value in ["-0.1", "-inf", "inf", "NaN", "2"] {
            assert!(parse_q(value) > 1.0);
        }
        for value in ["inf", "-inf", "NaN", "-1"] {
            assert_eq!(parse_intensity(value), 0.0);
        }
    }
    #[test]
    fn invalid_q_is_excluded_even_at_one_but_unfiltered_remains_explicit() {
        let data = vec![SampleData {
            sample_name: "S".into(),
            n_peptides: 2,
            n_groups: 1,
            observations: vec![
                observation("VALID", 10.0, parse_q("1.0")),
                observation("INVALID", 1000.0, parse_q("NaN")),
            ],
        }];
        let names = vec!["S".into()];
        let (_, matrix) = build_protein_group_matrix(&data, &names, Some(1.0), "protein_q");
        assert_eq!(matrix["P1"][0], Some(10.0));
        let (_, raw) = build_protein_group_matrix(&data, &names, None, "protein_q");
        assert_eq!(raw["P1"][0], Some(1010.0));
    }

    #[test]
    fn partial_search_or_lfq_rows_reject_the_sample() {
        let dir = std::env::temp_dir().join(format!("agg_partial_{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        let search =
            "peptide\tproteins\tlabel\tpeptide_q\tprotein_q\nPEPTIDEAK\tP1\t1\t0.001\t0.001\n";
        let lfq = "peptide\tcharge\tproteins\tq_value\tscore\tspectral_angle\tS.mzML\nPEPTIDEAK\t2\tP1\t0.5\t1\t0.9\t10\n";
        fs::write(
            dir.join("results.sage.tsv"),
            format!("{search}BROKEN\tP1\n"),
        )
        .unwrap();
        fs::write(dir.join("lfq.tsv"), lfq).unwrap();
        assert!(parse_lfq_sample("S", &dir, GroupKeyMode::MemberSet).is_none());
        fs::write(dir.join("results.sage.tsv"), search).unwrap();
        fs::write(dir.join("lfq.tsv"), format!("{lfq}BROKEN\t2\tP1\n")).unwrap();
        assert!(parse_lfq_sample("S", &dir, GroupKeyMode::MemberSet).is_none());
    }
}
