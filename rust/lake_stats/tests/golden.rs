//! Golden-file integration test for lake_stats.
//!
//! Builds a tiny FASTA with hand-enumerated tryptic peptides, runs the binary,
//! and asserts the output TSV matches expected D (unique peptide count) per length.
//!
//! Tryptic rule: cut after K/R unless followed by P. Peptides containing any of
//! B,J,O,U,X,Z,*,- are dropped (handled by `is_clean`).

use std::fs;
use std::process::Command;

/// Two synthetic proteins with hand-enumerated tryptic cuts.
///
/// p1: MAKPEPTIDERWAYR
///   cuts after K@3 (next char P → DO NOT CUT — "K|P" is suppressed!)
///   wait: rule is "cut after K/R unless followed by P". K@3 is followed by P → NO CUT.
///   So protein becomes one piece up to next K/R that is not followed by P.
///   M-A-K-P-E-P-T-I-D-E-R-W-A-Y-R
///       K3 (next P)  → no cut
///                   R10 (next W)  → CUT after R10 → "MAKPEPTIDER" (len 11)
///                            R14 (last)   → emit "WAYR" (len 4)
///   So p1 → ["MAKPEPTIDER" (11), "WAYR" (4)]
///
/// p2: GLAYERSEEK
///   G-L-A-Y-E-R-S-E-E-K
///       R5 (next S) → CUT after R5 → "GLAYER" (len 6)
///                S-E-E-K9 → emit "SEEK" (len 4)
///   So p2 → ["GLAYER" (6), "SEEK" (4)]
///
/// Combined unique peptides per length:
///   L=4:  WAYR, SEEK         → D=2
///   L=6:  GLAYER             → D=1
///   L=11: MAKPEPTIDER        → D=1
///   L=5,7,8,9,10:            → D=0
const FASTA: &str = ">p1\nMAKPEPTIDERWAYR\n>p2\nGLAYERSEEK\n";

#[test]
fn golden_two_protein_fasta() {
    let tmp = tempdir();
    let fasta_path = tmp.join("test.fasta");
    let out_path = tmp.join("out.tsv");

    fs::write(&fasta_path, FASTA).expect("write fasta");

    let bin = env!("CARGO_BIN_EXE_lake_stats");
    let status = Command::new(bin)
        .args([
            "--fasta",
            fasta_path.to_str().unwrap(),
            "--output",
            out_path.to_str().unwrap(),
            "--lengths",
            "4,5,6,7,8,9,10,11",
            "--name",
            "golden_test",
            "--threads",
            "1",
        ])
        .status()
        .expect("run lake_stats");
    assert!(status.success(), "lake_stats binary failed");

    let out = fs::read_to_string(&out_path).expect("read output");
    let mut by_len: std::collections::HashMap<usize, u64> = std::collections::HashMap::new();
    for line in out.lines().skip(1) {
        let cols: Vec<&str> = line.split('\t').collect();
        let len: usize = cols[1].parse().unwrap();
        let d: u64 = cols[3].parse().unwrap();
        by_len.insert(len, d);
    }
    assert_eq!(
        by_len.get(&4).copied(),
        Some(2),
        "L=4 should be {{WAYR, SEEK}}"
    );
    assert_eq!(by_len.get(&5).copied(), Some(0));
    assert_eq!(by_len.get(&6).copied(), Some(1), "L=6 should be {{GLAYER}}");
    assert_eq!(by_len.get(&7).copied(), Some(0));
    assert_eq!(by_len.get(&8).copied(), Some(0));
    assert_eq!(by_len.get(&9).copied(), Some(0));
    assert_eq!(by_len.get(&10).copied(), Some(0));
    assert_eq!(
        by_len.get(&11).copied(),
        Some(1),
        "L=11 should be {{MAKPEPTIDER}}"
    );

    // n_proteins must be 2
    let nprot: u64 = out
        .lines()
        .nth(1)
        .unwrap()
        .split('\t')
        .nth(5)
        .unwrap()
        .parse()
        .unwrap();
    assert_eq!(nprot, 2, "n_proteins should be 2");
}

#[test]
fn golden_drops_ambiguous_residues() {
    // Protein with X must be skipped; others retained.
    let tmp = tempdir();
    let fasta_path = tmp.join("amb.fasta");
    let out_path = tmp.join("amb.tsv");
    // q1: GLXAYER → tryptic at R6 (next end-of-prot) → "GLXAYER" (len 7) — has X, must be dropped
    // q2: GLAYER  → "GLAYER" (len 6) — clean
    fs::write(&fasta_path, ">q1\nGLXAYER\n>q2\nGLAYER\n").expect("write");

    let bin = env!("CARGO_BIN_EXE_lake_stats");
    let status = Command::new(bin)
        .args([
            "--fasta",
            fasta_path.to_str().unwrap(),
            "--output",
            out_path.to_str().unwrap(),
            "--lengths",
            "6,7",
            "--threads",
            "1",
        ])
        .status()
        .expect("run");
    assert!(status.success());

    let out = fs::read_to_string(&out_path).unwrap();
    let mut by_len: std::collections::HashMap<usize, u64> = std::collections::HashMap::new();
    for line in out.lines().skip(1) {
        let cols: Vec<&str> = line.split('\t').collect();
        by_len.insert(cols[1].parse().unwrap(), cols[3].parse().unwrap());
    }
    assert_eq!(by_len.get(&6).copied(), Some(1), "L=6 GLAYER kept");
    assert_eq!(
        by_len.get(&7).copied(),
        Some(0),
        "L=7 GLXAYER dropped (contains X)"
    );
}

#[test]
fn golden_drops_peptides_spanning_stop_or_gap() {
    // A peptide must NOT be concatenated across a `*` (stop codon) or `-` (gap):
    // it must reach `is_clean` (which drops it) rather than being stripped at
    // read time. Regression test for the reader that used to keep only
    // is_ascii_alphabetic bytes.
    let tmp = tempdir();
    let fasta_path = tmp.join("stopgap.fasta");
    let out_path = tmp.join("stopgap.tsv");
    // s1: GLAY*ER  -> single tryptic peptide "GLAY*ER" (len 7) — contains `*`, must be DROPPED.
    //     If `*` were stripped, it would become "GLAYER" (len 6) and be wrongly counted.
    // s2: GLAY-ER  -> single tryptic peptide "GLAY-ER" (len 7) — contains `-`, must be DROPPED.
    // s3: GLAYER   -> clean "GLAYER" (len 6), the only survivor.
    fs::write(&fasta_path, ">s1\nGLAY*ER\n>s2\nGLAY-ER\n>s3\nGLAYER\n").expect("write");

    let bin = env!("CARGO_BIN_EXE_lake_stats");
    let status = Command::new(bin)
        .args([
            "--fasta",
            fasta_path.to_str().unwrap(),
            "--output",
            out_path.to_str().unwrap(),
            "--lengths",
            "6,7",
            "--threads",
            "1",
        ])
        .status()
        .expect("run");
    assert!(status.success());

    let out = fs::read_to_string(&out_path).unwrap();
    let mut by_len: std::collections::HashMap<usize, u64> = std::collections::HashMap::new();
    for line in out.lines().skip(1) {
        let cols: Vec<&str> = line.split('\t').collect();
        by_len.insert(cols[1].parse().unwrap(), cols[3].parse().unwrap());
    }
    // Only the clean GLAYER (len 6) survives; both starred/gapped peptides are dropped.
    assert_eq!(
        by_len.get(&6).copied(),
        Some(1),
        "L=6 only GLAYER (from s3) — NOT the stripped GLAY*ER/GLAY-ER"
    );
    assert_eq!(
        by_len.get(&7).copied(),
        Some(0),
        "L=7 GLAY*ER and GLAY-ER both dropped (contain */-)"
    );
}

/// Minimal tempdir helper to avoid pulling in a tempfile dep.
fn tempdir() -> std::path::PathBuf {
    let pid = std::process::id();
    let nanos = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let p = std::env::temp_dir().join(format!("lake_stats_test_{pid}_{nanos}"));
    fs::create_dir_all(&p).expect("mkdir tempdir");
    p
}
