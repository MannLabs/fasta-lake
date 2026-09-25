//! `fasta_export entrapment` must list samples in a filesystem-independent order.
//!
//! Eight sample directories are created in reverse name order; the per-sample
//! table on stdout must come out in lexicographic order regardless of how the
//! filesystem enumerates them, and the run must be byte-identical when repeated.

use std::fs;
use std::process::Command;

const HEADER: &str = "peptide\tproteins\tspectrum_q\tlabel\n";

fn write_sample(root: &std::path::Path, name: &str, n_target: usize, n_entrap: usize) {
    let dir = root.join(name);
    fs::create_dir_all(&dir).unwrap();
    let mut tsv = String::from(HEADER);
    for i in 0..n_target {
        tsv.push_str(&format!("TARGET{i}K\tprot{i}\t0.001\t1\n"));
    }
    for i in 0..n_entrap {
        tsv.push_str(&format!("ENTRAP{i}K\tENTRAP_SHUF_{i}\t0.001\t1\n"));
    }
    tsv.push_str("DECOYK\trev_prot0\t0.001\t-1\n");
    fs::write(dir.join("results.sage.tsv"), tsv).unwrap();
}

fn run(root: &std::path::Path) -> String {
    let out = Command::new(env!("CARGO_BIN_EXE_fasta_export"))
        .args([
            "entrapment",
            "--input",
            root.to_str().unwrap(),
            "--prefix",
            "ENTRAP_SHUF_",
            "--target-db-size",
            "100",
            "--entrap-db-size",
            "10",
        ])
        .output()
        .expect("binary runs");
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8(out.stdout).unwrap()
}

#[test]
fn samples_are_listed_in_sorted_order_and_runs_are_identical() {
    let tmp = tempfile::tempdir().unwrap();
    let root = tmp.path();
    // Created in reverse order with unequal hit counts, so a wrong order changes
    // both the table and the float sum behind the summary mean.
    for i in (0..8).rev() {
        write_sample(root, &format!("s{i}"), 20 + i, i);
    }
    let first = run(root);
    let second = run(root);
    assert_eq!(
        first, second,
        "two runs on the same directory must be byte-identical"
    );

    let samples: Vec<&str> = first
        .lines()
        .skip(1)
        .map(|l| l.split('\t').next().unwrap())
        .collect();
    let mut sorted = samples.clone();
    sorted.sort();
    assert_eq!(
        samples, sorted,
        "table order must not depend on read_dir order"
    );
    assert_eq!(samples, ["s0", "s1", "s2", "s3", "s4", "s5", "s6", "s7"]);
}
