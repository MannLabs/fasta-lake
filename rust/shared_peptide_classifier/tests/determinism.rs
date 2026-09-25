//! Determinism regression test (D5).
//!
//! Runs the binary twice on the same fixture with different thread counts and
//! asserts byte-identical output. Before the D5 fix, row order came from
//! rayon thread scheduling and could differ between runs / thread counts.

use std::fs;
use std::process::Command;

const LAKE: &str = ">member1\nAAACCCDDDEEEPEPTIDEAAA\nGGGFFFFKKKKLLLLM\n\
>member2\nAAACCCDDDEEEXXXAAA\nGGGFFFFKKKKLLLLM\n\
>member3\nAAACCCDDDEEEYYYAAA\nGGGFFFFKKKKLLLLM\n\
>member4\nNNNNHHHHWWWWQQQQRARERAREPEPSEQK\n\
>member5\nZZZZZZZZZZ\n";

const CLUSTER: &str = "clusterA\tmember1\nclusterA\tmember2\nclusterA\tmember3\n\
clusterB\tmember4\nclusterB\tmember5\n";

// Many pairs across two reps so parallel scheduling has room to reorder rows.
fn rescued() -> String {
    let mut s = String::from("peptide\tcluster_rep\tscore\n");
    for i in 0..200 {
        let (pep, rep) = match i % 4 {
            0 => ("GGGFFFFK", "clusterA"),
            1 => ("PEPTIDE", "clusterA"),
            2 => ("PEPSEQK", "clusterB"),
            _ => ("GHOSTPEP", "clusterC"),
        };
        s.push_str(&format!("{pep}\t{rep}\t0.9\n"));
    }
    s
}

fn run(dir: &std::path::Path, threads: &str, out_name: &str) -> String {
    let out = dir.join(out_name);
    let status = Command::new(env!("CARGO_BIN_EXE_shared_peptide_classifier"))
        .args([
            "--rescued",
            dir.join("rescued.tsv").to_str().unwrap(),
            "--cluster-tsv",
            dir.join("cluster.tsv").to_str().unwrap(),
            "--lake",
            dir.join("lake.fasta").to_str().unwrap(),
            "--output",
            out.to_str().unwrap(),
            "--shared-threshold",
            "0.8",
            "--threads",
            threads,
        ])
        .status()
        .expect("failed to run binary");
    assert!(status.success(), "binary exited nonfatally");
    fs::read_to_string(out).unwrap()
}

#[test]
fn output_is_deterministic_across_thread_counts() {
    let dir = std::env::temp_dir().join(format!("spc_det_{}", std::process::id()));
    fs::create_dir_all(&dir).unwrap();
    fs::write(dir.join("lake.fasta"), LAKE).unwrap();
    fs::write(dir.join("cluster.tsv"), CLUSTER).unwrap();
    fs::write(dir.join("rescued.tsv"), rescued()).unwrap();

    let a = run(&dir, "1", "out1.tsv");
    let b = run(&dir, "8", "out8.tsv");
    let c = run(&dir, "4", "out4.tsv");

    let _ = &dir; // Keep test artifacts for release review.

    assert_eq!(a, b, "output differs between 1 and 8 threads");
    assert_eq!(a, c, "output differs between 1 and 4 threads");

    // Sanity: header + 150 non-ghost rows (GHOSTPEP/clusterC dropped).
    assert_eq!(a.lines().count(), 1 + 150);
    // Rows are sorted by (peptide, cluster_rep): first data row is GGGFFFFK.
    assert!(a
        .lines()
        .nth(1)
        .unwrap()
        .starts_with("GGGFFFFK\tclusterA\t"));
}
