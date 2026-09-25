use assert_cmd::Command;
use predicates::prelude::*;
use std::fs;
use tempfile::TempDir;

/// Helper: create a minimal AlphaNovo-format CSV with known peptides and scores.
fn write_test_csv(path: &std::path::Path, rows: &[(&str, f64)]) {
    let mut content = String::from("peptide_prediction_detokenized_unmodified,score\n");
    for (seq, score) in rows {
        content.push_str(&format!("{},{}\n", seq, score));
    }
    fs::write(path, content).expect("Failed to write test CSV");
}

/// Helper: create a minimal FASTA database with known protein sequences.
fn write_test_fasta(path: &std::path::Path, proteins: &[(&str, &str)]) {
    let mut content = String::new();
    for (id, seq) in proteins {
        content.push_str(&format!(">{}\n{}\n", id, seq));
    }
    fs::write(path, content).expect("Failed to write test FASTA");
}

/// Helper: count ">" headers in a FASTA file.
fn count_fasta_proteins(path: &std::path::Path) -> usize {
    let content = fs::read_to_string(path).unwrap_or_default();
    content.lines().filter(|l| l.starts_with('>')).count()
}

// `Command::cargo_bin` is deprecated in favor of the `cargo_bin_cmd!` macro, but the
// function still works and keeps this test compatible with the pinned assert_cmd 2.x.
#[allow(deprecated)]
fn binary() -> Command {
    Command::cargo_bin("fasta_extractor").expect("binary not found")
}

// ============================================================================
// Test 1: Single CSV file produces correct peptide matching
// ============================================================================
#[test]
fn test_single_csv_file_matches_proteins() {
    let tmp = TempDir::new().unwrap();

    // Peptides: ACDEFGHIK (9 aa), LMNPQRSTV (9 aa)
    let csv_path = tmp.path().join("sample_predictions.csv");
    write_test_csv(
        &csv_path,
        &[
            ("ACDEFGHIK", 0.95),
            ("LMNPQRSTV", 0.88),
            ("WYDKAMCFP", 0.72),
            ("SHORT", 0.99), // Too short (5 aa < 8 min)
        ],
    );

    // Proteins: one contains ACDEFGHIK, one contains LMNPQRSTV, one has no match
    let db_path = tmp.path().join("test.fasta");
    write_test_fasta(
        &db_path,
        &[
            ("prot1 Protein with first peptide", "MAAACDEFGHIKWWWWW"),
            ("prot2 Protein with second peptide", "MAALMNPQRSTVWWWWW"),
            ("prot3 No matching peptide", "MAAGGGGGGGGGGWWWWW"),
            ("prot4 Protein with third peptide", "MAAWYDKAMCFPWWWWW"),
        ],
    );

    let output_path = tmp.path().join("output.fasta");

    binary()
        .args([
            "-p",
            csv_path.to_str().unwrap(),
            "-d",
            db_path.to_str().unwrap(),
            "-o",
            output_path.to_str().unwrap(),
            "-s",
            "TEST",
            "--min-length",
            "8",
            "--max-length",
            "50",
            "--normalize-il",
            "false",
            "-t",
            "1",
        ])
        .assert()
        .success();

    // Should match 3 proteins (prot1, prot2, prot4), not prot3
    assert_eq!(count_fasta_proteins(&output_path), 3);

    // Verify output content
    let output = fs::read_to_string(&output_path).unwrap();
    assert!(output.contains("prot1"), "Should contain prot1");
    assert!(output.contains("prot2"), "Should contain prot2");
    assert!(output.contains("prot4"), "Should contain prot4");
    assert!(!output.contains("prot3"), "Should NOT contain prot3");
}

// ============================================================================
// Test 2: Directory mode produces identical results to single CSV
// ============================================================================
#[test]
fn test_directory_mode_matches_single_csv() {
    let tmp = TempDir::new().unwrap();

    let peptides = &[("ACDEFGHIK", 0.95), ("LMNPQRSTV", 0.88)];

    let proteins = &[
        ("prot1 Has peptide1", "MAAACDEFGHIKWWWWW"),
        ("prot2 Has peptide2", "MAALMNPQRSTVWWWWW"),
        ("prot3 No match", "MAAGGGGGGGGGGWWWWW"),
    ];

    let db_path = tmp.path().join("test.fasta");
    write_test_fasta(&db_path, proteins);

    // --- Single CSV mode ---
    let csv_path = tmp.path().join("sample_predictions.csv");
    write_test_csv(&csv_path, peptides);

    let output_single = tmp.path().join("output_single.fasta");
    binary()
        .args([
            "-p",
            csv_path.to_str().unwrap(),
            "-d",
            db_path.to_str().unwrap(),
            "-o",
            output_single.to_str().unwrap(),
            "-s",
            "TEST",
            "--min-length",
            "8",
            "--normalize-il",
            "false",
            "-t",
            "1",
        ])
        .assert()
        .success();

    // --- Directory mode ---
    let csv_dir = tmp.path().join("csv_dir");
    fs::create_dir(&csv_dir).unwrap();
    let csv_in_dir = csv_dir.join("sample_predictions.csv");
    write_test_csv(&csv_in_dir, peptides);

    let output_dir = tmp.path().join("output_dir.fasta");
    binary()
        .args([
            "-p",
            csv_dir.to_str().unwrap(),
            "-d",
            db_path.to_str().unwrap(),
            "-o",
            output_dir.to_str().unwrap(),
            "-s",
            "TEST",
            "--min-length",
            "8",
            "--normalize-il",
            "false",
            "-t",
            "1",
        ])
        .assert()
        .success();

    // Both modes should find the same proteins
    assert_eq!(
        count_fasta_proteins(&output_single),
        count_fasta_proteins(&output_dir),
        "Single CSV and directory mode should yield identical protein counts"
    );
    assert_eq!(count_fasta_proteins(&output_single), 2);
}

// ============================================================================
// Test 3: Plain text file with one peptide per line
// ============================================================================
#[test]
fn test_plain_text_peptide_file() {
    let tmp = TempDir::new().unwrap();

    let txt_path = tmp.path().join("peptides.txt");
    fs::write(&txt_path, "ACDEFGHIK\nLMNPQRSTV\nSHORT\n").unwrap();

    let db_path = tmp.path().join("test.fasta");
    write_test_fasta(
        &db_path,
        &[
            ("prot1 Match", "MAAACDEFGHIKWWWWW"),
            ("prot2 No match", "MAAGGGGGGGGGGWWWWW"),
        ],
    );

    let output_path = tmp.path().join("output.fasta");
    binary()
        .args([
            "-p",
            txt_path.to_str().unwrap(),
            "-d",
            db_path.to_str().unwrap(),
            "-o",
            output_path.to_str().unwrap(),
            "-s",
            "TEST",
            "--min-length",
            "8",
            "--normalize-il",
            "false",
            "-t",
            "1",
        ])
        .assert()
        .success();

    assert_eq!(count_fasta_proteins(&output_path), 1);
    let output = fs::read_to_string(&output_path).unwrap();
    assert!(output.contains("prot1"));
}

// ============================================================================
// Test 4: Score filtering with top-percent threshold
// ============================================================================
#[test]
fn test_score_filtering_top_percent() {
    let tmp = TempDir::new().unwrap();

    // 10 peptides with scores from 0.1 to 1.0
    // Top 30% = top 3 peptides: scores 1.0, 0.9, 0.8
    let csv_path = tmp.path().join("predictions.csv");
    write_test_csv(
        &csv_path,
        &[
            ("PEPTIDEAA", 1.0),
            ("PEPTIDEBB", 0.9),
            ("PEPTIDECC", 0.8),
            ("PEPTIDEDD", 0.7),
            ("PEPTIDEEE", 0.6),
            ("PEPTIDEFFA", 0.5),
            ("PEPTIDEGGA", 0.4),
            ("PEPTIDEHHA", 0.3),
            ("PEPTIDEIIA", 0.2),
            ("PEPTIDEJJA", 0.1),
        ],
    );

    // Database: one protein contains the top-scoring peptide, one contains low-scoring
    let db_path = tmp.path().join("test.fasta");
    write_test_fasta(
        &db_path,
        &[
            ("prot_top Top scoring match", "MAAPEPTIDEAAWWWWW"),
            ("prot_low Low scoring match", "MAAPEPTIDEJJAWWWWW"),
            ("prot_mid Mid scoring match", "MAAPEPTIDEDDWWWWW"),
        ],
    );

    let output_path = tmp.path().join("output.fasta");
    binary()
        .args([
            "-p",
            csv_path.to_str().unwrap(),
            "-d",
            db_path.to_str().unwrap(),
            "-o",
            output_path.to_str().unwrap(),
            "-s",
            "TEST",
            "--min-length",
            "8",
            "--max-length",
            "50",
            "--top-percent",
            "0.30",
            "--normalize-il",
            "false",
            "-t",
            "1",
        ])
        .assert()
        .success();

    // Only the top-scoring peptide's protein should match
    let output = fs::read_to_string(&output_path).unwrap();
    assert!(
        output.contains("prot_top"),
        "Top-scoring peptide should match"
    );
    assert!(
        !output.contains("prot_low"),
        "Low-scoring peptide should be filtered out"
    );
    // prot_mid has score 0.7 which is not in top 30% (threshold is 0.8)
    assert!(
        !output.contains("prot_mid"),
        "Mid-scoring peptide should be filtered out"
    );
}

// ============================================================================
// Test 5: I/L normalization matches I↔L variants
// ============================================================================
#[test]
fn test_il_normalization_proper() {
    let tmp = TempDir::new().unwrap();

    // Peptide: ACDEFGHIK (contains I)
    let csv_path = tmp.path().join("predictions.csv");
    write_test_csv(&csv_path, &[("ACDEFGHIK", 0.95)]);

    // Protein: same but with L instead of I: ACDEFGHLK
    let db_path = tmp.path().join("test.fasta");
    write_test_fasta(
        &db_path,
        &[
            ("prot_L L-variant protein", "MAAACDEFGHLKWWWWW"),
            ("prot_I I-variant protein", "MAAACDEFGHIKWWWWW"),
            ("prot_none No match", "MAAGGGGGGGGGGWWWWW"),
        ],
    );

    // With normalization ON: I→L in both peptide and protein, so ACDEFGHLK matches
    let output_on = tmp.path().join("output_on.fasta");
    binary()
        .args([
            "-p",
            csv_path.to_str().unwrap(),
            "-d",
            db_path.to_str().unwrap(),
            "-o",
            output_on.to_str().unwrap(),
            "-s",
            "TEST",
            "--min-length",
            "8",
            "--normalize-il",
            "true",
            "-t",
            "1",
        ])
        .assert()
        .success();

    let output = fs::read_to_string(&output_on).unwrap();
    assert!(
        output.contains("prot_L"),
        "L-variant should match with normalization ON"
    );
    assert!(
        output.contains("prot_I"),
        "I-variant should match with normalization ON"
    );
    assert_eq!(count_fasta_proteins(&output_on), 2);

    // With normalization OFF: only exact I match
    let output_off = tmp.path().join("output_off.fasta");
    binary()
        .args([
            "-p",
            csv_path.to_str().unwrap(),
            "-d",
            db_path.to_str().unwrap(),
            "-o",
            output_off.to_str().unwrap(),
            "-s",
            "TEST",
            "--min-length",
            "8",
            "--normalize-il",
            "false",
            "-t",
            "1",
        ])
        .assert()
        .success();

    let output = fs::read_to_string(&output_off).unwrap();
    assert!(
        output.contains("prot_I"),
        "I-variant should match without normalization"
    );
    assert!(
        !output.contains("prot_L"),
        "L-variant should NOT match without normalization"
    );
    assert_eq!(count_fasta_proteins(&output_off), 1);
}

// ============================================================================
// Test 7: Empty input produces error
// ============================================================================
#[test]
fn test_empty_csv_fails() {
    let tmp = TempDir::new().unwrap();

    let csv_path = tmp.path().join("empty.csv");
    fs::write(
        &csv_path,
        "peptide_prediction_detokenized_unmodified,score\n",
    )
    .unwrap();

    let db_path = tmp.path().join("test.fasta");
    write_test_fasta(&db_path, &[("prot1 Dummy", "MAAACDEFGHIKWWWWW")]);

    let output_path = tmp.path().join("output.fasta");
    binary()
        .args([
            "-p",
            csv_path.to_str().unwrap(),
            "-d",
            db_path.to_str().unwrap(),
            "-o",
            output_path.to_str().unwrap(),
            "-s",
            "TEST",
            "--normalize-il",
            "false",
            "-t",
            "1",
        ])
        .assert()
        .failure()
        .stderr(predicate::str::contains("No peptides loaded"));
}

// ============================================================================
// Test 8: Min/max length filtering works correctly
// ============================================================================
#[test]
fn test_length_filtering() {
    let tmp = TempDir::new().unwrap();

    let csv_path = tmp.path().join("predictions.csv");
    write_test_csv(
        &csv_path,
        &[
            ("ACDEFG", 0.95),                                               // 6 aa - below min 8
            ("ACDEFGHIK", 0.95),                                            // 9 aa - within range
            ("ACDEFGHIKLMNPQRST", 0.95),                                    // 17 aa - within range
            ("ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWYACDEFGHIKLMN", 0.95), // 51 aa - above max 50
        ],
    );

    let db_path = tmp.path().join("test.fasta");
    write_test_fasta(
        &db_path,
        &[
            ("prot_short Contains 6aa", "MAAACDEFGWWWWWWWW"),
            ("prot_ok Contains 9aa", "MAAACDEFGHIKWWWWW"),
            ("prot_long Contains 17aa", "MAAACDEFGHIKLMNPQRSTWWWWW"),
        ],
    );

    let output_path = tmp.path().join("output.fasta");
    binary()
        .args([
            "-p",
            csv_path.to_str().unwrap(),
            "-d",
            db_path.to_str().unwrap(),
            "-o",
            output_path.to_str().unwrap(),
            "-s",
            "TEST",
            "--min-length",
            "8",
            "--max-length",
            "50",
            "--normalize-il",
            "false",
            "-t",
            "1",
        ])
        .assert()
        .success();

    let output = fs::read_to_string(&output_path).unwrap();
    assert!(
        !output.contains("prot_short"),
        "6aa peptide should be filtered by min-length"
    );
    assert!(output.contains("prot_ok"), "9aa peptide should pass");
    assert!(output.contains("prot_long"), "17aa peptide should pass");
}

// ============================================================================
// Test 9: Stats JSON is produced alongside FASTA output
// ============================================================================
#[test]
fn test_stats_json_output() {
    let tmp = TempDir::new().unwrap();

    let csv_path = tmp.path().join("predictions.csv");
    write_test_csv(&csv_path, &[("ACDEFGHIK", 0.95)]);

    let db_path = tmp.path().join("test.fasta");
    write_test_fasta(
        &db_path,
        &[
            ("prot1 Match", "MAAACDEFGHIKWWWWW"),
            ("prot2 No match", "MAAGGGGGGGGGGWWWWW"),
        ],
    );

    let output_path = tmp.path().join("output.fasta");
    binary()
        .args([
            "-p",
            csv_path.to_str().unwrap(),
            "-d",
            db_path.to_str().unwrap(),
            "-o",
            output_path.to_str().unwrap(),
            "-s",
            "TEST",
            "--normalize-il",
            "false",
            "-t",
            "1",
        ])
        .assert()
        .success();

    // Check stats JSON exists and has correct values
    let stats_path = tmp.path().join("output.stats.json");
    assert!(stats_path.exists(), "Stats JSON should be created");

    let stats: serde_json::Value =
        serde_json::from_str(&fs::read_to_string(&stats_path).unwrap()).unwrap();

    assert_eq!(stats["proteins_matched"], 1);
    assert_eq!(stats["proteins_searched"], 2);
    assert_eq!(stats["total_peptides"], 1);
    assert_eq!(stats["source_tag"], "TEST");
}

// ============================================================================
// Test 10: Hits TSV output
// ============================================================================
#[test]
fn test_hits_tsv_output() {
    let tmp = TempDir::new().unwrap();

    let csv_path = tmp.path().join("predictions.csv");
    write_test_csv(&csv_path, &[("ACDEFGHIK", 0.95), ("LMNPQRSTV", 0.88)]);

    let db_path = tmp.path().join("test.fasta");
    // prot1 contains both peptides
    write_test_fasta(&db_path, &[("prot1 Both peptides", "ACDEFGHIKLMNPQRSTV")]);

    let output_path = tmp.path().join("output.fasta");
    let hits_path = tmp.path().join("hits.tsv");
    binary()
        .args([
            "-p",
            csv_path.to_str().unwrap(),
            "-d",
            db_path.to_str().unwrap(),
            "-o",
            output_path.to_str().unwrap(),
            "-s",
            "TEST",
            "--normalize-il",
            "false",
            "--output-hits",
            hits_path.to_str().unwrap(),
            "-t",
            "1",
        ])
        .assert()
        .success();

    let hits = fs::read_to_string(&hits_path).unwrap();
    assert!(
        hits.contains("peptide\tprotein_id\tposition"),
        "TSV should have header"
    );
    assert!(hits.contains("ACDEFGHIK"), "Should list first peptide hit");
    assert!(hits.contains("LMNPQRSTV"), "Should list second peptide hit");
    // Two peptide hits for one protein = 2 data lines + 1 header = 3 lines
    assert_eq!(hits.lines().count(), 3);
}
