//! FastaLake input checks for DIANovo; this is not model inference.
//! Upstream: https://github.com/hearthewind/dianovo at cc413f25140e47a5f6a4234bae15da0ef525ac3a.
//! See denovo_from_feature_detection/rnova/task/{sequence_generation_inference_utils,
//! optimal_path_inference_classic}.py and docs/how-to/de-novo-inputs.md.
//! Keep semantics aligned with fasta_lake/dianovo.py; cross-language tests compare
//! extraction and razor selection with independently scored generic predictions.

/// Parse a complete canonical peptide and its mean confidence, retaining zeros.
/// The space-delimited path vector alone may have one extra leading zero marker.
/// Never erase numeric mass gaps or silently discard malformed probabilities.
pub fn parse_prediction(sequence: &str, probabilities: &str) -> Option<(String, f64)> {
    let peptide: String = sequence
        .chars()
        .filter(|c| !c.is_ascii_whitespace())
        .map(|c| c.to_ascii_uppercase())
        .collect();
    if peptide.is_empty()
        || !peptide
            .bytes()
            .all(|aa| b"ACDEFGHIKLMNPQRSTVWY".contains(&aa))
    {
        return None;
    }
    let raw = probabilities.trim();
    let bracketed = raw.starts_with('[') && raw.ends_with(']');
    let raw = if bracketed {
        &raw[1..raw.len() - 1]
    } else {
        raw
    };
    let comma_format = bracketed || raw.contains(',');
    let tokens: Vec<&str> = if comma_format {
        raw.split(',').collect()
    } else {
        raw.split_whitespace().collect()
    };
    let values: Vec<f64> = tokens
        .iter()
        .map(|token| token.trim().parse::<f64>().ok())
        .collect::<Option<_>>()?;
    if values.is_empty()
        || values
            .iter()
            .any(|p| !p.is_finite() || !(0.0..=1.0).contains(p))
    {
        return None;
    }
    let residues = if !comma_format && values.len() == peptide.len() + 1 && values[0] == 0.0 {
        &values[1..]
    } else {
        &values[..]
    };
    if residues.len() != peptide.len() {
        return None;
    }
    let score = residues.iter().sum::<f64>() / residues.len() as f64;
    Some((peptide, score))
}

#[cfg(test)]
mod tests {
    use super::parse_prediction;

    #[test]
    fn zeros_and_path_markers_have_different_meanings() {
        assert_eq!(
            parse_prediction("AcD", "[0, 0.9, 0.9]"),
            Some(("ACD".into(), 0.6))
        );
        assert_eq!(
            parse_prediction("A c D", "0 0 0.9 0.9"),
            Some(("ACD".into(), 0.6))
        );
        assert_eq!(
            parse_prediction("AcD", "0 0.9 0.9"),
            Some(("ACD".into(), 0.6))
        );
        assert_eq!(
            parse_prediction("AcD", "[0, 0, 0]"),
            Some(("ACD".into(), 0.0))
        );
    }

    #[test]
    fn no_gap_joining_or_partial_vector_recovery() {
        for (seq, prob) in [
            ("AC113.084D", "[0.9, 0.9, 0.9]"),
            ("AXB", "[1, 1, 1]"),
            ("ACD", "[0, 1, 1, 1]"),
            ("ACD", "0.1 1 1 1"),
            ("ACD", "[1, bogus, 1]"),
            ("ACD", "[1, NaN, 1]"),
            ("ACD", "[1, inf, 1]"),
            ("ACD", "[1, -0.1, 1]"),
            ("ACD", "[1, 1.1, 1]"),
            ("ACD", "[1, 1, 1"),
            ("ACD", "[1, 1, 1,]"),
            ("ACD", "[]"),
        ] {
            assert_eq!(parse_prediction(seq, prob), None, "{seq}: {prob}");
        }
    }
}
