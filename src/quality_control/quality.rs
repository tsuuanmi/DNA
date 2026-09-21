//! Deterministic relative quality-score conversion.

/// Converts non-negative penalties into an uncalibrated bounded score.
pub(crate) fn relative_scores(penalties: &[i32], maximum: u8) -> Vec<u8> {
    let max_penalty = penalties.iter().copied().max().unwrap_or(0);
    if max_penalty <= 0 {
        return vec![maximum; penalties.len()];
    }
    penalties
        .iter()
        .map(|penalty| {
            let fraction = 1.0 - f64::from(*penalty) / f64::from(max_penalty);
            (f64::from(maximum) * fraction.clamp(0.0, 1.0)).floor() as u8
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn zero_penalty_receives_maximum() {
        assert_eq!(relative_scores(&[0, 0], 60), vec![60, 60]);
    }

    #[test]
    fn worst_penalty_receives_zero() {
        assert_eq!(relative_scores(&[0, 5], 60), vec![60, 0]);
    }
}
