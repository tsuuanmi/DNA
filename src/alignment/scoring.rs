//! Fixed-point profile substitution scoring and deterministic state ordering.

use crate::config::AlignmentConfig;
use crate::model::locus_evidence::EvidenceProfile;

pub(crate) const NEGATIVE_INFINITY: i64 = i64::MIN / 4;

/// Fixed-point denominator used by all alignment score deltas.
///
/// Multiplying every configured score by the same scale preserves the existing
/// clean one-hot alignment ordering while allowing mixed evidence to contribute
/// fractional expected match/mismatch support without floating-point DP state.
pub(crate) const SCORE_SCALE: i64 = 1024;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum State {
    Match = 0,
    Insertion = 1,
    Deletion = 2,
}

impl State {
    pub(crate) const fn from_bits(bits: u8) -> Self {
        match bits {
            1 => Self::Insertion,
            2 => Self::Deletion,
            _ => Self::Match,
        }
    }
}

/// Returns a fixed-point expected match/mismatch score for one evidence profile.
///
/// Canonical reference bases use the profile mass assigned to that base.
/// Missing profiles or non-canonical reference symbols receive the configured
/// ambiguous score. The DP itself never compares floating-point values.
pub(crate) fn substitution(
    profile: Option<EvidenceProfile>,
    reference: u8,
    config: &AlignmentConfig,
) -> i64 {
    let Some(reference_index) = canonical_index(reference) else {
        return scaled(config.ambiguous_score);
    };
    let Some(profile) = profile else {
        return scaled(config.ambiguous_score);
    };
    let support = quantize_weight(profile.weights[reference_index]);
    support * i64::from(config.match_score)
        + (SCORE_SCALE - support) * i64::from(config.mismatch_score)
}

pub(crate) const fn scaled(delta: i32) -> i64 {
    delta as i64 * SCORE_SCALE
}

pub(crate) const fn is_canonical(base: u8) -> bool {
    matches!(base, b'A' | b'C' | b'G' | b'T')
}

pub(crate) fn add(score: i64, delta: i64) -> i64 {
    if score <= NEGATIVE_INFINITY / 2 {
        NEGATIVE_INFINITY
    } else {
        score + delta
    }
}

const fn canonical_index(base: u8) -> Option<usize> {
    match base {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

/// Quantizes a validated normalized weight to 1/1024 units.
///
/// Rust round semantics use ties away from zero. Evidence weights are
/// non-negative, so exact half-unit ties round upward.
fn quantize_weight(weight: f64) -> i64 {
    (weight * SCORE_SCALE as f64)
        .round()
        .clamp(0.0, SCORE_SCALE as f64) as i64
}

#[cfg(test)]
mod tests {
    use super::*;

    fn config() -> AlignmentConfig {
        AlignmentConfig {
            match_score: 3,
            mismatch_score: -5,
            ambiguous_score: 0,
            gap_open_score: -10,
            gap_extension_score: -4,
            minimum_callable_bases: 1,
            minimum_identity: 0.8,
        }
    }

    #[test]
    fn one_hot_profiles_preserve_clean_match_and_mismatch_scores() {
        let profile = EvidenceProfile {
            weights: [1.0, 0.0, 0.0, 0.0],
        };
        assert_eq!(
            substitution(Some(profile), b'A', &config()),
            scaled(config().match_score)
        );
        assert_eq!(
            substitution(Some(profile), b'C', &config()),
            scaled(config().mismatch_score)
        );
    }

    #[test]
    fn mixed_profile_uses_quantized_expected_score() {
        let profile = EvidenceProfile {
            weights: [0.75, 0.25, 0.0, 0.0],
        };
        assert_eq!(
            substitution(Some(profile), b'A', &config()),
            768 * 3 + 256 * -5
        );
        assert_eq!(
            substitution(Some(profile), b'C', &config()),
            256 * 3 + 768 * -5
        );
    }

    #[test]
    fn missing_profile_or_ambiguous_reference_uses_ambiguous_score() {
        assert_eq!(substitution(None, b'A', &config()), 0);
        let profile = EvidenceProfile {
            weights: [1.0, 0.0, 0.0, 0.0],
        };
        assert_eq!(substitution(Some(profile), b'N', &config()), 0);
    }

    #[test]
    fn half_unit_quantization_rounds_up() {
        assert_eq!(quantize_weight(0.5 / SCORE_SCALE as f64), 1);
    }
}
