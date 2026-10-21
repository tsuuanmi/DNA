//! The "shortest dominant" reporting convention for a homopolymer length
//! heteroplasmy: an optional, labelled policy that turns a length
//! distribution into the length the lab records (ISFG, Parson et al. 2014,
//! 3.3.2: where two or more length variants are equally dominant, the
//! shortest dominant one is determined, as EMPOP does).
//!
//! It never decides a call: DNA's decided length stays the dominant one, and
//! the reading is a labelled suggestion beside it. The tract decoder stays
//! free of the convention (PROP-0012 non-goal).

/// The convention's name, recorded with every reading it gives.
pub const RULE: &str = "shortest_dominant";

/// Largest accepted ratio: a decode resolves two neighbouring lengths' share
/// ratio only to about 0.85-0.9, and forward reads under-read the shorter
/// neighbour (ratio 0.88-0.96 of the truth's), so a ratio nearer to 1 would
/// read noise (trace-simulator benchmark, 2026-10-21: 5,059 simulated reads
/// with known mixtures, noise fitted from development panel v5).
pub const MAX_RATIO: f64 = 0.85;

/// The policy: report the next-shorter length when its share is at least
/// `ratio` times the dominant's.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ShortestDominant {
    ratio: f64,
}

impl ShortestDominant {
    /// The policy with `ratio`, in `(0, MAX_RATIO]`; `None` outside it.
    #[must_use]
    pub fn new(ratio: f64) -> Option<Self> {
        (ratio > 0.0 && ratio <= MAX_RATIO).then_some(Self { ratio })
    }

    /// The ratio of the next-shorter length's share to the dominant's at
    /// which the shorter one is reported.
    #[must_use]
    pub const fn ratio(self) -> f64 {
        self.ratio
    }

    /// The length the convention records from one length distribution
    /// (`lengths` with their `shares`, in the same order): the dominant
    /// length, or the next-shorter one (dominant − 1) when its share is at
    /// least `ratio` times the dominant's, repeated while the shorter one
    /// keeps qualifying against the dominant. `None` without a positive
    /// share.
    #[must_use]
    pub fn reading(self, lengths: &[usize], shares: &[f64]) -> Option<usize> {
        let share = |length: usize| {
            lengths
                .iter()
                .position(|&kept| kept == length)
                .and_then(|index| shares.get(index))
                .copied()
                .unwrap_or(0.0)
        };
        let (dominant, top) = lengths
            .iter()
            .zip(shares)
            .filter(|(_, share)| **share > 0.0)
            .max_by(|left, right| left.1.total_cmp(right.1).then(right.0.cmp(left.0)))
            .map(|(&length, &share)| (length, share))?;
        let mut reading = dominant;
        while reading > 1 && share(reading - 1) >= self.ratio * top {
            reading -= 1;
        }
        Some(reading)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reads_the_dominant_or_the_qualifying_shorter_length() {
        let policy = ShortestDominant::new(0.6);
        assert!(ShortestDominant::new(0.0).is_none() && ShortestDominant::new(0.9).is_none());
        let Some(policy) = policy else {
            return;
        };
        assert!((policy.ratio() - 0.6).abs() < f64::EPSILON);
        // C9 dominant, C8 at 0.7 of it: the lab records C8.
        assert_eq!(policy.reading(&[8, 9, 10], &[0.35, 0.5, 0.15]), Some(8));
        // C8 below 0.6 of C9: the dominant.
        assert_eq!(policy.reading(&[8, 9, 10], &[0.25, 0.6, 0.15]), Some(9));
        // A longer neighbour never moves the reading.
        assert_eq!(policy.reading(&[9, 10], &[0.55, 0.45]), Some(9));
        // Equal shares: the shortest of the equally dominant lengths.
        assert_eq!(policy.reading(&[7, 8, 9], &[1.0 / 3.0; 3]), Some(7));
        // A gap stops the walk: C7 is not next to the dominant C9.
        assert_eq!(policy.reading(&[7, 9], &[0.45, 0.55]), Some(9));
        assert_eq!(policy.reading(&[], &[]), None);
    }
}
