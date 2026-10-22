//! The A/C/G/T channels of a call ranked by their signal: the one order
//! basecalling, signal processing, and callability rank channels by.

use std::cmp::Ordering;

/// The four channel indices by decreasing value under `order`, a tie going
/// to the lower index.
pub(crate) fn ranked<T>(values: &[T; 4], order: impl Fn(&T, &T) -> Ordering) -> [usize; 4] {
    let mut ranked = [0_usize, 1, 2, 3];
    ranked.sort_by(|left, right| {
        order(&values[*right], &values[*left]).then_with(|| left.cmp(right))
    });
    ranked
}

#[cfg(test)]
mod tests {
    use super::ranked;

    #[test]
    fn ranks_channels_by_decreasing_value_ties_to_the_lower_index() {
        assert_eq!(ranked(&[1, 5, 5, 0], Ord::cmp), [1, 2, 0, 3]);
        assert_eq!(ranked(&[0.5, 0.0, 2.0, 0.5], f64::total_cmp), [2, 0, 3, 1]);
    }
}
