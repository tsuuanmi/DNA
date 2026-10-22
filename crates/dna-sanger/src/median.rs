//! The medians the Sanger stages take, one function per definition: the
//! midpoint median, the upper median, and the lower median of a histogram.
//! Each stage keeps the definition it was specified with; unifying them
//! would change outputs.

use std::cmp::Ordering;

/// The midpoint median of `values` under `f64::total_cmp`: the middle value,
/// or the midpoint of the two middle values for an even count; `None` for
/// none. Found by selection, not a full sort; reorders `values`.
pub(crate) fn midpoint(values: &mut [f64]) -> Option<f64> {
    if values.is_empty() {
        return None;
    }
    let even = values.len().is_multiple_of(2);
    let middle = values.len() / 2;
    let (lower, upper, _) = values.select_nth_unstable_by(middle, f64::total_cmp);
    let upper = *upper;
    Some(match lower.iter().copied().max_by(f64::total_cmp) {
        Some(lower) if even => f64::midpoint(lower, upper),
        _ => upper,
    })
}

/// The upper median of `values` under `order`: the middle value, the upper
/// of the two middle values for an even count; `None` for none. Reorders
/// `values`.
pub(crate) fn upper<T: Copy>(values: &mut [T], order: impl FnMut(&T, &T) -> Ordering) -> Option<T> {
    if values.is_empty() {
        return None;
    }
    let middle = values.len() / 2;
    Some(*values.select_nth_unstable_by(middle, order).1)
}

/// The lower median of a histogram, given as (value, count) pairs in
/// increasing value order with `total` counts: the first value whose
/// cumulative count reaches half the total, rounded up; the last value when
/// none does.
pub(crate) fn lower_of_histogram(values: impl Iterator<Item = (u64, u64)>, total: u64) -> u64 {
    let middle = total.div_ceil(2).max(1);
    let mut seen = 0;
    let mut last = 0;
    for (value, count) in values {
        seen += count;
        last = value;
        if seen >= middle {
            return value;
        }
    }
    last
}

#[cfg(test)]
mod tests {
    use super::{lower_of_histogram, midpoint, upper};

    #[test]
    fn the_midpoint_median_averages_the_two_middle_values() {
        assert_eq!(midpoint(&mut []), None);
        assert_eq!(midpoint(&mut [3.0, 1.0, 2.0]), Some(2.0));
        assert_eq!(midpoint(&mut [4.0, 1.0, 3.0, 2.0]), Some(2.5));
    }

    #[test]
    fn the_upper_median_takes_the_upper_middle_value() {
        assert_eq!(upper::<usize>(&mut [], Ord::cmp), None);
        assert_eq!(upper(&mut [3, 1, 2], Ord::cmp), Some(2));
        assert_eq!(upper(&mut [4, 1, 3, 2], Ord::cmp), Some(3));
        assert_eq!(upper(&mut [4.0, 1.0, 3.0, 2.0], f64::total_cmp), Some(3.0));
    }

    #[test]
    fn the_histogram_median_takes_the_lower_middle_value() {
        let histogram = [(1, 2), (5, 2)];
        assert_eq!(lower_of_histogram(histogram.into_iter(), 4), 1);
        assert_eq!(lower_of_histogram([(1, 1), (5, 2)].into_iter(), 3), 5);
        assert_eq!(lower_of_histogram(std::iter::empty(), 0), 0);
    }
}
