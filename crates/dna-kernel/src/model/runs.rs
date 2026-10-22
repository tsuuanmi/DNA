//! Maximal runs of equal items: the one scan for runs of equal calls.

use std::ops::Range;

/// The maximal runs of consecutive `items` with equal `key`, in order, as
/// index ranges that together cover every item once. Every item starts or
/// extends a run, unresolved ones included: what a run of unresolved items
/// means is the caller's to decide.
pub fn maximal_runs<T, K: PartialEq>(
    items: &[T],
    key: impl Fn(&T) -> K,
) -> impl Iterator<Item = Range<usize>> {
    let mut start = 0;
    std::iter::from_fn(move || {
        let value = key(items.get(start)?);
        let length = items[start..]
            .iter()
            .take_while(|item| key(item) == value)
            .count()
            .max(1);
        let run = start..start + length;
        start += length;
        Some(run)
    })
}

#[cfg(test)]
mod tests {
    use super::maximal_runs;

    #[test]
    fn splits_items_into_maximal_runs() {
        let runs = maximal_runs(b"AACNNNT".as_slice(), |base| *base).collect::<Vec<_>>();
        assert_eq!(runs, [0..2, 2..3, 3..6, 6..7]);
        assert_eq!(maximal_runs(&[] as &[u8], |base| *base).count(), 0);
        let keyed = maximal_runs(&[1, 3, 2, 4], |value| value % 2).collect::<Vec<_>>();
        assert_eq!(keyed, [0..2, 2..4]);
    }

    #[test]
    fn an_item_unequal_to_itself_is_a_run_of_its_own() {
        let runs = maximal_runs(&[f64::NAN, f64::NAN], |value| *value).collect::<Vec<_>>();
        assert_eq!(runs, [0..1, 1..2]);
    }
}
