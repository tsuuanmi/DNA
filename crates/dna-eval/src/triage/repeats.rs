//! The reference geometry the one-strand exception needs: the repeat windows
//! and the reference's homopolymer runs (one-strand research, R6m and
//! R6m-h4).

use crate::compare::Region;

/// Positions within this many bases of a repeat window count as inside it:
/// the post-homopolymer exits of the reverse reads (299–301, 16176–16179).
pub(crate) const WINDOW_MARGIN: i64 = 3;

/// Per reference position, the length of the homopolymer run that holds it.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub(crate) struct Repeats {
    runs: Vec<usize>,
}

impl Repeats {
    /// The runs of a reference sequence.
    pub(crate) fn new(reference: &str) -> Self {
        let bases = reference.as_bytes();
        let mut runs = vec![0; bases.len()];
        let mut start = 0;
        while start < bases.len() {
            let mut end = start;
            while end + 1 < bases.len() && bases[end + 1].eq_ignore_ascii_case(&bases[start]) {
                end += 1;
            }
            runs[start..=end].fill(end - start + 1);
            start = end + 1;
        }
        Self { runs }
    }

    /// Whether the 1-based `position` lies in a reference run of at least
    /// `length` bases, or next to one.
    pub(crate) fn near_run(&self, position: usize, length: usize) -> bool {
        let Some(index) = position.checked_sub(1) else {
            return false;
        };
        [index.checked_sub(1), Some(index), index.checked_add(1)]
            .into_iter()
            .flatten()
            .filter_map(|index| self.runs.get(index))
            .any(|run| *run >= length)
    }
}

/// Whether the 1-based `position` lies in a repeat window or within
/// [`WINDOW_MARGIN`] of one.
pub(crate) fn in_window(position: usize) -> bool {
    let Ok(position) = i64::try_from(position) else {
        return false;
    };
    Region::WINDOWS
        .iter()
        .any(|(_, first, last)| (first - WINDOW_MARGIN..=last + WINDOW_MARGIN).contains(&position))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn runs_of_three_four_and_five_hold_their_own_and_neighbouring_positions() {
        // 1-based: C at 1, AAA at 2-4, G at 5, TTTT at 6-9, C at 10,
        // GGGGG at 11-15, A at 16.
        let repeats = Repeats::new("CAAAGTTTTCGGGGGA");
        for (position, three, four, five) in [
            (1, true, false, false),
            (3, true, false, false),
            (5, true, true, false),
            (7, true, true, false),
            (10, true, true, true),
            (13, true, true, true),
            (16, true, true, true),
        ] {
            assert_eq!(repeats.near_run(position, 3), three, "{position} >= 3");
            assert_eq!(repeats.near_run(position, 4), four, "{position} >= 4");
            assert_eq!(repeats.near_run(position, 5), five, "{position} >= 5");
        }
        assert!(!repeats.near_run(0, 1) && !repeats.near_run(99, 1));
    }

    #[test]
    fn windows_are_the_measures_windows_with_a_margin_of_three() {
        // The policy and the measures read the same table.
        for (_, first, last) in Region::WINDOWS {
            for (position, inside) in [
                (first - WINDOW_MARGIN - 1, false),
                (first - WINDOW_MARGIN, true),
                (last + WINDOW_MARGIN, true),
                (last + WINDOW_MARGIN + 1, false),
            ] {
                let position = usize::try_from(position).unwrap_or(0);
                assert_eq!(in_window(position), inside, "{position}");
            }
        }
        assert!(in_window(299) && in_window(16177) && !in_window(400));
    }
}
