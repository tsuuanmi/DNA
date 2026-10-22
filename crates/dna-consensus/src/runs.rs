//! Runs of one base: a sequence's runs, the reference run holding a
//! position, and the anchored pair of runs in a list of runs. The consensus
//! and the review read runs through these alone.

use crate::ReadRun;
use crate::model::consensus::ConsensusRun;

/// A run of one base with its length, whoever reads it.
pub trait Run {
    /// The run's base.
    fn base(&self) -> char;
    /// The run's length.
    fn length(&self) -> usize;
}

impl Run for (char, usize) {
    fn base(&self) -> char {
        self.0
    }

    fn length(&self) -> usize {
        self.1
    }
}

impl Run for ReadRun {
    fn base(&self) -> char {
        self.base
    }

    fn length(&self) -> usize {
        self.length
    }
}

impl Run for ConsensusRun {
    fn base(&self) -> char {
        self.base
    }

    fn length(&self) -> usize {
        self.length
    }
}

/// The index of the first run of `repeat` in `runs` that a single `anchor`
/// base and another run of `repeat` follow: an anchored pair, as an anchored
/// repeat window holds it. `None` without one.
#[must_use]
pub fn anchored_pair<R: Run>(runs: &[R], repeat: char, anchor: char) -> Option<usize> {
    runs.windows(3).position(|triple| {
        triple[0].base() == repeat
            && triple[1].base() == anchor
            && triple[1].length() == 1
            && triple[2].base() == repeat
    })
}

/// Runs of equal bases as `(base, length)`, in order.
#[must_use]
pub fn runs_of(bases: impl IntoIterator<Item = char>) -> Vec<(char, usize)> {
    let mut runs: Vec<(char, usize)> = Vec::new();
    for base in bases {
        match runs.last_mut() {
            Some((previous, count)) if *previous == base => *count += 1,
            _ => runs.push((base, 1)),
        }
    }
    runs
}

/// The run of one base of `bases` holding the 0-based `position`, as its
/// first and last positions (0-based, inclusive); `None` past the end. The
/// run does not wrap the reference origin.
#[must_use]
pub fn reference_run(bases: &[u8], position: usize) -> Option<(usize, usize)> {
    let &base = bases.get(position)?;
    let mut first = position;
    while first > 0 && bases[first - 1] == base {
        first -= 1;
    }
    let mut last = position;
    while last + 1 < bases.len() && bases[last + 1] == base {
        last += 1;
    }
    Some((first, last))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn splits_a_sequence_into_runs() {
        assert_eq!(runs_of("CCCTCC".chars()), [('C', 3), ('T', 1), ('C', 2)]);
        assert!(runs_of("".chars()).is_empty());
    }

    #[test]
    fn finds_the_first_anchored_pair() {
        let runs = runs_of("CCTCCCTTCCTCC".chars());
        assert_eq!(anchored_pair(&runs, 'C', 'T'), Some(0));
        assert_eq!(anchored_pair(&runs[1..], 'C', 'T'), Some(3));
        assert_eq!(anchored_pair(&runs_of("CCTTCC".chars()), 'C', 'T'), None);
        assert_eq!(anchored_pair(&runs_of("CC".chars()), 'C', 'T'), None);
    }

    #[test]
    fn finds_the_reference_run_holding_a_position() {
        let bases = b"ACCCTA";
        assert_eq!(reference_run(bases, 2), Some((1, 3)));
        assert_eq!(reference_run(bases, 0), Some((0, 0)));
        assert_eq!(reference_run(bases, 5), Some((5, 5)));
        assert_eq!(reference_run(bases, 6), None);
    }
}
