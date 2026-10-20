//! Where a tract opens, from the read's own primary calls (PROP-0012
//! increment 2).
//!
//! A read loses phase after a homopolymer run whose length varies between
//! its molecules. Given the call where the modality sees the read leave phase,
//! the tract is the latest run of at least `min_run` equal primaries starting
//! before it, extended back over each single call of another base preceded
//! by at least `min_before` calls of the run base: the mark between two runs
//! of one tract. It reads no reference.

/// The method's constants.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct LocateParams {
    /// The fewest equal primaries that open a tract.
    pub min_run: usize,
    /// The fewest run-base calls before a single-call mark that join the run
    /// before it to the tract.
    pub min_before: usize,
}

/// The tract that opens before call `onset` of `primaries` (channel
/// indices): its first call and its run base. `None` when no run of
/// `min_run` equal primaries starts before `onset`.
#[must_use]
pub fn locate(primaries: &[usize], onset: usize, params: &LocateParams) -> Option<(usize, usize)> {
    let mut latest = None;
    let mut index = 0;
    while index < primaries.len() && index < onset {
        let base = primaries[index];
        let length = primaries[index..]
            .iter()
            .take_while(|&&next| next == base)
            .count();
        if length >= params.min_run.max(1) {
            latest = Some(index);
        }
        index += length;
    }
    let mut start = latest?;
    let run = primaries[start];
    while start > params.min_before
        && primaries[start - 1] != run
        && (1..=params.min_before).all(|back| primaries[start - 1 - back] == run)
    {
        start -= 1;
        while start > 0 && primaries[start - 1] == run {
            start -= 1;
        }
    }
    Some((start, run))
}

#[cfg(test)]
mod tests {
    use super::{LocateParams, locate};

    const PARAMS: LocateParams = LocateParams {
        min_run: 4,
        min_before: 2,
    };

    fn primaries(sequence: &str) -> Vec<usize> {
        sequence
            .chars()
            .map(|base| "ACGT".find(base).unwrap_or_default())
            .collect()
    }

    #[test]
    fn takes_the_latest_run_before_the_onset() {
        let calls = primaries("GGGGATCCCCCAGT");
        assert_eq!(locate(&calls, 11, &PARAMS), Some((6, 1)));
        // A run starting at the onset is not before it.
        assert_eq!(locate(&calls, 6, &PARAMS), Some((0, 2)));
    }

    #[test]
    fn extends_back_over_a_single_call_mark() {
        // CCCCC T CCCCCC: the mark joins the two runs of one tract.
        let calls = primaries("AGCCCCCTCCCCCCGA");
        assert_eq!(locate(&calls, 14, &PARAMS), Some((2, 1)));
        // One run-base call before the mark is too few.
        let calls = primaries("AGACTCCCCCCGA");
        assert_eq!(locate(&calls, 11, &PARAMS), Some((5, 1)));
        // Two marks in a row end the tract.
        let calls = primaries("CCCCTTCCCCCGA");
        assert_eq!(locate(&calls, 11, &PARAMS), Some((6, 1)));
    }

    #[test]
    fn finds_nothing_without_a_run() {
        let calls = primaries("ACGTACGTCCCACGT");
        assert_eq!(locate(&calls, 12, &PARAMS), None);
        assert_eq!(locate(&[], 3, &PARAMS), None);
    }
}
