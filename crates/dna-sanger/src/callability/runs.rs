//! Repeat runs in the read's own primary calls.
//!
//! A run is a prior for the phase segmentation: polymerase slippage starts in
//! long homopolymers and dinucleotide tandem repeats, so the onset threshold is
//! lowered right after one. Runs never mask by themselves.

use crate::model::callability::{RepeatRun, RepeatUnit};

/// Finds every maximal homopolymer and dinucleotide run of at least
/// `minimum_length` calls; an unresolved call breaks a run.
pub(super) fn find(primary: &[Option<usize>], minimum_length: usize) -> Vec<RepeatRun> {
    let mut runs = homopolymers(primary, minimum_length);
    runs.extend(dinucleotides(primary, minimum_length));
    runs.sort_by_key(|run| (run.call_start_0based, run.call_end_0based_exclusive));
    runs
}

fn homopolymers(primary: &[Option<usize>], minimum_length: usize) -> Vec<RepeatRun> {
    let mut runs = Vec::new();
    let mut start = 0;
    while start < primary.len() {
        let Some(channel) = primary[start] else {
            start += 1;
            continue;
        };
        let length = primary[start..]
            .iter()
            .take_while(|&&next| next == Some(channel))
            .count();
        if length >= minimum_length {
            runs.push(RepeatRun {
                call_start_0based: start,
                call_end_0based_exclusive: start + length,
                unit: RepeatUnit::Homopolymer(channel),
            });
        }
        start += length;
    }
    runs
}

fn dinucleotides(primary: &[Option<usize>], minimum_length: usize) -> Vec<RepeatRun> {
    let mut runs = Vec::new();
    let mut start = 0;
    while start + 1 < primary.len() {
        let (Some(first), Some(second)) = (primary[start], primary[start + 1]) else {
            start += 1;
            continue;
        };
        if first == second {
            start += 1;
            continue;
        }
        let mut length = 2;
        while start + length < primary.len()
            && primary[start + length] == primary[start + length - 2]
        {
            length += 1;
        }
        if length >= minimum_length {
            runs.push(RepeatRun {
                call_start_0based: start,
                call_end_0based_exclusive: start + length,
                unit: RepeatUnit::Dinucleotide(first, second),
            });
            // The run's last call may open the next alternating run.
            start += length - 1;
        } else {
            start += 1;
        }
    }
    runs
}

#[cfg(test)]
mod tests {
    use super::*;

    fn channels(sequence: &str) -> Vec<Option<usize>> {
        sequence
            .chars()
            .map(|base| match base {
                'A' => Some(0),
                'C' => Some(1),
                'G' => Some(2),
                'T' => Some(3),
                _ => None,
            })
            .collect()
    }

    #[test]
    fn finds_long_homopolymers_and_ignores_short_or_broken_runs() {
        let runs = find(&channels("ACCCCCCCCTAAAAAAAACCCCNCCCC"), 8);
        assert_eq!(
            runs,
            [
                RepeatRun {
                    call_start_0based: 1,
                    call_end_0based_exclusive: 9,
                    unit: RepeatUnit::Homopolymer(1),
                },
                RepeatRun {
                    call_start_0based: 10,
                    call_end_0based_exclusive: 18,
                    unit: RepeatUnit::Homopolymer(0),
                },
            ]
        );
    }

    #[test]
    fn finds_dinucleotide_tandem_repeats_in_call_order() {
        let runs = find(&channels("GACACACACACGTGTGTGTGTT"), 8);
        assert_eq!(
            runs,
            [
                RepeatRun {
                    call_start_0based: 1,
                    call_end_0based_exclusive: 11,
                    unit: RepeatUnit::Dinucleotide(0, 1),
                },
                RepeatRun {
                    call_start_0based: 11,
                    call_end_0based_exclusive: 21,
                    unit: RepeatUnit::Dinucleotide(2, 3),
                },
            ]
        );
    }

    #[test]
    fn returns_no_runs_for_short_reads_and_unresolved_calls() {
        assert!(find(&channels("ACGTNNNN"), 8).is_empty());
        assert!(find(&[], 8).is_empty());
    }
}
