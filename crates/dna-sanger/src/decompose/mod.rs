//! Length mixtures after a repeat run (SRS-CALL-023, ADR-0082).
//!
//! A run that carries length heteroplasmy is a population of molecules whose
//! runs differ in length. After it, a read shows their sum: its own dominant
//! sequence `S` (its primary calls) superposed with copies of itself offset by
//! whole calls. [`mixtures`] fits that model to every segment callability marks
//! dephased, over the segment and the run that opens it: the longest
//! homopolymer of the read's own primary calls, of two calls or more, whose
//! end lies within `callability.repeat_window_calls` calls of the segment's
//! start and which starts at or after the trim interval's start; a stretch in
//! the trimmed start is start noise, too short to place. The
//! stretch ends at the segment's end, or where a run of the read's primary
//! calls at least as long as the opening one starts, a new source of mixture;
//! shorter runs inside it are part of the tract the offsets explain:
//!
//! `y_i = Σ_d f_d · onehot(S[i − d])`, `f_d ≥ 0`,
//!
//! where `y_i` is call `i`'s normalised A/C/G/T profile from signal
//! processing. The offsets are bounded by the run's own length. The fit starts
//! from offset 0 alone and adds, one at a time, the offset that lowers the
//! residual most, while it lowers it by more than the read's own noise level:
//! the root mean square of `profile − onehot(primary)` over the read's in-phase
//! single-channel calls. A call of the segment whose own residual exceeds the
//! `signal_processing.decompose_unexplained_quantile` nearest-rank quantile of
//! those calls' residuals is unexplained. The mixture is explained when its
//! own residual is at most that quantile: the components fit the stretch as
//! well as the read's in-phase calls fit their own bases. Otherwise they are
//! only the fit's best, and no sequence is read from them. A read with fewer
//! such calls than `callability.minimum_callable_calls` asserts no mixture.
//!
//! The dominant component is the read's own primary calls over the masked
//! segment; the fit asserts that they are the dominant sequence of a
//! self-shifted mixture and gives the minor offsets and fractions. Nothing
//! reads a reference.

use dna_kernel::model::nucleotide::Nucleotide;
use dna_kernel::read_evidence::{MixtureComponent, PhaseMixture};

use crate::model::basecalls::BaseCalls;
use crate::model::callability::{PhaseState, ReadCallability};
use crate::model::signal::SignalAnalysis;

/// The phase mixtures of one read, one per dephased segment that a repeat
/// run starting at or after `trim_start` opens, in call order. `None` when the read has too
/// few in-phase single-channel calls to measure its own noise.
#[must_use]
pub fn mixtures(
    calls: &BaseCalls,
    signal: &SignalAnalysis,
    callability: &ReadCallability,
    trim_start: usize,
    unexplained_quantile: f64,
    minimum_calls: usize,
    onset_window: usize,
) -> Option<Vec<PhaseMixture>> {
    let profile = |index: usize| {
        signal
            .loci
            .get(index)
            .filter(|locus| locus.call_index_0based == index)
            .and_then(|locus| locus.profile)
            .map(|profile| profile.weights)
    };
    let onehot = |index: usize| -> [f64; 4] {
        let mut vector = [0.0; 4];
        if let Some(base) = calls
            .calls
            .get(index)
            .and_then(|call| Nucleotide::from_char(call.primary))
        {
            vector[base.channel_index()] = 1.0;
        }
        vector
    };
    // The read's own residuals where it has no mixture.
    let mut in_phase = (0..calls.calls.len())
        .filter(|&index| {
            callability.mask.get(index).is_some_and(Option::is_none)
                && calls.calls[index].qualifying_channels.len() == 1
        })
        .filter_map(|index| Some(distance(&profile(index)?, &onehot(index))))
        .collect::<Vec<_>>();
    if in_phase.len() < minimum_calls.max(1) {
        tracing::debug!(
            event = "phase_mixture_skipped",
            in_phase_calls = in_phase.len(),
            minimum = minimum_calls,
        );
        return None;
    }
    let noise =
        (in_phase.iter().map(|value| value * value).sum::<f64>() / in_phase.len() as f64).sqrt();
    in_phase.sort_by(f64::total_cmp);
    // The quantile lies in (0, 1), so the nearest rank lies in 1..=len.
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let rank =
        ((unexplained_quantile * in_phase.len() as f64).ceil() as usize).clamp(1, in_phase.len());
    let unexplained_level = in_phase[rank - 1];

    let mut found = Vec::new();
    for segment in &callability.segments {
        if segment.state != PhaseState::Dephased {
            continue;
        }
        let Some((start, run_end)) = opening_run(calls, segment.call_start_0based, onset_window)
        else {
            continue;
        };
        if start < trim_start {
            continue;
        }
        let run_calls = run_end - start;
        // The mixture ends at the segment's end, or where a run at least as
        // long as the opening one starts: a new source of mixture.
        let end = next_run(calls, run_end, run_calls)
            .unwrap_or(usize::MAX)
            .min(segment.call_end_0based_exclusive)
            .min(calls.calls.len());
        if end <= run_end {
            continue;
        }
        let rows = (start..end)
            .filter_map(|index| Some((index, profile(index)?)))
            .collect::<Vec<_>>();
        if rows.is_empty() {
            continue;
        }
        // The model column of offset `d` at call `i`: the dominant call `d`
        // calls before it.
        let column = |offset: i32, index: usize| {
            index
                .checked_add_signed(-(offset as isize))
                .filter(|&source| source < calls.calls.len())
                .map_or([0.0; 4], &onehot)
        };
        let fit = |offsets: &[i32]| -> (Vec<f64>, f64) {
            let design = rows
                .iter()
                .flat_map(|&(index, _)| (0..4).map(move |channel| (index, channel)))
                .map(|(index, channel)| {
                    offsets
                        .iter()
                        .map(|&offset| column(offset, index)[channel])
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>();
            let observed = rows
                .iter()
                .flat_map(|(_, profile)| profile.iter().copied())
                .collect::<Vec<_>>();
            let weights = nnls(&design, &observed);
            let squared = design
                .iter()
                .zip(&observed)
                .map(|(row, y)| {
                    let fitted = row.iter().zip(&weights).map(|(a, w)| a * w).sum::<f64>();
                    (y - fitted).powi(2)
                })
                .sum::<f64>();
            (weights, (squared / rows.len() as f64).sqrt())
        };
        let bound = i32::try_from(run_calls).unwrap_or(i32::MAX);
        let mut offsets = vec![0];
        let (mut weights, mut residual) = fit(&offsets);
        loop {
            let best = (-bound..=bound)
                .filter(|offset| !offsets.contains(offset))
                .map(|offset| {
                    let mut trial = offsets.clone();
                    trial.push(offset);
                    let (trial_weights, trial_residual) = fit(&trial);
                    (offset, trial_weights, trial_residual)
                })
                .min_by(|left, right| {
                    left.2
                        .total_cmp(&right.2)
                        .then(left.0.abs().cmp(&right.0.abs()))
                });
            match best {
                Some((offset, trial_weights, trial_residual))
                    if residual - trial_residual > noise =>
                {
                    offsets.push(offset);
                    weights = trial_weights;
                    residual = trial_residual;
                }
                _ => break,
            }
        }
        let total = weights.iter().sum::<f64>();
        if total <= 0.0 {
            continue;
        }
        let dominant = MixtureComponent {
            offset: 0,
            fraction: weights[0] / total,
        };
        let mut minors = offsets
            .iter()
            .zip(&weights)
            .skip(1)
            .filter(|(_, weight)| **weight > 0.0)
            .map(|(&offset, &weight)| MixtureComponent {
                offset,
                fraction: weight / total,
            })
            .collect::<Vec<_>>();
        minors.sort_by(|left, right| {
            right
                .fraction
                .total_cmp(&left.fraction)
                .then(left.offset.cmp(&right.offset))
        });
        let components = std::iter::once(dominant).chain(minors).collect::<Vec<_>>();
        let unexplained = rows
            .iter()
            .filter(|&&(index, _)| index >= segment.call_start_0based)
            .filter(|&&(index, profile)| {
                let mut fitted = [0.0; 4];
                for component in &components {
                    let source = column(component.offset, index);
                    for channel in 0..4 {
                        fitted[channel] += component.fraction * source[channel];
                    }
                }
                distance(&profile, &fitted) > unexplained_level
            })
            .map(|&(index, _)| index)
            .collect();
        found.push(PhaseMixture {
            call_start: start,
            call_end: end,
            run_calls,
            components,
            residual,
            unexplained,
            explained: residual <= unexplained_level,
        });
    }
    Some(found)
}

/// The run that opens a segment starting at call `start`: the longest maximal
/// run of two or more equal canonical primary calls whose end lies within
/// `window` calls of `start`, the nearer on a tie, as `(first, end)`. The
/// longest run there is the one whose length varies most.
fn opening_run(calls: &BaseCalls, start: usize, window: usize) -> Option<(usize, usize)> {
    let primary = |index: usize| {
        calls
            .calls
            .get(index)
            .map(|call| call.primary)
            .filter(|&base| Nucleotide::from_char(base).is_some())
    };
    let mut runs = Vec::new();
    let mut index = start.saturating_sub(window);
    // Begin at the start of the run holding the window's first call.
    while index > 0 && primary(index).is_some() && primary(index - 1) == primary(index) {
        index -= 1;
    }
    let limit = (start + window).min(calls.calls.len());
    while index < limit {
        let Some(base) = primary(index) else {
            index += 1;
            continue;
        };
        let mut end = index + 1;
        while primary(end) == Some(base) {
            end += 1;
        }
        if end - index >= 2 && end.abs_diff(start) <= window {
            runs.push((index, end));
        }
        index = end;
    }
    runs.into_iter()
        .min_by_key(|&(first, end)| (std::cmp::Reverse(end - first), end.abs_diff(start)))
}

/// The first call, at or after `from`, of a run of equal canonical primary
/// calls at least `length` long.
fn next_run(calls: &BaseCalls, from: usize, length: usize) -> Option<usize> {
    let primary = |index: usize| {
        calls
            .calls
            .get(index)
            .map(|call| call.primary)
            .filter(|&base| Nucleotide::from_char(base).is_some())
    };
    let mut index = from;
    while index < calls.calls.len() {
        let Some(base) = primary(index) else {
            index += 1;
            continue;
        };
        let mut end = index + 1;
        while primary(end) == Some(base) {
            end += 1;
        }
        if end - index >= length {
            return Some(index);
        }
        index = end;
    }
    None
}

/// Euclidean distance of two A/C/G/T vectors.
fn distance(left: &[f64; 4], right: &[f64; 4]) -> f64 {
    left.iter()
        .zip(right)
        .map(|(a, b)| (a - b).powi(2))
        .sum::<f64>()
        .sqrt()
}

/// Gradient below which no column enters the NNLS active set.
const TOLERANCE: f64 = 1e-10;

/// Non-negative least squares, `min ‖A·x − y‖` with `x ≥ 0`, by the
/// Lawson–Hanson active-set method; `design` holds the rows of `A`.
fn nnls(design: &[Vec<f64>], observed: &[f64]) -> Vec<f64> {
    let columns = design.first().map_or(0, Vec::len);
    let mut solution = vec![0.0; columns];
    let mut passive = vec![false; columns];
    let gradient = |solution: &[f64]| -> Vec<f64> {
        (0..columns)
            .map(|column| {
                design
                    .iter()
                    .zip(observed)
                    .map(|(row, y)| {
                        let fitted = row.iter().zip(solution).map(|(a, x)| a * x).sum::<f64>();
                        row[column] * (y - fitted)
                    })
                    .sum::<f64>()
            })
            .collect()
    };
    for _ in 0..(3 * columns + 3) {
        let descent = gradient(&solution);
        let Some(entering) = (0..columns)
            .filter(|&column| !passive[column] && descent[column] > TOLERANCE)
            .max_by(|&left, &right| descent[left].total_cmp(&descent[right]))
        else {
            break;
        };
        passive[entering] = true;
        loop {
            let unconstrained = least_squares(design, observed, &passive);
            if (0..columns)
                .filter(|&column| passive[column])
                .all(|column| unconstrained[column] > 0.0)
            {
                solution = unconstrained;
                break;
            }
            let step = (0..columns)
                .filter(|&column| passive[column] && unconstrained[column] <= 0.0)
                .map(|column| solution[column] / (solution[column] - unconstrained[column]))
                .fold(f64::INFINITY, f64::min);
            for column in 0..columns {
                if passive[column] {
                    solution[column] += step * (unconstrained[column] - solution[column]);
                    if solution[column] <= TOLERANCE {
                        solution[column] = 0.0;
                        passive[column] = false;
                    }
                }
            }
            if !passive.iter().any(|&flag| flag) {
                break;
            }
        }
    }
    solution
}

/// Unconstrained least squares over the `passive` columns by the normal
/// equations; the other columns are zero.
fn least_squares(design: &[Vec<f64>], observed: &[f64], passive: &[bool]) -> Vec<f64> {
    let chosen = (0..passive.len())
        .filter(|&column| passive[column])
        .collect::<Vec<_>>();
    let size = chosen.len();
    let mut matrix = vec![vec![0.0; size + 1]; size];
    for (row, y) in design.iter().zip(observed) {
        for (i, &left) in chosen.iter().enumerate() {
            for (j, &right) in chosen.iter().enumerate() {
                matrix[i][j] += row[left] * row[right];
            }
            matrix[i][size] += row[left] * y;
        }
    }
    // Gaussian elimination with partial pivoting; a singular column stays 0.
    for pivot in 0..size {
        let Some(best) = (pivot..size).max_by(|&left, &right| {
            matrix[left][pivot]
                .abs()
                .total_cmp(&matrix[right][pivot].abs())
        }) else {
            continue;
        };
        matrix.swap(pivot, best);
        let value = matrix[pivot][pivot];
        if value.abs() < 1e-12 {
            continue;
        }
        let pivot_row = matrix[pivot].clone();
        for (index, row) in matrix.iter_mut().enumerate() {
            if index != pivot {
                let factor = row[pivot] / value;
                for (cell, above) in row.iter_mut().zip(&pivot_row).skip(pivot) {
                    *cell -= factor * above;
                }
            }
        }
    }
    let mut solution = vec![0.0; passive.len()];
    for (i, &column) in chosen.iter().enumerate() {
        if matrix[i][i].abs() >= 1e-12 {
            solution[column] = matrix[i][size] / matrix[i][i];
        }
    }
    solution
}

// SRS-CALL-023
#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::basecalls::{BaseCall, ChannelPeak, PeakSource};
    use crate::model::callability::PhaseSegment;
    use crate::model::locus_evidence::LocusEvidence;
    use crate::model::signal::SangerIntegrity;
    use dna_kernel::read_evidence::EvidenceProfile;

    /// In-phase calls before the run, the run, and the tail after it.
    const PREFIX: &str = "ACGTTGCAAGCTTGACATGCAGTCAGTACGATGCATCAGT";
    const RUN: &str = "CCCCCCC";
    const TAIL: &str = "TACGGATCAATGGTACATCGGATACTGACT";

    fn sequence() -> String {
        format!("{PREFIX}{RUN}{TAIL}")
    }

    fn onehot(base: char) -> [f64; 4] {
        let mut vector = [0.0; 4];
        if let Some(base) = Nucleotide::from_char(base) {
            vector[base.channel_index()] = 1.0;
        }
        vector
    }

    /// A read of `sequence()` whose calls from the run on carry `profile`,
    /// and in-phase calls before it carry a little of the next channel.
    fn read(
        profile: impl Fn(usize, &[char]) -> [f64; 4],
        dephased: bool,
    ) -> (BaseCalls, SignalAnalysis, ReadCallability) {
        let bases = sequence().chars().collect::<Vec<_>>();
        let count = bases.len();
        let run_start = PREFIX.len();
        let run_end = run_start + RUN.len();
        let calls = bases
            .iter()
            .enumerate()
            .map(|(index, &base)| BaseCall {
                index_0based: index,
                locus_position_0based: index * 10,
                window_start_0based: index * 10,
                window_end_0based_exclusive: index * 10 + 10,
                peaks: std::array::from_fn(|channel| ChannelPeak {
                    base: Nucleotide::ALL[channel],
                    height: 100,
                    position_0based: index * 10,
                    source: PeakSource::LocalMaximum,
                }),
                primary_peak_evidence: None,
                primary: base,
                ambiguity: base,
                qualifying_channels: vec![Nucleotide::from_char(base).unwrap_or(Nucleotide::A)],
                vendor_agrees: None,
            })
            .collect::<Vec<_>>();
        let loci = (0..count)
            .map(|index| {
                let weights = if index < run_start {
                    // In phase: 0.96 on the call's base, 0.04 on the next.
                    let mut weights = onehot(bases[index]).map(|value| value * 0.96);
                    let next = (Nucleotide::from_char(bases[index])
                        .map_or(0, Nucleotide::channel_index)
                        + 1)
                        % 4;
                    weights[next] += 0.04;
                    weights
                } else {
                    profile(index, &bases)
                };
                LocusEvidence {
                    call_index_0based: index,
                    locus_position_0based: index * 10,
                    window_start_0based: index * 10,
                    window_end_0based_exclusive: index * 10 + 10,
                    context_call_start_0based: index,
                    context_call_end_0based_exclusive: index + 1,
                    context_sample_start_0based: index * 10,
                    context_sample_end_0based_exclusive: index * 10 + 10,
                    event_position_0based: index * 10,
                    channel_heights: [0; 4],
                    channel_baselines: [0.0; 4],
                    channel_noise_sigmas: [1.0; 4],
                    corrected_amplitudes: [0.0; 4],
                    snrs: [0.0; 4],
                    profile: Some(EvidenceProfile { weights }),
                    channel_widths: [4; 4],
                }
            })
            .collect();
        let state = if dephased {
            PhaseState::Dephased
        } else {
            PhaseState::InPhase
        };
        let callability = ReadCallability {
            repeats: Vec::new(),
            segments: vec![
                PhaseSegment {
                    call_start_0based: 0,
                    call_end_0based_exclusive: run_end,
                    state: PhaseState::InPhase,
                    after_repeat: false,
                    shadow: None,
                },
                PhaseSegment {
                    call_start_0based: run_end,
                    call_end_0based_exclusive: count,
                    state,
                    after_repeat: false,
                    shadow: None,
                },
            ],
            mask: (0..count)
                .map(|index| (dephased && index >= run_end).then_some(state))
                .collect(),
            callable_start_0based: 0,
            callable_end_0based_exclusive: if dephased { run_end } else { count },
        };
        (
            BaseCalls {
                primary_sequence: sequence(),
                calls,
            },
            SignalAnalysis {
                integrity: SangerIntegrity {
                    locus_count: count,
                    duplicate_locus_count: 0,
                    vendor_primary_count: None,
                    vendor_quality_count: None,
                    minimum_locus_spacing: None,
                    median_locus_spacing: None,
                    maximum_locus_spacing: None,
                    clipped_channel_samples: 0,
                    maximum_to_median_event_signal_ratio: None,
                },
                loci,
                windows: Vec::new(),
                noisy_regions: Vec::new(),
            },
            callability,
        )
    }

    /// The profile of a mixture of `components` (offset, fraction) at call
    /// `index`.
    fn mixed(components: &[(i32, f64)]) -> impl Fn(usize, &[char]) -> [f64; 4] + '_ {
        move |index, bases| {
            let mut weights = [0.0; 4];
            for &(offset, fraction) in components {
                if let Some(&base) = index
                    .checked_add_signed(-(offset as isize))
                    .and_then(|source| bases.get(source))
                {
                    let vector = onehot(base);
                    for channel in 0..4 {
                        weights[channel] += fraction * vector[channel];
                    }
                }
            }
            weights
        }
    }

    fn fit(
        profile: impl Fn(usize, &[char]) -> [f64; 4],
        dephased: bool,
    ) -> Option<Vec<PhaseMixture>> {
        let (calls, signal, callability) = read(profile, dephased);
        mixtures(&calls, &signal, &callability, 0, 0.99, 20, 16)
    }

    fn offsets(mixture: &PhaseMixture) -> Vec<(i32, f64)> {
        mixture
            .components
            .iter()
            .map(|component| {
                (
                    component.offset,
                    (component.fraction * 100.0).round() / 100.0,
                )
            })
            .collect()
    }

    #[test]
    fn a_single_sequence_is_one_component() {
        let found = fit(mixed(&[(0, 1.0)]), true).unwrap_or_default();
        assert_eq!(found.len(), 1);
        assert_eq!(offsets(&found[0]), [(0, 1.0)]);
        assert_eq!(found[0].call_start, PREFIX.len());
        assert_eq!(found[0].run_calls, RUN.len());
        assert!(found[0].unexplained.is_empty());
        assert!(found[0].explained);
    }

    #[test]
    fn a_longer_run_in_a_fifth_of_the_molecules_is_a_minor_component() {
        let found = fit(mixed(&[(0, 0.8), (1, 0.2)]), true).unwrap_or_default();
        assert_eq!(offsets(&found[0]), [(0, 0.8), (1, 0.2)]);
        assert!(found[0].explained);
        let shorter = fit(mixed(&[(0, 0.8), (-1, 0.2)]), true).unwrap_or_default();
        assert_eq!(offsets(&shorter[0]), [(0, 0.8), (-1, 0.2)]);
    }

    #[test]
    fn three_lengths_are_three_components() {
        let found = fit(mixed(&[(0, 0.6), (-1, 0.25), (1, 0.15)]), true).unwrap_or_default();
        assert_eq!(offsets(&found[0]), [(0, 0.6), (-1, 0.25), (1, 0.15)]);
    }

    #[test]
    fn a_tail_call_no_component_explains_is_unexplained() {
        let substituted = PREFIX.len() + RUN.len() + 10;
        let profile = |index: usize, bases: &[char]| {
            if index == substituted {
                onehot('G')
                    .map(|value| value * 0.5)
                    .iter()
                    .zip(onehot('T'))
                    .map(|(a, b)| a + 0.5 * b)
                    .collect::<Vec<_>>()
                    .try_into()
                    .unwrap_or([0.0; 4])
            } else {
                mixed(&[(0, 0.8), (1, 0.2)])(index, bases)
            }
        };
        let found = fit(profile, true).unwrap_or_default();
        assert!(found[0].unexplained.contains(&substituted));
    }

    #[test]
    fn a_tail_no_offsets_explain_is_not_explained() {
        // Every channel equal: no sum of shifted copies fits it as well as
        // the read's in-phase calls fit their own bases.
        let tail = PREFIX.len() + RUN.len();
        let profile = |index: usize, bases: &[char]| {
            if index < tail {
                onehot(bases[index])
            } else {
                [0.25; 4]
            }
        };
        let found = fit(profile, true).unwrap_or_default();
        assert_eq!(found.len(), 1);
        assert!(!found[0].explained);
    }

    #[test]
    fn a_run_in_the_trimmed_start_opens_no_mixture() {
        let (calls, signal, callability) = read(mixed(&[(0, 0.8), (1, 0.2)]), true);
        assert_eq!(
            mixtures(
                &calls,
                &signal,
                &callability,
                PREFIX.len() + 1,
                0.99,
                20,
                16
            ),
            Some(Vec::new())
        );
    }

    #[test]
    fn an_in_phase_tail_after_a_run_is_no_mixture() {
        assert_eq!(fit(mixed(&[(0, 1.0)]), false), Some(Vec::new()));
    }

    #[test]
    fn a_read_too_short_to_measure_its_noise_asserts_none() {
        let (calls, signal, callability) = read(mixed(&[(0, 0.8), (1, 0.2)]), true);
        assert_eq!(
            mixtures(&calls, &signal, &callability, 0, 0.99, 1000, 16),
            None
        );
    }

    #[test]
    fn a_mixture_ends_where_a_run_as_long_as_the_opening_one_starts() {
        let (calls, _, _) = read(mixed(&[(0, 1.0)]), true);
        // The run of seven C opens at the prefix's end; none as long follows.
        assert_eq!(next_run(&calls, 0, RUN.len()), Some(PREFIX.len()));
        assert_eq!(next_run(&calls, PREFIX.len() + RUN.len(), RUN.len()), None);
        // The first run of two is the prefix's TT at 3.
        assert_eq!(next_run(&calls, 0, 2), Some(3));
    }
}
