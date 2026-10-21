//! The second-channel signal a run's reads share (PROP-0009, SRS-CALL-019).
//!
//! A Sanger run shows the same shadow of a second base in the same sequence
//! context in most of its reads, at a level that changes from run to run. The
//! run signal records, per key, how often each second-channel share occurs
//! among a run's clean calls, so a call can later be judged against its own
//! run rather than a fixed share.
//!
//! The key is the read-strand context of the read's own primary calls around
//! a call, `run_signal_flank` calls on each side (a 5-mer by default), the
//! call's base, and one other channel: no reference position, so
//! the modality stays reference-free (ADR-0069). Every selected call counts
//! under each of its three other channels, so a key holds every call of its
//! context, not only those whose strongest other channel it is. The shares are kept as sparse
//! histograms in `run_signal_resolution` steps rather than quantiles, so the estimator can
//! change without rebuilding a run signal and a read's own contributions can
//! be taken out exactly.
//!
//! The selection is one function, [`contributions`], shared by building a run
//! signal and by taking a read out of it; [`SELECTION`] versions it.
//!
//! [`threshold`] is the one second-channel threshold every rule that judges a
//! call's second peak share applies: the configured floor, raised to the
//! key's run level when a run signal holds the key.

use std::collections::{BTreeMap, BTreeSet};

use dna_kernel::model::nucleotide::is_canonical;
use dna_kernel::read_evidence::{Placement, ReadEvidence};

use crate::error::{Result, RunSignalError};
use crate::model::sanger::Chromatogram;
use crate::model::second_cause::{SecondCause, SecondPeak};
use crate::read_processing::evidence::{PEAK_BELOW_MINIMUM, RELATIVE_QUALITY_NOT_ABOVE_THRESHOLD};
use crate::read_processing::{self, Admission, SangerConfig};
use crate::signal_processing::{RunSignalParameters, SignalProcessingConfig};
use dna_kernel::model::nucleotide::Nucleotide;

/// Version of the call selection [`contributions`] applies.
pub const SELECTION: u32 = 2;

/// What a call's share is counted under: the read-strand context of the
/// read's primary calls centred on the call, the call's base, and one of its
/// other channels.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Key {
    /// The read's primary calls from `flank` before to `flank` after the
    /// call.
    pub context: String,
    /// The call's base.
    pub base: char,
    /// The other channel whose share is counted.
    pub second: char,
}

impl Key {
    /// The key of call `index` of a read and its channel `second`, with
    /// `flank` calls of context on each side: `None` within `flank` calls of
    /// either end or without a canonical context.
    #[must_use]
    pub fn of(evidence: &ReadEvidence, index: usize, second: char, flank: usize) -> Option<Self> {
        let calls = evidence.calls();
        let context = calls.get(index.checked_sub(flank)?..=index + flank)?;
        context
            .iter()
            .all(|call| is_canonical(call.base))
            .then(|| Self {
                context: context.iter().map(|call| call.base).collect(),
                base: context[flank].base,
                second,
            })
    }
}

/// The second-channel threshold of one call (PROP-0009).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Threshold {
    /// The share, as a fraction, the call's second channel must exceed.
    pub share: f64,
    /// The key's calls in the run signal, when its run level was usable.
    pub calls: Option<u64>,
}

/// The canonical second-channel threshold (SRS-CALL-020): the configured
/// `second_floor`, raised to the run level of the call's `key` when `signal`
/// holds it with enough reads and calls ([`RunSignal::level`]). Without a
/// key or a run signal the floor alone decides.
#[must_use]
pub fn threshold(
    key: Option<&Key>,
    signal: Option<&RunSignal>,
    config: &SignalProcessingConfig,
) -> Threshold {
    let floor = config.second_floor;
    key.zip(signal)
        .and_then(|(key, signal)| {
            let level = signal.level(key, &config.run_signal)?;
            Some(Threshold {
                share: level.max(floor),
                calls: Some(signal.keys.get(key)?.values().sum()),
            })
        })
        .unwrap_or(Threshold {
            share: floor,
            calls: None,
        })
}

/// The cause of a call's second `peak` with the run signal (SRS-CALL-021): a
/// `second_template` peak whose share does not exceed the call's threshold
/// (SRS-CALL-020) at a usable run level of its `key` is a shadow the run's
/// reads share, `context`. Any other cause stands, and without a run signal
/// or a usable run level the shape alone decides.
#[must_use]
pub fn refine(
    peak: &SecondPeak,
    key: Option<&Key>,
    signal: Option<&RunSignal>,
    config: &SignalProcessingConfig,
) -> SecondCause {
    let threshold = threshold(key, signal, config);
    if peak.cause == SecondCause::SecondTemplate
        && threshold.calls.is_some()
        && peak.share <= threshold.share
    {
        SecondCause::Context
    } else {
        peak.cause
    }
}

/// Counts of calls by a channel's share in `resolution` steps (0 to
/// `resolution`): its peak
/// against the call's primary peak and itself, the measure of a call's
/// second base (SRS-CALL-015).
pub type Histogram = BTreeMap<u16, u64>;

/// The calls one read contributes to a run signal, by key.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Contributions {
    /// Selected calls.
    pub calls: u64,
    /// Their shares by key.
    pub keys: BTreeMap<Key, Histogram>,
}

/// The selected calls of one read (SRS-CALL-019): scored, informative, clear
/// of a signal edge for their base, without the peak-floor or
/// relative-quality reason (the mixed reason is kept: its shadow is what a
/// run signal measures), with a canonical base and canonical calls on both
/// sides, `parameters.flank` deep. A read's first and last `flank` calls have
/// no context. Each counts under each of its three other channels, with that
/// channel's share of its peak and the primary peak, in
/// `parameters.resolution` steps, from the read's `heights`: each call's
/// selected A/C/G/T peak heights in its window ([`BaseCalls::heights`](crate::model::basecalls::BaseCalls::heights)).
#[must_use]
pub fn contributions(
    heights: &[[i32; 4]],
    evidence: &ReadEvidence,
    parameters: &RunSignalParameters,
) -> Contributions {
    let (flank, resolution) = (parameters.flank, f64::from(parameters.resolution));
    let verdicts = evidence.calls();
    let mut selected = Contributions::default();
    for index in flank..verdicts.len().saturating_sub(flank) {
        let call = &verdicts[index];
        let context = verdicts[index - flank..=index + flank]
            .iter()
            .map(|call| call.base)
            .collect::<String>();
        let Some(peaks) = heights.get(index) else {
            continue;
        };
        let Some(primary) = Nucleotide::from_char(call.base) else {
            continue;
        };
        let height = peaks[primary.channel_index()];
        if !evidence.informative_call(index)
            || call.placement != Placement::Scored
            || call.identity.edge
            || call.identity.withheld.contains(PEAK_BELOW_MINIMUM)
            || call
                .identity
                .withheld
                .contains(RELATIVE_QUALITY_NOT_ABOVE_THRESHOLD)
            || !context.chars().all(is_canonical)
            || height <= 0
        {
            continue;
        }
        for other in Nucleotide::ALL {
            if other == primary {
                continue;
            }
            let peak = peaks[other.channel_index()].max(0);
            let total = f64::from(height) + f64::from(peak);
            // A share in [0, 1], so the rounded steps fit the resolution.
            #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
            let share = (resolution * f64::from(peak) / total)
                .round()
                .clamp(0.0, resolution) as u16;
            *selected
                .keys
                .entry(Key {
                    context: context.clone(),
                    base: call.base,
                    second: other.as_char(),
                })
                .or_default()
                .entry(share)
                .or_default() += 1;
        }
        selected.calls += 1;
    }
    selected
}

/// The contributions of one trace, processed as a sample command processes
/// it up to its core evidence; `None` for a read too short to process or with
/// too few callable calls, which a sample command rejects.
///
/// # Errors
///
/// Returns the failure of a Sanger stage other than those rejections.
pub fn read(trace: &Chromatogram, config: &SangerConfig) -> Result<Option<Contributions>> {
    let Admission::Processed(processed) = read_processing::admit(trace, config)? else {
        return Ok(None);
    };
    let evidence = read_processing::read_evidence(
        &processed.calls,
        &processed.signal,
        &processed.callability,
        &processed.quality,
        &processed.second,
        &processed.mixtures,
        &processed.tracts,
        &config.sanger_evidence,
    )?;
    Ok(Some(contributions(
        &processed.calls.heights(),
        &evidence,
        &config.signal_processing.run_signal,
    )))
}

/// A run's second-channel signal: per key, the shares of the selected calls
/// of every read of the run.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RunSignal {
    /// The run's identifier.
    pub run_id: String,
    /// The selection the contributions were made with ([`SELECTION`]).
    pub selection: u32,
    /// Calls of context on each side of a key's call.
    pub flank: usize,
    /// Steps of a share from 0 to 1.
    pub resolution: u16,
    /// Content SHA-256 of every trace of the run, rejected ones included.
    pub traces: BTreeSet<String>,
    /// Reads that contributed.
    pub reads: usize,
    /// Calls that contributed.
    pub calls: u64,
    /// Shares by key.
    pub keys: BTreeMap<Key, Histogram>,
}

impl RunSignal {
    /// The run signal of `run_id` over its `traces` (content SHA-256) and the
    /// contributions of those reads that could be processed, made with
    /// `parameters`' flank and resolution.
    #[must_use]
    pub fn build(
        run_id: String,
        traces: BTreeSet<String>,
        reads: impl IntoIterator<Item = Contributions>,
        parameters: &RunSignalParameters,
    ) -> Self {
        let mut signal = Self {
            run_id,
            selection: SELECTION,
            flank: parameters.flank,
            resolution: parameters.resolution,
            traces,
            reads: 0,
            calls: 0,
            keys: BTreeMap::new(),
        };
        for read in reads {
            signal.reads += 1;
            signal.calls += read.calls;
            for (key, histogram) in read.keys {
                let into = signal.keys.entry(key).or_default();
                for (share, count) in histogram {
                    *into.entry(share).or_default() += count;
                }
            }
        }
        signal
    }

    /// Adds `reads` reads, of content `traces`, and their `contributions`.
    pub fn include(
        &mut self,
        traces: impl IntoIterator<Item = String>,
        reads: usize,
        contributions: &Contributions,
    ) {
        self.traces.extend(traces);
        self.reads += reads;
        self.calls += contributions.calls;
        for (key, histogram) in &contributions.keys {
            let into = self.keys.entry(key.clone()).or_default();
            for (share, count) in histogram {
                *into.entry(*share).or_default() += count;
            }
        }
    }

    /// Refuses a run signal made with another flank or resolution than
    /// `parameters`': its keys or shares would not be this configuration's.
    ///
    /// # Errors
    ///
    /// Returns [`RunSignalError::Parameters`] on a mismatch.
    pub fn check(
        &self,
        parameters: &RunSignalParameters,
    ) -> std::result::Result<(), RunSignalError> {
        if self.flank == parameters.flank && self.resolution == parameters.resolution {
            return Ok(());
        }
        Err(RunSignalError::Parameters {
            flank: self.flank,
            resolution: self.resolution,
            expected_flank: parameters.flank,
            expected_resolution: parameters.resolution,
        })
    }

    /// The level above which a call of `key` stands out from its run:
    /// `median + max(margin, mad_factor * MAD)` of the key's shares, as a
    /// fraction. `None` when the run has fewer reads than `minimum_reads` or
    /// the key fewer calls than `minimum_calls`.
    #[must_use]
    pub fn level(&self, key: &Key, parameters: &RunSignalParameters) -> Option<f64> {
        if self.reads < parameters.minimum_reads {
            return None;
        }
        let histogram = self.keys.get(key)?;
        let calls = histogram.values().sum::<u64>();
        if calls < parameters.minimum_calls {
            return None;
        }
        let median = median(
            histogram
                .iter()
                .map(|(&share, &count)| (u64::from(share), count)),
            calls,
        );
        let mut deviations = BTreeMap::<u64, u64>::new();
        for (&share, &count) in histogram {
            *deviations
                .entry(u64::from(share).abs_diff(median))
                .or_default() += count;
        }
        let mad = median_of(&deviations, calls);
        let resolution = f64::from(self.resolution);
        let median = median as f64 / resolution;
        let spread = (parameters.mad_factor * mad as f64 / resolution).max(parameters.margin);
        Some(median + spread)
    }

    /// The same run signal without one read's contributions, so that a read
    /// is never judged against itself. Counts never fall below zero, and
    /// empty bins and keys are dropped; a read whose contributions exceed the
    /// run signal's is not of this run, which is logged once.
    #[must_use]
    pub fn without(&self, read: &Contributions) -> Self {
        let mut signal = self.clone();
        let mut saturated = read.calls > signal.calls;
        signal.reads = signal.reads.saturating_sub(1);
        signal.calls = signal.calls.saturating_sub(read.calls);
        for (key, histogram) in &read.keys {
            let Some(into) = signal.keys.get_mut(key) else {
                continue;
            };
            for (share, count) in histogram {
                if let Some(bin) = into.get_mut(share) {
                    saturated |= *count > *bin;
                    *bin = bin.saturating_sub(*count);
                    if *bin == 0 {
                        into.remove(share);
                    }
                } else {
                    saturated = true;
                }
            }
            if into.is_empty() {
                signal.keys.remove(key);
            }
        }
        if saturated {
            tracing::warn!(
                event = "run_signal_read_not_in_run",
                run_id = %signal.run_id,
                calls = read.calls,
            );
        }
        signal
    }
}

/// The lower median of `values`, given as (value, count) pairs in increasing
/// value order with `total` counts.
fn median(values: impl Iterator<Item = (u64, u64)>, total: u64) -> u64 {
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

/// The lower median of a histogram with `total` counts.
fn median_of(histogram: &BTreeMap<u64, u64>, total: u64) -> u64 {
    median(
        histogram.iter().map(|(&value, &count)| (value, count)),
        total,
    )
}

#[cfg(test)]
mod tests {
    use dna_kernel::read_evidence::{CallEvidence, EvidenceProfile};

    use super::*;

    type TestResult = std::result::Result<(), Box<dyn std::error::Error>>;

    /// The configured parameters: a 5-mer, shares in thousandths.
    const PARAMETERS: RunSignalParameters = RunSignalParameters {
        flank: 2,
        resolution: 1000,
        minimum_reads: 40,
        minimum_calls: 40,
        margin: 0.05,
        mad_factor: 4.0,
    };

    /// A read's window peak heights and its evidence.
    type Read = (Vec<[i32; 4]>, ReadEvidence);

    /// A read of `bases`, each call's primary peak 1000 and its other
    /// channels at `others` (A/C/G/T order, the primary's entry ignored).
    fn read(
        bases: &str,
        others: [i32; 4],
    ) -> std::result::Result<Read, Box<dyn std::error::Error>> {
        let heights = bases
            .chars()
            .map(|base| {
                std::array::from_fn(|channel| {
                    if Nucleotide::ALL[channel].as_char() == base {
                        1000
                    } else {
                        others[channel]
                    }
                })
            })
            .collect::<Vec<_>>();
        let evidence = bases
            .chars()
            .map(|base| {
                CallEvidence::scored(
                    base,
                    EvidenceProfile::from_corrected_amplitudes([1.0, 1.0, 1.0, 1.0]),
                )
            })
            .collect::<Vec<_>>();
        let count = evidence.len();
        Ok((heights, ReadEvidence::new(evidence, 0..count, Vec::new())?))
    }

    // SRS-CALL-019
    #[test]
    fn counts_every_selected_call_under_each_other_channel() -> TestResult {
        let (calls, evidence) = read("CGTNCGTAC", [0, 250, 0, 100])?;
        let selected = contributions(&calls, &evidence, &PARAMETERS);
        // Calls 2-5 have the N in their context, and the first and last two
        // have no context: call 6 remains, a T.
        assert_eq!(selected.calls, 1);
        let key = |second| Key {
            context: "CGTAC".into(),
            base: 'T',
            second,
        };
        assert_eq!(selected.keys[&key('A')][&0], 1);
        // 250 against 1000 and itself: 200 thousandths.
        assert_eq!(selected.keys[&key('C')][&200], 1);
        assert_eq!(selected.keys[&key('G')][&0], 1);
        assert!(!selected.keys.contains_key(&key('T')));
        // With one call of context on each side and shares in hundredths,
        // calls 1, 5, and 7 too have a context, and the share is 20.
        let coarse = RunSignalParameters {
            flank: 1,
            resolution: 100,
            ..PARAMETERS
        };
        let selected = contributions(&calls, &evidence, &coarse);
        assert_eq!(selected.calls, 4);
        let key = Key {
            context: "GTA".into(),
            base: 'T',
            second: 'C',
        };
        assert_eq!(selected.keys[&key][&20], 1);
        Ok(())
    }

    // SRS-CALL-019
    #[test]
    fn taking_a_read_out_equals_building_without_it() -> TestResult {
        let (calls, evidence) = read("ACGTACGTAC", [100, 0, 50, 0])?;
        let first = contributions(&calls, &evidence, &PARAMETERS);
        let (calls, evidence) = read("ACGTACGTAC", [300, 0, 0, 20])?;
        let second = contributions(&calls, &evidence, &PARAMETERS);
        let traces = BTreeSet::from(["a".to_owned(), "b".to_owned()]);
        let both = RunSignal::build(
            "run".into(),
            traces.clone(),
            [first.clone(), second.clone()],
            &PARAMETERS,
        );
        let alone = RunSignal::build("run".into(), traces, [first], &PARAMETERS);
        assert_eq!(both.without(&second), alone);
        Ok(())
    }

    // SRS-CALL-020
    #[test]
    fn a_level_needs_enough_reads_and_calls_and_stands_above_the_spread() {
        let key = Key {
            context: "ACGTA".into(),
            base: 'G',
            second: 'C',
        };
        // 40 calls: 30 at 20 thousandths, 10 at 40.
        let histogram = Histogram::from([(20, 30), (40, 10)]);
        let signal = RunSignal {
            run_id: "run".into(),
            selection: SELECTION,
            flank: 2,
            resolution: 1000,
            traces: BTreeSet::new(),
            reads: 50,
            calls: 40,
            keys: BTreeMap::from([(key.clone(), histogram)]),
        };
        let parameters = PARAMETERS;
        // Median 0.020, MAD 0: the margin sets the spread.
        assert!((signal.level(&key, &parameters).unwrap_or(0.0) - 0.07).abs() < 1e-9);
        let narrow = RunSignalParameters {
            margin: 0.0001,
            ..parameters
        };
        assert!((signal.level(&key, &narrow).unwrap_or(0.0) - 0.0201).abs() < 1e-9);
        let few = RunSignalParameters {
            minimum_calls: 41,
            ..parameters
        };
        assert_eq!(signal.level(&key, &few), None);
        let small = RunSignalParameters {
            minimum_reads: 51,
            ..parameters
        };
        assert_eq!(signal.level(&key, &small), None);
    }

    // SRS-CALL-020
    #[test]
    fn the_threshold_is_the_floor_raised_to_a_usable_run_level() {
        let key = Key {
            context: "ACGTA".into(),
            base: 'G',
            second: 'C',
        };
        // 40 calls at 300 thousandths: level 0.35 above a 0.21 floor.
        let signal = RunSignal {
            run_id: "run".into(),
            selection: SELECTION,
            flank: 2,
            resolution: 1000,
            traces: BTreeSet::new(),
            reads: 50,
            calls: 40,
            keys: BTreeMap::from([(key.clone(), Histogram::from([(300, 40)]))]),
        };
        let config = SignalProcessingConfig {
            window_size_bases: 10,
            minimum_primary_snr: 3.0,
            minimum_noisy_windows: 2,
            minimum_noise_sigma: 1.0,
            second_floor: 0.21,
            run_signal: PARAMETERS,
            second_cause: crate::signal_processing::SecondCauseParameters {
                noise_ratio: 3.0,
                tail_offset: 0.4,
                blob_width: 1.8,
                width_calls: 6,
            },
            decompose_unexplained_quantile: 0.99,
            spike: crate::signal_processing::SpikeParameters {
                height: 3.0,
                span: 1.5,
                level_calls: 6,
            },
        };
        let floor = Threshold {
            share: 0.21,
            calls: None,
        };
        assert_eq!(threshold(None, Some(&signal), &config), floor);
        assert_eq!(threshold(Some(&key), None, &config), floor);
        let raised = threshold(Some(&key), Some(&signal), &config);
        assert!((raised.share - 0.35).abs() < 1e-9);
        assert_eq!(raised.calls, Some(40));
        let low = RunSignal {
            keys: BTreeMap::from([(key.clone(), Histogram::from([(20, 40)]))]),
            ..signal.clone()
        };
        assert_eq!(
            threshold(Some(&key), Some(&low), &config),
            Threshold {
                share: 0.21,
                calls: Some(40)
            }
        );
        let sparse = RunSignalParameters {
            minimum_calls: 41,
            ..config.run_signal
        };
        let config = SignalProcessingConfig {
            run_signal: sparse,
            ..config
        };
        assert_eq!(threshold(Some(&key), Some(&signal), &config), floor);
    }

    // SRS-CALL-019
    #[test]
    fn a_run_signal_of_another_flank_or_resolution_is_refused() {
        let signal = RunSignal::build("run".into(), BTreeSet::new(), [], &PARAMETERS);
        assert!(signal.check(&PARAMETERS).is_ok());
        for other in [
            RunSignalParameters {
                flank: 3,
                ..PARAMETERS
            },
            RunSignalParameters {
                resolution: 100,
                ..PARAMETERS
            },
        ] {
            assert!(matches!(
                signal.check(&other),
                Err(RunSignalError::Parameters { .. })
            ));
        }
    }

    // SRS-CALL-020, SRS-CALL-021
    #[test]
    fn a_second_template_within_its_run_level_is_context() {
        let key = Key {
            context: "ACGTA".into(),
            base: 'G',
            second: 'C',
        };
        // 40 calls at 300 thousandths: run level 0.35.
        let signal = RunSignal::build("run".into(), BTreeSet::new(), [], &PARAMETERS);
        let signal = RunSignal {
            reads: 50,
            calls: 40,
            keys: BTreeMap::from([(key.clone(), Histogram::from([(300, 40)]))]),
            ..signal
        };
        let config = SignalProcessingConfig {
            window_size_bases: 10,
            minimum_primary_snr: 3.0,
            minimum_noisy_windows: 2,
            minimum_noise_sigma: 1.0,
            second_floor: 0.20,
            run_signal: PARAMETERS,
            second_cause: crate::signal_processing::SecondCauseParameters {
                noise_ratio: 3.0,
                tail_offset: 0.4,
                blob_width: 1.8,
                width_calls: 6,
            },
            decompose_unexplained_quantile: 0.99,
            spike: crate::signal_processing::SpikeParameters {
                height: 3.0,
                span: 1.5,
                level_calls: 6,
            },
        };
        let peak = |share, cause| SecondPeak {
            base: 'C',
            share,
            stands_out: true,
            cause,
        };
        let template = peak(0.30, SecondCause::SecondTemplate);
        assert_eq!(
            refine(&template, Some(&key), Some(&signal), &config),
            SecondCause::Context
        );
        // Above the run level, or without a run signal, it stays.
        let above = peak(0.40, SecondCause::SecondTemplate);
        assert_eq!(
            refine(&above, Some(&key), Some(&signal), &config),
            SecondCause::SecondTemplate
        );
        assert_eq!(
            refine(&template, Some(&key), None, &config),
            SecondCause::SecondTemplate
        );
        // A cause the shape gives stands.
        let tail = peak(0.30, SecondCause::NeighbourTail);
        assert_eq!(
            refine(&tail, Some(&key), Some(&signal), &config),
            SecondCause::NeighbourTail
        );
    }
}
