//! The second-channel signal a run's reads share (PROP-0009, SRS-CALL-019).
//!
//! A Sanger run shows the same shadow of a second base in the same sequence
//! context in most of its reads, at a level that changes from run to run. The
//! run signal records, per key, how often each second-channel share occurs
//! among a run's clean calls, so a call can later be judged against its own
//! run rather than a fixed share.
//!
//! The key is the read-strand 5-mer of the read's own primary calls around a
//! call, the call's base, and one other channel: no reference position, so
//! the modality stays reference-free (ADR-0069). Every selected call counts
//! under each of its three other channels, so a key holds every call of its
//! context, not only those whose strongest other channel it is. The shares are kept as sparse
//! histograms in thousandths rather than quantiles, so the estimator can
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

use crate::error::Result;
use crate::model::basecalls::BaseCalls;
use crate::model::sanger::Chromatogram;
use crate::model::signal::ShortRead;
use crate::read_processing::evidence::{PEAK_BELOW_MINIMUM, RELATIVE_QUALITY_NOT_ABOVE_THRESHOLD};
use crate::read_processing::{self, SangerConfig};
use crate::signal_processing::{RunSignalParameters, SignalProcessingConfig};
use dna_kernel::model::nucleotide::Nucleotide;

/// Version of the call selection [`contributions`] applies.
pub const SELECTION: u32 = 2;

/// Calls on each side of a call that form its context.
const FLANK: usize = 2;

/// What a call's share is counted under: the read-strand 5-mer of the read's
/// primary calls centred on the call, the call's base, and one of its other
/// channels.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Key {
    /// The read's primary calls from two before to two after the call.
    pub context: [char; 5],
    /// The call's base.
    pub base: char,
    /// The other channel whose share is counted.
    pub second: char,
}

impl Key {
    /// The key of call `index` of a read and its channel `second`: `None`
    /// within two calls of either end or without a canonical context.
    #[must_use]
    pub fn of(evidence: &ReadEvidence, index: usize, second: char) -> Option<Self> {
        let calls = evidence.calls();
        let context = calls.get(index.checked_sub(FLANK)?..=index + FLANK)?;
        let context: [char; 5] = std::array::from_fn(|offset| context[offset].base);
        context
            .iter()
            .all(|&base| is_canonical(base))
            .then(|| Self {
                context,
                base: context[FLANK],
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

/// Counts of calls by a channel's share in thousandths (0 to 1000): its peak
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
/// sides, two deep. A read's first and last two calls have no context. Each
/// counts under each of its three other channels, with that channel's share
/// of its peak and the primary peak, from the read's `calls`.
#[must_use]
pub fn contributions(calls: &BaseCalls, evidence: &ReadEvidence) -> Contributions {
    let verdicts = evidence.calls();
    let mut selected = Contributions::default();
    for index in FLANK..verdicts.len().saturating_sub(FLANK) {
        let call = &verdicts[index];
        let context = std::array::from_fn(|offset| verdicts[index - FLANK + offset].base);
        let Some(peaks) = calls.calls.get(index).map(|call| call.peaks) else {
            continue;
        };
        let Some(primary) = Nucleotide::from_char(call.base) else {
            continue;
        };
        let height = peaks[primary.channel_index()].height;
        if !evidence.informative_call(index)
            || call.placement != Placement::Scored
            || call.identity.edge
            || call.identity.withheld.contains(PEAK_BELOW_MINIMUM)
            || call
                .identity
                .withheld
                .contains(RELATIVE_QUALITY_NOT_ABOVE_THRESHOLD)
            || !context.iter().all(|&base| is_canonical(base))
            || height <= 0
        {
            continue;
        }
        for other in Nucleotide::ALL {
            if other == primary {
                continue;
            }
            let peak = peaks[other.channel_index()].height.max(0);
            let total = f64::from(height) + f64::from(peak);
            // A share in [0, 1], so the rounded thousandths fit.
            #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
            let share = (1000.0 * f64::from(peak) / total)
                .round()
                .clamp(0.0, 1000.0) as u16;
            *selected
                .keys
                .entry(Key {
                    context,
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
pub fn read(
    trace: &Chromatogram,
    config: &SangerConfig,
    edge_margin: usize,
) -> Result<Option<Contributions>> {
    let prepared = match read_processing::prepare(trace, config) {
        Ok(prepared) => prepared,
        Err(error) if ShortRead::of(&error).is_some() => return Ok(None),
        Err(error) => return Err(error),
    };
    if prepared.rejection(config).is_some() {
        return Ok(None);
    }
    let processed = read_processing::finish(trace, prepared, config, edge_margin)?;
    let evidence = read_processing::read_evidence(
        &processed.calls,
        &processed.signal,
        &processed.callability,
        &processed.quality,
        &config.sanger_evidence,
        edge_margin,
    )?;
    Ok(Some(contributions(&processed.calls, &evidence)))
}

/// A run's second-channel signal: per key, the shares of the selected calls
/// of every read of the run.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RunSignal {
    /// The run's identifier.
    pub run_id: String,
    /// The selection the contributions were made with ([`SELECTION`]).
    pub selection: u32,
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
    /// contributions of those reads that could be processed.
    #[must_use]
    pub fn build(
        run_id: String,
        traces: BTreeSet<String>,
        reads: impl IntoIterator<Item = Contributions>,
    ) -> Self {
        let mut signal = Self {
            run_id,
            selection: SELECTION,
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
            let into = self.keys.entry(*key).or_default();
            for (share, count) in histogram {
                *into.entry(*share).or_default() += count;
            }
        }
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
        let median = median as f64 / 1000.0;
        let spread = (parameters.mad_factor * mad as f64 / 1000.0).max(parameters.margin);
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
    use crate::model::basecalls::{BaseCall, ChannelPeak, PeakSource};

    type TestResult = std::result::Result<(), Box<dyn std::error::Error>>;

    /// A read of `bases`, each call's primary peak 1000 and its other
    /// channels at `others` (A/C/G/T order, the primary's entry ignored).
    fn read(
        bases: &str,
        others: [i32; 4],
    ) -> std::result::Result<(BaseCalls, ReadEvidence), Box<dyn std::error::Error>> {
        let calls = bases
            .chars()
            .enumerate()
            .map(|(index, base)| BaseCall {
                index_0based: index,
                locus_position_0based: index * 4,
                window_start_0based: index * 4,
                window_end_0based_exclusive: index * 4 + 1,
                peaks: std::array::from_fn(|channel| ChannelPeak {
                    base: Nucleotide::ALL[channel],
                    height: if Nucleotide::ALL[channel].as_char() == base {
                        1000
                    } else {
                        others[channel]
                    },
                    position_0based: index * 4,
                    source: PeakSource::LocalMaximum,
                }),
                primary_peak_evidence: None,
                primary: base,
                ambiguity: base,
                qualifying_channels: Vec::new(),
                vendor_agrees: None,
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
        Ok((
            BaseCalls {
                primary_sequence: bases.into(),
                calls,
            },
            ReadEvidence::new(evidence, 0..count, Vec::new())?,
        ))
    }

    #[test]
    fn counts_every_selected_call_under_each_other_channel() -> TestResult {
        let (calls, evidence) = read("CGTNCGTAC", [0, 250, 0, 100])?;
        let selected = contributions(&calls, &evidence);
        // Calls 2-5 have the N in their context, and the first and last two
        // have no context: call 6 remains, a T.
        assert_eq!(selected.calls, 1);
        let key = |second| Key {
            context: ['C', 'G', 'T', 'A', 'C'],
            base: 'T',
            second,
        };
        assert_eq!(selected.keys[&key('A')][&0], 1);
        // 250 against 1000 and itself: 200 thousandths.
        assert_eq!(selected.keys[&key('C')][&200], 1);
        assert_eq!(selected.keys[&key('G')][&0], 1);
        assert!(!selected.keys.contains_key(&key('T')));
        Ok(())
    }

    #[test]
    fn taking_a_read_out_equals_building_without_it() -> TestResult {
        let (calls, evidence) = read("ACGTACGTAC", [100, 0, 50, 0])?;
        let first = contributions(&calls, &evidence);
        let (calls, evidence) = read("ACGTACGTAC", [300, 0, 0, 20])?;
        let second = contributions(&calls, &evidence);
        let traces = BTreeSet::from(["a".to_owned(), "b".to_owned()]);
        let both = RunSignal::build(
            "run".into(),
            traces.clone(),
            [first.clone(), second.clone()],
        );
        let alone = RunSignal::build("run".into(), traces, [first]);
        assert_eq!(both.without(&second), alone);
        Ok(())
    }

    #[test]
    fn a_level_needs_enough_reads_and_calls_and_stands_above_the_spread() {
        let key = Key {
            context: ['A', 'C', 'G', 'T', 'A'],
            base: 'G',
            second: 'C',
        };
        // 40 calls: 30 at 20 thousandths, 10 at 40.
        let histogram = Histogram::from([(20, 30), (40, 10)]);
        let signal = RunSignal {
            run_id: "run".into(),
            selection: SELECTION,
            traces: BTreeSet::new(),
            reads: 50,
            calls: 40,
            keys: BTreeMap::from([(key, histogram)]),
        };
        let parameters = RunSignalParameters {
            minimum_reads: 40,
            minimum_calls: 40,
            margin: 0.05,
            mad_factor: 4.0,
        };
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

    #[test]
    fn the_threshold_is_the_floor_raised_to_a_usable_run_level() {
        let key = Key {
            context: ['A', 'C', 'G', 'T', 'A'],
            base: 'G',
            second: 'C',
        };
        // 40 calls at 300 thousandths: level 0.35 above a 0.21 floor.
        let signal = RunSignal {
            run_id: "run".into(),
            selection: SELECTION,
            traces: BTreeSet::new(),
            reads: 50,
            calls: 40,
            keys: BTreeMap::from([(key, Histogram::from([(300, 40)]))]),
        };
        let config = SignalProcessingConfig {
            window_size_bases: 10,
            minimum_primary_snr: 3.0,
            minimum_noisy_windows: 2,
            second_floor: 0.21,
            run_signal: RunSignalParameters {
                minimum_reads: 40,
                minimum_calls: 40,
                margin: 0.05,
                mad_factor: 4.0,
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
            keys: BTreeMap::from([(key, Histogram::from([(20, 40)]))]),
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
}
