//! Observation-only rolling signal-quality analysis.

mod config;

pub use config::{RawSignalProcessingConfig, SignalProcessingConfig};

mod features;
mod integrity;
mod locus_evidence;
mod regions;
mod statistics;

pub(crate) use statistics::round_metric;

use crate::error::{Result, SignalError};
use crate::model::basecalls::BaseCalls;
use crate::model::sanger::Chromatogram;
use crate::model::signal::SignalAnalysis;

/// Fails a read with fewer calls than one window (SRS-SIG-007).
///
/// This is the one length rule of signal processing: every rolling window and
/// locus context needs a full window of calls, and the window of 5 to 10 calls
/// also covers the two loci that locus-window geometry needs. Read processing
/// applies it before basecalling.
///
/// # Errors
///
/// Returns `SignalError::TooFewCalls` for a read shorter than one window.
pub(crate) fn check_length(call_count: usize, config: &SignalProcessingConfig) -> Result<()> {
    if call_count < config.window_size_bases {
        return Err(SignalError::TooFewCalls {
            calls: call_count,
            window: config.window_size_bases,
        }
        .into());
    }
    Ok(())
}

/// Calculates rolling SNR features, basecall-independent locus evidence, and merged noisy regions.
///
/// Both feature windows and locus contexts estimate channel statistics over
/// runs of whole locus windows; each distinct span is estimated once. The
/// read must pass [`check_length`]; a shorter one fails typed.
pub(crate) fn analyze(
    trace: &Chromatogram,
    calls: &BaseCalls,
    config: &SignalProcessingConfig,
) -> Result<SignalAnalysis> {
    let mut estimates = statistics::SpanEstimates::new(trace);
    let windows = features::calculate(calls, config, &mut estimates)?;
    let loci = locus_evidence::calculate(config, &mut estimates)?;
    let noisy_regions = regions::merge(&windows, config.minimum_noisy_windows);
    let integrity = integrity::assess(trace, &loci)?;
    Ok(SignalAnalysis {
        integrity,
        loci,
        windows,
        noisy_regions,
    })
}

#[cfg(test)]
mod tests {
    use crate::error::{CallabilityError, SangerError as Error};
    use crate::model::signal::ShortRead;

    use super::*;

    fn config(window_size_bases: usize) -> SignalProcessingConfig {
        SignalProcessingConfig {
            window_size_bases,
            minimum_primary_snr: 3.0,
            minimum_noisy_windows: 2,
        }
    }

    /// A read of exactly one window passes; one call fewer fails with the
    /// error that `ShortRead::of` reads back (SRS-SIG-007).
    #[test]
    fn fails_a_read_shorter_than_one_window() {
        assert!(check_length(10, &config(10)).is_ok());
        let error = check_length(9, &config(10)).err();
        assert!(matches!(
            error,
            Some(Error::Signal(SignalError::TooFewCalls {
                calls: 9,
                window: 10
            }))
        ));
        let short = error.as_ref().and_then(ShortRead::of);
        assert_eq!(
            short.map(|short| (short.calls, short.minimum_calls)),
            Some((9, 10))
        );
        let other = Error::from(CallabilityError::Inconsistent("other"));
        assert_eq!(ShortRead::of(&other), None);
    }
}
