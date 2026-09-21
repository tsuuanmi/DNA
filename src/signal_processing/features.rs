//! Rolling sample-domain baseline, noise, and peak-SNR features.

use crate::config::DNAProcessingConfig;
use crate::error::{Error, Result};
use crate::model::basecalls::BaseCalls;
use crate::model::signal::DNAWindow;
use crate::model::trace::Chromatogram;

use super::statistics;

/// Computes one feature record per full-width, stride-one base window.
pub(super) fn calculate(
    trace: &Chromatogram,
    calls: &BaseCalls,
    config: &DNAProcessingConfig,
) -> Result<Vec<DNAWindow>> {
    if calls.len() < config.window_size_bases {
        return Err(Error::DNAProcessing(format!(
            "{} calls are fewer than window_size_bases {}",
            calls.len(),
            config.window_size_bases
        )));
    }

    let mut windows = Vec::with_capacity(calls.len() - config.window_size_bases + 1);
    for call_start in 0..=calls.len() - config.window_size_bases {
        let call_end = call_start + config.window_size_bases;
        let selected = &calls.calls[call_start..call_end];
        let sample_start = selected[0].window_start_0based;
        let sample_end = selected[selected.len() - 1].window_end_0based_exclusive;
        if sample_start >= sample_end || sample_end > trace.sample_count() {
            return Err(Error::DNAProcessing(format!(
                "invalid sample interval {sample_start}..{sample_end} for call window {call_start}..{call_end}"
            )));
        }

        let mut baselines = [0.0; 4];
        let mut noise_sigmas = [0.0; 4];
        for channel in 0..4 {
            let statistics =
                statistics::estimate(&trace.channels[channel][sample_start..sample_end])?;
            baselines[channel] = statistics.baseline;
            noise_sigmas[channel] = statistics.noise_sigma;
        }

        let mut minimum_primary_snr = f64::INFINITY;
        let mut maximum_secondary_snr = 0.0_f64;
        for call in selected {
            let mut ranked = [0_usize, 1, 2, 3];
            let corrected = std::array::from_fn::<_, 4, _>(|channel| {
                statistics::corrected_amplitude(call.peaks[channel].height, baselines[channel])
            });
            ranked.sort_by(|left, right| {
                corrected[*right]
                    .total_cmp(&corrected[*left])
                    .then_with(|| left.cmp(right))
            });
            let primary = ranked[0];
            let secondary = ranked[1];
            minimum_primary_snr =
                minimum_primary_snr.min(statistics::snr(corrected[primary], noise_sigmas[primary]));
            maximum_secondary_snr = maximum_secondary_snr.max(statistics::snr(
                corrected[secondary],
                noise_sigmas[secondary],
            ));
        }

        let minimum_primary_snr = statistics::round_metric(minimum_primary_snr);
        let maximum_secondary_snr = statistics::round_metric(maximum_secondary_snr);
        windows.push(DNAWindow {
            call_start_0based: call_start,
            call_end_0based_exclusive: call_end,
            sample_start_0based: sample_start,
            sample_end_0based_exclusive: sample_end,
            minimum_primary_snr,
            maximum_secondary_snr,
            candidate_noisy: minimum_primary_snr < config.minimum_primary_snr,
        });
    }
    Ok(windows)
}

#[cfg(test)]
mod tests {
    use crate::model::basecalls::{BaseCall, ChannelPeak, PeakSource};
    use crate::model::nucleotide::Nucleotide;
    use crate::model::trace::VendorEvidence;

    use super::*;

    fn evidence(count: usize) -> (Chromatogram, BaseCalls) {
        let sample_count = count * 4 + 4;
        let mut channels: [Vec<i32>; 4] = std::array::from_fn(|_| vec![0; sample_count]);
        let mut records = Vec::with_capacity(count);
        for index in 0..count {
            let ploc = 2 + index * 4;
            let primary_channel = index % 4;
            let primary = Nucleotide::ALL[primary_channel];
            channels[primary_channel][ploc] = 1_000;
            records.push(BaseCall {
                index_0based: index,
                ploc_0based: ploc,
                window_start_0based: ploc - 2,
                window_end_0based_exclusive: ploc + 2,
                peaks: std::array::from_fn(|channel| ChannelPeak {
                    base: Nucleotide::ALL[channel],
                    height: if channel == primary_channel { 1_000 } else { 0 },
                    position_0based: ploc,
                    source: if channel == primary_channel {
                        PeakSource::LocalMaximum
                    } else {
                        PeakSource::PlocFallback
                    },
                }),
                primary_peak_evidence: Some(crate::model::basecalls::PrimaryPeakEvidence {
                    position_0based: ploc,
                    channel_heights: std::array::from_fn(|channel| channels[channel][ploc]),
                }),
                primary: primary.as_char(),
                ambiguity: primary.as_char(),
                qualifying_channels: vec![primary],
                vendor_agrees: None,
            });
        }
        let sequence = records.iter().map(|call| call.primary).collect::<String>();
        (
            Chromatogram {
                source_name: "synthetic.ab1".into(),
                source_sha256: String::new(),
                channels,
                base_locations: records.iter().map(|call| call.ploc_0based).collect(),
                vendor: VendorEvidence::default(),
            },
            BaseCalls {
                calls: records,
                primary_sequence: sequence,
            },
        )
    }

    #[test]
    fn computes_full_windows_and_keeps_exact_threshold_clean() -> Result<()> {
        let (trace, calls) = evidence(7);
        let windows = calculate(
            &trace,
            &calls,
            &DNAProcessingConfig {
                window_size_bases: 5,
                minimum_primary_snr: 1_000.0,
                minimum_noisy_windows: 2,
            },
        )?;
        assert_eq!(windows.len(), 3);
        assert_eq!(windows[0].call_start_0based, 0);
        assert_eq!(windows[0].call_end_0based_exclusive, 5);
        assert_eq!(windows[0].minimum_primary_snr, 1_000.0);
        assert_eq!(windows[0].maximum_secondary_snr, 0.0);
        assert!(!windows[0].candidate_noisy);
        assert!(
            windows
                .iter()
                .all(|window| window.minimum_primary_snr.is_finite())
        );
        Ok(())
    }

    #[test]
    fn rejects_a_trace_shorter_than_the_configured_window() {
        let (trace, calls) = evidence(4);
        assert!(matches!(
            calculate(
                &trace,
                &calls,
                &DNAProcessingConfig {
                    window_size_bases: 5,
                    minimum_primary_snr: 3.0,
                    minimum_noisy_windows: 2,
                },
            ),
            Err(Error::DNAProcessing(message)) if message.contains("fewer than")
        ));
    }
}
