//! Rolling sample-domain baseline, noise, and peak-SNR features.

use crate::error::{Result, SignalError};
use crate::model::basecalls::BaseCalls;
use crate::model::signal::SignalWindow;
use crate::signal_processing::SignalProcessingConfig;

use super::statistics::{self, SpanEstimates};

/// Computes one feature record per full-width, stride-one base window of the
/// trace that `estimates` holds.
pub(super) fn calculate(
    calls: &BaseCalls,
    config: &SignalProcessingConfig,
    estimates: &mut SpanEstimates<'_>,
) -> Result<Vec<SignalWindow>> {
    let trace = estimates.trace();
    let width = config.window_size_bases;
    let mut windows = Vec::with_capacity((calls.len() + 1).saturating_sub(width));
    for (call_start, selected) in calls.calls.windows(width).enumerate() {
        let call_end = call_start + width;
        let sample_start = selected[0].window_start_0based;
        let sample_end = selected[selected.len() - 1].window_end_0based_exclusive;
        if sample_start >= sample_end || sample_end > trace.sample_count() {
            return Err(SignalError::InvalidSampleInterval {
                sample_start,
                sample_end,
                call_start,
                call_end,
            }
            .into());
        }

        let span = estimates.span(sample_start, sample_end)?;
        let (baselines, noise_sigmas) = (span.baselines, span.noise_sigmas);

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
        windows.push(SignalWindow {
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
    use crate::model::sanger::{Chromatogram, VendorEvidence};
    use dna_kernel::model::nucleotide::Nucleotide;

    use super::*;

    fn evidence(count: usize) -> (Chromatogram, BaseCalls) {
        let sample_count = count * 4 + 4;
        let mut channels: [Vec<i32>; 4] = std::array::from_fn(|_| vec![0; sample_count]);
        let mut records = Vec::with_capacity(count);
        for index in 0..count {
            let locus_position = 2 + index * 4;
            let primary_channel = index % 4;
            let primary = Nucleotide::ALL[primary_channel];
            channels[primary_channel][locus_position] = 1_000;
            records.push(BaseCall {
                index_0based: index,
                locus_position_0based: locus_position,
                window_start_0based: locus_position - 2,
                window_end_0based_exclusive: locus_position + 2,
                peaks: std::array::from_fn(|channel| ChannelPeak {
                    base: Nucleotide::ALL[channel],
                    height: if channel == primary_channel { 1_000 } else { 0 },
                    position_0based: locus_position,
                    source: if channel == primary_channel {
                        PeakSource::LocalMaximum
                    } else {
                        PeakSource::LocusFallback
                    },
                }),
                primary_peak_evidence: Some(crate::model::basecalls::PrimaryPeakEvidence {
                    position_0based: locus_position,
                    channel_heights: std::array::from_fn(|channel| {
                        channels[channel][locus_position]
                    }),
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
                locus_positions: records
                    .iter()
                    .map(|call| call.locus_position_0based)
                    .collect(),
                duplicate_loci: 0,
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
            &calls,
            &SignalProcessingConfig {
                window_size_bases: 5,
                minimum_primary_snr: 1_000.0,
                minimum_noisy_windows: 2,
                minimum_noise_sigma: 1.0,
                second_floor: 0.21,
                run_signal: crate::signal_processing::RunSignalParameters {
                    flank: 2,
                    resolution: 1000,
                    minimum_reads: 40,
                    minimum_calls: 40,
                    margin: 0.05,
                    mad_factor: 4.0,
                },
            },
            &mut SpanEstimates::new(&trace, 1.0),
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
}
