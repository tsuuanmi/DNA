//! Basecall-independent signal evidence at each PLOC-defined locus.

use crate::config::DNAProcessingConfig;
use crate::error::{Error, Result};
use crate::locus::{self, LocusWindow};
use crate::model::locus_evidence::{EvidenceProfile, LocusEvidence};
use crate::model::trace::Chromatogram;

use super::statistics;

/// Calculates one immutable signal-evidence record per vendor-defined locus.
///
/// Event refinement, local statistics, and profile construction use analyzed
/// channel evidence plus PLOC geometry directly. They do not consume basecall
/// records, primary/ambiguity calls, selected basecall peaks, or qualifying
/// channels.
pub(super) fn calculate(
    trace: &Chromatogram,
    config: &DNAProcessingConfig,
) -> Result<Vec<LocusEvidence>> {
    let locus_count = trace.call_count();
    if locus_count < config.window_size_bases {
        return Err(Error::DNAProcessing(format!(
            "{locus_count} loci are fewer than window_size_bases {}",
            config.window_size_bases
        )));
    }

    let locus_windows = locus::windows(trace).map_err(Error::DNAProcessing)?;
    if locus_windows.len() != locus_count {
        return Err(Error::DNAProcessing(format!(
            "expected {locus_count} locus windows, found {}",
            locus_windows.len()
        )));
    }

    let mut evidence = Vec::with_capacity(locus_count);
    for (locus_index, (&ploc, &locus_window)) in trace
        .base_locations
        .iter()
        .zip(locus_windows.iter())
        .enumerate()
    {
        let context_start = context_start(locus_index, locus_count, config.window_size_bases);
        let context_end = context_start + config.window_size_bases;
        let context_sample_start = locus_windows[context_start].start;
        let context_sample_end = locus_windows[context_end - 1].end;
        let (channel_baselines, channel_noise_sigmas) =
            local_statistics(trace, context_sample_start, context_sample_end)?;
        let event_position = select_event_position(trace, locus_window, ploc, channel_baselines)?;
        let channel_heights =
            std::array::from_fn(|channel| trace.channels[channel][event_position]);
        let corrected_amplitudes = std::array::from_fn(|channel| {
            statistics::corrected_amplitude(channel_heights[channel], channel_baselines[channel])
        });
        let snrs = std::array::from_fn(|channel| {
            statistics::round_metric(statistics::snr(
                corrected_amplitudes[channel],
                channel_noise_sigmas[channel],
            ))
        });

        let record = LocusEvidence {
            call_index_0based: locus_index,
            ploc_0based: ploc,
            window_start_0based: locus_window.start,
            window_end_0based_exclusive: locus_window.end,
            context_call_start_0based: context_start,
            context_call_end_0based_exclusive: context_end,
            context_sample_start_0based: context_sample_start,
            context_sample_end_0based_exclusive: context_sample_end,
            event_position_0based: event_position,
            channel_heights,
            channel_baselines,
            channel_noise_sigmas,
            corrected_amplitudes,
            snrs,
            profile: EvidenceProfile::from_corrected_amplitudes(corrected_amplitudes),
        };
        validate_evidence(&record, config.window_size_bases)?;
        evidence.push(record);
    }
    Ok(evidence)
}

fn validate_evidence(evidence: &LocusEvidence, context_width: usize) -> Result<()> {
    let valid_coordinates = evidence.window_start_0based < evidence.window_end_0based_exclusive
        && evidence.window_start_0based <= evidence.ploc_0based
        && evidence.ploc_0based < evidence.window_end_0based_exclusive
        && evidence.window_start_0based <= evidence.event_position_0based
        && evidence.event_position_0based < evidence.window_end_0based_exclusive
        && evidence.context_call_start_0based <= evidence.call_index_0based
        && evidence.call_index_0based < evidence.context_call_end_0based_exclusive
        && evidence.context_call_end_0based_exclusive - evidence.context_call_start_0based
            == context_width
        && evidence.context_sample_start_0based <= evidence.window_start_0based
        && evidence.window_end_0based_exclusive <= evidence.context_sample_end_0based_exclusive;
    let valid_metrics = evidence
        .channel_baselines
        .iter()
        .all(|value| value.is_finite())
        && evidence
            .channel_noise_sigmas
            .iter()
            .all(|value| value.is_finite() && *value >= 1.0)
        && evidence
            .corrected_amplitudes
            .iter()
            .all(|value| value.is_finite() && *value >= 0.0)
        && evidence
            .snrs
            .iter()
            .all(|value| value.is_finite() && *value >= 0.0)
        && evidence
            .channel_heights
            .iter()
            .zip(evidence.channel_baselines)
            .zip(evidence.corrected_amplitudes)
            .all(|((&height, baseline), corrected)| {
                corrected
                    .total_cmp(&statistics::corrected_amplitude(height, baseline))
                    .is_eq()
            });
    let total = evidence.corrected_amplitudes.iter().sum::<f64>();
    let valid_profile = match (total > 0.0, evidence.profile) {
        (false, None) => true,
        (true, Some(profile)) => profile
            .weights
            .iter()
            .zip(evidence.corrected_amplitudes)
            .all(|(&weight, amplitude)| weight.total_cmp(&(amplitude / total)).is_eq()),
        _ => false,
    };
    if valid_coordinates && valid_metrics && valid_profile {
        Ok(())
    } else {
        Err(Error::DNAProcessing(format!(
            "inconsistent locus evidence at call {}",
            evidence.call_index_0based
        )))
    }
}

fn local_statistics(
    trace: &Chromatogram,
    sample_start: usize,
    sample_end: usize,
) -> Result<([f64; 4], [f64; 4])> {
    if sample_start >= sample_end || sample_end > trace.sample_count() {
        return Err(Error::DNAProcessing(format!(
            "invalid locus context sample interval {sample_start}..{sample_end}"
        )));
    }
    let mut baselines = [0.0; 4];
    let mut noise_sigmas = [0.0; 4];
    for channel in 0..4 {
        let local = statistics::estimate(&trace.channels[channel][sample_start..sample_end])?;
        baselines[channel] = local.baseline;
        noise_sigmas[channel] = local.noise_sigma;
    }
    Ok((baselines, noise_sigmas))
}

/// Refines the locus event without consulting a basecall verdict.
///
/// Candidate events are positive local maxima of total non-negative
/// baseline-corrected A/C/G/T amplitude. The nearest candidate to PLOC wins;
/// equal-distance candidates prefer greater total signal, then the lower sample
/// coordinate. If the locus window has no positive total-signal local maximum,
/// the validated PLOC sample is used directly.
fn select_event_position(
    trace: &Chromatogram,
    window: LocusWindow,
    ploc: usize,
    baselines: [f64; 4],
) -> Result<usize> {
    if window.start >= window.end
        || window.end > trace.sample_count()
        || ploc < window.start
        || ploc >= window.end
    {
        return Err(Error::DNAProcessing(format!(
            "invalid locus event window {}..{} at PLOC {ploc}",
            window.start, window.end
        )));
    }

    let search_start = window.start.max(1);
    let search_end = window.end.min(trace.sample_count().saturating_sub(1));
    let mut best: Option<(usize, f64)> = None;

    for position in search_start..search_end {
        let left = corrected_total(trace, position - 1, baselines);
        let current = corrected_total(trace, position, baselines);
        let right = corrected_total(trace, position + 1, baselines);
        let local = (left <= current && current > right) || (left < current && current >= right);
        if !local || current <= 0.0 {
            continue;
        }

        let replace = best.is_none_or(|(best_position, best_total)| {
            let distance = position.abs_diff(ploc);
            let best_distance = best_position.abs_diff(ploc);
            distance < best_distance
                || (distance == best_distance
                    && (current.total_cmp(&best_total).is_gt()
                        || (current.total_cmp(&best_total).is_eq() && position < best_position)))
        });
        if replace {
            best = Some((position, current));
        }
    }

    Ok(best.map_or(ploc, |(position, _)| position))
}

fn corrected_total(trace: &Chromatogram, position: usize, baselines: [f64; 4]) -> f64 {
    (0..4)
        .map(|channel| {
            statistics::corrected_amplitude(trace.channels[channel][position], baselines[channel])
        })
        .sum()
}

fn context_start(locus_index: usize, locus_count: usize, window_size_bases: usize) -> usize {
    let left_context = (window_size_bases - 1) / 2;
    locus_index
        .saturating_sub(left_context)
        .min(locus_count - window_size_bases)
}

#[cfg(test)]
mod tests {
    use crate::config::DNAProcessingConfig;
    use crate::model::trace::{Chromatogram, VendorEvidence};

    use super::*;

    fn config(window_size_bases: usize) -> DNAProcessingConfig {
        DNAProcessingConfig {
            window_size_bases,
            minimum_primary_snr: 3.0,
            minimum_noisy_windows: 2,
        }
    }

    fn trace() -> Chromatogram {
        let mut channels = std::array::from_fn(|_| vec![0; 12]);
        for (index, position) in [1_usize, 3, 5, 7, 9].into_iter().enumerate() {
            channels[index % 4][position] = 100;
        }
        Chromatogram {
            source_name: "synthetic.ab1".into(),
            source_sha256: String::new(),
            channels,
            base_locations: vec![1, 3, 5, 7, 9],
            vendor: VendorEvidence::default(),
        }
    }

    #[test]
    fn selects_centered_and_edge_contexts_deterministically() {
        assert_eq!(
            (0..10)
                .map(|index| context_start(index, 10, 5))
                .collect::<Vec<_>>(),
            vec![0, 0, 0, 1, 2, 3, 4, 5, 5, 5]
        );
        assert_eq!(
            (0..10)
                .map(|index| context_start(index, 10, 6))
                .collect::<Vec<_>>(),
            vec![0, 0, 0, 1, 2, 3, 4, 4, 4, 4]
        );
    }

    #[test]
    fn produces_profile_without_any_basecall_input() -> Result<()> {
        let evidence = calculate(&trace(), &config(5))?;
        assert_eq!(evidence.len(), 5);
        let locus = &evidence[2];
        assert_eq!(locus.call_index_0based, 2);
        assert_eq!(locus.ploc_0based, 5);
        assert_eq!(locus.event_position_0based, 5);
        assert_eq!(locus.channel_heights, [0, 0, 100, 0]);
        assert_eq!(locus.corrected_amplitudes, [0.0, 0.0, 100.0, 0.0]);
        assert_eq!(locus.snrs, [0.0, 0.0, 100.0, 0.0]);
        assert_eq!(
            locus.profile.as_ref().map(|profile| profile.weights),
            Some([0.0, 0.0, 1.0, 0.0])
        );
        Ok(())
    }

    #[test]
    fn profile_retains_mixed_channel_mass_without_threshold_membership() -> Result<()> {
        let mut trace = trace();
        trace.channels[1][5] = 40;
        let evidence = calculate(&trace, &config(5))?;
        let profile = evidence[2]
            .profile
            .as_ref()
            .ok_or_else(|| Error::DNAProcessing("missing profile".into()))?;
        assert_eq!(evidence[2].event_position_0based, 5);
        assert_eq!(evidence[2].corrected_amplitudes, [0.0, 40.0, 100.0, 0.0]);
        assert_eq!(profile.weights, [0.0, 2.0 / 7.0, 5.0 / 7.0, 0.0]);
        Ok(())
    }

    #[test]
    fn event_refinement_prefers_nearest_local_event_over_stronger_neighbor() -> Result<()> {
        let mut channels = std::array::from_fn(|_| vec![0; 40]);
        channels[0][15] = 220;
        channels[2][18] = 120;
        let trace = Chromatogram {
            source_name: "synthetic.ab1".into(),
            source_sha256: String::new(),
            channels,
            base_locations: vec![2, 10, 18, 26, 34],
            vendor: VendorEvidence::default(),
        };

        let evidence = calculate(&trace, &config(5))?;
        assert_eq!(evidence[2].ploc_0based, 18);
        assert_eq!(evidence[2].event_position_0based, 18);
        assert_eq!(evidence[2].channel_heights, [0, 0, 120, 0]);
        Ok(())
    }

    #[test]
    fn event_refinement_breaks_equal_distance_by_signal_then_coordinate() -> Result<()> {
        let mut channels = std::array::from_fn(|_| vec![0; 20]);
        channels[0][5] = 100;
        channels[2][7] = 120;
        let trace = Chromatogram {
            source_name: "synthetic.ab1".into(),
            source_sha256: String::new(),
            channels,
            base_locations: vec![2, 6, 10, 14, 18],
            vendor: VendorEvidence::default(),
        };
        let window = LocusWindow { start: 3, end: 9 };

        assert_eq!(select_event_position(&trace, window, 6, [0.0; 4])?, 7);

        let mut tied = trace.clone();
        tied.channels[2][7] = 100;
        assert_eq!(select_event_position(&tied, window, 6, [0.0; 4])?, 5);
        Ok(())
    }

    #[test]
    fn event_refinement_falls_back_to_ploc_without_positive_local_maximum() -> Result<()> {
        let trace = Chromatogram {
            source_name: "synthetic.ab1".into(),
            source_sha256: String::new(),
            channels: std::array::from_fn(|_| vec![0; 12]),
            base_locations: vec![1, 3, 5, 7, 9],
            vendor: VendorEvidence::default(),
        };

        assert_eq!(
            select_event_position(&trace, LocusWindow { start: 4, end: 7 }, 5, [0.0; 4])?,
            5
        );
        Ok(())
    }

    #[test]
    fn zero_corrected_signal_has_no_profile() -> Result<()> {
        let mut trace = trace();
        for channel in &mut trace.channels {
            channel.fill(0);
        }
        let evidence = calculate(&trace, &config(5))?;
        assert!(evidence.iter().all(|locus| locus.profile.is_none()));
        Ok(())
    }
}
