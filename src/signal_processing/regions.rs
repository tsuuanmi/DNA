//! Deterministic union of candidate-noisy rolling windows.

use crate::model::signal::{DNAWindow, NoisyRegion};

/// Merges candidate-noisy runs that contain enough consecutive windows.
pub(super) fn merge(windows: &[DNAWindow], minimum_noisy_windows: usize) -> Vec<NoisyRegion> {
    let mut regions = Vec::new();
    let mut current = None;
    let mut run_length = 0;

    for window in windows {
        if !window.candidate_noisy {
            append_if_supported(
                &mut regions,
                current.take(),
                run_length,
                minimum_noisy_windows,
            );
            run_length = 0;
            continue;
        }

        if let Some(region) = current.as_mut()
            && window.call_start_0based > region.call_end_0based_exclusive
        {
            append_if_supported(
                &mut regions,
                current.take(),
                run_length,
                minimum_noisy_windows,
            );
            run_length = 0;
        }

        if let Some(region) = current.as_mut() {
            region.call_end_0based_exclusive = region
                .call_end_0based_exclusive
                .max(window.call_end_0based_exclusive);
            region.sample_end_0based_exclusive = region
                .sample_end_0based_exclusive
                .max(window.sample_end_0based_exclusive);
            region.minimum_primary_snr = region.minimum_primary_snr.min(window.minimum_primary_snr);
        } else {
            current = Some(NoisyRegion {
                call_start_0based: window.call_start_0based,
                call_end_0based_exclusive: window.call_end_0based_exclusive,
                sample_start_0based: window.sample_start_0based,
                sample_end_0based_exclusive: window.sample_end_0based_exclusive,
                minimum_primary_snr: window.minimum_primary_snr,
            });
        }
        run_length += 1;
    }

    append_if_supported(&mut regions, current, run_length, minimum_noisy_windows);
    regions
}

fn append_if_supported(
    regions: &mut Vec<NoisyRegion>,
    region: Option<NoisyRegion>,
    run_length: usize,
    minimum_noisy_windows: usize,
) {
    if run_length >= minimum_noisy_windows
        && let Some(region) = region
    {
        regions.push(region);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn window(start: usize, end: usize, noisy: bool) -> DNAWindow {
        DNAWindow {
            call_start_0based: start,
            call_end_0based_exclusive: end,
            sample_start_0based: start * 4,
            sample_end_0based_exclusive: end * 4,
            minimum_primary_snr: start as f64,
            maximum_secondary_snr: 0.0,
            candidate_noisy: noisy,
        }
    }

    #[test]
    fn requires_two_candidate_windows_for_a_region() {
        let regions = merge(&[window(0, 5, true), window(1, 6, false)], 2);
        assert!(regions.is_empty());

        let regions = merge(&[window(0, 5, true), window(1, 6, true)], 2);
        assert_eq!(regions.len(), 1);
        assert_eq!(regions[0].call_start_0based, 0);
        assert_eq!(regions[0].call_end_0based_exclusive, 6);
    }

    #[test]
    fn merges_overlapping_and_adjacent_windows() {
        let regions = merge(
            &[window(0, 5, true), window(1, 6, true), window(6, 11, true)],
            2,
        );
        assert_eq!(regions.len(), 1);
        assert_eq!(regions[0].call_start_0based, 0);
        assert_eq!(regions[0].call_end_0based_exclusive, 11);
        assert_eq!(regions[0].sample_end_0based_exclusive, 44);
        assert_eq!(regions[0].minimum_primary_snr, 0.0);
    }

    #[test]
    fn keeps_regions_separate_across_a_clean_gap() {
        let regions = merge(
            &[
                window(0, 5, true),
                window(1, 6, true),
                window(2, 7, false),
                window(7, 12, true),
                window(8, 13, true),
            ],
            2,
        );
        assert_eq!(regions.len(), 2);
        assert_eq!(regions[0].call_end_0based_exclusive, 6);
        assert_eq!(regions[1].call_start_0based, 7);
    }
}
