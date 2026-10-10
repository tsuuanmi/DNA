//! Observation-only Sanger integrity evidence derived from canonical loci and analyzed channels.

use crate::model::locus_evidence::LocusEvidence;
use crate::model::sanger::Chromatogram;
use crate::model::signal::SangerIntegrity;
use dna_kernel::error::{Result, SignalError};

use super::statistics;

/// Derives structural and signal-scale evidence without changing calls or alignment.
pub(super) fn assess(trace: &Chromatogram, loci: &[LocusEvidence]) -> Result<SangerIntegrity> {
    if loci.len() != trace.call_count() {
        return Err(SignalError::LocusCountMismatch {
            expected: trace.call_count(),
            found: loci.len(),
        }
        .into());
    }

    let spacings = trace
        .locus_positions
        .windows(2)
        .map(|pair| pair[1] - pair[0])
        .collect::<Vec<_>>();
    let (minimum_locus_spacing, median_locus_spacing, maximum_locus_spacing) =
        if spacings.is_empty() {
            (None, None, None)
        } else {
            let minimum = spacings.iter().copied().min();
            let maximum = spacings.iter().copied().max();
            let median = Some(statistics::median_usize(&spacings)?);
            (minimum, median, maximum)
        };

    let clipped_channel_samples = trace
        .channels
        .iter()
        .flat_map(|channel| channel.iter())
        .filter(|&&value| value == i32::from(i16::MIN) || value == i32::from(i16::MAX))
        .count();

    let positive_event_signal = loci
        .iter()
        .map(|locus| locus.corrected_amplitudes.iter().sum::<f64>())
        .filter(|total| *total > 0.0)
        .collect::<Vec<_>>();
    let maximum_to_median_event_signal_ratio = if positive_event_signal.is_empty() {
        None
    } else {
        let maximum = positive_event_signal
            .iter()
            .copied()
            .fold(0.0_f64, f64::max);
        let median = statistics::median_f64(&positive_event_signal)?;
        (median > 0.0).then(|| statistics::round_metric(maximum / median))
    };

    Ok(SangerIntegrity {
        locus_count: trace.call_count(),
        duplicate_locus_count: trace.duplicate_loci,
        vendor_primary_count: trace.vendor.primary.as_ref().map(String::len),
        vendor_quality_count: trace.vendor.qualities.as_ref().map(Vec::len),
        minimum_locus_spacing,
        median_locus_spacing,
        maximum_locus_spacing,
        clipped_channel_samples,
        maximum_to_median_event_signal_ratio,
    })
}

#[cfg(test)]
mod tests {
    use crate::model::sanger::VendorEvidence;
    use dna_kernel::read_evidence::EvidenceProfile;

    use super::*;

    fn locus(index: usize, height: i32) -> LocusEvidence {
        let total = f64::from(height);
        LocusEvidence {
            call_index_0based: index,
            locus_position_0based: index * 4 + 2,
            window_start_0based: index * 4,
            window_end_0based_exclusive: index * 4 + 4,
            context_call_start_0based: 0,
            context_call_end_0based_exclusive: 3,
            context_sample_start_0based: 0,
            context_sample_end_0based_exclusive: 12,
            event_position_0based: index * 4 + 2,
            channel_heights: [height, 0, 0, 0],
            channel_baselines: [0.0; 4],
            channel_noise_sigmas: [1.0; 4],
            corrected_amplitudes: [total, 0.0, 0.0, 0.0],
            snrs: [total, 0.0, 0.0, 0.0],
            profile: EvidenceProfile::from_corrected_amplitudes([total, 0.0, 0.0, 0.0]),
        }
    }

    #[test]
    fn preserves_vendor_length_mismatch_as_integrity_evidence() -> Result<()> {
        let trace = Chromatogram {
            source_name: "synthetic.ab1".into(),
            source_sha256: String::new(),
            channels: std::array::from_fn(|_| vec![0; 12]),
            locus_positions: vec![2, 6, 10],
            duplicate_loci: 0,
            vendor: VendorEvidence {
                primary: Some("AAAA".into()),
                qualities: Some(vec![40; 2]),
            },
        };

        let integrity = assess(&trace, &[locus(0, 10), locus(1, 20), locus(2, 40)])?;
        assert_eq!(integrity.locus_count, 3);
        assert_eq!(integrity.vendor_primary_count, Some(4));
        assert_eq!(integrity.vendor_quality_count, Some(2));
        assert_eq!(integrity.vendor_length_mismatch_count(), 2);
        assert_eq!(integrity.minimum_locus_spacing, Some(4));
        assert_eq!(integrity.median_locus_spacing, Some(4.0));
        assert_eq!(integrity.maximum_locus_spacing, Some(4));
        assert_eq!(integrity.maximum_to_median_event_signal_ratio, Some(2.0));
        Ok(())
    }

    #[test]
    fn reports_exact_channel_clipping_without_classifying_other_amplitude_outliers() -> Result<()> {
        let mut channels = std::array::from_fn(|_| vec![0; 12]);
        channels[0][3] = i32::from(i16::MAX);
        channels[1][7] = i32::from(i16::MIN);
        let trace = Chromatogram {
            source_name: "synthetic.ab1".into(),
            source_sha256: String::new(),
            channels,
            locus_positions: vec![2, 6, 10],
            duplicate_loci: 0,
            vendor: VendorEvidence::default(),
        };

        let integrity = assess(&trace, &[locus(0, 10), locus(1, 10), locus(2, 1000)])?;
        assert_eq!(integrity.clipped_channel_samples, 2);
        assert_eq!(integrity.maximum_to_median_event_signal_ratio, Some(100.0));
        Ok(())
    }
}
