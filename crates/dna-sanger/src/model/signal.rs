//! Observation-only signal-quality windows, locus evidence, and noisy regions.

use crate::model::locus_evidence::LocusEvidence;

/// Observation-only structural and amplitude integrity evidence for one trace.
#[derive(Debug, Clone, PartialEq)]
pub struct SangerIntegrity {
    /// Canonical locus count from `PLOC.2`, after merging repeated positions.
    pub locus_count: usize,
    /// Source peak locations merged into an equal predecessor.
    pub duplicate_locus_count: usize,
    /// Vendor base count, when present.
    pub vendor_primary_count: Option<usize>,
    /// Vendor quality count, when present.
    pub vendor_quality_count: Option<usize>,
    /// Smallest spacing of adjacent loci, with at least two loci.
    pub minimum_locus_spacing: Option<usize>,
    /// Median spacing of adjacent loci, with at least two loci.
    pub median_locus_spacing: Option<f64>,
    /// Largest spacing of adjacent loci, with at least two loci.
    pub maximum_locus_spacing: Option<usize>,
    /// Channel samples at the signed 16-bit limits.
    pub clipped_channel_samples: usize,
    /// Largest to median corrected event signal, when defined.
    pub maximum_to_median_event_signal_ratio: Option<f64>,
}

impl SangerIntegrity {
    /// Number of present vendor series whose length differs from the canonical locus series.
    #[must_use]
    pub fn vendor_length_mismatch_count(&self) -> usize {
        [self.vendor_primary_count, self.vendor_quality_count]
            .into_iter()
            .flatten()
            .filter(|&count| count != self.locus_count)
            .count()
    }
}

/// Signal-quality features for one rolling base-call window.
#[derive(Debug, Clone)]
pub(crate) struct SignalWindow {
    pub(crate) call_start_0based: usize,
    pub(crate) call_end_0based_exclusive: usize,
    pub(crate) sample_start_0based: usize,
    pub(crate) sample_end_0based_exclusive: usize,
    pub(crate) minimum_primary_snr: f64,
    pub(crate) maximum_secondary_snr: f64,
    pub(crate) candidate_noisy: bool,
}

/// Union of overlapping or adjacent candidate-noisy windows.
#[derive(Debug, Clone)]
pub struct NoisyRegion {
    /// First call of the region (0-based, inclusive).
    pub call_start_0based: usize,
    /// Call after the region (0-based, exclusive).
    pub call_end_0based_exclusive: usize,
    /// First trace sample of the region (0-based, inclusive).
    pub sample_start_0based: usize,
    /// Trace sample after the region (0-based, exclusive).
    pub sample_end_0based_exclusive: usize,
    /// Smallest primary-channel SNR of the region's windows.
    pub minimum_primary_snr: f64,
}

/// Complete observation-only signal analysis.
#[derive(Debug, Clone)]
pub struct SignalAnalysis {
    /// Trace-integrity evidence.
    pub integrity: SangerIntegrity,
    pub(crate) loci: Vec<LocusEvidence>,
    pub(crate) windows: Vec<SignalWindow>,
    /// Merged candidate-noisy regions, in call order.
    pub noisy_regions: Vec<NoisyRegion>,
}

impl SignalAnalysis {
    /// Number of rolling windows classified as candidate-noisy.
    pub(crate) fn noisy_window_count(&self) -> usize {
        self.windows
            .iter()
            .filter(|window| window.candidate_noisy)
            .count()
    }

    /// Number of distinct calls covered by merged candidate-noisy regions.
    pub(crate) fn noisy_call_count(&self) -> usize {
        self.noisy_regions
            .iter()
            .map(|region| region.call_end_0based_exclusive - region.call_start_0based)
            .sum()
    }
}
