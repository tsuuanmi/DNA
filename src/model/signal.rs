//! Observation-only signal-quality windows, locus evidence, and noisy regions.

use crate::model::locus_evidence::LocusEvidence;

/// Observation-only structural and amplitude integrity evidence for one trace.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct SangerIntegrity {
    pub(crate) locus_count: usize,
    pub(crate) vendor_primary_count: Option<usize>,
    pub(crate) vendor_quality_count: Option<usize>,
    pub(crate) minimum_locus_spacing: Option<usize>,
    pub(crate) median_locus_spacing: Option<f64>,
    pub(crate) maximum_locus_spacing: Option<usize>,
    pub(crate) clipped_channel_samples: usize,
    pub(crate) maximum_to_median_event_signal_ratio: Option<f64>,
}

impl SangerIntegrity {
    /// Number of present vendor series whose length differs from the canonical locus series.
    pub(crate) fn vendor_length_mismatch_count(&self) -> usize {
        [self.vendor_primary_count, self.vendor_quality_count]
            .into_iter()
            .flatten()
            .filter(|&count| count != self.locus_count)
            .count()
    }
}

/// DNA-quality features for one rolling base-call window.
#[derive(Debug, Clone)]
pub(crate) struct DNAWindow {
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
pub(crate) struct NoisyRegion {
    pub(crate) call_start_0based: usize,
    pub(crate) call_end_0based_exclusive: usize,
    pub(crate) sample_start_0based: usize,
    pub(crate) sample_end_0based_exclusive: usize,
    pub(crate) minimum_primary_snr: f64,
}

/// Complete observation-only signal analysis.
#[derive(Debug, Clone)]
pub(crate) struct DNAAnalysis {
    pub(crate) integrity: SangerIntegrity,
    pub(crate) loci: Vec<LocusEvidence>,
    pub(crate) windows: Vec<DNAWindow>,
    pub(crate) noisy_regions: Vec<NoisyRegion>,
}

impl DNAAnalysis {
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
