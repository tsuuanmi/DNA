//! Basecall-independent signal evidence at one vendor-defined locus.

use dna_kernel::read_evidence::EvidenceProfile;

/// Immutable signal evidence at one PLOC-defined locus.
///
/// Event refinement and profile construction are derived from analyzed channel
/// values directly. They do not depend on the primary call, ambiguity code,
/// selected basecall peaks, or qualifying-channel set.
#[derive(Debug, Clone)]
pub struct LocusEvidence {
    pub(crate) call_index_0based: usize,
    pub(crate) locus_position_0based: usize,
    pub(crate) window_start_0based: usize,
    pub(crate) window_end_0based_exclusive: usize,
    pub(crate) context_call_start_0based: usize,
    pub(crate) context_call_end_0based_exclusive: usize,
    pub(crate) context_sample_start_0based: usize,
    pub(crate) context_sample_end_0based_exclusive: usize,
    pub(crate) event_position_0based: usize,
    pub(crate) channel_heights: [i32; 4],
    pub(crate) channel_baselines: [f64; 4],
    pub(crate) channel_noise_sigmas: [f64; 4],
    /// Baseline-corrected A/C/G/T amplitudes at the event.
    pub corrected_amplitudes: [f64; 4],
    /// A/C/G/T signal-to-noise ratios at the event.
    pub snrs: [f64; 4],
    /// Normalized A/C/G/T profile; absent without positive signal.
    pub profile: Option<EvidenceProfile>,
    /// A/C/G/T peak widths at the event, in samples: the contiguous samples
    /// around it, within the context, where the channel stays at or above
    /// half its corrected amplitude over its baseline; zero without positive
    /// signal.
    pub(crate) channel_widths: [usize; 4],
}

/// The width of channel `channel`'s peak at call `call`'s event, from the
/// locus that records that call; `None` when none does or the channel has no
/// positive signal there.
pub(crate) fn channel_width(loci: &[LocusEvidence], call: usize, channel: usize) -> Option<usize> {
    loci.get(call)
        .filter(|locus| locus.call_index_0based == call)
        .and_then(|locus| locus.channel_widths.get(channel).copied())
        .filter(|&width| width > 0)
}

#[cfg(test)]
mod tests {
    use super::*;

    // SRS-SIG-009
    #[test]
    fn normalizes_positive_corrected_signal() {
        let Some(profile) = EvidenceProfile::from_corrected_amplitudes([0.0, 40.0, 100.0, 0.0])
        else {
            panic!("positive signal should produce a profile");
        };
        assert_eq!(profile.weights, [0.0, 2.0 / 7.0, 5.0 / 7.0, 0.0]);
    }

    #[test]
    fn complements_channel_weights() {
        let profile = EvidenceProfile {
            weights: [0.1, 0.2, 0.3, 0.4],
        };
        assert_eq!(profile.complemented().weights, [0.4, 0.3, 0.2, 0.1]);
    }

    // SRS-SIG-009
    #[test]
    fn zero_signal_has_no_profile() {
        assert!(EvidenceProfile::from_corrected_amplitudes([0.0; 4]).is_none());
    }
}
