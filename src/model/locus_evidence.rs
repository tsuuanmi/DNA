//! Basecall-independent signal evidence at one vendor-defined locus.

/// Normalized non-negative A/C/G/T evidence derived from corrected channel amplitudes.
///
/// Channel order follows `Nucleotide::ALL`: A, C, G, T. A profile exists only
/// when the locus contains positive corrected signal.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct EvidenceProfile {
    pub(crate) weights: [f64; 4],
}

impl EvidenceProfile {
    pub(crate) fn from_corrected_amplitudes(amplitudes: [f64; 4]) -> Option<Self> {
        let total = amplitudes.iter().sum::<f64>();
        (total > 0.0).then(|| Self {
            weights: amplitudes.map(|amplitude| amplitude / total),
        })
    }

    /// Complements A/C/G/T evidence while preserving its total mass.
    pub(crate) const fn complemented(self) -> Self {
        Self {
            weights: [
                self.weights[3],
                self.weights[2],
                self.weights[1],
                self.weights[0],
            ],
        }
    }
}

/// Immutable signal evidence at one PLOC-defined locus.
///
/// Event refinement and profile construction are derived from analyzed channel
/// values directly. They do not depend on the primary call, ambiguity code,
/// selected basecall peaks, or qualifying-channel set.
#[derive(Debug, Clone)]
pub(crate) struct LocusEvidence {
    pub(crate) call_index_0based: usize,
    pub(crate) ploc_0based: usize,
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
    pub(crate) corrected_amplitudes: [f64; 4],
    pub(crate) snrs: [f64; 4],
    pub(crate) profile: Option<EvidenceProfile>,
}

#[cfg(test)]
mod tests {
    use super::*;

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

    #[test]
    fn zero_signal_has_no_profile() {
        assert!(EvidenceProfile::from_corrected_amplitudes([0.0; 4]).is_none());
    }
}
