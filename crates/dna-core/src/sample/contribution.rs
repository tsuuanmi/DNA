//! Structural nucleotide-contribution eligibility for retained sample-locus observations.

use crate::model::sample_evidence::{LocusState, NucleotideContribution};
use dna_kernel::read_evidence::EvidenceProfile;

/// Classifies whether one retained locus observation has nucleotide-profile evidence.
pub(super) fn classify(
    state: LocusState,
    profile: Option<EvidenceProfile>,
) -> NucleotideContribution {
    if state == LocusState::Deletion {
        NucleotideContribution::DeletionEvent
    } else if state == LocusState::Masked {
        NucleotideContribution::MaskedCall
    } else if profile.is_some() {
        NucleotideContribution::Eligible
    } else {
        NucleotideContribution::MissingProfile
    }
}
