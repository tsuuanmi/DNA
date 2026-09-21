//! Structural nucleotide-contribution eligibility for retained sample-locus observations.

use crate::model::sample_evidence::{CallDNAEvidence, LocusState, NucleotideContribution};

/// Classifies whether one retained locus observation has nucleotide-profile evidence.
pub(super) fn classify(
    state: LocusState,
    signal: Option<CallDNAEvidence>,
) -> NucleotideContribution {
    if state == LocusState::Deletion {
        NucleotideContribution::DeletionEvent
    } else if signal.and_then(|evidence| evidence.profile).is_some() {
        NucleotideContribution::Eligible
    } else {
        NucleotideContribution::MissingProfile
    }
}
