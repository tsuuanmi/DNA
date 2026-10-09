//! Reference-oriented evidence profile of one source call.

use crate::error::{Result, SampleError};
use crate::model::alignment::Orientation;
use crate::model::called_read::CalledRead;
use crate::read_evidence::EvidenceProfile;

/// Resolves one source call's evidence profile and projects its A/C/G/T
/// channels onto the selected reference strand.
pub(super) fn profile(
    read: &CalledRead,
    call_index_0based: usize,
) -> Result<Option<EvidenceProfile>> {
    let call =
        read.evidence
            .calls()
            .get(call_index_0based)
            .ok_or(SampleError::MissingLocusEvidence {
                call: call_index_0based,
            })?;
    Ok(call
        .profile
        .map(|profile| match read.alignment.orientation {
            Orientation::Forward => profile,
            Orientation::Reverse => profile.complemented(),
        }))
}
