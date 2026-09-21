//! Reference-oriented projection of call-backed signal evidence.

use crate::error::{Error, Result};
use crate::model::alignment::Orientation;
use crate::model::read_observation::ReadObservation;
use crate::model::sample_evidence::CallDNAEvidence;

/// Resolves one source call's quantitative signal evidence and projects A/C/G/T
/// channels onto the selected reference strand.
pub(super) fn for_call(
    read: &ReadObservation,
    call_index_0based: usize,
) -> Result<CallDNAEvidence> {
    let locus = read
        .signal
        .loci
        .get(call_index_0based)
        .filter(|locus| locus.call_index_0based == call_index_0based)
        .ok_or_else(|| {
            Error::Sample(format!(
                "call index {call_index_0based} lacks matching locus evidence"
            ))
        })?;

    let orientation = read.alignment.orientation;
    let profile = locus.profile.map(|profile| match orientation {
        Orientation::Forward => profile,
        Orientation::Reverse => profile.complemented(),
    });
    let in_noisy_region = read.signal.noisy_regions.iter().any(|region| {
        region.call_start_0based <= call_index_0based
            && call_index_0based < region.call_end_0based_exclusive
    });

    Ok(CallDNAEvidence {
        corrected_amplitudes: orientation.reference_signal_values(locus.corrected_amplitudes),
        snrs: orientation.reference_signal_values(locus.snrs),
        profile,
        in_noisy_region,
    })
}
