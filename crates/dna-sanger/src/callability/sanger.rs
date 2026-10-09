//! The Sanger adapter: per-position evidence from the chromatogram's loci,
//! basecall-independent locus evidence, and the primary calls.
//!
//! This is the only file of the module that knows the Sanger types; it owns no
//! thresholds.

use crate::model::basecalls::BaseCalls;
use crate::model::callability::PositionEvidence;
use crate::model::sanger::Chromatogram;
use crate::model::signal::SignalAnalysis;
use dna_kernel::error::{CallabilityError, Result};
use dna_kernel::model::nucleotide::Nucleotide;

/// Builds one evidence record per call locus.
pub(super) fn evidence(
    trace: &Chromatogram,
    calls: &BaseCalls,
    signal: &SignalAnalysis,
) -> Result<Vec<PositionEvidence>> {
    for found in [signal.loci.len(), trace.locus_positions.len()] {
        if found != calls.len() {
            return Err(CallabilityError::PositionCountMismatch {
                expected: calls.len(),
                found,
            }
            .into());
        }
    }
    Ok(calls
        .calls
        .iter()
        .zip(&signal.loci)
        .zip(&trace.locus_positions)
        .map(|((call, locus), &coordinate)| PositionEvidence {
            amplitudes: locus.corrected_amplitudes,
            coordinate,
            primary: Nucleotide::from_char(call.primary).map(Nucleotide::channel_index),
        })
        .collect())
}
