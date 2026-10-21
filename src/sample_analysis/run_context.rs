//! A run context in memory (`dna.run_context/v1`, PROP-0014): what a
//! laboratory batch's own samples give the analysis of each of them, and how
//! one sample's consensus reads it.
//!
//! Increment 1 holds the length section: the batch's cross-strand pairs of
//! the two-strand rule's views over the profile's recorded anchored windows
//! (SRS-CONS-014). A sample analysed with a run context takes its own pairs
//! out when its reads are among the context's (leave-one-out), and builds its
//! consensus with τ rescaled by the rest.

use dna_call::model::called_read::CalledRead;

use super::windows::{recorded_windows, sample_pairs};
use crate::config::Config;
use crate::error::{Result, RunContextError};
pub use crate::model::run_context::{ContextUse, LENGTH_METHOD, RunContext};
use dna_kernel::model::reference::Reference;
use dna_kernel::profile::Profile;

impl RunContext {
    /// Reads and validates a `dna.run_context/v1` document.
    ///
    /// # Errors
    ///
    /// As `dna consensus --run-context` for an unreadable, oversized, or
    /// invalid document.
    pub fn read(path: &std::path::Path) -> Result<Self> {
        crate::input::run_context::load(path)
    }
}

/// The consensus configuration a run context calibrates for one sample whose
/// admitted reads are `reads`, with how it was used.
///
/// # Errors
///
/// [`RunContextError::Reference`] when the context was made against another
/// reference, [`RunContextError::LengthParameters`] when it was made with
/// another τ, other bin edges, or over other windows than the
/// configuration's and the profile's, and
/// [`RunContextError::PartialReads`] when only some of the sample's reads are
/// among the context's.
pub(crate) fn calibrate(
    context: &RunContext,
    profile: &Profile,
    reference: &Reference,
    config: &Config,
    reads: &[&CalledRead],
) -> Result<(dna_consensus::ConsensusConfig, ContextUse)> {
    let consensus = &config.consensus;
    let calibration = consensus.length_calibration();
    if context.reference_sha256 != reference.sequence_sha256 {
        return Err(RunContextError::Reference {
            found: context.reference_sha256.clone(),
            expected: reference.sequence_sha256.clone(),
        }
        .into());
    }
    if !context
        .length
        .made_with(consensus.length_rule(), &calibration.edges)
        || context.windows != recorded_windows(profile)
    {
        return Err(RunContextError::LengthParameters.into());
    }
    let inside = reads
        .iter()
        .filter(|read| context.reads.contains(&read.input_sha256))
        .count();
    if inside != 0 && inside != reads.len() {
        return Err(RunContextError::PartialReads.into());
    }
    let left_out = inside != 0;
    let statistics = if left_out {
        context
            .length
            .clone()
            .without(&sample_pairs(profile, reference, config, reads))
    } else {
        context.length.clone()
    };
    let scale = statistics.scale(calibration);
    Ok((
        consensus.with_length_rule(consensus.length_rule().calibrated(scale)),
        ContextUse {
            run_id: context.run_id.clone(),
            sha256: context.sha256.clone(),
            left_out,
            scale,
        },
    ))
}
