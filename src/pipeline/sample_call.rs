//! The call of a sample analysis's consensus, in process (SRS-API-012).
//!
//! It is the public [`SampleAnalysis::call`]: it lives here because it
//! composes the `consensus` report, the consensus FASTA, and the `call`
//! pipeline, which depend on `sample_analysis`.

use crate::input::sequence::{self, CallInputs};
use crate::pipeline::plugins as plugin;
use crate::report::{self, CompletedConsensus};
use crate::sample_analysis::{ConsensusCall, NotationCall, ReferenceSegment, SampleAnalysis};
use dna_kernel::error::Result;

use super::call;
use super::sample::validate_sample_id;

impl SampleAnalysis {
    /// Calls the consensus sequences of this analysis exactly as `dna call`
    /// reads them from the consensus FASTA `dna consensus` publishes for
    /// `sample_id`, with the configuration, profile, and reference of the
    /// analysis, without CLI logging or publication side effects.
    ///
    /// # Errors
    ///
    /// Returns [`Error`](dna_kernel::error::Error) when the sample identifier
    /// is invalid or a stage of `call` fails, as the command would; for
    /// example when the consensus has no segment.
    pub fn call(&self, sample_id: &str) -> Result<ConsensusCall> {
        validate_sample_id(sample_id)?;
        let document = report::build_consensus(CompletedConsensus {
            sample_id: sample_id.to_owned(),
            reference: self.inputs.reference.clone(),
            profile: self.profile.clone(),
            evidence: self.evidence.clone(),
            rejected: self.rejected.clone(),
            consensus: self.consensus.clone(),
            plugins: plugin::CONSENSUS,
        })?;
        let reads = sequence::from_fasta(&report::consensus_fasta(&document))?;
        let completed = call::complete(
            sample_id,
            CallInputs {
                config: self.inputs.config.clone(),
                profile: self.inputs.profile.clone(),
                reads,
                reference: self.inputs.reference.clone(),
            },
        )?;
        let notation = completed
            .notation
            .as_ref()
            .map(|notation| report::notation_calls(&completed.reference, notation))
            .transpose()?;
        let variants = report::build_call(completed)?;
        Ok(ConsensusCall {
            callable_segments: variants
                .reads
                .iter()
                .flat_map(|read| &read.alignment.callable_reference_segments)
                .map(|segment| ReferenceSegment {
                    start_0based: segment.start,
                    end_0based_exclusive: segment.end,
                })
                .collect(),
            notation: notation.map(|calls| {
                calls
                    .into_iter()
                    .map(|call| NotationCall {
                        position: call.position,
                        ordinal: call.ordinal,
                        reference: call.reference,
                        alternate: call.alternate,
                        text: call.text,
                    })
                    .collect()
            }),
        })
    }
}
