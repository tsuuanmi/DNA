//! The call of a sample analysis's consensus, in process (SRS-API-012).
//!
//! It is the public [`SampleAnalysis::call`]: it lives here because it
//! composes the consensus FASTA of `report` and the scientific step of the
//! `call` pipeline, which depend on `sample_analysis`.

use crate::input::sequence;
use crate::report;
use crate::sample_analysis::{ConsensusCall, NotationCall, ReferenceSegment, SampleAnalysis};
use dna_kernel::error::Result;

use super::call::{self, Called};
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
        let loaded = &self.loaded;
        let sequences = sequence::from_fasta(&report::consensus_segments_fasta(
            sample_id,
            &self.consensus,
        ))?;
        let Called { reads, notation } = call::call_sequences(
            &sequences,
            &loaded.config,
            &loaded.profile,
            &loaded.reference,
        )?;
        let calls = notation
            .as_ref()
            .map(|notation| report::notation_calls(&loaded.reference, notation))
            .transpose()?;
        Ok(ConsensusCall {
            callable_segments: reads
                .iter()
                .flat_map(|read| &read.alignment.callable_segments)
                .map(|segment| ReferenceSegment {
                    start_0based: segment.start_0based,
                    end_0based_exclusive: segment.end_0based_exclusive,
                })
                .collect(),
            notation: calls.map(|calls| {
                calls
                    .into_iter()
                    .map(|(call, reads)| NotationCall {
                        position: call.position,
                        ordinal: call.ordinal,
                        reference: call.reference,
                        alternate: call.alternate,
                        text: call.text,
                        convention: reads.convention.cloned(),
                    })
                    .collect()
            }),
        })
    }
}
