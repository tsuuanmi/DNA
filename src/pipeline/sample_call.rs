//! The call of a sample analysis's consensus, in process (SRS-API-012).
//!
//! It is the public [`SampleAnalysis::call`]: it lives here because it
//! composes the consensus FASTA of `report` and the scientific step of the
//! `call` pipeline, which depend on `sample_analysis`.

use crate::error::Result;
use crate::input::sequence;
use crate::report;
use crate::sample_analysis::{ConsensusCall, NotationCall, ReferenceSegment, SampleAnalysis};
use dna_empop::export::{self, Record};
use dna_kernel::variant::Variant;

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
    /// Returns [`Error`](crate::error::Error) when the sample identifier
    /// is invalid or a stage of `call` fails, as the command would; for
    /// example when the consensus has no segment.
    pub fn call(&self, sample_id: &str) -> Result<ConsensusCall> {
        validate_sample_id(sample_id)?;
        let loaded = &self.loaded;
        let sequences = sequence::from_fasta(&report::consensus_segments_fasta(
            sample_id,
            &self.consensus,
        ))?;
        let Called { reads, notation } = call::call_consensus(
            &self.consensus,
            &sequences,
            &loaded.config,
            &loaded.profile,
            &loaded.reference,
        )?;
        let calls = notation
            .as_ref()
            .map(|notation| report::notation_calls(&loaded.reference, notation))
            .transpose()?;
        let empop = calls
            .as_ref()
            .filter(|_| loaded.profile.empop)
            .map(|calls| {
                let records = reads
                    .iter()
                    .map(|read| {
                        let segments =
                            |segments: &[dna_core::model::alignment::ReferenceSegment]| {
                                segments
                                    .iter()
                                    .map(|segment| {
                                        segment.start_0based..segment.end_0based_exclusive
                                    })
                                    .collect()
                            };
                        Record {
                            reference: segments(&read.alignment.reference_segments),
                            callable: segments(&read.alignment.callable_segments),
                            excluded: read
                                .variants
                                .observed
                                .iter()
                                .filter(|observed| !observed.exclusion_reasons.is_empty())
                                .map(|observed| Variant::from(&observed.variant))
                                .collect(),
                        }
                    })
                    .collect::<Vec<_>>();
                let calls = calls
                    .keys()
                    .map(|call| export::Call {
                        position: call.position,
                        ordinal: call.ordinal,
                        text: &call.text,
                    })
                    .collect::<Vec<_>>();
                export::export(&loaded.profile.regions, &records, &calls)
            });
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
            empop,
        })
    }
}
