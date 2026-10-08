//! Sample-level validation, deterministic read ordering, and evidence assembly.

use std::collections::BTreeSet;

use crate::config::SampleReconciliationConfig;
use crate::error::{Result, SampleError};
use crate::model::read_observation::ReadObservation;
use crate::model::sample_evidence::{
    SampleEvidence, SampleReadAlignmentEvidence, SampleReadEvidence,
};

use super::{coverage, loci, overlap, variants};

/// Aggregates independently processed reads without using filenames or pair labels as merge keys.
pub(crate) fn aggregate(
    reads: &[ReadObservation],
    config: &SampleReconciliationConfig,
) -> Result<SampleEvidence> {
    let ordered = validated_ordered_reads(reads)?;
    let first = ordered.first().ok_or(SampleError::NoReads)?;
    let reference_sha256 = first.reference_sha256.clone();
    let configuration_sha256 = first.configuration_sha256.clone();

    let read_evidence: Vec<SampleReadEvidence> = ordered
        .iter()
        .map(|read| SampleReadEvidence {
            input_name: read.input_name.clone(),
            input_sha256: read.input_sha256.clone(),
            integrity: read.signal.integrity.clone(),
            alignment: SampleReadAlignmentEvidence {
                orientation: read.alignment.orientation,
                callable_bases: read.alignment.metrics.callable_columns,
                identity: read.alignment.metrics.callable_identity,
                gap_opens: read.alignment.metrics.gap_opens,
                unresolved_bases: read.alignment.metrics.unresolved_query_bases,
                reference_segments: read.alignment.reference_segments.clone(),
                wraps_origin: read.alignment.wraps_origin,
            },
        })
        .collect();

    let coverage = coverage::summarize(&read_evidence)?;

    Ok(SampleEvidence {
        reference_sha256,
        configuration_sha256,
        reads: read_evidence,
        coverage,
        overlaps: overlap::assess(&ordered, config)?,
        locus_differences: loci::aggregate(&ordered)?,
        variants: variants::aggregate(&ordered)?,
    })
}

pub(crate) fn validated_ordered_reads(reads: &[ReadObservation]) -> Result<Vec<&ReadObservation>> {
    let first = reads.first().ok_or(SampleError::NoReads)?;

    let mut identities = BTreeSet::new();
    for read in reads {
        if read.reference_sha256 != first.reference_sha256 {
            return Err(SampleError::MixedReference.into());
        }
        if read.configuration_sha256 != first.configuration_sha256 {
            return Err(SampleError::MixedConfiguration.into());
        }
        if !identities.insert(read.input_sha256.as_str()) {
            return Err(SampleError::DuplicateTrace.into());
        }
    }

    let mut ordered = reads.iter().collect::<Vec<_>>();
    ordered.sort_by(|left, right| left.input_sha256.cmp(&right.input_sha256));
    Ok(ordered)
}

#[cfg(test)]
mod tests;
