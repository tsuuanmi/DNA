//! Sample-level validation, deterministic read ordering, and evidence assembly.

use std::collections::BTreeSet;

use crate::error::{Result, SampleError};
use crate::model::called_read::CalledRead;
use crate::model::sample_evidence::{
    RejectedSampleRead, SampleEvidence, SampleReadAlignmentEvidence, SampleReadEvidence,
};
use crate::sample::SampleReconciliationConfig;

use super::{coverage, loci, overlap, variants};

/// Aggregates independently processed reads without using filenames or pair
/// labels as merge keys; rejected reads are recorded but contribute nothing.
pub(crate) fn aggregate(
    reads: &[&CalledRead],
    rejected: &[RejectedSampleRead],
    config: &SampleReconciliationConfig,
) -> Result<SampleEvidence> {
    if reads.is_empty() && !rejected.is_empty() {
        return Err(SampleError::NoAdmittedReads {
            rejected: rejected.len(),
        }
        .into());
    }
    let ordered = validated_ordered_reads(reads)?;
    let first = ordered.first().ok_or(SampleError::NoReads)?;
    let mut identities = ordered
        .iter()
        .map(|read| read.input_sha256.as_str())
        .collect::<BTreeSet<_>>();
    if !rejected
        .iter()
        .all(|read| identities.insert(read.input_sha256.as_str()))
    {
        return Err(SampleError::DuplicateTrace.into());
    }
    let mut rejected_reads = rejected.to_vec();
    rejected_reads.sort_by(|left, right| left.input_sha256.cmp(&right.input_sha256));
    let reference_sha256 = first.reference_sha256.clone();
    let configuration_sha256 = first.configuration_sha256.clone();

    let read_evidence: Vec<SampleReadEvidence> = ordered
        .iter()
        .map(|read| SampleReadEvidence {
            input_name: read.input_name.clone(),
            input_sha256: read.input_sha256.clone(),
            alignment: SampleReadAlignmentEvidence {
                orientation: read.alignment.orientation,
                callable_bases: read.alignment.metrics.callable_columns,
                identity: read.alignment.metrics.callable_identity,
                gap_opens: read.alignment.metrics.gap_opens,
                unresolved_bases: read.alignment.metrics.unresolved_query_bases,
                masked_bases: read.alignment.metrics.masked_query_bases,
                reference_segments: read.alignment.reference_segments.clone(),
                callable_reference_segments: read.alignment.callable_segments.clone(),
                wraps_origin: read.alignment.wraps_origin,
            },
        })
        .collect();

    let coverage = coverage::summarize(&read_evidence)?;

    Ok(SampleEvidence {
        reference_sha256,
        configuration_sha256,
        reads: read_evidence,
        rejected_reads,
        coverage,
        overlaps: overlap::assess(&ordered, config)?,
        locus_differences: loci::aggregate(&ordered)?,
        variants: variants::aggregate(&ordered)?,
    })
}

pub(crate) fn validated_ordered_reads<'a>(reads: &[&'a CalledRead]) -> Result<Vec<&'a CalledRead>> {
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

    let mut ordered = reads.to_vec();
    ordered.sort_by(|left, right| left.input_sha256.cmp(&right.input_sha256));
    Ok(ordered)
}

#[cfg(test)]
mod tests;
