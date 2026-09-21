//! Reference-coordinate production sample-locus evidence.

use std::collections::{BTreeMap, BTreeSet};

use crate::error::{Error, Result};
use crate::model::alignment::AlignmentColumn;
use crate::model::read_observation::ReadObservation;
use crate::model::sample_evidence::{
    LocusState, LocusSupportTopology, SampleLocusEvidence, SampleLocusObservation,
};

use super::{call_evidence, contribution, nucleotide_support};

struct LocusBuilder {
    reference_base: char,
    observations: Vec<SampleLocusObservation>,
}

pub(crate) fn aggregate(reads: &[&ReadObservation]) -> Result<Vec<SampleLocusEvidence>> {
    let mut loci: BTreeMap<usize, LocusBuilder> = BTreeMap::new();

    for read in reads {
        for column in &read.alignment.columns {
            let Some(reference_index_0based) = column.reference_index_0based else {
                continue;
            };
            if column.reference_base == '-' {
                return Err(Error::Sample(
                    "reference-coordinate locus cannot contain an insertion column".into(),
                ));
            }
            let state = classify(column);
            if state == LocusState::Reference {
                continue;
            }
            let position_1based = reference_index_0based
                .checked_add(1)
                .ok_or_else(|| Error::Sample("reference coordinate overflow".into()))?;
            let entry = loci.entry(position_1based).or_insert_with(|| LocusBuilder {
                reference_base: column.reference_base,
                observations: Vec::new(),
            });
            if entry.reference_base != column.reference_base {
                return Err(Error::Sample(format!(
                    "reference base disagrees at position {position_1based}"
                )));
            }
        }
    }

    if loci.is_empty() {
        return Ok(Vec::new());
    }

    for (read_index, read) in reads.iter().enumerate() {
        let mut seen = BTreeSet::new();
        for column in &read.alignment.columns {
            let Some(reference_index_0based) = column.reference_index_0based else {
                continue;
            };
            let position_1based = reference_index_0based
                .checked_add(1)
                .ok_or_else(|| Error::Sample("reference coordinate overflow".into()))?;
            let Some(entry) = loci.get_mut(&position_1based) else {
                continue;
            };
            if !seen.insert(position_1based) {
                return Err(Error::Sample(format!(
                    "read {} contains duplicate reference coordinate {position_1based}",
                    read.input_sha256
                )));
            }
            if entry.reference_base != column.reference_base {
                return Err(Error::Sample(format!(
                    "reference base disagrees at position {position_1based}"
                )));
            }
            entry
                .observations
                .push(observation(read_index, read, column)?);
        }
    }

    loci.into_iter()
        .map(|(position_1based, built)| {
            let support_topology = support_topology(&built.observations, reads)?;
            let nucleotide_support = nucleotide_support::aggregate(&built.observations, reads)?;
            Ok(SampleLocusEvidence {
                position_1based,
                reference_base: built.reference_base,
                support_topology,
                nucleotide_support,
                observations: built.observations,
            })
        })
        .collect()
}

fn support_topology(
    observations: &[SampleLocusObservation],
    reads: &[&ReadObservation],
) -> Result<LocusSupportTopology> {
    let mut topology = LocusSupportTopology {
        reads: observations.len(),
        forward_reads: 0,
        reverse_reads: 0,
        reference_reads: 0,
        alternate_reads: 0,
        unresolved_reads: 0,
        deletion_reads: 0,
        profile_reads: 0,
        profile_forward_reads: 0,
        profile_reverse_reads: 0,
    };

    for observation in observations {
        let read = reads.get(observation.read_index).ok_or_else(|| {
            Error::Sample(format!(
                "locus observation references missing read {}",
                observation.read_index
            ))
        })?;
        let has_profile = observation
            .signal
            .as_ref()
            .and_then(|signal| signal.profile)
            .is_some();
        match read.alignment.orientation {
            crate::model::alignment::Orientation::Forward => {
                topology.forward_reads += 1;
                if has_profile {
                    topology.profile_forward_reads += 1;
                }
            }
            crate::model::alignment::Orientation::Reverse => {
                topology.reverse_reads += 1;
                if has_profile {
                    topology.profile_reverse_reads += 1;
                }
            }
        }
        if has_profile {
            topology.profile_reads += 1;
        }
        match observation.state {
            LocusState::Reference => topology.reference_reads += 1,
            LocusState::Alternate => topology.alternate_reads += 1,
            LocusState::Unresolved => topology.unresolved_reads += 1,
            LocusState::Deletion => topology.deletion_reads += 1,
        }
    }

    if topology.reads != topology.forward_reads + topology.reverse_reads
        || topology.reads
            != topology.reference_reads
                + topology.alternate_reads
                + topology.unresolved_reads
                + topology.deletion_reads
        || topology.profile_reads != topology.profile_forward_reads + topology.profile_reverse_reads
        || topology.profile_reads > topology.reads - topology.deletion_reads
    {
        return Err(Error::Sample(
            "locus support topology counts are inconsistent".into(),
        ));
    }

    Ok(topology)
}

fn observation(
    read_index: usize,
    read: &ReadObservation,
    column: &AlignmentColumn,
) -> Result<SampleLocusObservation> {
    let state = classify(column);
    if state == LocusState::Deletion {
        let signal = None;
        return Ok(SampleLocusObservation {
            read_index,
            state,
            base: None,
            quality: None,
            signal,
            nucleotide_contribution: contribution::classify(state, signal),
        });
    }

    let call_index_0based = column
        .original_call_index_0based
        .ok_or_else(|| Error::Sample("aligned query base lacks original call index".into()))?;
    let quality = read
        .quality
        .per_call
        .get(call_index_0based)
        .filter(|quality| quality.index_0based == call_index_0based)
        .ok_or_else(|| Error::Sample("aligned call lacks matching quality evidence".into()))?;

    let signal = Some(call_evidence::for_call(read, call_index_0based)?);
    Ok(SampleLocusObservation {
        read_index,
        state,
        base: Some(column.query_base),
        quality: Some(quality.relative_quality_score),
        signal,
        nucleotide_contribution: contribution::classify(state, signal),
    })
}

fn classify(column: &AlignmentColumn) -> LocusState {
    if column.query_base == '-' {
        LocusState::Deletion
    } else if !is_canonical(column.query_base) || !is_canonical(column.reference_base) {
        LocusState::Unresolved
    } else if column.query_base == column.reference_base {
        LocusState::Reference
    } else {
        LocusState::Alternate
    }
}

const fn is_canonical(base: char) -> bool {
    matches!(base, 'A' | 'C' | 'G' | 'T')
}
