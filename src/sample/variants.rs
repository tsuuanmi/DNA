//! Deterministic aggregation of normalized variant observations across reads.

use std::collections::{BTreeMap, BTreeSet};

use crate::error::{Result, SampleError};
use crate::model::alignment::Orientation;
use crate::model::read_observation::ReadObservation;
use crate::model::reference_call;
use crate::model::sample_evidence::{
    VariantCallEvidence, VariantEvidence, VariantSupport, VariantSupportTopology,
};
use crate::model::variant::{VariantCallMapping, VariantKind};

use super::call_evidence;

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
struct VariantKey {
    position_1based: usize,
    reference: String,
    alternate: String,
    kind: VariantKind,
}

pub(super) fn aggregate(reads: &[&ReadObservation]) -> Result<Vec<VariantEvidence>> {
    let mut variants: BTreeMap<VariantKey, Vec<VariantSupport>> = BTreeMap::new();

    for (read_index, read) in reads.iter().enumerate() {
        let mut seen = BTreeSet::new();
        for observed in &read.variants.observed {
            let variant = &observed.variant;
            let key = VariantKey {
                position_1based: variant.position_1based,
                reference: variant.reference.clone(),
                alternate: variant.alternate.clone(),
                kind: variant.kind,
            };
            if !seen.insert(key.clone()) {
                return Err(SampleError::DuplicateVariant {
                    read: read.input_sha256.clone(),
                }
                .into());
            }
            variants.entry(key).or_default().push(VariantSupport {
                read_index,
                eligible: observed.eligible(),
                exclusion_reasons: observed.exclusion_reasons.clone(),
                calls: variant_calls(read, &variant.calls)?,
            });
        }
    }

    variants
        .into_iter()
        .map(|(key, support)| {
            let support_topology = support_topology(&support, reads)?;
            Ok(VariantEvidence {
                position_1based: key.position_1based,
                reference: key.reference,
                alternate: key.alternate,
                kind: key.kind,
                support_topology,
                support,
            })
        })
        .collect()
}

fn support_topology(
    support: &[VariantSupport],
    reads: &[&ReadObservation],
) -> Result<VariantSupportTopology> {
    let mut topology = VariantSupportTopology {
        reads: support.len(),
        eligible_reads: 0,
        forward_reads: 0,
        reverse_reads: 0,
        eligible_forward_reads: 0,
        eligible_reverse_reads: 0,
    };
    for item in support {
        let read = reads.get(item.read_index).ok_or(SampleError::MissingRead {
            context: "variant support",
            index: item.read_index,
        })?;
        match read.alignment.orientation {
            Orientation::Forward => {
                topology.forward_reads += 1;
                if item.eligible {
                    topology.eligible_forward_reads += 1;
                }
            }
            Orientation::Reverse => {
                topology.reverse_reads += 1;
                if item.eligible {
                    topology.eligible_reverse_reads += 1;
                }
            }
        }
        if item.eligible {
            topology.eligible_reads += 1;
        }
    }
    if topology.reads != topology.forward_reads + topology.reverse_reads
        || topology.eligible_reads
            != topology.eligible_forward_reads + topology.eligible_reverse_reads
    {
        return Err(
            SampleError::Inconsistent("variant support topology counts are inconsistent").into(),
        );
    }
    Ok(topology)
}

fn variant_calls(
    read: &ReadObservation,
    mappings: &[VariantCallMapping],
) -> Result<Vec<VariantCallEvidence>> {
    reference_call::resolve_public_calls(
        &read.calls,
        &read.quality,
        read.alignment.orientation,
        mappings,
    )
    .map_err(SampleError::CallEvidence)?
    .into_iter()
    .map(|call| {
        Ok(VariantCallEvidence {
            role: call.mapping.role,
            base: call.evidence.base,
            peak_heights: call.evidence.peak_heights,
            quality: call.evidence.quality,
            signal: call_evidence::for_call(read, call.mapping.call_index_0based)?,
        })
    })
    .collect()
}
