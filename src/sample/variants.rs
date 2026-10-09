//! Deterministic aggregation of normalized variant observations across reads.

use std::collections::{BTreeMap, BTreeSet};

use crate::error::{Result, SampleError};
use crate::model::alignment::Orientation;
use crate::model::read_observation::ReadObservation;
use crate::model::reference_call;
use crate::model::sample_evidence::{
    VariantCallEvidence, VariantEvidence, VariantOpposition, VariantSupport, VariantSupportTopology,
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
            let opposition = opposition(&key, &support, reads);
            Ok(VariantEvidence {
                position_1based: key.position_1based,
                reference: key.reference,
                alternate: key.alternate,
                kind: key.kind,
                support_topology,
                support,
                opposition,
            })
        })
        .collect()
}

/// 0-based reference positions a read must observe callably to oppose a
/// variant: the anchored reference allele, plus the base after an insertion so
/// that both sides of the inserted bases are seen.
fn evidence_span(key: &VariantKey) -> std::ops::Range<usize> {
    let start = key.position_1based.saturating_sub(1);
    let length = key.reference.len() + usize::from(key.kind == VariantKind::Ins);
    start..start.saturating_add(length)
}

/// Reads that callably observe the variant's evidence span without supporting
/// the variant.
fn opposition(
    key: &VariantKey,
    support: &[VariantSupport],
    reads: &[&ReadObservation],
) -> VariantOpposition {
    let supporting = support
        .iter()
        .map(|item| item.read_index)
        .collect::<BTreeSet<_>>();
    let span = evidence_span(key);
    let mut opposition = VariantOpposition {
        read_indices: Vec::new(),
        forward_reads: 0,
        reverse_reads: 0,
    };
    for (index, read) in reads.iter().enumerate() {
        if supporting.contains(&index) {
            continue;
        }
        let callable = span.clone().all(|position| {
            read.alignment.callable_segments.iter().any(|segment| {
                (segment.start_0based..segment.end_0based_exclusive).contains(&position)
            })
        });
        if !callable {
            continue;
        }
        opposition.read_indices.push(index);
        match read.alignment.orientation {
            Orientation::Forward => opposition.forward_reads += 1,
            Orientation::Reverse => opposition.reverse_reads += 1,
        }
    }
    opposition
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

#[cfg(test)]
mod tests {
    use super::*;

    fn key(position_1based: usize, reference: &str, kind: VariantKind) -> VariantKey {
        VariantKey {
            position_1based,
            reference: reference.into(),
            alternate: String::new(),
            kind,
        }
    }

    #[test]
    fn spans_the_anchored_reference_allele_and_both_sides_of_an_insertion() {
        assert_eq!(evidence_span(&key(73, "A", VariantKind::Snv)), 72..73);
        assert_eq!(evidence_span(&key(566, "CAC", VariantKind::Del)), 565..568);
        assert_eq!(evidence_span(&key(309, "C", VariantKind::Ins)), 308..310);
    }
}
