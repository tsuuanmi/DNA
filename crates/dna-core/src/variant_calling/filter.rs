//! Profile-region and modality support-veto eligibility filters.

use super::eligibility::ReadEligibility;
use crate::model::alignment::Orientation;
use crate::model::variant::{
    ExcludedVariant, ObservedVariant, Variant, VariantCallRole, VariantCallingResult,
    VariantExclusionReason, VariantKind,
};
use dna_kernel::error::{Result, VariantError};
use dna_kernel::read_evidence::{ReadEvidence, VetoScope, VetoSet};

/// Removes normalized candidates outside the profile regions, vetoed by their
/// supporting calls' modality evidence, or ineligible for their read.
///
/// Reasons are ordered: region, support vetoes in the read's vocabulary order,
/// then `read_end`, then `run_boundary`, then mask reasons.
pub(super) fn apply(
    extracted: VariantCallingResult,
    evidence: &ReadEvidence,
    eligibility: &ReadEligibility<'_>,
    orientation: Orientation,
    regions: &[[usize; 2]],
) -> Result<VariantCallingResult> {
    let mut reported = Vec::with_capacity(extracted.reported.len());
    let mut observed = Vec::with_capacity(extracted.reported.len());
    let mut excluded = extracted.excluded;
    for variant in extracted.reported {
        let mut reasons = Vec::new();
        if !in_region(variant.position_1based, regions) {
            reasons.push(VariantExclusionReason::OutsideTargetRegion);
        }
        reasons.extend(support_vetoes(&variant, evidence)?);
        reasons.extend(eligibility.reasons(&variant, orientation));
        observed.push(ObservedVariant {
            variant: variant.clone(),
            exclusion_reasons: reasons.clone(),
        });
        if reasons.is_empty() {
            reported.push(variant);
        } else {
            excluded.push(ExcludedVariant {
                contig: variant.contig,
                position_1based: Some(variant.position_1based),
                kind: variant.kind,
                reasons,
            });
        }
    }
    Ok(VariantCallingResult {
        reported,
        observed,
        excluded,
    })
}

fn in_region(position_1based: usize, regions: &[[usize; 2]]) -> bool {
    regions
        .iter()
        .any(|[start, end]| *start <= position_1based && position_1based <= *end)
}

/// The vetoes raised by any supporting call that apply to the variant's kind,
/// in vocabulary order. Deletions have no supporting calls and no vetoes.
fn support_vetoes(
    variant: &Variant,
    evidence: &ReadEvidence,
) -> Result<Vec<VariantExclusionReason>> {
    if variant.kind == VariantKind::Del {
        return Ok(Vec::new());
    }
    let supporting = variant
        .calls
        .iter()
        .filter(|mapping| mapping.role == VariantCallRole::Supporting)
        .collect::<Vec<_>>();
    if supporting.is_empty() {
        return Err(VariantError::NoSupportingCalls {
            kind: variant.kind.label(),
            position: variant.position_1based,
        }
        .into());
    }
    let mut raised = VetoSet::default();
    for mapping in supporting {
        let index = mapping.call_index_0based;
        let call = evidence
            .calls()
            .get(index)
            .ok_or(VariantError::MissingCall { index })?;
        raised = raised.union(call.vetoes);
    }
    Ok(evidence
        .support_vetoes()
        .iter()
        .enumerate()
        .filter(|(index, veto)| {
            raised.contains(*index)
                && (veto.scope == VetoScope::SubstitutionsAndInsertions
                    || variant.kind == VariantKind::Snv)
        })
        .map(|(_, veto)| VariantExclusionReason::Evidence(veto.reason))
        .collect())
}

#[cfg(test)]
mod tests {
    use crate::model::variant::VariantCallMapping;
    use crate::variant_calling::VariantCallingConfig;
    use dna_kernel::error::{Error, VariantError};
    use dna_kernel::read_evidence::{CallEvidence, EvidenceReason, SupportVeto};

    use super::*;

    const PEAK: usize = 0;
    const QUALITY: usize = 1;
    const MIXED: usize = 2;
    const VOCABULARY: [SupportVeto; 3] = [
        SupportVeto {
            reason: EvidenceReason::new("peak_below_minimum"),
            scope: VetoScope::SubstitutionsAndInsertions,
        },
        SupportVeto {
            reason: EvidenceReason::new("relative_quality_not_above_threshold"),
            scope: VetoScope::SubstitutionsAndInsertions,
        },
        SupportVeto {
            reason: EvidenceReason::new("mixed_supporting_dna"),
            scope: VetoScope::Substitutions,
        },
    ];

    fn reason(index: usize) -> VariantExclusionReason {
        VariantExclusionReason::Evidence(VOCABULARY[index].reason)
    }

    fn config(read_end_margin: usize) -> VariantCallingConfig {
        VariantCallingConfig {
            max_indel_length: 50,
            read_end_margin,
        }
    }

    /// Evidence whose call `i` raises the vetoes listed at `vetoes[i]`; its
    /// ends are vouched, so they bound every run.
    fn evidence(vetoes: &[&[usize]]) -> Result<ReadEvidence> {
        Ok(unvouched(vetoes)?.with_vouched_ends())
    }

    /// Like [`evidence`], with ends that bound nothing.
    fn unvouched(vetoes: &[&[usize]]) -> Result<ReadEvidence> {
        let calls = vetoes
            .iter()
            .map(|raised| {
                let mut set = VetoSet::default();
                for &index in *raised {
                    set.insert(index);
                }
                CallEvidence {
                    base: 'A',
                    profile: None,
                    mask: None,
                    vetoes: set,
                }
            })
            .collect::<Vec<_>>();
        let informative = 0..calls.len();
        ReadEvidence::new(calls, informative, VOCABULARY.to_vec())
    }

    fn filter(
        extracted: VariantCallingResult,
        evidence: &ReadEvidence,
        read_end_margin: usize,
        regions: &[[usize; 2]],
    ) -> Result<VariantCallingResult> {
        let config = config(read_end_margin);
        let eligibility = ReadEligibility::new(evidence, &config);
        apply(
            extracted,
            evidence,
            &eligibility,
            Orientation::Forward,
            regions,
        )
    }

    fn mapping(role: VariantCallRole, index: usize) -> VariantCallMapping {
        VariantCallMapping {
            role,
            call_index_0based: index,
            reference_position_0based: (role == VariantCallRole::Flanking).then_some(index),
        }
    }

    fn variant(kind: VariantKind, position: usize, calls: Vec<VariantCallMapping>) -> Variant {
        Variant {
            contig: "ref".into(),
            position_1based: position,
            reference: "A".into(),
            alternate: "T".into(),
            kind,
            calls,
        }
    }

    fn snvs(count: usize) -> VariantCallingResult {
        VariantCallingResult {
            reported: (0..count)
                .map(|index| {
                    variant(
                        VariantKind::Snv,
                        index + 1,
                        vec![mapping(VariantCallRole::Supporting, index)],
                    )
                })
                .collect(),
            observed: Vec::new(),
            excluded: Vec::new(),
        }
    }

    #[test]
    fn keeps_inclusive_region_endpoints_and_orders_region_before_vetoes() -> Result<()> {
        let read = evidence(&[&[PEAK], &[], &[], &[]])?;
        let mut extracted = snvs(4);
        for (variant, position) in extracted.reported.iter_mut().zip([9, 10, 20, 21]) {
            variant.position_1based = position;
        }
        extracted.excluded.push(ExcludedVariant {
            contig: "ref".into(),
            position_1based: None,
            kind: VariantKind::Ins,
            reasons: vec![VariantExclusionReason::IndelLengthExceeded],
        });

        let result = filter(extracted, &read, 0, &[[10, 20]])?;

        assert_eq!(result.reported.len(), 2);
        assert_eq!(result.reported[0].position_1based, 10);
        assert_eq!(result.reported[1].position_1based, 20);
        assert_eq!(result.excluded_count(), 3);
        assert_eq!(
            result.excluded[1].reasons,
            vec![VariantExclusionReason::OutsideTargetRegion, reason(PEAK)]
        );
        assert_eq!(result.excluded[2].position_1based, Some(21));
        Ok(())
    }

    #[test]
    fn reports_vetoes_of_all_supporting_calls_in_vocabulary_order() -> Result<()> {
        let read = evidence(&[&[QUALITY], &[PEAK]])?;
        let extracted = VariantCallingResult {
            reported: vec![variant(
                VariantKind::Ins,
                1,
                vec![
                    mapping(VariantCallRole::Supporting, 0),
                    mapping(VariantCallRole::Supporting, 1),
                ],
            )],
            observed: Vec::new(),
            excluded: Vec::new(),
        };

        let result = filter(extracted, &read, 0, &[[1, 1]])?;

        assert_eq!(
            result.observed[0].exclusion_reasons,
            vec![reason(PEAK), reason(QUALITY)]
        );
        Ok(())
    }

    #[test]
    fn applies_substitution_only_vetoes_to_snvs() -> Result<()> {
        let read = evidence(&[&[MIXED], &[MIXED]])?;
        let extracted = VariantCallingResult {
            reported: vec![
                variant(
                    VariantKind::Snv,
                    1,
                    vec![mapping(VariantCallRole::Supporting, 0)],
                ),
                variant(
                    VariantKind::Ins,
                    2,
                    vec![mapping(VariantCallRole::Supporting, 1)],
                ),
            ],
            observed: Vec::new(),
            excluded: Vec::new(),
        };

        let result = filter(extracted, &read, 0, &[[1, 2]])?;

        assert_eq!(result.reported.len(), 1);
        assert_eq!(result.reported[0].kind, VariantKind::Ins);
        assert_eq!(result.observed[0].exclusion_reasons, vec![reason(MIXED)]);
        assert!(result.observed[1].eligible());
        Ok(())
    }

    #[test]
    fn ignores_insertion_flanks_and_deletion_flanks() -> Result<()> {
        let read = evidence(&[&[], &[PEAK, QUALITY]])?;
        let extracted = VariantCallingResult {
            reported: vec![
                variant(
                    VariantKind::Ins,
                    1,
                    vec![
                        mapping(VariantCallRole::Supporting, 0),
                        mapping(VariantCallRole::Flanking, 1),
                    ],
                ),
                variant(
                    VariantKind::Del,
                    2,
                    vec![
                        mapping(VariantCallRole::Flanking, 0),
                        mapping(VariantCallRole::Flanking, 1),
                    ],
                ),
            ],
            observed: Vec::new(),
            excluded: Vec::new(),
        };

        let result = filter(extracted, &read, 0, &[[1, 2]])?;

        assert_eq!(result.reported.len(), 2);
        assert_eq!(result.excluded_count(), 0);
        Ok(())
    }

    #[test]
    fn rejects_missing_supporting_call_mapping() -> Result<()> {
        let read = evidence(&[&[]])?;
        let extracted = VariantCallingResult {
            reported: vec![variant(
                VariantKind::Snv,
                1,
                vec![mapping(VariantCallRole::Supporting, 2)],
            )],
            observed: Vec::new(),
            excluded: Vec::new(),
        };

        assert!(matches!(
            filter(extracted, &read, 0, &[[1, 1]]),
            Err(Error::Variant(VariantError::MissingCall { index: 2 }))
        ));
        Ok(())
    }

    #[test]
    fn marks_a_variant_supported_at_the_read_end_ineligible() -> Result<()> {
        let none: &[usize] = &[];
        let read = unvouched(&[none; 6])?;

        let result = filter(snvs(3), &read, 2, &[[1, 100]])?;

        assert_eq!(result.reported.len(), 1);
        assert_eq!(result.reported[0].position_1based, 3);
        assert_eq!(
            result.observed[0].exclusion_reasons,
            [VariantExclusionReason::ReadEnd]
        );
        Ok(())
    }
}
