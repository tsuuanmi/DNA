//! Profile-region and modality support-veto eligibility filters.

use super::eligibility::ReadEligibility;
use crate::model::alignment::Orientation;
use crate::model::variant::{
    ExcludedVariant, ObservedVariant, Variant, VariantCallRole, VariantCallingResult,
    VariantExclusionReason, VariantKind,
};
use dna_kernel::error::{Result, VariantError};
use dna_kernel::model::reference::{Reference, ReferenceTopology};
use dna_kernel::read_evidence::ReadEvidence;

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
    reference: &Reference,
) -> Result<VariantCallingResult> {
    let mut reported = Vec::with_capacity(extracted.reported.len());
    let mut observed = Vec::with_capacity(extracted.reported.len());
    let mut excluded = extracted.excluded;
    for variant in extracted.reported {
        let mut reasons = Vec::new();
        if !in_region(variant.position_1based, regions) {
            reasons.push(VariantExclusionReason::OutsideTargetRegion);
        }
        let eligibility_reasons = eligibility.reasons(&variant, orientation);
        // An insertion into a repeat the read bounds is judged on its
        // placements along the repeat, not on the calls its 3'-most placement
        // happens to take as inserted.
        let repeat = variant.kind == VariantKind::Ins
            && !eligibility_reasons.contains(&VariantExclusionReason::RunBoundary)
            && eligibility.repeat_unvetoed(&supporting_calls(&variant));
        let vetoes = support_vetoes(&variant, evidence)?;
        if !repeat {
            reasons.extend(vetoes);
        }
        reasons.extend(eligibility_reasons);
        if variant.kind == VariantKind::Snv
            && !reasons.contains(&VariantExclusionReason::UnresolvedFlank)
            && unanchored(&variant, eligibility, reference)
        {
            reasons.push(VariantExclusionReason::UnresolvedFlank);
        }
        observed.push(ObservedVariant {
            variant: variant.clone(),
            exclusion_reasons: reasons.clone(),
        });
        if reasons.is_empty() {
            reported.push(variant);
        } else {
            excluded.push(ExcludedVariant::of(&variant, reasons));
        }
    }
    Ok(VariantCallingResult {
        reported,
        observed,
        excluded,
    })
}

/// Whether a substitution rests on a call vouched for by its base only that
/// may sit one position off (see `ReadEligibility::unanchored`).
fn unanchored(variant: &Variant, eligibility: &ReadEligibility<'_>, reference: &Reference) -> bool {
    let bases = reference.sequence.as_bytes();
    let length = bases.len();
    let Some(position) = variant.position_1based.checked_sub(1) else {
        return false;
    };
    let circular = reference.topology == ReferenceTopology::Circular;
    let at = |offset: isize| {
        let shifted = position.checked_add_signed(offset)?;
        let index = if circular && length > 0 {
            (shifted + length) % length
        } else {
            shifted
        };
        bases.get(index).map(|&base| char::from(base))
    };
    let adjacent = [
        if position == 0 && circular {
            bases.last().map(|&base| char::from(base))
        } else {
            at(-1)
        },
        at(1),
    ];
    let Some(base) = variant.alternate.chars().next() else {
        return false;
    };
    supporting_calls(variant)
        .iter()
        .any(|&index| eligibility.unanchored(index, base, adjacent))
}

fn in_region(position_1based: usize, regions: &[[usize; 2]]) -> bool {
    regions
        .iter()
        .any(|[start, end]| *start <= position_1based && position_1based <= *end)
}

/// The variant's supporting calls, in mapping order.
fn supporting_calls(variant: &Variant) -> Vec<usize> {
    variant
        .calls
        .iter()
        .filter(|mapping| mapping.role == VariantCallRole::Supporting)
        .map(|mapping| mapping.call_index_0based)
        .collect()
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
    // A substitution rests on its supporting calls' base, an insertion on
    // their being inserted.
    let mut raised = Vec::new();
    for mapping in supporting {
        let index = mapping.call_index_0based;
        let call = evidence
            .calls()
            .get(index)
            .ok_or(VariantError::MissingCall { index })?;
        let verdict = if variant.kind == VariantKind::Snv {
            call.identity
        } else {
            call.indel
        };
        raised.extend(verdict.withheld.iter());
    }
    Ok(evidence
        .reason_order()
        .iter()
        .filter(|reason| raised.contains(reason))
        .map(|&reason| VariantExclusionReason::Evidence(reason))
        .collect())
}

#[cfg(test)]
mod tests {
    use crate::model::variant::VariantCallMapping;
    use crate::variant_calling::VariantCallingConfig;
    use dna_kernel::error::{Error, VariantError};
    use dna_kernel::read_evidence::{CallEvidence, EvidenceReason};

    use super::*;

    const PEAK: usize = 0;
    const QUALITY: usize = 1;
    const MIXED: usize = 2;
    /// The reasons in reporting order; the first two also withhold an
    /// insertion, the mixed one only a substitution.
    const VOCABULARY: [EvidenceReason; 3] = [
        EvidenceReason::new("peak_below_minimum"),
        EvidenceReason::new("relative_quality_not_above_threshold"),
        EvidenceReason::new("mixed_supporting_dna"),
    ];

    fn reason(index: usize) -> VariantExclusionReason {
        VariantExclusionReason::Evidence(VOCABULARY[index])
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
                let mut call = CallEvidence::scored('A', None);
                for &index in *raised {
                    assert!(call.identity.withheld.insert(VOCABULARY[index]).is_ok());
                    if index != MIXED {
                        assert!(call.indel.withheld.insert(VOCABULARY[index]).is_ok());
                    }
                }
                call
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
            &Reference {
                name: "test".into(),
                sequence: "ACGT".repeat(8),
                topology: ReferenceTopology::Linear,
                sequence_sha256: String::new(),
            },
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
            undecided: Vec::new(),
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
