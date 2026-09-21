//! Configured biological-region and supporting-signal eligibility filters.

use crate::config::VariantCallingConfig;
use crate::error::{Error, Result};
use crate::model::basecalls::BaseCalls;
use crate::model::quality::QualityControlResult;
use crate::model::variant::{
    ExcludedVariant, ObservedVariant, Variant, VariantCallMapping, VariantCallRole,
    VariantCallingResult, VariantExclusionReason, VariantKind,
};

/// Removes normalized candidates outside configured regions or below supporting evidence floors.
pub(super) fn apply(
    extracted: VariantCallingResult,
    calls: &BaseCalls,
    quality: &QualityControlResult,
    config: &VariantCallingConfig,
) -> Result<VariantCallingResult> {
    let mut reported = Vec::with_capacity(extracted.reported.len());
    let mut observed = Vec::with_capacity(extracted.reported.len());
    let mut excluded = extracted.excluded;
    for variant in extracted.reported {
        let mut reasons = Vec::new();
        if !in_configured_region(variant.position_1based, &config.regions) {
            reasons.push(VariantExclusionReason::OutsideConfiguredRegion);
        }
        reasons.extend(supporting_evidence_reasons(
            &variant, calls, quality, config,
        )?);
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

fn in_configured_region(position_1based: usize, regions: &[[usize; 2]]) -> bool {
    regions
        .iter()
        .any(|[start, end]| *start <= position_1based && position_1based <= *end)
}

fn supporting_evidence_reasons(
    variant: &Variant,
    calls: &BaseCalls,
    quality: &QualityControlResult,
    config: &VariantCallingConfig,
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
        return Err(Error::Variant(format!(
            "{:?} at position {} has no supporting calls",
            variant.kind, variant.position_1based
        )));
    }
    let mut peak_failed = false;
    let mut quality_failed = false;
    let mut mixed_supporting_dna = false;
    for mapping in supporting {
        let assessment = assess_call(mapping, calls, quality, config)?;
        peak_failed |= !assessment.peak_passes;
        quality_failed |= !assessment.quality_passes;
        mixed_supporting_dna |= variant.kind == VariantKind::Snv && assessment.mixed_signal;
    }
    let mut reasons = Vec::new();
    if peak_failed {
        reasons.push(VariantExclusionReason::PeakBelowMinimum);
    }
    if quality_failed {
        reasons.push(VariantExclusionReason::RelativeQualityNotAboveThreshold);
    }
    if mixed_supporting_dna {
        reasons.push(VariantExclusionReason::MixedSupportingDNA);
    }
    Ok(reasons)
}

struct SupportingCallAssessment {
    peak_passes: bool,
    quality_passes: bool,
    mixed_signal: bool,
}

fn assess_call(
    mapping: &VariantCallMapping,
    calls: &BaseCalls,
    quality: &QualityControlResult,
    config: &VariantCallingConfig,
) -> Result<SupportingCallAssessment> {
    let index = mapping.call_index_0based;
    let call = calls.calls.get(index).ok_or_else(|| {
        Error::Variant(format!(
            "variant filter references missing call index {index}"
        ))
    })?;
    let score = quality.per_call.get(index).ok_or_else(|| {
        Error::Variant(format!(
            "variant filter references missing quality index {index}"
        ))
    })?;
    if call.index_0based != index || score.index_0based != index {
        return Err(Error::Variant(format!(
            "variant filter call index {index} does not match call/quality records"
        )));
    }
    let highest_peak = call
        .peaks
        .iter()
        .map(|peak| peak.height)
        .max()
        .ok_or_else(|| Error::Variant(format!("call index {index} has no channel peaks")))?;
    Ok(SupportingCallAssessment {
        peak_passes: highest_peak >= config.minimum_peak_height,
        quality_passes: score.relative_quality_score > config.relative_quality_threshold,
        mixed_signal: call.qualifying_channels.len() > 1,
    })
}

#[cfg(test)]
mod tests {
    use crate::model::basecalls::{BaseCall, ChannelPeak, PeakSource};
    use crate::model::nucleotide::Nucleotide;
    use crate::model::quality::CallQuality;

    use super::*;

    fn config(regions: Vec<[usize; 2]>) -> VariantCallingConfig {
        VariantCallingConfig {
            max_indel_length: 50,
            minimum_peak_height: 150,
            relative_quality_threshold: 30,
            regions,
        }
    }

    fn evidence(peaks: &[i32], scores: &[u8]) -> (BaseCalls, QualityControlResult) {
        let calls = peaks
            .iter()
            .enumerate()
            .map(|(index, &height)| BaseCall {
                index_0based: index,
                ploc_0based: index * 4,
                window_start_0based: index * 4,
                window_end_0based_exclusive: index * 4 + 1,
                peaks: std::array::from_fn(|channel| ChannelPeak {
                    base: Nucleotide::ALL[channel],
                    height,
                    position_0based: index * 4,
                    source: PeakSource::LocalMaximum,
                }),
                primary_peak_evidence: None,
                primary: 'A',
                ambiguity: 'A',
                qualifying_channels: vec![Nucleotide::A],
                vendor_agrees: None,
            })
            .collect::<Vec<_>>();
        let per_call = scores
            .iter()
            .enumerate()
            .map(|(index, &score)| CallQuality {
                index_0based: index,
                penalty: 0,
                relative_quality_score: score,
                vendor_quality_applies: false,
            })
            .collect::<Vec<_>>();
        (
            BaseCalls {
                primary_sequence: "A".repeat(calls.len()),
                calls,
            },
            QualityControlResult {
                retained_sequence: "A".repeat(per_call.len()),
                trim_start_0based: 0,
                trim_end_0based_exclusive: per_call.len(),
                per_call,
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

    fn prior_exclusion() -> ExcludedVariant {
        ExcludedVariant {
            contig: "ref".into(),
            position_1based: None,
            kind: VariantKind::Ins,
            reasons: vec![VariantExclusionReason::IndelLengthExceeded],
        }
    }

    #[test]
    fn keeps_inclusive_region_endpoints_and_counts_each_rejection_once() -> Result<()> {
        let (calls, quality) = evidence(&[149, 150, 150, 150], &[31; 4]);
        let extracted = VariantCallingResult {
            reported: vec![
                variant(
                    VariantKind::Snv,
                    9,
                    vec![mapping(VariantCallRole::Supporting, 0)],
                ),
                variant(
                    VariantKind::Snv,
                    10,
                    vec![mapping(VariantCallRole::Supporting, 1)],
                ),
                variant(
                    VariantKind::Snv,
                    20,
                    vec![mapping(VariantCallRole::Supporting, 2)],
                ),
                variant(
                    VariantKind::Snv,
                    21,
                    vec![mapping(VariantCallRole::Supporting, 3)],
                ),
            ],
            observed: Vec::new(),
            excluded: vec![prior_exclusion()],
        };

        let result = apply(extracted, &calls, &quality, &config(vec![[10, 20]]))?;

        assert_eq!(result.reported.len(), 2);
        assert_eq!(result.reported[0].position_1based, 10);
        assert_eq!(result.reported[1].position_1based, 20);
        assert_eq!(result.excluded_count(), 3);
        assert_eq!(
            result.excluded[1].reasons,
            vec![
                VariantExclusionReason::OutsideConfiguredRegion,
                VariantExclusionReason::PeakBelowMinimum,
            ]
        );
        assert_eq!(result.excluded[1].position_1based, Some(9));
        assert_eq!(result.excluded[2].position_1based, Some(21));
        Ok(())
    }

    #[test]
    fn applies_peak_and_strict_quality_boundaries() -> Result<()> {
        let (calls, quality) = evidence(&[149, 150, 150], &[31, 30, 31]);
        let extracted = VariantCallingResult {
            reported: (0..3)
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
        };

        let result = apply(extracted, &calls, &quality, &config(vec![[1, 3]]))?;

        assert_eq!(result.reported.len(), 1);
        assert_eq!(result.reported[0].position_1based, 3);
        assert_eq!(result.observed.len(), 3);
        assert_eq!(
            result.observed[0].exclusion_reasons,
            vec![VariantExclusionReason::PeakBelowMinimum]
        );
        assert_eq!(
            result.observed[1].exclusion_reasons,
            vec![VariantExclusionReason::RelativeQualityNotAboveThreshold]
        );
        assert!(result.observed[2].eligible());
        assert_eq!(result.excluded_count(), 2);
        assert_eq!(
            result.excluded[0].reasons,
            vec![VariantExclusionReason::PeakBelowMinimum]
        );
        assert_eq!(
            result.excluded[1].reasons,
            vec![VariantExclusionReason::RelativeQualityNotAboveThreshold]
        );
        Ok(())
    }

    #[test]
    fn excludes_mixed_snv_but_does_not_apply_the_gate_to_insertions() -> Result<()> {
        let (mut calls, quality) = evidence(&[200, 200], &[31, 31]);
        for call in &mut calls.calls {
            call.qualifying_channels = vec![Nucleotide::A, Nucleotide::C];
            call.ambiguity = 'M';
        }
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

        let result = apply(extracted, &calls, &quality, &config(vec![[1, 2]]))?;

        assert_eq!(result.reported.len(), 1);
        assert_eq!(result.reported[0].kind, VariantKind::Ins);
        assert_eq!(
            result.observed[0].exclusion_reasons,
            vec![VariantExclusionReason::MixedSupportingDNA]
        );
        assert!(result.observed[1].eligible());
        assert_eq!(
            result.excluded[0].reasons,
            vec![VariantExclusionReason::MixedSupportingDNA]
        );
        Ok(())
    }

    #[test]
    fn requires_every_inserted_base_but_ignores_insertion_flanks() -> Result<()> {
        let (calls, quality) = evidence(&[150, 149, 1], &[31, 31, 0]);
        let extracted = VariantCallingResult {
            reported: vec![
                variant(
                    VariantKind::Ins,
                    1,
                    vec![
                        mapping(VariantCallRole::Supporting, 0),
                        mapping(VariantCallRole::Supporting, 1),
                        mapping(VariantCallRole::Flanking, 2),
                    ],
                ),
                variant(
                    VariantKind::Ins,
                    2,
                    vec![
                        mapping(VariantCallRole::Supporting, 0),
                        mapping(VariantCallRole::Flanking, 2),
                    ],
                ),
            ],
            observed: Vec::new(),
            excluded: Vec::new(),
        };

        let result = apply(extracted, &calls, &quality, &config(vec![[1, 2]]))?;

        assert_eq!(result.reported.len(), 1);
        assert_eq!(result.reported[0].position_1based, 2);
        assert_eq!(result.excluded_count(), 1);
        Ok(())
    }

    #[test]
    fn deletion_flanks_are_exempt_from_signal_thresholds() -> Result<()> {
        let (calls, quality) = evidence(&[1, 1], &[0, 0]);
        let extracted = VariantCallingResult {
            reported: vec![variant(
                VariantKind::Del,
                5,
                vec![
                    mapping(VariantCallRole::Flanking, 0),
                    mapping(VariantCallRole::Flanking, 1),
                ],
            )],
            observed: Vec::new(),
            excluded: Vec::new(),
        };

        let result = apply(extracted, &calls, &quality, &config(vec![[5, 5]]))?;

        assert_eq!(result.reported.len(), 1);
        assert_eq!(result.excluded_count(), 0);
        Ok(())
    }

    #[test]
    fn rejects_missing_supporting_call_mapping() {
        let (calls, quality) = evidence(&[150], &[31]);
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
            apply(extracted, &calls, &quality, &config(vec![[1, 1]])),
            Err(Error::Variant(message)) if message.contains("missing call index 2")
        ));
    }
}
