//! Alignment-difference extraction with original-call/reference mappings.

use crate::model::alignment::{Alignment, AlignmentColumn};
use crate::model::variant::{
    ExcludedVariant, Variant, VariantCallMapping, VariantCallingResult, VariantExclusionReason,
    VariantKind,
};
use crate::variant_calling::VariantCallingConfig;
use crate::variant_calling::{anchor, mapping};
use dna_kernel::error::{Result, VariantError};
use dna_kernel::model::nucleotide::is_canonical;
use dna_kernel::model::reference::Reference;

use super::eligibility::ReadEligibility;

/// Extracts normalized primary-sequence differences.
///
/// A masked call aligns as unresolved: a column on a masked call is no SNV
/// candidate at all, and an insertion with a masked inserted call is excluded
/// with the reasons of its masks.
pub(crate) fn call(
    alignment: &Alignment,
    reference: &Reference,
    eligibility: &ReadEligibility<'_>,
    config: &VariantCallingConfig,
) -> Result<VariantCallingResult> {
    let mut reported = Vec::new();
    let mut excluded = Vec::new();
    let mut index = 0;
    while index < alignment.columns.len() {
        let column = &alignment.columns[index];
        if column.query_base == '-' {
            let start = index;
            let previous_reference = previous_reference(&alignment.columns, start);
            let first_deleted_reference =
                column
                    .reference_index_0based
                    .ok_or(VariantError::Inconsistent(
                        "deletion column lacks a reference coordinate",
                    ))?;
            let previous_flank = previous_flank(&alignment.columns, start)?;
            let mut deleted = String::new();
            while index < alignment.columns.len() && alignment.columns[index].query_base == '-' {
                if alignment.columns[index].reference_base != '-' {
                    deleted.push(alignment.columns[index].reference_base);
                }
                index += 1;
            }
            let next_reference = next_reference(&alignment.columns, index);
            let next_flank = next_flank(&alignment.columns, index)?;
            let reasons = allele_exclusion_reasons(&deleted, config.max_indel_length);
            if reasons.is_empty() {
                reported.push(anchor::deletion(
                    reference,
                    previous_reference,
                    first_deleted_reference,
                    next_reference,
                    &deleted,
                    mapping::sort_dedup(optional_pair(previous_flank, next_flank)),
                )?);
            } else {
                excluded.push(ExcludedVariant {
                    contig: reference.name.clone(),
                    position_1based: None,
                    kind: VariantKind::Del,
                    reasons,
                });
            }
            continue;
        }
        if column.reference_base == '-' {
            let start = index;
            let previous_reference = previous_reference(&alignment.columns, start);
            let previous_flank = previous_flank(&alignment.columns, start)?;
            let mut inserted = String::new();
            let mut calls = Vec::new();
            while index < alignment.columns.len() && alignment.columns[index].reference_base == '-'
            {
                inserted.push(alignment.columns[index].query_base);
                calls.push(mapping::supporting(&alignment.columns[index])?);
                index += 1;
            }
            let next_reference = next_reference(&alignment.columns, index);
            let next_flank = next_flank(&alignment.columns, index)?;
            let mut reasons = Vec::new();
            for reason in calls
                .iter()
                .filter_map(|call| eligibility.mask_reason(call.call_index_0based))
            {
                if !reasons.contains(&reason) {
                    reasons.push(reason);
                }
            }
            if reasons.is_empty() {
                reasons = allele_exclusion_reasons(&inserted, config.max_indel_length);
            }
            if reasons.is_empty() {
                calls.extend(optional_pair(previous_flank, next_flank));
                reported.push(anchor::insertion(
                    reference,
                    previous_reference,
                    next_reference,
                    &inserted,
                    mapping::sort_dedup(calls),
                )?);
            } else {
                excluded.push(ExcludedVariant {
                    contig: reference.name.clone(),
                    position_1based: None,
                    kind: VariantKind::Ins,
                    reasons,
                });
            }
            continue;
        }
        let masked = column
            .original_call_index_0based
            .is_some_and(|call| eligibility.mask_reason(call).is_some());
        if column.query_base != column.reference_base && !masked {
            if is_canonical(column.query_base) && is_canonical(column.reference_base) {
                let reference_position =
                    column
                        .reference_index_0based
                        .ok_or(VariantError::Inconsistent(
                            "SNV column lacks a reference coordinate",
                        ))?;
                reported.push(anchor::snv(
                    reference,
                    reference_position,
                    column.query_base,
                    vec![mapping::supporting(column)?],
                )?);
            } else {
                excluded.push(ExcludedVariant {
                    contig: reference.name.clone(),
                    position_1based: column
                        .reference_index_0based
                        .and_then(|position| position.checked_add(1)),
                    kind: VariantKind::Snv,
                    reasons: vec![VariantExclusionReason::NonCanonicalAllele],
                });
            }
        }
        index += 1;
    }
    reported.sort_by(|left, right| {
        (
            &left.contig,
            left.position_1based,
            &left.reference,
            &left.alternate,
        )
            .cmp(&(
                &right.contig,
                right.position_1based,
                &right.reference,
                &right.alternate,
            ))
    });
    let mut merged: Vec<Variant> = Vec::with_capacity(reported.len());
    for variant in reported {
        if let Some(previous) = merged.last_mut()
            && previous.contig == variant.contig
            && previous.position_1based == variant.position_1based
            && previous.reference == variant.reference
            && previous.alternate == variant.alternate
        {
            previous.calls.extend(variant.calls);
            previous.calls = mapping::sort_dedup(std::mem::take(&mut previous.calls));
            continue;
        }
        merged.push(variant);
    }
    Ok(VariantCallingResult {
        reported: merged,
        observed: Vec::new(),
        excluded,
    })
}

fn allele_exclusion_reasons(allele: &str, max_indel_length: usize) -> Vec<VariantExclusionReason> {
    let mut reasons = Vec::new();
    if allele.len() > max_indel_length {
        reasons.push(VariantExclusionReason::IndelLengthExceeded);
    }
    if !allele.chars().all(is_canonical) {
        reasons.push(VariantExclusionReason::NonCanonicalAllele);
    }
    reasons
}

fn previous_reference(columns: &[AlignmentColumn], index: usize) -> Option<usize> {
    columns[..index]
        .iter()
        .rev()
        .find_map(|column| column.reference_index_0based)
}

fn next_reference(columns: &[AlignmentColumn], index: usize) -> Option<usize> {
    columns[index..]
        .iter()
        .find_map(|column| column.reference_index_0based)
}

fn previous_flank(columns: &[AlignmentColumn], index: usize) -> Result<Option<VariantCallMapping>> {
    columns[..index]
        .iter()
        .rev()
        .find(|column| {
            column.original_call_index_0based.is_some() && column.reference_index_0based.is_some()
        })
        .map(mapping::flanking)
        .transpose()
}

fn next_flank(columns: &[AlignmentColumn], index: usize) -> Result<Option<VariantCallMapping>> {
    columns[index..]
        .iter()
        .find(|column| {
            column.original_call_index_0based.is_some() && column.reference_index_0based.is_some()
        })
        .map(mapping::flanking)
        .transpose()
}

fn optional_pair(
    left: Option<VariantCallMapping>,
    right: Option<VariantCallMapping>,
) -> Vec<VariantCallMapping> {
    left.into_iter().chain(right).collect()
}

#[cfg(test)]
mod tests {
    use crate::model::alignment::{AlignmentMetrics, Orientation};
    use crate::model::variant::VariantCallRole;
    use dna_kernel::model::reference::ReferenceTopology;
    use dna_kernel::read_evidence::{CallMask, EvidenceReason, MaskedAlignment, ReadEvidence};

    use super::*;

    fn alignment(columns: Vec<AlignmentColumn>) -> Alignment {
        Alignment {
            orientation: Orientation::Forward,
            score: 0,
            reference_segments: Vec::new(),
            callable_segments: Vec::new(),
            wraps_origin: false,
            metrics: AlignmentMetrics {
                exact_matches: 0,
                mismatches: 0,
                gap_opens: 1,
                callable_columns: 1,
                callable_identity: 1.0,
                unresolved_query_bases: 0,
                masked_query_bases: 0,
            },
            columns,
        }
    }

    fn reference(topology: ReferenceTopology) -> Reference {
        Reference {
            name: "ref".into(),
            sequence: "TTCG".into(),
            topology,
            sequence_sha256: String::new(),
        }
    }

    fn config() -> VariantCallingConfig {
        VariantCallingConfig {
            max_indel_length: 50,
            read_end_margin: 0,
        }
    }

    /// Extracts with every call of a 16-call read callable.
    fn extract(alignment: &Alignment, reference: &Reference) -> Result<VariantCallingResult> {
        extract_masked(
            alignment,
            reference,
            &ReadEvidence::clean(&"ACGT".repeat(4)),
        )
    }

    fn extract_masked(
        alignment: &Alignment,
        reference: &Reference,
        evidence: &ReadEvidence,
    ) -> Result<VariantCallingResult> {
        let eligibility = ReadEligibility::new(evidence, &config());
        call(alignment, reference, &eligibility, &config())
    }

    /// A 16-call read whose calls from `from` on are masked as dephased.
    fn masked_from(from: usize) -> ReadEvidence {
        let mut read = ReadEvidence::clean(&"ACGT".repeat(4));
        for index in from..16 {
            let mut call = read.calls()[index];
            call.mask = Some(CallMask {
                alignment: MaskedAlignment::Anchoring,
                reason: EvidenceReason::new("dephased_signal"),
            });
            read = read.with_call(index, call);
        }
        read
    }

    #[test]
    fn leading_alignment_deletion_uses_reference_predecessor() -> Result<()> {
        let result = extract(
            &alignment(vec![
                AlignmentColumn {
                    query_base: '-',
                    reference_base: 'C',
                    original_call_index_0based: None,
                    reference_index_0based: Some(2),
                },
                AlignmentColumn {
                    query_base: 'G',
                    reference_base: 'G',
                    original_call_index_0based: Some(7),
                    reference_index_0based: Some(3),
                },
            ]),
            &reference(ReferenceTopology::Linear),
        )?;
        let variant = &result.reported[0];
        assert_eq!(variant.position_1based, 2);
        assert_eq!(variant.reference, "TC");
        assert_eq!(variant.alternate, "T");
        assert_eq!(variant.calls.len(), 1);
        assert_eq!(variant.calls[0].role, VariantCallRole::Flanking);
        assert_eq!(variant.calls[0].call_index_0based, 7);
        assert_eq!(variant.calls[0].reference_position_0based, Some(3));
        Ok(())
    }

    #[test]
    fn unresolved_primary_difference_is_excluded() -> Result<()> {
        let result = extract(
            &alignment(vec![AlignmentColumn {
                query_base: 'N',
                reference_base: 'C',
                original_call_index_0based: Some(4),
                reference_index_0based: Some(2),
            }]),
            &reference(ReferenceTopology::Linear),
        )?;
        assert!(result.reported.is_empty());
        assert_eq!(result.excluded_count(), 1);
        assert_eq!(result.excluded[0].kind, VariantKind::Snv);
        assert_eq!(result.excluded[0].position_1based, Some(3));
        assert_eq!(
            result.excluded[0].reasons,
            vec![VariantExclusionReason::NonCanonicalAllele]
        );
        Ok(())
    }

    #[test]
    fn leading_alignment_insertion_keeps_inserted_and_flanking_mappings() -> Result<()> {
        let result = extract(
            &alignment(vec![
                AlignmentColumn {
                    query_base: 'A',
                    reference_base: '-',
                    original_call_index_0based: Some(4),
                    reference_index_0based: None,
                },
                AlignmentColumn {
                    query_base: 'G',
                    reference_base: 'G',
                    original_call_index_0based: Some(5),
                    reference_index_0based: Some(2),
                },
            ]),
            &reference(ReferenceTopology::Linear),
        )?;
        let variant = &result.reported[0];
        assert_eq!(variant.position_1based, 2);
        assert_eq!(variant.reference, "T");
        assert_eq!(variant.alternate, "TA");
        assert_eq!(variant.calls.len(), 2);
        assert_eq!(variant.calls[0].role, VariantCallRole::Supporting);
        assert_eq!(variant.calls[0].call_index_0based, 4);
        assert_eq!(variant.calls[0].reference_position_0based, None);
        assert_eq!(variant.calls[1].role, VariantCallRole::Flanking);
        assert_eq!(variant.calls[1].call_index_0based, 5);
        assert_eq!(variant.calls[1].reference_position_0based, Some(2));
        Ok(())
    }

    #[test]
    fn masked_calls_are_no_snv_candidates_and_exclude_insertions() -> Result<()> {
        let columns = vec![
            AlignmentColumn {
                query_base: 'N',
                reference_base: 'C',
                original_call_index_0based: Some(4),
                reference_index_0based: Some(2),
            },
            AlignmentColumn {
                query_base: 'N',
                reference_base: '-',
                original_call_index_0based: Some(5),
                reference_index_0based: None,
            },
            AlignmentColumn {
                query_base: 'N',
                reference_base: 'G',
                original_call_index_0based: Some(6),
                reference_index_0based: Some(3),
            },
        ];
        let result = extract_masked(
            &alignment(columns),
            &reference(ReferenceTopology::Linear),
            &masked_from(4),
        )?;
        assert!(result.reported.is_empty());
        assert_eq!(result.excluded_count(), 1);
        assert_eq!(result.excluded[0].kind, VariantKind::Ins);
        assert_eq!(
            result.excluded[0].reasons,
            vec![VariantExclusionReason::Evidence(EvidenceReason::new(
                "dephased_signal"
            ))]
        );
        Ok(())
    }
}
