//! Alignment-difference extraction with original-call/reference mappings.

use crate::VariantCallingConfig;
use crate::model::variant::{
    ExcludedVariant, Variant, VariantCallMapping, VariantCallingResult, VariantExclusionReason,
    VariantKind,
};
use crate::{anchor, mapping};
use dna_align::model::alignment::{Alignment, AlignmentColumn};
use dna_kernel::error::{Result, VariantError};
use dna_kernel::model::nucleotide::is_canonical;
use dna_kernel::model::reference::Reference;

use super::eligibility::Eligibility;

/// Extracts normalized primary-sequence differences.
///
/// A masked call aligns as unresolved: a column on a masked call is no SNV
/// candidate at all, and an insertion with a masked inserted call is excluded
/// with the reasons of its masks.
pub(crate) fn call(
    alignment: &Alignment,
    reference: &Reference,
    eligibility: &Eligibility<'_>,
    config: &VariantCallingConfig,
) -> Result<VariantCallingResult> {
    let mut reported = Vec::new();
    let mut excluded = Vec::new();
    let mut index = 0;
    while index < alignment.columns().len() {
        let column = &alignment.columns()[index];
        if column.query_base() == '-' {
            let start = index;
            let previous_reference = previous_reference(alignment.columns(), start);
            let first_deleted_reference =
                column
                    .reference_index_0based()
                    .ok_or(VariantError::Inconsistent(
                        "deletion column lacks a reference coordinate",
                    ))?;
            let previous_flank = previous_flank(alignment.columns(), start)?;
            let mut deleted = String::new();
            let mut deleted_positions = Vec::new();
            while index < alignment.columns().len()
                && alignment.columns()[index].query_base() == '-'
            {
                if alignment.columns()[index].reference_base() != '-' {
                    deleted.push(alignment.columns()[index].reference_base());
                    deleted_positions.extend(alignment.columns()[index].reference_index_0based());
                }
                index += 1;
            }
            let next_reference = next_reference(alignment.columns(), index);
            let next_flank = next_flank(alignment.columns(), index)?;
            let mut reasons = allele_exclusion_reasons(&deleted, config.max_indel_length);
            // A deletion neither of whose flanking calls reads a base observes
            // nothing: where its gap lies among unresolved calls is arbitrary.
            if reasons.is_empty()
                && ![previous_flank, next_flank]
                    .iter()
                    .flatten()
                    .any(|flank| eligibility.reads_base(flank.call_index_0based))
            {
                reasons.push(VariantExclusionReason::UnresolvedFlank);
            }
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
                    undecided: deleted_positions,
                });
            }
            continue;
        }
        if column.reference_base() == '-' {
            let start = index;
            let previous_reference = previous_reference(alignment.columns(), start);
            let previous_flank = previous_flank(alignment.columns(), start)?;
            let mut inserted = String::new();
            let mut calls = Vec::new();
            while index < alignment.columns().len()
                && alignment.columns()[index].reference_base() == '-'
            {
                inserted.push(alignment.columns()[index].query_base());
                calls.push(mapping::supporting(&alignment.columns()[index])?);
                index += 1;
            }
            let next_reference = next_reference(alignment.columns(), index);
            let next_flank = next_flank(alignment.columns(), index)?;
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
                    undecided: [previous_reference, next_reference]
                        .into_iter()
                        .flatten()
                        .collect(),
                });
            }
            continue;
        }
        let masked = column
            .original_call_index_0based()
            .is_some_and(|call| eligibility.mask_reason(call).is_some());
        if column.query_base() != column.reference_base() && !masked {
            if is_canonical(column.query_base()) && is_canonical(column.reference_base()) {
                let reference_position =
                    column
                        .reference_index_0based()
                        .ok_or(VariantError::Inconsistent(
                            "SNV column lacks a reference coordinate",
                        ))?;
                reported.push(anchor::snv(
                    reference,
                    reference_position,
                    column.query_base(),
                    vec![mapping::supporting(column)?],
                )?);
            } else {
                excluded.push(ExcludedVariant {
                    contig: reference.name.clone(),
                    position_1based: column
                        .reference_index_0based()
                        .and_then(|position| position.checked_add(1)),
                    kind: VariantKind::Snv,
                    reasons: vec![VariantExclusionReason::NonCanonicalAllele],
                    undecided: column.reference_index_0based().into_iter().collect(),
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
        .find_map(AlignmentColumn::reference_index_0based)
}

fn next_reference(columns: &[AlignmentColumn], index: usize) -> Option<usize> {
    columns[index..]
        .iter()
        .find_map(AlignmentColumn::reference_index_0based)
}

fn previous_flank(columns: &[AlignmentColumn], index: usize) -> Result<Option<VariantCallMapping>> {
    columns[..index]
        .iter()
        .rev()
        .find(|column| {
            column.original_call_index_0based().is_some()
                && column.reference_index_0based().is_some()
        })
        .map(mapping::flanking)
        .transpose()
}

fn next_flank(columns: &[AlignmentColumn], index: usize) -> Result<Option<VariantCallMapping>> {
    columns[index..]
        .iter()
        .find(|column| {
            column.original_call_index_0based().is_some()
                && column.reference_index_0based().is_some()
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
    use crate::model::variant::VariantCallRole;
    use dna_align::model::alignment::Orientation;
    use dna_kernel::model::reference::ReferenceTopology;
    use dna_kernel::read_evidence::{EvidenceReason, Placement, ReadEvidence};

    use super::*;

    fn alignment(columns: Vec<AlignmentColumn>) -> Alignment {
        Alignment::fixture(Orientation::Forward, Vec::new(), Vec::new(), columns)
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
        let eligibility = Eligibility::new(evidence, &config());
        call(alignment, reference, &eligibility, &config())
    }

    /// A 16-call read whose calls from `from` on are masked as dephased.
    fn masked_from(from: usize) -> ReadEvidence {
        let mut read = ReadEvidence::clean(&"ACGT".repeat(4));
        for index in from..16 {
            let mut call = read.calls()[index];
            call = call.masked(Placement::Anchoring(EvidenceReason::new("dephased_signal")));
            read = read.with_call(index, call);
        }
        read
    }

    // SRS-VAR-004, SRS-VAR-006
    #[test]
    fn leading_alignment_deletion_uses_reference_predecessor() -> Result<()> {
        let result = extract(
            &alignment(vec![
                AlignmentColumn::new('-', 'C', None, Some(2)),
                AlignmentColumn::new('G', 'G', Some(7), Some(3)),
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

    // SRS-VAR-018
    #[test]
    fn a_deletion_no_flank_reads_is_excluded() -> Result<()> {
        // The deleted C lies between two calls read as N.
        let column = |query, reference, call, position| {
            AlignmentColumn::new(query, reference, call, Some(position))
        };
        let result = extract_masked(
            &alignment(vec![
                column('N', 'T', Some(0), 1),
                column('-', 'C', None, 2),
                column('N', 'G', Some(1), 3),
            ]),
            &reference(ReferenceTopology::Linear),
            &ReadEvidence::clean("NN"),
        )?;
        assert!(result.reported.is_empty());
        let deletions = result
            .excluded
            .iter()
            .filter(|excluded| excluded.kind == VariantKind::Del)
            .map(|excluded| excluded.reasons.clone())
            .collect::<Vec<_>>();
        assert_eq!(deletions, [[VariantExclusionReason::UnresolvedFlank]]);
        Ok(())
    }

    // SRS-VAR-002
    #[test]
    fn unresolved_primary_difference_is_excluded() -> Result<()> {
        let result = extract(
            &alignment(vec![AlignmentColumn::new('N', 'C', Some(4), Some(2))]),
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

    // SRS-VAR-004
    #[test]
    fn leading_alignment_insertion_keeps_inserted_and_flanking_mappings() -> Result<()> {
        let result = extract(
            &alignment(vec![
                AlignmentColumn::new('A', '-', Some(4), None),
                AlignmentColumn::new('G', 'G', Some(5), Some(2)),
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

    // SRS-VAR-013
    #[test]
    fn masked_calls_are_no_snv_candidates_and_exclude_insertions() -> Result<()> {
        let columns = vec![
            AlignmentColumn::new('N', 'C', Some(4), Some(2)),
            AlignmentColumn::new('N', '-', Some(5), None),
            AlignmentColumn::new('N', 'G', Some(6), Some(3)),
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

    // SRS-VAR-001
    #[test]
    fn variants_derive_only_from_the_selected_alignment_columns() -> Result<()> {
        // The read's call 5 carries a profile favouring another base, but the
        // selected alignment pairs it with the reference base: no variant.
        let mut evidence = ReadEvidence::clean(&"ACGT".repeat(4));
        let mut call = evidence.calls()[5];
        call.profile = Some(dna_kernel::read_evidence::EvidenceProfile {
            weights: [0.0, 0.0, 0.0, 1.0],
        });
        evidence = evidence.with_call(5, call);
        let matching = alignment(vec![
            AlignmentColumn::new('T', 'T', Some(4), Some(0)),
            AlignmentColumn::new('C', 'C', Some(5), Some(2)),
        ]);
        let result = extract_masked(&matching, &reference(ReferenceTopology::Linear), &evidence)?;
        assert!(result.reported.is_empty());
        assert_eq!(result.excluded_count(), 0);
        // A column difference is the only source of a variant.
        let differing = alignment(vec![AlignmentColumn::new('A', 'C', Some(5), Some(2))]);
        let result = extract_masked(&differing, &reference(ReferenceTopology::Linear), &evidence)?;
        assert_eq!(result.reported.len(), 1);
        assert_eq!(result.reported[0].position_1based, 3);
        assert_eq!(result.reported[0].alternate, "A");
        Ok(())
    }

    // SRS-VAR-003
    #[test]
    fn indel_length_limit_excludes_the_anchor_and_labels_over_limit() -> Result<()> {
        let reference = Reference {
            name: "ref".into(),
            sequence: "TACGTACG".into(),
            topology: ReferenceTopology::Linear,
            sequence_sha256: String::new(),
        };
        let limited = VariantCallingConfig {
            max_indel_length: 2,
            read_end_margin: 0,
        };
        let run = |columns| -> Result<VariantCallingResult> {
            let evidence = ReadEvidence::clean(&"ACGT".repeat(4));
            let eligibility = Eligibility::new(&evidence, &limited);
            call(&alignment(columns), &reference, &eligibility, &limited)
        };
        let deletion = |length: usize| {
            let mut columns = vec![AlignmentColumn::new('T', 'T', Some(4), Some(0))];
            columns.extend((1..=length).map(|offset| {
                AlignmentColumn::new(
                    '-',
                    reference.sequence.as_bytes()[offset] as char,
                    None,
                    Some(offset),
                )
            }));
            columns.push(AlignmentColumn::new('G', 'G', Some(5), Some(length + 1)));
            columns
        };
        // Two deleted bases with the anchor: the reference allele is three
        // long, yet the change is within the limit.
        let at_limit = run(deletion(2))?;
        assert_eq!(at_limit.reported.len(), 1);
        assert_eq!(at_limit.reported[0].reference.len(), 3);
        // Three deleted bases exceed it: labelled, not reported.
        let over = run(deletion(3))?;
        assert!(over.reported.is_empty());
        assert_eq!(over.excluded.len(), 1);
        assert_eq!(over.excluded[0].kind, VariantKind::Del);
        assert_eq!(
            over.excluded[0].reasons,
            vec![VariantExclusionReason::IndelLengthExceeded]
        );
        let insertion = |length: usize| {
            let mut columns = vec![AlignmentColumn::new('T', 'T', Some(4), Some(0))];
            columns.extend(
                (0..length).map(|offset| AlignmentColumn::new('A', '-', Some(5 + offset), None)),
            );
            columns.push(AlignmentColumn::new('A', 'A', Some(5 + length), Some(1)));
            columns
        };
        let at_limit = run(insertion(2))?;
        assert_eq!(at_limit.reported.len(), 1);
        assert_eq!(at_limit.reported[0].alternate.len(), 3);
        let over = run(insertion(3))?;
        assert!(over.reported.is_empty());
        assert_eq!(over.excluded[0].kind, VariantKind::Ins);
        assert_eq!(
            over.excluded[0].reasons,
            vec![VariantExclusionReason::IndelLengthExceeded]
        );
        Ok(())
    }
}
