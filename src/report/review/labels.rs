//! Label order and the summary: how labels sort within an item, an
//! item's priority, and the counts over the items. The label table itself
//! is `LabelKind` in `model::review_result`.

use std::collections::BTreeSet;

use super::ranges::inside;
use crate::model::review_result::{
    CONTEXT_RANK, ItemKind, ItemResult, LabelCounts, LabelResult, PositionCounts, ReviewStatus,
    SummaryResult, ZoneCounts, ZoneResult,
};

/// Priority of an unlabelled variant, after every label's rank.
pub(super) const UNLABELLED: u8 = CONTEXT_RANK + 1;

/// Label order within an item: rank, name, then the position or read.
pub(super) fn key(label: &LabelResult) -> (u8, &'static str, usize, &str) {
    let (position, read) = match label {
        LabelResult::PossibleHeteroplasmy { position, .. } => (*position, ""),
        LabelResult::ReferenceSecond { position, read, .. } => (*position, read.as_str()),
        LabelResult::EdgeDeletion { reads, .. }
        | LabelResult::RunDeletion { reads, .. }
        | LabelResult::RunIndel { reads, .. } => (0, reads.first().map_or("", String::as_str)),
        LabelResult::OneRead { read }
        | LabelResult::ReadDisagrees { read, .. }
        | LabelResult::NoisySignal { read, .. }
        | LabelResult::PhaseShift { read, .. }
        | LabelResult::ReadEdge { read }
        | LabelResult::PoorSignal { read } => (0, read.as_str()),
        LabelResult::StrandsDisagree { .. }
        | LabelResult::SameStrandDissent { .. }
        | LabelResult::NotCalled
        | LabelResult::LongRun { .. }
        | LabelResult::OneStrandInsertion { .. }
        | LabelResult::EstimatedLength { .. }
        | LabelResult::LengthUnconfirmed { .. }
        | LabelResult::LengthNotRecorded { .. }
        | LabelResult::TractLength { .. }
        | LabelResult::Undecided
        | LabelResult::OneStrand { .. }
        | LabelResult::NotCovered
        | LabelResult::DephasedEnd { .. }
        | LabelResult::BesideNotAnalysed => (0, ""),
    };
    (label.rank(), label.kind().name(), position, read)
}

/// An item's priority: its best label rank, or [`UNLABELLED`].
pub(super) fn priority(labels: &[LabelResult]) -> u8 {
    labels
        .iter()
        .map(LabelResult::rank)
        .min()
        .unwrap_or(UNLABELLED)
}

/// Counts over the items, the reportable positions inside `regions` that
/// `check` items hold, the `reference_only` positions, and the `zones` the
/// checks touch; `clear` when the call succeeded and no item is a check.
pub(super) fn summarize(
    items: &[ItemResult],
    status: ReviewStatus,
    regions: &[[usize; 2]],
    zones: &[ZoneResult],
    reference_only: &[usize],
    not_analysed: usize,
) -> SummaryResult {
    let mut labels = LabelCounts::default();
    for item in items {
        let kinds = item
            .labels
            .iter()
            .map(LabelResult::kind)
            .collect::<BTreeSet<_>>();
        for kind in kinds {
            labels.add(kind);
        }
    }
    let checks = items
        .iter()
        .filter(|item| item.kind == ItemKind::Check)
        .count();
    let check = items
        .iter()
        .filter(|item| item.kind == ItemKind::Check)
        .flat_map(|item| item.first.position..=item.last.position)
        .filter(|&position| inside(regions, position))
        .collect::<BTreeSet<_>>()
        .len();
    SummaryResult {
        items: items.len(),
        checks,
        variants: items
            .iter()
            .filter(|item| item.kind == ItemKind::Variant)
            .count(),
        labels,
        positions: PositionCounts {
            reportable: regions.iter().map(|[start, end]| end - start + 1).sum(),
            check,
            reference_only: reference_only.len(),
            not_analysed,
        },
        zones: ZoneCounts {
            total: zones.len(),
            check: zones.iter().filter(|zone| zone.checks > 0).count(),
        },
        clear: status == ReviewStatus::Called && checks == 0,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::review_result::LabelKind;
    use crate::report::review::test_support::{F, R, item, read};

    // SRS-API-014
    #[test]
    fn labels_and_items_sort_by_priority() {
        let mut labels = vec![
            LabelResult::ReadEdge { read: read("b") },
            LabelResult::NotCovered,
            LabelResult::ReadEdge { read: read("a") },
            LabelResult::OneStrand {
                strand: F,
                reads: vec![read("a")],
            },
            LabelResult::Undecided,
            LabelResult::NotCalled,
            LabelResult::EstimatedLength { runs: Vec::new() },
            LabelResult::ReadDisagrees {
                read: read("a"),
                strand: F,
                sequence: "CC".to_owned(),
            },
            LabelResult::StrandsDisagree { reads: Vec::new() },
        ];
        labels.sort_by(|left, right| key(left).cmp(&key(right)));
        assert_eq!(
            labels
                .iter()
                .map(|label| label.kind().name())
                .collect::<Vec<_>>(),
            [
                "read_disagrees",
                "strands_disagree",
                "estimated_length",
                "not_called",
                "undecided",
                "one_strand",
                "not_covered",
                "read_edge",
                "read_edge",
            ]
        );
        assert_eq!(labels[8], LabelResult::ReadEdge { read: read("b") });
        assert_eq!(priority(&labels), 0);
        assert_eq!(priority(&labels[5..]), 3);
        assert_eq!(priority(&[]), UNLABELLED);
    }

    // SRS-API-014
    #[test]
    fn only_labels_that_leave_an_item_unassured_make_a_check() {
        assert!(LabelResult::Undecided.unassured());
        assert!(LabelResult::NotCovered.unassured());
        assert!(
            !LabelResult::OneRead {
                read: "r".to_owned()
            }
            .unassured()
        );
        assert!(
            !LabelResult::ReadEdge {
                read: "r".to_owned()
            }
            .unassured()
        );
        // SRS-CONS-015: one strand cannot assure an insertion.
        let insertion = LabelResult::OneStrandInsertion {
            strand: F,
            reads: vec!["r".to_owned()],
            other: "not_covered",
            from_start: 25,
        };
        assert!(insertion.unassured());
        assert_eq!(
            (insertion.rank(), insertion.kind().name()),
            (2, "one_strand_insertion")
        );
    }

    // SRS-API-014
    #[test]
    fn summary_counts_items_per_label_and_clears_only_called_samples_without_checks() {
        let variant = item(ItemKind::Variant, Vec::new());
        let one_strand = item(
            ItemKind::Variant,
            vec![
                LabelResult::OneStrand {
                    strand: R,
                    reads: vec![read("a")],
                },
                LabelResult::ReadEdge { read: read("a") },
                LabelResult::ReadEdge { read: read("b") },
            ],
        );
        let summary = summarize(
            &[variant, one_strand],
            ReviewStatus::Called,
            &[[1, 10]],
            &[],
            &[],
            0,
        );
        assert_eq!((summary.items, summary.checks, summary.variants), (2, 0, 2));
        assert_eq!(summary.labels.get(LabelKind::OneStrand), 1);
        assert_eq!(summary.labels.get(LabelKind::ReadEdge), 1);
        assert_eq!(
            summary.positions,
            PositionCounts {
                reportable: 10,
                check: 0,
                reference_only: 0,
                not_analysed: 0
            }
        );
        assert!(summary.clear);

        let check = item(ItemKind::Check, vec![LabelResult::NotCovered]);
        let summary = summarize(&[check], ReviewStatus::Called, &[[1, 10]], &[], &[], 0);
        assert_eq!(
            (summary.checks, summary.labels.get(LabelKind::NotCovered)),
            (1, 1)
        );
        // The check item holds position 1 of the ten reportable positions.
        assert_eq!(
            summary.positions,
            PositionCounts {
                reportable: 10,
                check: 1,
                reference_only: 0,
                not_analysed: 0
            }
        );
        assert!(!summary.clear);
        assert!(!summarize(&[], ReviewStatus::CallFailed, &[[1, 10]], &[], &[], 0).clear);
        assert!(!summarize(&[], ReviewStatus::Rejected, &[[1, 10]], &[], &[], 0).clear);
        assert!(summarize(&[], ReviewStatus::Called, &[[1, 10]], &[], &[], 0).clear);
        // Reference-only positions are counted, and never make a check.
        let summary = summarize(&[], ReviewStatus::Called, &[[1, 10]], &[], &[3, 4], 0);
        assert_eq!((summary.positions.reference_only, summary.clear), (2, true));
        // Positions not analysed are counted, and never make a check
        // (SRS-CONS-012).
        let summary = summarize(&[], ReviewStatus::Called, &[[1, 10]], &[], &[], 6);
        assert_eq!((summary.positions.not_analysed, summary.clear), (6, true));
    }
}
