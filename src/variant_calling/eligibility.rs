//! Read eligibility: which calls of a read can support a variant.
//!
//! A call is informative when it lies inside the read's informative interval
//! and is not masked as unresolved; anchoring masked calls still carry
//! alignment information. A mapped call within `read_end_margin` calls of an
//! uninformative call cannot support a variant (`read_end`): the alignment has
//! no information there, so a difference close to it may be an edge artifact.
//! The read's ends count as uninformative unless its modality vouches for them.
//! A masked evidence call — a supporting call, or a flanking call of a
//! deletion, which has none — contributes its modality's mask reason.

use crate::model::variant::{
    VariantCallMapping, VariantCallRole, VariantExclusionReason, VariantKind,
};
use crate::read_evidence::{MaskedAlignment, ReadEnds, ReadEvidence};
use crate::variant_calling::VariantCallingConfig;

/// Trusted calls of one read.
pub(super) struct ReadEligibility<'a> {
    /// Per call: at least `read_end_margin` informative calls separate it from
    /// the nearest uninformative call on both sides.
    trusted: Vec<bool>,
    evidence: &'a ReadEvidence,
}

impl<'a> ReadEligibility<'a> {
    /// Derives the trusted calls of one read from its evidence.
    pub(super) fn new(evidence: &'a ReadEvidence, config: &VariantCallingConfig) -> Self {
        let interval = evidence.informative();
        let informative = evidence
            .calls()
            .iter()
            .enumerate()
            .map(|(index, call)| {
                interval.contains(&index)
                    && call
                        .mask
                        .is_none_or(|mask| mask.alignment == MaskedAlignment::Anchoring)
            })
            .collect::<Vec<_>>();
        let calls = informative.len();
        let margin = config.read_end_margin;
        // Informative run lengths ending at and starting at every call; a
        // vouched read end counts as a full margin of informative calls.
        let edge = match evidence.ends() {
            ReadEnds::Unvouched => 0,
            ReadEnds::Vouched => margin,
        };
        let mut before = vec![0_usize; calls];
        let mut run = edge;
        for index in 0..calls {
            run = if informative[index] { run + 1 } else { 0 };
            before[index] = run;
        }
        let mut after = vec![0_usize; calls];
        run = edge;
        for index in (0..calls).rev() {
            run = if informative[index] { run + 1 } else { 0 };
            after[index] = run;
        }
        let trusted = (0..calls)
            .map(|index| before[index] > margin && after[index] > margin)
            .collect();
        Self { trusted, evidence }
    }

    /// The exclusion reason of a masked call, or `None` for an unmasked one.
    pub(super) fn mask_reason(&self, index: usize) -> Option<VariantExclusionReason> {
        let mask = self.evidence.calls().get(index)?.mask?;
        Some(VariantExclusionReason::Evidence(mask.reason))
    }

    /// Eligibility reasons against a variant's calls, each reported once:
    /// `read_end` for any mapped call, then the mask reasons of its evidence
    /// calls in call order.
    pub(super) fn reasons(
        &self,
        kind: VariantKind,
        mappings: &[VariantCallMapping],
    ) -> Vec<VariantExclusionReason> {
        let evidence_role = if kind == VariantKind::Del {
            VariantCallRole::Flanking
        } else {
            VariantCallRole::Supporting
        };
        let evidence = mappings
            .iter()
            .filter(|mapping| mapping.role == evidence_role)
            .map(|mapping| mapping.call_index_0based)
            .collect::<Vec<_>>();
        let mut reasons = Vec::new();
        if mappings.iter().any(|mapping| {
            !self
                .trusted
                .get(mapping.call_index_0based)
                .copied()
                .unwrap_or(false)
        }) {
            reasons.push(VariantExclusionReason::ReadEnd);
        }
        for reason in evidence.iter().filter_map(|&index| self.mask_reason(index)) {
            if !reasons.contains(&reason) {
                reasons.push(reason);
            }
        }
        reasons
    }
}

#[cfg(test)]
mod tests {
    use crate::read_evidence::{CallMask, EvidenceReason};

    use super::*;

    fn settings(read_end_margin: usize) -> VariantCallingConfig {
        VariantCallingConfig {
            max_indel_length: 50,
            read_end_margin,
        }
    }

    /// Evidence over `calls` positions with masks per `(start, end, alignment, reason)`
    /// and the given informative interval.
    fn evidence(
        calls: usize,
        informative: std::ops::Range<usize>,
        masks: &[(usize, usize, MaskedAlignment, &'static str)],
    ) -> ReadEvidence {
        let mut read = ReadEvidence::clean(&"ACGT".repeat(calls / 4 + 1)[..calls]);
        for &(start, end, alignment, reason) in masks {
            for index in start..end {
                let mut call = read.calls()[index];
                call.mask = Some(CallMask {
                    alignment,
                    reason: EvidenceReason::new(reason),
                });
                read = read.with_call(index, call);
            }
        }
        let calls = read.calls().to_vec();
        ReadEvidence::new(calls, informative, Vec::new()).unwrap_or(read)
    }

    fn mapping(role: VariantCallRole, index: usize) -> VariantCallMapping {
        VariantCallMapping {
            role,
            call_index_0based: index,
            reference_position_0based: Some(index),
        }
    }

    #[test]
    fn flags_calls_within_the_margin_of_the_informative_interval_ends() {
        let read = evidence(32, 5..25, &[]);
        let eligibility = ReadEligibility::new(&read, &settings(3));
        for (index, read_end) in [(7, true), (8, false), (21, false), (22, true)] {
            assert_eq!(
                eligibility.reasons(
                    VariantKind::Snv,
                    &[mapping(VariantCallRole::Supporting, index)]
                ) == [VariantExclusionReason::ReadEnd],
                read_end,
                "call {index}"
            );
        }
    }

    #[test]
    fn trusts_calls_at_vouched_read_ends() {
        let unvouched = evidence(16, 0..16, &[]);
        let vouched = unvouched.clone().with_vouched_ends();
        let masked = evidence(
            16,
            0..16,
            &[(4, 5, MaskedAlignment::Unresolved, "weak_signal")],
        )
        .with_vouched_ends();
        let snv = |read: &ReadEvidence, index| {
            ReadEligibility::new(read, &settings(3))
                .reasons(
                    VariantKind::Snv,
                    &[mapping(VariantCallRole::Supporting, index)],
                )
                .is_empty()
        };
        assert!(!snv(&unvouched, 0));
        assert!(!snv(&unvouched, 15));
        assert!(snv(&vouched, 0));
        assert!(snv(&vouched, 15));
        assert!(snv(&masked, 0));
        assert!(!snv(&masked, 7));
        assert!(snv(&masked, 8));
    }

    #[test]
    fn checks_supporting_calls_of_insertions_and_flanks_of_deletions() {
        let read = evidence(
            20,
            0..20,
            &[(10, 20, MaskedAlignment::Anchoring, "post_homopolymer")],
        );
        let eligibility = ReadEligibility::new(&read, &settings(0));
        let calls = [
            mapping(VariantCallRole::Supporting, 9),
            mapping(VariantCallRole::Flanking, 10),
        ];
        assert!(eligibility.reasons(VariantKind::Ins, &calls).is_empty());
        let flanks = [
            mapping(VariantCallRole::Flanking, 9),
            mapping(VariantCallRole::Flanking, 10),
        ];
        assert_eq!(
            eligibility.reasons(VariantKind::Del, &flanks),
            [VariantExclusionReason::Evidence(EvidenceReason::new(
                "post_homopolymer"
            ))]
        );
    }

    #[test]
    fn reports_the_modality_reason_of_a_masked_call() {
        let read = evidence(
            20,
            0..20,
            &[
                (0, 10, MaskedAlignment::Unresolved, "weak_signal"),
                (15, 20, MaskedAlignment::Anchoring, "dephased_signal"),
            ],
        );
        let eligibility = ReadEligibility::new(&read, &settings(0));
        assert_eq!(
            [0, 12, 17].map(|index| eligibility.mask_reason(index)),
            [
                Some(VariantExclusionReason::Evidence(EvidenceReason::new(
                    "weak_signal"
                ))),
                None,
                Some(VariantExclusionReason::Evidence(EvidenceReason::new(
                    "dephased_signal"
                ))),
            ]
        );
    }

    #[test]
    fn flags_evidence_near_an_unresolved_mask_but_not_an_anchoring_one() {
        let read = evidence(
            40,
            0..40,
            &[
                (15, 20, MaskedAlignment::Unresolved, "mixed_signal"),
                (30, 35, MaskedAlignment::Anchoring, "dephased_signal"),
            ],
        );
        let eligibility = ReadEligibility::new(&read, &settings(3));
        let read_end = |index| {
            eligibility
                .reasons(
                    VariantKind::Snv,
                    &[mapping(VariantCallRole::Supporting, index)],
                )
                .contains(&VariantExclusionReason::ReadEnd)
        };
        for (index, expected) in [
            (11, false),
            (12, true),
            (22, true),
            (23, false),
            (29, false),
        ] {
            assert_eq!(read_end(index), expected, "call {index}");
        }
    }

    #[test]
    fn reports_each_reason_once_per_variant() {
        let read = evidence(
            12,
            0..12,
            &[(8, 12, MaskedAlignment::Unresolved, "mixed_signal")],
        );
        let eligibility = ReadEligibility::new(&read, &settings(3));
        let calls = (8..12)
            .map(|index| mapping(VariantCallRole::Supporting, index))
            .collect::<Vec<_>>();
        assert_eq!(
            eligibility.reasons(VariantKind::Ins, &calls),
            [
                VariantExclusionReason::ReadEnd,
                VariantExclusionReason::Evidence(EvidenceReason::new("mixed_signal"))
            ]
        );
    }
}
