//! Read eligibility: which calls of a read can support a variant.
//!
//! A mapped call within `read_end_margin` calls of either end of the trim
//! interval, which is the read's callable span, cannot support a variant
//! (`read_end`): beyond the span the alignment has no information, so a
//! difference close to it may be an edge artifact. A masked evidence call — a
//! supporting call, or a flanking call of a deletion, which has none —
//! contributes the reason of its phase segment.

use crate::config::VariantCallingConfig;
use crate::model::callability::{PhaseState, ReadCallability};
use crate::model::quality::QualityControlResult;
use crate::model::variant::{
    VariantCallMapping, VariantCallRole, VariantExclusionReason, VariantKind,
};

/// Trusted calls of one read.
pub(super) struct ReadEligibility<'a> {
    /// First call index outside the leading end margin.
    trusted_start: usize,
    /// End (exclusive) of the calls outside the trailing end margin.
    trusted_end: usize,
    callability: &'a ReadCallability,
}

impl<'a> ReadEligibility<'a> {
    /// Derives the trusted calls of one read from its trim interval and
    /// callability.
    pub(super) fn new(
        quality: &QualityControlResult,
        callability: &'a ReadCallability,
        config: &VariantCallingConfig,
    ) -> Self {
        Self {
            trusted_start: quality
                .trim_start_0based
                .saturating_add(config.read_end_margin),
            trusted_end: quality
                .trim_end_0based_exclusive
                .saturating_sub(config.read_end_margin),
            callability,
        }
    }

    /// The exclusion reason of a masked call, or `None` for a callable one:
    /// `post_homopolymer` when its segment starts right after a repeat run,
    /// otherwise the segment's phase state.
    pub(super) fn mask_reason(&self, index: usize) -> Option<VariantExclusionReason> {
        let segment = self.callability.segment_at(index)?;
        if segment.after_repeat && segment.state != PhaseState::InPhase {
            return Some(VariantExclusionReason::PostHomopolymer);
        }
        match segment.state {
            PhaseState::InPhase => None,
            PhaseState::Dephased => Some(VariantExclusionReason::DephasedSignal),
            PhaseState::Mixed => Some(VariantExclusionReason::MixedSignal),
            PhaseState::Weak => Some(VariantExclusionReason::WeakSignal),
            PhaseState::Irregular => Some(VariantExclusionReason::IrregularSpacing),
        }
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
            mapping.call_index_0based < self.trusted_start
                || mapping.call_index_0based >= self.trusted_end
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
    use crate::model::callability::PhaseSegment;
    use crate::model::quality::{CallQuality, QualityControlResult};

    use super::*;

    fn quality(calls: usize, trim: (usize, usize)) -> QualityControlResult {
        QualityControlResult {
            per_call: (0..calls)
                .map(|index| CallQuality {
                    index_0based: index,
                    penalty: 0,
                    relative_quality_score: 40,
                    vendor_quality_applies: false,
                })
                .collect(),
            trim_start_0based: trim.0,
            trim_end_0based_exclusive: trim.1,
            retained_sequence: String::new(),
        }
    }

    fn settings(read_end_margin: usize) -> VariantCallingConfig {
        VariantCallingConfig {
            max_indel_length: 50,
            minimum_peak_height: 150,
            relative_quality_threshold: 30,
            read_end_margin,
        }
    }

    /// A read with the given segments over `calls` positions.
    fn callability(calls: usize, segments: &[(usize, usize, PhaseState, bool)]) -> ReadCallability {
        let mut read = ReadCallability::in_phase(calls);
        read.segments = segments
            .iter()
            .map(|&(start, end, state, after_repeat)| PhaseSegment {
                call_start_0based: start,
                call_end_0based_exclusive: end,
                state,
                after_repeat,
                shadow: None,
            })
            .collect();
        for &(start, end, state, _) in segments {
            for entry in &mut read.mask[start..end] {
                *entry = (state != PhaseState::InPhase).then_some(state);
            }
        }
        read
    }

    fn mapping(role: VariantCallRole, index: usize) -> VariantCallMapping {
        VariantCallMapping {
            role,
            call_index_0based: index,
            reference_position_0based: Some(index),
        }
    }

    #[test]
    fn flags_evidence_within_the_margin_of_the_trim_interval_ends() {
        let read = ReadCallability::in_phase(32);
        let eligibility = ReadEligibility::new(&quality(32, (5, 25)), &read, &settings(3));
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
    fn checks_supporting_calls_of_insertions_and_flanks_of_deletions() {
        let read = callability(
            20,
            &[
                (0, 10, PhaseState::InPhase, false),
                (10, 20, PhaseState::Dephased, true),
            ],
        );
        let eligibility = ReadEligibility::new(&quality(20, (0, 20)), &read, &settings(0));
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
            [VariantExclusionReason::PostHomopolymer]
        );
    }

    #[test]
    fn reports_the_phase_state_of_a_masked_segment() {
        let read = callability(
            40,
            &[
                (0, 10, PhaseState::Weak, false),
                (10, 20, PhaseState::InPhase, false),
                (20, 25, PhaseState::Dephased, false),
                (25, 30, PhaseState::Mixed, false),
                (30, 35, PhaseState::Irregular, false),
                (35, 40, PhaseState::Mixed, true),
            ],
        );
        let eligibility = ReadEligibility::new(&quality(40, (0, 40)), &read, &settings(0));
        assert_eq!(
            [0, 15, 22, 27, 32, 37].map(|index| eligibility.mask_reason(index)),
            [
                Some(VariantExclusionReason::WeakSignal),
                None,
                Some(VariantExclusionReason::DephasedSignal),
                Some(VariantExclusionReason::MixedSignal),
                Some(VariantExclusionReason::IrregularSpacing),
                Some(VariantExclusionReason::PostHomopolymer),
            ]
        );
    }

    #[test]
    fn reports_each_reason_once_per_variant() {
        let read = callability(
            12,
            &[
                (0, 8, PhaseState::InPhase, false),
                (8, 12, PhaseState::Mixed, false),
            ],
        );
        let eligibility = ReadEligibility::new(&quality(12, (0, 12)), &read, &settings(3));
        let calls = (8..12)
            .map(|index| mapping(VariantCallRole::Supporting, index))
            .collect::<Vec<_>>();
        assert_eq!(
            eligibility.reasons(VariantKind::Ins, &calls),
            [
                VariantExclusionReason::ReadEnd,
                VariantExclusionReason::MixedSignal
            ]
        );
    }
}
