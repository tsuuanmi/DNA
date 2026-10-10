//! Read eligibility: which calls of a read can support a variant.
//!
//! A call is informative when it lies inside the read's informative interval
//! and is not masked as unresolved; anchoring masked calls still carry
//! alignment information. A mapped call within `read_end_margin` calls of an
//! uninformative call cannot support a variant (`read_end`): the alignment has
//! no information there, so a difference close to it may be an edge artifact.
//! The read's ends count as uninformative unless its modality vouches for them.
//! An edit that changes the length of a run of one base can be judged only
//! where the read resolves that run's end (`run_boundary`): the call beyond it
//! is informative and reads another base, unmasked or anchoring. A masked
//! evidence call — a supporting call, or a flanking call of a deletion, which
//! has none — contributes its modality's mask reason.

use crate::model::alignment::Orientation;
use crate::model::variant::{Variant, VariantCallRole, VariantExclusionReason, VariantKind};
use crate::variant_calling::VariantCallingConfig;
use dna_kernel::read_evidence::{MaskedAlignment, ReadEnds, ReadEvidence};

/// Trusted calls of one read.
pub(crate) struct ReadEligibility<'a> {
    /// Per call: at least `read_end_margin` informative calls separate it from
    /// the nearest uninformative call on both sides.
    trusted: Vec<bool>,
    evidence: &'a ReadEvidence,
}

impl<'a> ReadEligibility<'a> {
    /// Derives the trusted calls of one read from its evidence.
    pub(crate) fn new(evidence: &'a ReadEvidence, config: &VariantCallingConfig) -> Self {
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

    /// Whether `variant` changes the length of a run of one base at an end of
    /// that run the read does not resolve.
    ///
    /// The run is the read's maximal stretch of the edited base around the
    /// edit: the substituted base, the inserted bases, or what a deletion
    /// leaves of a run. An insertion or deletion needs both ends resolved and
    /// every call of the run unmasked; a substitution needs the ends it lies
    /// on. Edits of more than one distinct base, a lone substituted or
    /// inserted base, and a deletion of a whole run change no run's length.
    fn unbounded_run(
        &self,
        variant: &Variant,
        orientation: Orientation,
        evidence: &[usize],
    ) -> bool {
        let calls = self.evidence.calls();
        let edited = match variant.kind {
            VariantKind::Snv => &variant.alternate[..],
            VariantKind::Ins => variant.alternate.get(1..).unwrap_or_default(),
            VariantKind::Del => variant.reference.get(1..).unwrap_or_default(),
        };
        let mut bases = edited.chars();
        let Some(reference_base) = bases.next() else {
            return false;
        };
        if bases.any(|other| other != reference_base) {
            return false;
        }
        // The calls read the trace strand; complementing is its own inverse.
        let base = orientation.reference_base(reference_base);
        let (Some(&first), Some(&last)) = (evidence.iter().min(), evidence.iter().max()) else {
            return false;
        };
        let reads = |index: usize| calls.get(index).is_some_and(|call| call.base == base);
        // A deletion's run is the rest of the run, read by a flank.
        let (mut start, mut end) = match variant.kind {
            VariantKind::Del if reads(first) => (first, first),
            VariantKind::Del if reads(last) => (last, last),
            VariantKind::Del => return false,
            VariantKind::Snv | VariantKind::Ins => (first, last),
        };
        while start > 0 && reads(start - 1) {
            start -= 1;
        }
        while reads(end + 1) {
            end += 1;
        }
        // A lone substituted or inserted base starts no run.
        if variant.kind != VariantKind::Del && start == end {
            return false;
        }
        let before = start
            .checked_sub(1)
            .map_or(self.vouched(), |index| self.resolves(index, base));
        let after = if end + 1 < calls.len() {
            self.resolves(end + 1, base)
        } else {
            self.vouched()
        };
        if variant.kind == VariantKind::Snv {
            return (first == start && !before) || (first == end && !after);
        }
        !(before && after && calls[start..=end].iter().all(|call| call.mask.is_none()))
    }

    /// Whether the read's physical ends bound its evidence (ADR-0069).
    fn vouched(&self) -> bool {
        self.evidence.ends() == ReadEnds::Vouched
    }

    /// Whether the call at `index` resolves the end of a run of `base`: it is
    /// informative and reads another resolved base, unmasked or anchoring.
    fn resolves(&self, index: usize, base: char) -> bool {
        let Some(call) = self.evidence.calls().get(index) else {
            return false;
        };
        self.evidence.informative().contains(&index)
            && call.base != base
            && call.base != 'N'
            && call
                .mask
                .is_none_or(|mask| mask.alignment == MaskedAlignment::Anchoring)
    }

    /// Whether the call at `index` has at least `read_end_margin` informative
    /// calls between it and the nearest uninformative call on both sides.
    pub(crate) fn trusted(&self, index: usize) -> bool {
        self.trusted.get(index).copied().unwrap_or(false)
    }

    /// The exclusion reason of a masked call, or `None` for an unmasked one.
    pub(super) fn mask_reason(&self, index: usize) -> Option<VariantExclusionReason> {
        let mask = self.evidence.calls().get(index)?.mask?;
        Some(VariantExclusionReason::Evidence(mask.reason))
    }

    /// Eligibility reasons against a variant's calls, each reported once:
    /// `read_end` for any mapped call, `run_boundary` for a run-length edit
    /// the read does not bound, then the mask reasons of its evidence calls in
    /// call order.
    ///
    /// `orientation` is the read's placement, which relates the variant's
    /// reference-strand alleles to the read's calls.
    pub(super) fn reasons(
        &self,
        variant: &Variant,
        orientation: Orientation,
    ) -> Vec<VariantExclusionReason> {
        let kind = variant.kind;
        let mappings = &variant.calls;
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
        if self.unbounded_run(variant, orientation, &evidence) {
            reasons.push(VariantExclusionReason::RunBoundary);
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
    use crate::model::variant::VariantCallMapping;
    use dna_kernel::read_evidence::{CallMask, EvidenceReason};

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

    /// A variant over `calls` whose alleles the read's own bases support; a
    /// deletion removes a base neither flank reads.
    fn variant(read: &ReadEvidence, kind: VariantKind, calls: &[VariantCallMapping]) -> Variant {
        let bases = calls
            .iter()
            .map(|mapping| read.calls()[mapping.call_index_0based].base)
            .collect::<String>();
        let supporting = calls
            .iter()
            .filter(|mapping| mapping.role == VariantCallRole::Supporting)
            .map(|mapping| read.calls()[mapping.call_index_0based].base)
            .collect::<String>();
        let (reference, alternate) = match kind {
            VariantKind::Snv => ("N".to_owned(), supporting),
            VariantKind::Ins => ("N".to_owned(), format!("N{supporting}")),
            VariantKind::Del => {
                let deleted = "ACGT"
                    .chars()
                    .find(|base| !bases.contains(*base))
                    .unwrap_or('A');
                (format!("N{deleted}"), "N".to_owned())
            }
        };
        Variant {
            contig: "test".to_owned(),
            position_1based: 1,
            reference,
            alternate,
            kind,
            calls: calls.to_vec(),
        }
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
                    &variant(
                        &read,
                        VariantKind::Snv,
                        &[mapping(VariantCallRole::Supporting, index)]
                    ),
                    Orientation::Forward
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
                    &variant(
                        read,
                        VariantKind::Snv,
                        &[mapping(VariantCallRole::Supporting, index)],
                    ),
                    Orientation::Forward,
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
        assert!(
            eligibility
                .reasons(
                    &variant(&read, VariantKind::Ins, &calls),
                    Orientation::Forward
                )
                .is_empty()
        );
        let flanks = [
            mapping(VariantCallRole::Flanking, 9),
            mapping(VariantCallRole::Flanking, 10),
        ];
        assert_eq!(
            eligibility.reasons(
                &variant(&read, VariantKind::Del, &flanks),
                Orientation::Forward
            ),
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
                    &variant(
                        &read,
                        VariantKind::Snv,
                        &[mapping(VariantCallRole::Supporting, index)],
                    ),
                    Orientation::Forward,
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
            eligibility.reasons(
                &variant(&read, VariantKind::Ins, &calls),
                Orientation::Forward
            ),
            [
                VariantExclusionReason::ReadEnd,
                VariantExclusionReason::Evidence(EvidenceReason::new("mixed_signal"))
            ]
        );
    }

    /// A read of `sequence` with every call informative and the given masks.
    fn sequence(sequence: &str, masks: &[(usize, MaskedAlignment)]) -> ReadEvidence {
        let mut read = ReadEvidence::clean(sequence);
        for &(index, alignment) in masks {
            let mut call = read.calls()[index];
            call.mask = Some(CallMask {
                alignment,
                reason: EvidenceReason::new("dephased_signal"),
            });
            read = read.with_call(index, call);
        }
        read
    }

    fn run_boundary(read: &ReadEvidence, kind: VariantKind, calls: &[VariantCallMapping]) -> bool {
        ReadEligibility::new(read, &settings(0))
            .reasons(&variant(read, kind, calls), Orientation::Forward)
            .contains(&VariantExclusionReason::RunBoundary)
    }

    #[test]
    fn needs_both_ends_of_a_run_an_insertion_lengthens() {
        let inserted = [mapping(VariantCallRole::Supporting, 4)];
        // The run CCCCC ends at an anchoring call that still reads T.
        let resolved = sequence("AGTACCCCCTGA", &[(9, MaskedAlignment::Anchoring)]);
        assert!(!run_boundary(&resolved, VariantKind::Ins, &inserted));
        // A masked call that reads the run's base hides where the run ends.
        let hidden = sequence("AGTACCCCCCGA", &[(9, MaskedAlignment::Anchoring)]);
        assert!(run_boundary(&hidden, VariantKind::Ins, &inserted));
        // An unresolved end bounds nothing.
        let unresolved = sequence("AGTACCCCCTGA", &[(9, MaskedAlignment::Unresolved)]);
        assert!(run_boundary(&unresolved, VariantKind::Ins, &inserted));
        // A masked call inside the run leaves its length unread.
        let masked = sequence("AGTACCCCCTGA", &[(6, MaskedAlignment::Anchoring)]);
        assert!(run_boundary(&masked, VariantKind::Ins, &inserted));
    }

    #[test]
    fn reads_reference_alleles_on_the_trace_strand_of_a_reverse_read() {
        // The trace run GGGGG is a reference run of C on a reverse read; its end
        // hides behind a masked G.
        let read = sequence("AGTAGGGGGGCT", &[(9, MaskedAlignment::Anchoring)]);
        let inserted = [mapping(VariantCallRole::Supporting, 4)];
        let insertion = Variant {
            reference: "N".to_owned(),
            alternate: "NC".to_owned(),
            ..variant(&read, VariantKind::Ins, &inserted)
        };
        let flagged = |orientation| {
            ReadEligibility::new(&read, &settings(0))
                .reasons(&insertion, orientation)
                .contains(&VariantExclusionReason::RunBoundary)
        };
        assert!(flagged(Orientation::Reverse));
        assert!(!flagged(Orientation::Forward));
    }

    #[test]
    fn needs_only_the_run_end_a_substitution_lies_on() {
        // Calls 4..=9 read CCCCCC; the masked call 10 resolves no base.
        let read = sequence("AGTACCCCCCNGA", &[(10, MaskedAlignment::Anchoring)]);
        let snv = |index| [mapping(VariantCallRole::Supporting, index)];
        assert!(!run_boundary(&read, VariantKind::Snv, &snv(4)));
        assert!(!run_boundary(&read, VariantKind::Snv, &snv(7)));
        let extended = sequence("AGTACCCCCCAGA", &[(10, MaskedAlignment::Anchoring)]);
        assert!(!run_boundary(&extended, VariantKind::Snv, &snv(9)));
        assert!(run_boundary(&read, VariantKind::Snv, &snv(9)));
    }

    #[test]
    fn bounds_what_a_deletion_leaves_of_a_run() {
        let flanks = [
            mapping(VariantCallRole::Flanking, 5),
            mapping(VariantCallRole::Flanking, 6),
        ];
        // A deleted C whose remaining run CC ends at a masked C.
        let mut read = sequence("AGTACCCGA", &[(6, MaskedAlignment::Anchoring)]);
        let deletion = |read: &ReadEvidence, deleted: char| Variant {
            reference: format!("N{deleted}"),
            alternate: "N".to_owned(),
            ..variant(read, VariantKind::Del, &flanks)
        };
        let reasons = |read: &ReadEvidence, deleted| {
            ReadEligibility::new(read, &settings(0))
                .reasons(&deletion(read, deleted), Orientation::Forward)
        };
        assert!(reasons(&read, 'C').contains(&VariantExclusionReason::RunBoundary));
        read = sequence("AGTACCTGA", &[(6, MaskedAlignment::Anchoring)]);
        assert!(!reasons(&read, 'C').contains(&VariantExclusionReason::RunBoundary));
        // A whole run deleted leaves no run to measure.
        read = sequence("AGTACTGA", &[(5, MaskedAlignment::Unresolved)]);
        assert!(!reasons(&read, 'A').contains(&VariantExclusionReason::RunBoundary));
    }

    #[test]
    fn leaves_a_lone_edited_base_to_the_other_rules() {
        // A substituted T between C and an unresolved call lengthens no run.
        let read = sequence("AGTACTNGA", &[(6, MaskedAlignment::Unresolved)]);
        let snv = [mapping(VariantCallRole::Supporting, 5)];
        assert!(!run_boundary(&read, VariantKind::Snv, &snv));
        let inserted = [mapping(VariantCallRole::Supporting, 5)];
        assert!(!run_boundary(&read, VariantKind::Ins, &inserted));
    }

    #[test]
    fn lets_vouched_ends_bound_a_run_and_ignores_mixed_edits() {
        let read = sequence("CCCCAGT", &[]);
        let inserted = [mapping(VariantCallRole::Supporting, 1)];
        assert!(run_boundary(&read, VariantKind::Ins, &inserted));
        assert!(!run_boundary(
            &read.clone().with_vouched_ends(),
            VariantKind::Ins,
            &inserted
        ));
        let mixed = [
            mapping(VariantCallRole::Supporting, 3),
            mapping(VariantCallRole::Supporting, 4),
        ];
        assert!(!run_boundary(&read, VariantKind::Ins, &mixed));
        // `run_boundary` follows `read_end`.
        assert_eq!(
            ReadEligibility::new(&read, &settings(3)).reasons(
                &variant(&read, VariantKind::Ins, &inserted),
                Orientation::Forward
            ),
            [
                VariantExclusionReason::ReadEnd,
                VariantExclusionReason::RunBoundary
            ]
        );
    }
}
