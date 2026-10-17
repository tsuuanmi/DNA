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
//! has none — contributes its modality's mask reason. A mapped call read as
//! `N` places nothing (`unresolved_flank`): an `N` beside a deletion leaves
//! which bases it removes undecided.
//!
//! A run the read bounds is exempt from `read_end` where the margin's edge is
//! the read losing phase right after a repeat run, not a read end (ADR-0076):
//! an insertion or deletion of that run is not `read_end`, and the consensus
//! trusts the run's matches and insertions; a substituted call or a bound
//! never gains trust.

use crate::model::alignment::Orientation;
use crate::model::trust::{CallTrust, Withholder};
use crate::model::variant::{Variant, VariantCallRole, VariantExclusionReason, VariantKind};
use crate::variant_calling::VariantCallingConfig;
use dna_kernel::read_evidence::{MaskOnset, MaskedAlignment, ReadEnds, ReadEvidence, VetoScope};

/// Informative calls between a run's bound and a post-repeat edge for the
/// run to be exempt from `read_end` (ADR-0076).
const EXEMPTION_MINIMUM_CALLS: usize = 2;

/// Trusted calls of one read.
pub(crate) struct ReadEligibility<'a> {
    /// Per call: at least `read_end_margin` informative calls separate it from
    /// the nearest uninformative call on both sides.
    trusted: Vec<bool>,
    /// Per call: informative calls ending at it, and starting at it (a
    /// vouched read end counts as a full margin).
    before: Vec<usize>,
    after: Vec<usize>,
    margin: usize,
    /// Per call: it lies on a bounded run of one base that the exemption lifts
    /// from `read_end` (the run's bounds are not included).
    run_exempt: Vec<bool>,
    /// Per call: it bounds such a run.
    bound_exempt: Vec<bool>,
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
        let mut eligibility = Self {
            trusted,
            before,
            after,
            margin,
            run_exempt: vec![false; calls],
            bound_exempt: vec![false; calls],
            evidence,
        };
        {
            let calls = evidence.calls();
            let mut start = 0;
            while start < calls.len() {
                let base = calls[start].base;
                let mut end = start;
                while end + 1 < calls.len() && calls[end + 1].base == base {
                    end += 1;
                }
                if end > start
                    && dna_kernel::model::nucleotide::is_canonical(base)
                    && eligibility.bounded_run_exempt(start, end, base)
                {
                    // The run's own calls; its bounds are recorded apart, since
                    // the consensus trusts them only where they read the
                    // reference base.
                    eligibility.run_exempt[start..=end].fill(true);
                    eligibility.bound_exempt[start - 1] = true;
                    eligibility.bound_exempt[end + 1] = true;
                }
                start = end + 1;
            }
        }
        eligibility
    }

    /// Whether a call lies on a bounded run the exemption lifts from
    /// `read_end`. The consensus trusts such a call only where it matches
    /// the reference or is inserted: the run's length may change, its
    /// composition may not.
    pub(crate) fn run_exempt(&self, index: usize) -> bool {
        self.run_exempt.get(index).copied().unwrap_or(false)
    }

    /// Whether a call bounds a run the exemption lifts. The consensus trusts
    /// such a bound only where it reads the reference base, so it can never
    /// make a variant.
    pub(crate) fn bound_exempt(&self, index: usize) -> bool {
        self.bound_exempt.get(index).copied().unwrap_or(false)
    }

    /// Whether the run `start..=end` of `base` is bounded by resolving,
    /// unmasked calls with every run call unmasked, and each side's nearest
    /// uninformative call within the margin starts a masked segment right
    /// after a repeat run, at least [`EXEMPTION_MINIMUM_CALLS`] informative
    /// calls from the bound.
    fn bounded_run_exempt(&self, start: usize, end: usize, base: char) -> bool {
        let calls = self.evidence.calls();
        let (Some(left), right) = (start.checked_sub(1), end + 1) else {
            return false;
        };
        if right >= calls.len()
            || !self.resolves(left, base)
            || !self.resolves(right, base)
            || calls[left..=right].iter().any(|call| call.mask.is_some())
        {
            return false;
        }
        let allowed = |index: Option<usize>| {
            index
                .and_then(|index| calls.get(index))
                .and_then(|call| call.mask)
                .is_some_and(|mask| mask.onset == MaskOnset::AfterRun)
        };
        // The nearest uninformative call on each side; `None` is the read's
        // physical end, which never exempts.
        let side = |run: usize, edge: Option<usize>| {
            run > self.margin || (run >= EXEMPTION_MINIMUM_CALLS && allowed(edge))
        };
        let exempt = side(self.before[left], left.checked_sub(self.before[left]))
            && side(
                self.after[right],
                Some(right + self.after[right]).filter(|index| *index < calls.len()),
            );
        // Only a run the margin actually distrusts needs the exemption.
        exempt && (left..=right).any(|index| !self.trusted(index))
    }

    /// The read's maximal run of the edited base around an edit, as
    /// [`Self::unbounded_run`] defines it: its first and last call, and the
    /// trace-strand base. `None` when the edit changes no run's length.
    fn run(
        &self,
        kind: VariantKind,
        edited: &str,
        orientation: Orientation,
        evidence: &[usize],
    ) -> Option<(usize, usize, char)> {
        let calls = self.evidence.calls();
        let mut bases = edited.chars();
        let reference_base = bases.next()?;
        if bases.any(|other| other != reference_base) {
            return None;
        }
        // The calls read the trace strand; complementing is its own inverse.
        let base = orientation.reference_base(reference_base);
        let (first, last) = (*evidence.iter().min()?, *evidence.iter().max()?);
        let reads = |index: usize| calls.get(index).is_some_and(|call| call.base == base);
        // A deletion's run is the rest of the run, read by a flank.
        let (mut start, mut end) = match kind {
            VariantKind::Del if reads(first) => (first, first),
            VariantKind::Del if reads(last) => (last, last),
            VariantKind::Del => return None,
            VariantKind::Snv | VariantKind::Ins => (first, last),
        };
        while start > 0 && reads(start - 1) {
            start -= 1;
        }
        while reads(end + 1) {
            end += 1;
        }
        // A lone substituted or inserted base starts no run.
        (kind == VariantKind::Del || start != end).then_some((start, end, base))
    }

    /// Whether the exemption lifts `read_end` from an edit: an
    /// insertion or deletion of a bounded run (see `bounded_run_exempt`)
    /// whose mapped calls lie on the run and its bounds.
    fn exempt(
        &self,
        kind: VariantKind,
        edited: &str,
        orientation: Orientation,
        evidence: &[usize],
        mapped: &[usize],
    ) -> bool {
        if kind == VariantKind::Snv {
            return false;
        }
        let Some((start, end, base)) = self.run(kind, edited, orientation, evidence) else {
            return false;
        };
        let Some(left) = start.checked_sub(1) else {
            return false;
        };
        mapped.iter().all(|index| (left..=end + 1).contains(index))
            && self.bounded_run_exempt(start, end, base)
    }

    /// Whether an edit of `kind` changes the length of a run of one base at an
    /// end of that run the read does not resolve.
    ///
    /// `edited` holds the edit's reference-strand bases: the substituted base,
    /// the inserted bases, or the deleted ones. `evidence` holds the edit's
    /// calls: the supporting calls, or a deletion's flanking calls. The run is
    /// the read's maximal stretch of the edited base around the edit: the
    /// substituted base, the inserted bases, or what a deletion leaves of a
    /// run. An insertion or deletion needs both ends resolved and every call
    /// of the run unmasked; a substitution needs the ends it lies on. Edits of
    /// more than one distinct base, a lone substituted or inserted base, and a
    /// deletion of a whole run change no run's length.
    pub(crate) fn unbounded_run(
        &self,
        kind: VariantKind,
        edited: &str,
        orientation: Orientation,
        evidence: &[usize],
    ) -> bool {
        let calls = self.evidence.calls();
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
        let (mut start, mut end) = match kind {
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
        if kind != VariantKind::Del && start == end {
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
        if kind == VariantKind::Snv {
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

    /// Whether an insertion whose inserted calls are `inserted` has a
    /// placement along the read's repeat of its inserted unit whose calls
    /// raise no support veto that applies to insertions. In a repeat, which
    /// calls the alignment takes as the inserted ones is arbitrary: they are
    /// those of its 3'-most placement (ADR-0047), and any window of as many
    /// calls of the repeat describes the same sequence. So a veto on the calls
    /// that happen to be taken must not decide the insertion's support. `false`
    /// when the inserted calls are not contiguous or extend no repeat.
    pub(crate) fn repeat_unvetoed(&self, inserted: &[usize]) -> bool {
        let calls = self.evidence.calls();
        let (Some(&first), Some(&last)) = (inserted.iter().min(), inserted.iter().max()) else {
            return false;
        };
        let length = last - first + 1;
        if length != inserted.len() {
            return false;
        }
        let base = |index: usize| {
            calls
                .get(index)
                .map(|call| call.base)
                .filter(|base| matches!(base, 'A' | 'C' | 'G' | 'T'))
        };
        let Some(period) = (1..=length).find(|period| {
            length % period == 0
                && (first..=last).all(|index| {
                    base(index).is_some() && base(index) == base(first + (index - first) % period)
                })
        }) else {
            return false;
        };
        let mut start = first;
        while start > 0 && base(start - 1).is_some() && base(start - 1) == base(start - 1 + period)
        {
            start -= 1;
        }
        let mut end = last + 1;
        while end < calls.len() && base(end).is_some() && base(end) == base(end - period) {
            end += 1;
        }
        if end - start <= length {
            return false;
        }
        let vetoes = self
            .evidence
            .support_vetoes()
            .iter()
            .enumerate()
            .filter(|(_, veto)| veto.scope == VetoScope::SubstitutionsAndInsertions)
            .map(|(index, _)| index)
            .collect::<Vec<_>>();
        (start..=end - length).any(|offset| {
            (offset..offset + length).all(|index| {
                calls[index].mask.is_none()
                    && !vetoes
                        .iter()
                        .any(|&veto| calls[index].vetoes.contains(veto))
            })
        })
    }

    /// What withholds trust from the call at `index`, before then after it
    /// in read order: on each side whose informative run is within the
    /// margin, the nearest uninformative call or the read's own end. Empty
    /// when the call is trusted; an uninformative call is withheld by itself.
    pub(crate) fn withholders(&self, index: usize) -> Vec<Withholder> {
        let calls = self.evidence.calls();
        let interval = self.evidence.informative();
        let withholder = |call: usize| {
            if interval.contains(&call) {
                calls[call]
                    .mask
                    .map_or(Withholder::Untrimmed, Withholder::Masked)
            } else {
                Withholder::Untrimmed
            }
        };
        let Some(&before) = self.before.get(index) else {
            return Vec::new();
        };
        if before == 0 {
            return vec![withholder(index)];
        }
        let mut withheld = Vec::new();
        if before <= self.margin {
            withheld.push(
                index
                    .checked_sub(before)
                    .map_or(Withholder::ReadEnd, withholder),
            );
        }
        let after = self.after[index];
        if after <= self.margin {
            let beyond = index + after;
            withheld.push(if beyond < calls.len() {
                withholder(beyond)
            } else {
                Withholder::ReadEnd
            });
        }
        withheld
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
    /// the read does not bound, `unresolved_flank` for a mapped call read as
    /// `N`, then the mask reasons of its evidence calls in call order.
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
        let edited = match kind {
            VariantKind::Snv => &variant.alternate[..],
            VariantKind::Ins => variant.alternate.get(1..).unwrap_or_default(),
            VariantKind::Del => variant.reference.get(1..).unwrap_or_default(),
        };
        let mapped = mappings
            .iter()
            .map(|mapping| mapping.call_index_0based)
            .collect::<Vec<_>>();
        if mapped.iter().any(|index| !self.trusted(*index))
            && !self.exempt(kind, edited, orientation, &evidence, &mapped)
        {
            reasons.push(VariantExclusionReason::ReadEnd);
        }
        if self.unbounded_run(kind, edited, orientation, &evidence) {
            reasons.push(VariantExclusionReason::RunBoundary);
        }
        if mapped.iter().any(|&index| {
            self.evidence
                .calls()
                .get(index)
                .is_some_and(|call| call.base == 'N')
        }) {
            reasons.push(VariantExclusionReason::UnresolvedFlank);
        }
        for reason in evidence.iter().filter_map(|&index| self.mask_reason(index)) {
            if !reasons.contains(&reason) {
                reasons.push(reason);
            }
        }
        reasons
    }
}

/// The core's trust in every call of a read, as its variant calling and
/// consensus judge it (diagnostic evidence, not a calling contract).
#[must_use]
pub fn read_trust(evidence: &ReadEvidence, config: &VariantCallingConfig) -> Vec<CallTrust> {
    let eligibility = ReadEligibility::new(evidence, config);
    let vocabulary = evidence.support_vetoes();
    evidence
        .calls()
        .iter()
        .enumerate()
        .map(|(index, call)| CallTrust {
            mask: call.mask,
            unresolved: call.base == 'N',
            vetoes: vocabulary
                .iter()
                .enumerate()
                .filter(|(veto, _)| call.vetoes.contains(*veto))
                .map(|(_, veto)| veto.reason)
                .collect(),
            trusted: eligibility.trusted(index),
            exempt: eligibility.run_exempt(index) || eligibility.bound_exempt(index),
            withheld_by: if eligibility.trusted(index) {
                Vec::new()
            } else {
                eligibility.withholders(index)
            },
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use crate::model::variant::VariantCallMapping;
    use dna_kernel::read_evidence::{CallMask, EvidenceReason, MaskOnset};

    use super::*;

    /// A read of `sequence` whose calls in `vetoed` carry the peak veto,
    /// which applies to insertions, or in `mixed` the mixed veto, which does
    /// not.
    fn vetoed_read(sequence: &str, vetoed: &[usize], mixed: &[usize]) -> ReadEvidence {
        use dna_kernel::read_evidence::{EvidenceReason, SupportVeto, VetoScope};
        let vocabulary = vec![
            SupportVeto {
                reason: EvidenceReason::new("peak_below_minimum"),
                scope: VetoScope::SubstitutionsAndInsertions,
            },
            SupportVeto {
                reason: EvidenceReason::new("mixed_supporting_dna"),
                scope: VetoScope::Substitutions,
            },
        ];
        let mut calls = ReadEvidence::clean(sequence).calls().to_vec();
        for &index in vetoed {
            calls[index].vetoes.insert(0);
        }
        for &index in mixed {
            calls[index].vetoes.insert(1);
        }
        let length = calls.len();
        ReadEvidence::new(calls, 0..length, vocabulary)
            .unwrap_or_else(|_| ReadEvidence::clean(sequence))
    }

    #[test]
    fn names_what_withholds_trust_on_each_side() {
        use crate::model::trust::Withholder;
        // Calls 2-19 informative, 10-11 masked as unresolved; a margin of 3.
        let read = evidence(
            20,
            2..20,
            &[(10, 12, MaskedAlignment::Unresolved, "dephased_signal")],
        );
        let eligibility = ReadEligibility::new(&read, &settings(3));
        assert_eq!(eligibility.withholders(4), [Withholder::Untrimmed]);
        let masked = |index: usize| match eligibility.withholders(index).as_slice() {
            [Withholder::Masked(mask)] => Some(mask.reason.label()),
            _ => None,
        };
        assert_eq!(masked(8), Some("dephased_signal"));
        assert_eq!(masked(11), Some("dephased_signal"));
        assert_eq!(eligibility.withholders(18), [Withholder::ReadEnd]);
        assert!(eligibility.withholders(6).is_empty() || eligibility.trusted(6));
        let trust = read_trust(&read, &settings(3));
        assert!(trust[6].trusted && trust[6].withheld_by.is_empty());
        assert!(!trust[18].trusted);
    }

    #[test]
    fn an_insertion_into_a_repeat_is_unvetoed_if_some_placement_is() {
        let settings = settings(0);
        // GCCCCCCT: one C inserted into a run of five, taken as call 6.
        let read = vetoed_read("GCCCCCCT", &[6], &[]);
        assert!(ReadEligibility::new(&read, &settings).repeat_unvetoed(&[6]));
        // Every C vetoed: no placement is clean.
        let read = vetoed_read("GCCCCCCT", &[1, 2, 3, 4, 5, 6], &[]);
        assert!(!ReadEligibility::new(&read, &settings).repeat_unvetoed(&[6]));
        // A mixed veto does not apply to insertions.
        let read = vetoed_read("GCCCCCCT", &[], &[6]);
        assert!(ReadEligibility::new(&read, &settings).repeat_unvetoed(&[6]));
        // An AC unit inserted into GACACACT, taken as calls 5-6: a clean
        // window of two calls of the repeat exists.
        let read = vetoed_read("GACACACT", &[5, 6], &[]);
        assert!(ReadEligibility::new(&read, &settings).repeat_unvetoed(&[5, 6]));
        // An inserted base that extends no repeat is judged on itself.
        let read = vetoed_read("GACT", &[2], &[]);
        assert!(!ReadEligibility::new(&read, &settings).repeat_unvetoed(&[2]));
    }

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
                    onset: MaskOnset::Other,
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

    /// A read `ACGTACGTAC G CCCCC T AG`, then a masked segment of `onset` from
    /// call 19 (or the read's end there when `onset` is `None`), then 13
    /// clean calls.
    fn run_read(onset: Option<MaskOnset>) -> ReadEvidence {
        let tail = if onset.is_some() {
            "AAAAAACGTACGTACGTACG"
        } else {
            ""
        };
        let mut read = ReadEvidence::clean(&format!("ACGTACGTACGCCCCCTAG{tail}"));
        if let Some(onset) = onset {
            for index in 19..26 {
                let mut call = read.calls()[index];
                call.mask = Some(CallMask {
                    alignment: MaskedAlignment::Unresolved,
                    reason: EvidenceReason::new("post_homopolymer"),
                    onset,
                });
                read = read.with_call(index, call);
            }
        }
        read
    }

    #[test]
    fn a_bounded_run_before_a_post_repeat_edge_is_exempt_but_not_its_bounds() {
        let read = run_read(Some(MaskOnset::AfterRun));
        let eligibility = ReadEligibility::new(&read, &settings(8));
        // The run lies within the margin of the edge, but is exempt.
        assert!((11..=15).all(|index| !eligibility.trusted(index)));
        assert!((11..=15).all(|index| eligibility.run_exempt(index)));
        // The bound right of the run keeps its own trust: never a clean base.
        assert!(!eligibility.run_exempt(16) && !eligibility.trusted(16));
        // Too few informative calls between the bound and the edge: the
        // masked segment reaches the call right after the bound.
        let mut close = read.clone();
        for index in [17, 18] {
            let mut call = close.calls()[index];
            call.mask = read.calls()[19].mask;
            close = close.with_call(index, call);
        }
        assert!(!ReadEligibility::new(&close, &settings(8)).run_exempt(13));
    }

    #[test]
    fn other_edges_never_exempt_a_run() {
        // Another onset, or the read's physical end.
        for onset in [Some(MaskOnset::Other), None] {
            let read = run_read(onset);
            assert!(
                !ReadEligibility::new(&read, &settings(8)).run_exempt(13),
                "{onset:?}"
            );
        }
    }

    #[test]
    fn the_exemption_lifts_read_end_from_run_indels_only() {
        let read = run_read(Some(MaskOnset::AfterRun));
        let eligibility = ReadEligibility::new(&read, &settings(8));
        assert!(eligibility.exempt(VariantKind::Ins, "C", Orientation::Forward, &[15], &[15]));
        assert!(!eligibility.exempt(VariantKind::Snv, "C", Orientation::Forward, &[15], &[15]));
        // A mapped call beyond the run and its bounds.
        assert!(!eligibility.exempt(
            VariantKind::Ins,
            "C",
            Orientation::Forward,
            &[15],
            &[15, 20]
        ));
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
                onset: MaskOnset::Other,
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
