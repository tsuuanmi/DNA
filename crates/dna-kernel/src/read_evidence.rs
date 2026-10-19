//! Modality-neutral per-read evidence: the contract between a sequencing
//! modality and the core caller.
//!
//! A modality (Sanger today) turns its raw data into one [`ReadEvidence`] per
//! read: the called bases, an optional basecall-independent A/C/G/T profile per
//! call, how each call takes part in alignment, the modality's verdicts on
//! what each call shows (its base, an insertion or deletion, a run's end), and
//! the interval of calls that carries information. The core aligns and calls
//! variants from this record alone and never interprets a modality's reason
//! labels; it only reports them, in the read's reason order.

use std::ops::Range;

use crate::error::{EvidenceError, Result};

/// Normalized non-negative A/C/G/T evidence at one call.
///
/// Channel order is A, C, G, T. A profile exists only when the call has
/// positive signal.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EvidenceProfile {
    /// Normalized A/C/G/T weights summing to one.
    pub weights: [f64; 4],
}

impl EvidenceProfile {
    /// The normalized profile of corrected A/C/G/T amplitudes, or `None` when no channel has positive signal.
    #[must_use]
    pub fn from_corrected_amplitudes(amplitudes: [f64; 4]) -> Option<Self> {
        let total = amplitudes.iter().sum::<f64>();
        (total > 0.0).then(|| Self {
            weights: amplitudes.map(|amplitude| amplitude / total),
        })
    }

    /// Complements A/C/G/T evidence while preserving its total mass.
    #[must_use]
    pub const fn complemented(self) -> Self {
        Self {
            weights: [
                self.weights[3],
                self.weights[2],
                self.weights[1],
                self.weights[0],
            ],
        }
    }
}

/// A modality's stable `snake_case` reason label, reported verbatim.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct EvidenceReason(&'static str);

impl EvidenceReason {
    /// Wraps a label; [`ReadEvidence::new`] validates it.
    #[must_use]
    pub const fn new(label: &'static str) -> Self {
        Self(label)
    }

    /// The published label.
    #[must_use]
    pub const fn label(self) -> &'static str {
        self.0
    }
}

/// How the core aligns a call, and why a masked call is masked.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Placement {
    /// The call aligns with its base and profile.
    Scored,
    /// Masked, but anchoring: the call still places the read with its base
    /// and profile, and shows no base. The modality's reason, verbatim.
    Anchoring(EvidenceReason),
    /// Masked and unresolved: the call aligns as `N` without a profile. The
    /// modality's reason, verbatim.
    Unresolved(EvidenceReason),
}

impl Placement {
    /// The modality's reason for a masked call; `None` for a scored one.
    #[must_use]
    pub const fn reason(self) -> Option<EvidenceReason> {
        match self {
            Self::Scored => None,
            Self::Anchoring(reason) | Self::Unresolved(reason) => Some(reason),
        }
    }
}

/// Largest number of reasons one verdict holds.
const REASON_CAPACITY: usize = 4;

/// The modality's reasons against one use of a call, each at most once.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct Reasons([Option<EvidenceReason>; REASON_CAPACITY]);

impl Reasons {
    /// Adds `reason`.
    ///
    /// # Errors
    ///
    /// Returns `EvidenceError::TooManyReasons` when the verdict already holds
    /// as many reasons as it can: a reason is never dropped silently, since a
    /// dropped reason would raise the call's trust.
    pub fn insert(&mut self, reason: EvidenceReason) -> Result<()> {
        if self.contains(reason) {
            return Ok(());
        }
        let slot =
            self.0
                .iter_mut()
                .find(|slot| slot.is_none())
                .ok_or(EvidenceError::TooManyReasons {
                    maximum: REASON_CAPACITY,
                })?;
        *slot = Some(reason);
        Ok(())
    }

    /// Whether `reason` is held.
    #[must_use]
    pub fn contains(&self, reason: EvidenceReason) -> bool {
        self.0.contains(&Some(reason))
    }

    /// Whether no reason is held.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.0.iter().all(Option::is_none)
    }

    /// The held reasons, in insertion order.
    pub fn iter(&self) -> impl Iterator<Item = EvidenceReason> + '_ {
        self.0.iter().filter_map(|slot| *slot)
    }
}

/// The modality's verdict on one use of a call.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct Verdict {
    /// The call lies at an edge of the read's usable signal: the core reports
    /// what rests on it as `read_end`.
    pub edge: bool,
    /// The modality's own reasons against this use, reported verbatim.
    pub withheld: Reasons,
}

impl Verdict {
    /// Neither at an edge nor withheld.
    #[must_use]
    pub fn trusted(&self) -> bool {
        !self.edge && self.withheld.is_empty()
    }
}

/// How a call shows where a run of another base ends (ADR-0071, ADR-0074).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Bound {
    /// In phase: the call shows the end, and counts towards its own run.
    InPhase,
    /// Anchoring: the call shows at most the dominant end of a mixture.
    Anchored,
    /// The call shows no run's end and counts towards no run's length.
    Withheld,
}

/// A call's part in a bounded run before a post-repeat phase loss: a run of
/// one base the read bounds right where it loses phase after a repeat run
/// (ADR-0076), as the modality finds it on the read's own calls. One call can
/// be both: a call of one such run that bounds the next.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct PostRepeatRun {
    /// The call belongs to such a run.
    pub call: bool,
    /// The call bounds such a run.
    pub bound: bool,
}

/// The second base of a call's signal, as the modality sees it.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SecondBase {
    /// The second base on the read strand: A, C, G, or T.
    pub base: char,
    /// Its share of the call's two strongest bases, in `[0, 0.5]`.
    pub share: f64,
    /// The modality asserts that the share stands out: it exceeds the
    /// modality's own floor for a second base.
    pub stands_out: bool,
}

/// Evidence for one called position, in call order.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CallEvidence {
    /// The called base on the read strand: A, C, G, T, or N when unresolved.
    pub base: char,
    /// Basecall-independent A/C/G/T profile, when the modality has one.
    pub profile: Option<EvidenceProfile>,
    /// How the call takes part in alignment, and why it is masked.
    pub placement: Placement,
    /// Whether the call shows its base: a substitution it supports, and its
    /// observation of a position.
    pub identity: Verdict,
    /// Whether the call places an insertion or deletion: an inserted call, or
    /// a call flanking an insertion or deletion.
    pub indel: Verdict,
    /// Whether the call shows where a run ends, and counts towards a run.
    pub bound: Bound,
    /// The call's part in a run the read bounds before a post-repeat phase
    /// loss (ADR-0076).
    pub post_repeat: PostRepeatRun,
    /// The second base of the call's signal, when the modality sees one.
    /// Evidence for review; the core reads only its base and whether it
    /// stands out, to keep a reference-only call from hiding a second
    /// sequence (ADR-0080).
    pub second: Option<SecondBase>,
    /// The modality asserts that the call shows a second peak standing out
    /// that it does not explain by the signal's own shape: possibly a second
    /// sequence (ADR-0080).
    pub second_unexplained: bool,
    /// The modality withholds the call's base, yet asserts that its signal
    /// shows that base where it is the reference's: the call may confirm the
    /// reference at its own position, and nothing else (ADR-0080).
    pub reference_only: bool,
}

impl CallEvidence {
    /// A scored call with `base` and `profile` that the modality trusts for
    /// everything it shows, without a second base.
    #[must_use]
    pub const fn scored(base: char, profile: Option<EvidenceProfile>) -> Self {
        Self {
            base,
            profile,
            placement: Placement::Scored,
            identity: Verdict {
                edge: false,
                withheld: Reasons([None; REASON_CAPACITY]),
            },
            indel: Verdict {
                edge: false,
                withheld: Reasons([None; REASON_CAPACITY]),
            },
            bound: Bound::InPhase,
            post_repeat: PostRepeatRun {
                call: false,
                bound: false,
            },
            second: None,
            second_unexplained: false,
            reference_only: false,
        }
    }

    /// The same call masked: anchoring calls still show where a run ends,
    /// unresolved ones show nothing; neither shows its base or an indel,
    /// whatever its verdicts hold, and the placement carries the reason.
    #[must_use]
    pub fn masked(mut self, placement: Placement) -> Self {
        self.bound = match placement {
            Placement::Scored => self.bound,
            Placement::Anchoring(_) => Bound::Anchored,
            Placement::Unresolved(_) => Bound::Withheld,
        };
        self.placement = placement;
        self
    }

    /// Whether the call is vouched for by its base only: scored, but showing
    /// no run's end and placing no indel.
    #[must_use]
    pub fn base_only(&self) -> bool {
        self.placement == Placement::Scored && self.bound == Bound::Withheld
    }
}

/// Labels the core owns; a modality may not reuse them.
const CORE_REASONS: [&str; 8] = [
    "non_canonical_allele",
    "indel_length_exceeded",
    "outside_target_region",
    "read_end",
    "base_only",
    "run_boundary",
    "unresolved_flank",
    "placement_conflict",
];

/// Whether the modality vouches for a read's calls up to its physical ends.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReadEnds {
    /// The read's ends bound its evidence: like an uninformative call, an end
    /// keeps the calls within `read_end_margin` of it from supporting a variant.
    Unvouched,
    /// The modality vouches for every call up to the read's ends, as for a
    /// reviewed consensus sequence; only uninformative calls bound the evidence.
    Vouched,
}

/// One read's evidence for the core.
#[derive(Debug, Clone)]
pub struct ReadEvidence {
    calls: Vec<CallEvidence>,
    informative: Range<usize>,
    reason_order: Vec<EvidenceReason>,
    ends: ReadEnds,
    mixtures: Vec<PhaseMixture>,
}

/// One sequence a phase mixture holds: the read's dominant sequence offset by
/// `offset` calls, with its share of the signal (ADR-0082).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MixtureComponent {
    /// Calls by which this sequence lags the dominant one: a positive offset
    /// is a longer run before the mixture, a negative one a shorter run.
    pub offset: i32,
    /// Its share of the fitted signal, in `[0, 1]`; the shares sum to one.
    pub fraction: f64,
}

/// The modality's assertion that, over calls `call_start..call_end`, the
/// read's signal is its own dominant sequence (its calls there) superposed
/// with copies of itself offset by whole calls, as a length mixture of the
/// run that opens the stretch makes (ADR-0082). It decides nothing: the
/// calls keep their verdicts.
#[derive(Debug, Clone, PartialEq)]
pub struct PhaseMixture {
    /// First call of the stretch: the start of the run whose length varies
    /// (0-based, inclusive).
    pub call_start: usize,
    /// Call after the stretch (0-based, exclusive).
    pub call_end: usize,
    /// Calls of the run whose length varies, from `call_start`.
    pub run_calls: usize,
    /// The sequences the stretch holds, offset 0 (the dominant one) first,
    /// the rest by decreasing fraction.
    pub components: Vec<MixtureComponent>,
    /// What the components leave of the signal, as the root mean square per
    /// call on the modality's normalised scale.
    pub residual: f64,
    /// Calls of the stretch whose signal the components do not explain.
    pub unexplained: Vec<usize>,
    /// Whether the components explain the stretch: its residual is at most
    /// the residual the read's own in-phase calls reach at the modality's
    /// unexplained quantile. The components of an unexplained mixture are
    /// the fit's best, not sequences the read is shown to hold.
    pub explained: bool,
}

impl ReadEvidence {
    /// Validates and wraps one read's evidence.
    ///
    /// `informative` is the interval of calls the modality retains.
    /// `reason_order` lists the modality's verdict reasons in the order the
    /// core reports them; every reason a verdict holds must be listed there.
    /// Reason labels, those of masked placements included, must be non-empty
    /// `snake_case` and not core labels.
    ///
    /// # Errors
    ///
    /// Returns `EvidenceError` when the interval exceeds the calls, a verdict
    /// holds a reason the order does not list, or a reason label is invalid,
    /// repeated in the order, or a core label.
    pub fn new(
        calls: Vec<CallEvidence>,
        informative: Range<usize>,
        reason_order: Vec<EvidenceReason>,
    ) -> Result<Self> {
        if informative.start > informative.end || informative.end > calls.len() {
            return Err(EvidenceError::InvalidInformativeInterval {
                start: informative.start,
                end: informative.end,
                calls: calls.len(),
            }
            .into());
        }
        let mut seen = Vec::with_capacity(reason_order.len());
        for reason in &reason_order {
            let label = reason.label();
            if !valid_label(label) || seen.contains(&label) {
                return Err(EvidenceError::InvalidReason(label).into());
            }
            seen.push(label);
        }
        for (index, call) in calls.iter().enumerate() {
            if let Some(reason) = call.placement.reason()
                && !valid_label(reason.label())
            {
                return Err(EvidenceError::InvalidReason(reason.label()).into());
            }
            if call
                .identity
                .withheld
                .iter()
                .chain(call.indel.withheld.iter())
                .any(|reason| !reason_order.contains(&reason))
            {
                return Err(EvidenceError::UndeclaredReason { index }.into());
            }
        }
        Ok(Self {
            calls,
            informative,
            reason_order,
            ends: ReadEnds::Unvouched,
            mixtures: Vec::new(),
        })
    }

    /// The same evidence with the modality's phase mixtures (ADR-0082).
    ///
    /// # Errors
    ///
    /// Returns `EvidenceError::InvalidMixture` when a mixture's stretch does
    /// not fit the read's calls or its run, it has no offset-0 component
    /// first, a fraction is not finite in `[0, 1]`, an offset is longer than
    /// its run, or an unexplained call lies outside its stretch.
    pub fn with_mixtures(mut self, mixtures: Vec<PhaseMixture>) -> Result<Self> {
        for (index, mixture) in mixtures.iter().enumerate() {
            let fits = mixture.call_start < mixture.call_end
                && mixture.call_end <= self.calls.len()
                && mixture.run_calls > 0
                && mixture.call_start + mixture.run_calls <= mixture.call_end
                && mixture
                    .components
                    .first()
                    .is_some_and(|first| first.offset == 0)
                && mixture.residual.is_finite()
                && mixture.residual >= 0.0
                && mixture.components.iter().all(|component| {
                    component.fraction.is_finite()
                        && (0.0..=1.0).contains(&component.fraction)
                        && component.offset.unsigned_abs() as usize <= mixture.run_calls
                })
                && mixture
                    .unexplained
                    .iter()
                    .all(|call| (mixture.call_start..mixture.call_end).contains(call));
            if !fits {
                return Err(EvidenceError::InvalidMixture { index }.into());
            }
        }
        self.mixtures = mixtures;
        Ok(self)
    }

    /// The modality's phase mixtures, in call order (ADR-0082).
    #[must_use]
    pub fn mixtures(&self) -> &[PhaseMixture] {
        &self.mixtures
    }

    /// The same evidence with every call within `margin` calls of an
    /// uninformative call placed at a signal edge, for its base and for
    /// indels. An uninformative call lies outside the informative interval or
    /// is masked as unresolved; the read's ends count as uninformative unless
    /// the modality vouches for them. A call is clear of the edges when more
    /// than `margin` informative calls, itself included, separate it from the
    /// nearest uninformative call on each side.
    #[must_use]
    pub fn with_signal_edges(mut self, margin: usize) -> Self {
        let informative = (0..self.calls.len())
            .map(|index| self.informative_call(index))
            .collect::<Vec<_>>();
        let end = match self.ends {
            ReadEnds::Unvouched => 0,
            ReadEnds::Vouched => margin,
        };
        let mut before = vec![0_usize; informative.len()];
        let mut run = end;
        for (index, &flag) in informative.iter().enumerate() {
            run = if flag { run + 1 } else { 0 };
            before[index] = run;
        }
        run = end;
        for index in (0..informative.len()).rev() {
            run = if informative[index] { run + 1 } else { 0 };
            let edge = before[index] <= margin || run <= margin;
            self.calls[index].identity.edge = edge;
            self.calls[index].indel.edge = edge;
        }
        self
    }

    /// Whether the call at `index` carries information: inside the
    /// informative interval and not masked as unresolved.
    #[must_use]
    pub fn informative_call(&self, index: usize) -> bool {
        self.informative.contains(&index)
            && self
                .calls
                .get(index)
                .is_some_and(|call| !matches!(call.placement, Placement::Unresolved(_)))
    }

    /// The same evidence with the modality vouching for the read's ends.
    #[must_use]
    pub fn with_vouched_ends(mut self) -> Self {
        self.ends = ReadEnds::Vouched;
        self
    }

    /// Whether the modality vouches for the read's calls up to its ends.
    #[must_use]
    pub fn ends(&self) -> ReadEnds {
        self.ends
    }

    /// Every call, in call order.
    #[must_use]
    pub fn calls(&self) -> &[CallEvidence] {
        &self.calls
    }

    /// The interval of calls the modality retains.
    #[must_use]
    pub fn informative(&self) -> Range<usize> {
        self.informative.clone()
    }

    /// The modality's verdict reasons, in reporting order.
    #[must_use]
    pub fn reason_order(&self) -> &[EvidenceReason] {
        &self.reason_order
    }

    /// Clean evidence for `sequence`: one-hot profiles, no masks, no vetoes,
    /// every call informative.
    #[cfg(any(test, feature = "test-support"))]
    #[must_use]
    pub fn clean(sequence: &str) -> Self {
        let calls = sequence
            .chars()
            .map(|base| {
                let mut weights = [0.0; 4];
                let profile = "ACGT".find(base).map(|channel| {
                    weights[channel] = 1.0;
                    EvidenceProfile { weights }
                });
                CallEvidence::scored(base, profile)
            })
            .collect::<Vec<_>>();
        let informative = 0..calls.len();
        Self {
            calls,
            informative,
            reason_order: Vec::new(),
            ends: ReadEnds::Unvouched,
            mixtures: Vec::new(),
        }
    }

    /// Replaces one call, for tests of consuming stages.
    #[cfg(any(test, feature = "test-support"))]
    #[must_use]
    pub fn with_call(mut self, index: usize, call: CallEvidence) -> Self {
        self.calls[index] = call;
        self
    }
}

fn valid_label(label: &str) -> bool {
    !label.is_empty()
        && label
            .bytes()
            .all(|byte| byte.is_ascii_lowercase() || byte.is_ascii_digit() || byte == b'_')
        && !CORE_REASONS.contains(&label)
}

#[cfg(test)]
mod tests {
    use super::*;

    const PEAK: EvidenceReason = EvidenceReason::new("peak_below_minimum");

    #[test]
    fn rejects_inconsistent_evidence() {
        let call = CallEvidence::scored('A', None);
        assert!(ReadEvidence::new(vec![call; 2], 1..3, Vec::new()).is_err());
        let mut vetoed = call;
        assert!(vetoed.identity.withheld.insert(PEAK).is_ok());
        // A verdict reason the order does not list.
        assert!(ReadEvidence::new(vec![vetoed], 0..1, Vec::new()).is_err());
        assert!(ReadEvidence::new(vec![call], 0..1, vec![PEAK, PEAK]).is_err());
        // Every label the core owns is refused to a modality.
        for label in CORE_REASONS {
            assert!(
                ReadEvidence::new(vec![call], 0..1, vec![EvidenceReason::new(label)]).is_err(),
                "{label}"
            );
        }
        let masked = call.masked(Placement::Unresolved(EvidenceReason::new("Weak Signal")));
        assert!(ReadEvidence::new(vec![masked], 0..1, Vec::new()).is_err());
        assert!(ReadEvidence::new(vec![vetoed, call], 0..2, vec![PEAK]).is_ok());
    }

    #[test]
    fn a_masked_call_carries_its_reason_and_keeps_an_anchoring_bound() {
        let reason = EvidenceReason::new("dephased_signal");
        let anchoring = CallEvidence::scored('A', None).masked(Placement::Anchoring(reason));
        assert_eq!(anchoring.placement.reason(), Some(reason));
        assert_eq!(anchoring.bound, Bound::Anchored);
        let unresolved = CallEvidence::scored('A', None).masked(Placement::Unresolved(reason));
        assert_eq!(unresolved.bound, Bound::Withheld);
        assert!(!unresolved.base_only());
    }

    #[test]
    fn a_verdict_refuses_a_reason_it_cannot_hold() {
        let mut reasons = Reasons::default();
        for label in ["a", "b", "c", "d"] {
            assert!(reasons.insert(EvidenceReason::new(label)).is_ok());
        }
        // A reason it holds already is no new one.
        assert!(reasons.insert(EvidenceReason::new("a")).is_ok());
        assert!(reasons.insert(EvidenceReason::new("e")).is_err());
    }

    #[test]
    fn marks_calls_within_the_margin_of_an_uninformative_call() -> Result<()> {
        let call = CallEvidence::scored('A', None);
        let masked = call.masked(Placement::Unresolved(EvidenceReason::new("weak_signal")));
        let mut calls = vec![call; 12];
        calls[7] = masked;
        let edges = |evidence: ReadEvidence| {
            evidence
                .with_signal_edges(2)
                .calls()
                .iter()
                .map(|call| match (call.identity.edge, call.indel.edge) {
                    (true, true) => 'e',
                    (false, false) => '.',
                    _ => '?',
                })
                .collect::<String>()
        };
        // Informative calls 1..11; the masked call at 7 and both ends of the
        // interval are edges.
        let evidence = ReadEvidence::new(calls.clone(), 1..11, Vec::new())?;
        assert_eq!(edges(evidence), "eee..eeeeeee");
        let vouched = ReadEvidence::new(calls, 0..12, Vec::new())?.with_vouched_ends();
        assert_eq!(edges(vouched), ".....eeeee..");
        Ok(())
    }

    #[test]
    fn builds_clean_evidence_with_one_hot_profiles() {
        let evidence = ReadEvidence::clean("ACN");
        assert_eq!(evidence.informative(), 0..3);
        assert_eq!(
            evidence.calls()[1].profile,
            Some(EvidenceProfile {
                weights: [0.0, 1.0, 0.0, 0.0]
            })
        );
        assert_eq!(evidence.calls()[2].profile, None);
    }
}
