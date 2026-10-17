//! Modality-neutral per-read evidence: the contract between a sequencing
//! modality and the core caller.
//!
//! A modality (Sanger today) turns its raw data into one [`ReadEvidence`] per
//! read: the called bases, an optional basecall-independent A/C/G/T profile per
//! call, which calls are masked and why, the support vetoes its own evidence
//! raises against each call, and the interval of calls that carries
//! information. The core aligns and calls variants from this record alone and
//! never interprets a modality's reason labels; it only reports them.

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

/// How the core aligns a masked call.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MaskedAlignment {
    /// The call carries no alignment information: it aligns as `N` without a
    /// profile.
    Unresolved,
    /// The call still anchors the alignment with its base and profile, but
    /// supports no variant.
    Anchoring,
}

/// Where a masked call's segment begins, as the modality places it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MaskOnset {
    /// The segment starts in the window right after a repeat run of the
    /// read's own calls: the read loses phase leaving the run (ADR-0076).
    AfterRun,
    /// Any other start, or one the modality does not place.
    Other,
}

/// A masked call: how it aligns and why it is masked.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CallMask {
    /// How the masked call takes part in alignment.
    pub alignment: MaskedAlignment,
    /// The modality's reason label, reported verbatim.
    pub reason: EvidenceReason,
    /// Where the call's masked segment begins.
    pub onset: MaskOnset,
}

/// Which variant kinds a support veto applies to.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VetoScope {
    /// Substitutions only.
    Substitutions,
    /// Substitutions and the inserted bases of insertions.
    SubstitutionsAndInsertions,
}

/// One entry of a read's support-veto vocabulary.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SupportVeto {
    /// The modality's reason label, reported verbatim.
    pub reason: EvidenceReason,
    /// The variant kinds the veto applies to.
    pub scope: VetoScope,
}

/// The support vetoes raised against one call, as bits into its read's
/// vocabulary.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct VetoSet(u32);

impl VetoSet {
    /// Largest supported vocabulary.
    pub(crate) const CAPACITY: usize = 32;

    /// Adds vocabulary entry `index`.
    pub fn insert(&mut self, index: usize) {
        if let Some(bit) = u32::try_from(index)
            .ok()
            .and_then(|index| 1_u32.checked_shl(index))
        {
            self.0 |= bit;
        }
    }

    /// Whether vocabulary entry `index` is set.
    #[must_use]
    pub fn contains(self, index: usize) -> bool {
        u32::try_from(index)
            .ok()
            .and_then(|index| 1_u32.checked_shl(index))
            .is_some_and(|bit| self.0 & bit != 0)
    }

    /// Union with another set.
    #[must_use]
    pub const fn union(self, other: Self) -> Self {
        Self(self.0 | other.0)
    }

    fn fits(self, entries: usize) -> bool {
        entries >= Self::CAPACITY || self.0 >> entries == 0
    }
}

/// The second base of a call's signal, as the modality sees it.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SecondBase {
    /// The second base on the read strand: A, C, G, or T.
    pub base: char,
    /// Its share of the call's two strongest bases, in `[0, 0.5]`.
    pub share: f64,
}

/// Evidence for one called position, in call order.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CallEvidence {
    /// The called base on the read strand: A, C, G, T, or N when unresolved.
    pub base: char,
    /// Basecall-independent A/C/G/T profile, when the modality has one.
    pub profile: Option<EvidenceProfile>,
    /// Why and how the call is masked; `None` for a callable call.
    pub mask: Option<CallMask>,
    /// Support vetoes the modality raises against the call.
    pub vetoes: VetoSet,
    /// The modality withholds the call as too close to an edge of its usable
    /// signal: the core reports such a call's variants as `read_end`.
    pub edge: bool,
    /// The modality vouches for the call's base only: it neither counts
    /// towards a run's length nor shows where a run ends, and it supports no
    /// insertion and bounds no deletion.
    pub identity_only: bool,
    /// The second base of the call's signal, when the modality sees one;
    /// evidence for review, read by no calling rule.
    pub second: Option<SecondBase>,
}

impl CallEvidence {
    /// A call with `base`, `profile`, `mask`, and `vetoes`, outside every
    /// signal edge, trusted for runs, and without a second base.
    #[must_use]
    pub const fn new(
        base: char,
        profile: Option<EvidenceProfile>,
        mask: Option<CallMask>,
        vetoes: VetoSet,
    ) -> Self {
        Self {
            base,
            profile,
            mask,
            vetoes,
            edge: false,
            identity_only: false,
            second: None,
        }
    }
}

/// Labels the core owns; a modality may not reuse them.
const CORE_REASONS: [&str; 7] = [
    "non_canonical_allele",
    "indel_length_exceeded",
    "outside_target_region",
    "read_end",
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
    support_vetoes: Vec<SupportVeto>,
    ends: ReadEnds,
}

impl ReadEvidence {
    /// Validates and wraps one read's evidence.
    ///
    /// `informative` is the interval of calls the modality retains; every
    /// veto bit must name an entry of `support_vetoes`; reason labels must be
    /// non-empty `snake_case`, distinct across vetoes, and not core labels.
    ///
    /// # Errors
    ///
    /// Returns `EvidenceError` when the interval exceeds the calls, there are
    /// more than `VetoSet::CAPACITY` vetoes, a veto bit names no vocabulary
    /// entry, or a reason label is invalid or a core label.
    pub fn new(
        calls: Vec<CallEvidence>,
        informative: Range<usize>,
        support_vetoes: Vec<SupportVeto>,
    ) -> Result<Self> {
        if informative.start > informative.end || informative.end > calls.len() {
            return Err(EvidenceError::InvalidInformativeInterval {
                start: informative.start,
                end: informative.end,
                calls: calls.len(),
            }
            .into());
        }
        if support_vetoes.len() > VetoSet::CAPACITY {
            return Err(EvidenceError::TooManyVetoes {
                count: support_vetoes.len(),
                maximum: VetoSet::CAPACITY,
            }
            .into());
        }
        let mut seen = Vec::with_capacity(support_vetoes.len());
        for veto in &support_vetoes {
            let label = veto.reason.label();
            if !valid_label(label) || seen.contains(&label) {
                return Err(EvidenceError::InvalidReason(label).into());
            }
            seen.push(label);
        }
        for (index, call) in calls.iter().enumerate() {
            if !call.vetoes.fits(support_vetoes.len()) {
                return Err(EvidenceError::UndeclaredVeto { index }.into());
            }
            if let Some(mask) = call.mask
                && !valid_label(mask.reason.label())
            {
                return Err(EvidenceError::InvalidReason(mask.reason.label()).into());
            }
        }
        Ok(Self {
            calls,
            informative,
            support_vetoes,
            ends: ReadEnds::Unvouched,
        })
    }

    /// The same evidence with every call within `margin` calls of an
    /// uninformative call marked as at a signal edge (`edge`). An
    /// uninformative call lies outside the informative interval or is masked
    /// without anchoring; the read's ends count as uninformative unless the
    /// modality vouches for them. A call is clear of the edges when more than
    /// `margin` informative calls, itself included, separate it from the
    /// nearest uninformative call on each side.
    #[must_use]
    pub fn with_signal_edges(mut self, margin: usize) -> Self {
        let informative = self
            .calls
            .iter()
            .enumerate()
            .map(|(index, call)| {
                self.informative.contains(&index)
                    && call
                        .mask
                        .is_none_or(|mask| mask.alignment == MaskedAlignment::Anchoring)
            })
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
            self.calls[index].edge = before[index] <= margin || run <= margin;
        }
        self
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

    /// The support-veto vocabulary, in reporting order.
    #[must_use]
    pub fn support_vetoes(&self) -> &[SupportVeto] {
        &self.support_vetoes
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
                CallEvidence::new(base, profile, None, VetoSet::default())
            })
            .collect::<Vec<_>>();
        let informative = 0..calls.len();
        Self {
            calls,
            informative,
            support_vetoes: Vec::new(),
            ends: ReadEnds::Unvouched,
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

    const PEAK: SupportVeto = SupportVeto {
        reason: EvidenceReason::new("peak_below_minimum"),
        scope: VetoScope::SubstitutionsAndInsertions,
    };

    #[test]
    fn rejects_inconsistent_evidence() {
        let call = CallEvidence::new('A', None, None, VetoSet::default());
        assert!(ReadEvidence::new(vec![call; 2], 1..3, Vec::new()).is_err());
        let mut vetoed = call;
        vetoed.vetoes.insert(1);
        assert!(ReadEvidence::new(vec![vetoed], 0..1, vec![PEAK]).is_err());
        assert!(ReadEvidence::new(vec![call], 0..1, vec![PEAK, PEAK]).is_err());
        // Every label the core owns is refused to a modality.
        for label in CORE_REASONS {
            let core = SupportVeto {
                reason: EvidenceReason::new(label),
                ..PEAK
            };
            assert!(
                ReadEvidence::new(vec![call], 0..1, vec![core]).is_err(),
                "{label}"
            );
        }
        let masked = CallEvidence {
            mask: Some(CallMask {
                alignment: MaskedAlignment::Unresolved,
                reason: EvidenceReason::new("Weak Signal"),
                onset: MaskOnset::Other,
            }),
            ..call
        };
        assert!(ReadEvidence::new(vec![masked], 0..1, Vec::new()).is_err());
        vetoed = call;
        vetoed.vetoes.insert(0);
        assert!(ReadEvidence::new(vec![vetoed, call], 0..2, vec![PEAK]).is_ok());
    }

    #[test]
    fn marks_calls_within_the_margin_of_an_uninformative_call() -> Result<()> {
        let call = CallEvidence::new('A', None, None, VetoSet::default());
        let masked = CallEvidence {
            mask: Some(CallMask {
                alignment: MaskedAlignment::Unresolved,
                reason: EvidenceReason::new("weak_signal"),
                onset: MaskOnset::Other,
            }),
            ..call
        };
        let mut calls = vec![call; 12];
        calls[7] = masked;
        let edges = |evidence: ReadEvidence| {
            evidence
                .with_signal_edges(2)
                .calls()
                .iter()
                .map(|call| if call.edge { 'e' } else { '.' })
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
