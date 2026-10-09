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
pub(crate) struct EvidenceProfile {
    pub(crate) weights: [f64; 4],
}

impl EvidenceProfile {
    pub(crate) fn from_corrected_amplitudes(amplitudes: [f64; 4]) -> Option<Self> {
        let total = amplitudes.iter().sum::<f64>();
        (total > 0.0).then(|| Self {
            weights: amplitudes.map(|amplitude| amplitude / total),
        })
    }

    /// Complements A/C/G/T evidence while preserving its total mass.
    pub(crate) const fn complemented(self) -> Self {
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
pub(crate) struct EvidenceReason(&'static str);

impl EvidenceReason {
    /// Wraps a label; [`ReadEvidence::new`] validates it.
    pub(crate) const fn new(label: &'static str) -> Self {
        Self(label)
    }

    /// The published label.
    pub(crate) const fn label(self) -> &'static str {
        self.0
    }
}

/// How the core aligns a masked call.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum MaskedAlignment {
    /// The call carries no alignment information: it aligns as `N` without a
    /// profile.
    Unresolved,
    /// The call still anchors the alignment with its base and profile, but
    /// supports no variant.
    Anchoring,
}

/// A masked call: how it aligns and why it is masked.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct CallMask {
    pub(crate) alignment: MaskedAlignment,
    pub(crate) reason: EvidenceReason,
}

/// Which variant kinds a support veto applies to.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum VetoScope {
    /// Substitutions only.
    Substitutions,
    /// Substitutions and the inserted bases of insertions.
    SubstitutionsAndInsertions,
}

/// One entry of a read's support-veto vocabulary.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct SupportVeto {
    pub(crate) reason: EvidenceReason,
    pub(crate) scope: VetoScope,
}

/// The support vetoes raised against one call, as bits into its read's
/// vocabulary.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub(crate) struct VetoSet(u32);

impl VetoSet {
    /// Largest supported vocabulary.
    pub(crate) const CAPACITY: usize = 32;

    /// Adds vocabulary entry `index`.
    pub(crate) fn insert(&mut self, index: usize) {
        if let Some(bit) = u32::try_from(index)
            .ok()
            .and_then(|index| 1_u32.checked_shl(index))
        {
            self.0 |= bit;
        }
    }

    /// Whether vocabulary entry `index` is set.
    pub(crate) fn contains(self, index: usize) -> bool {
        u32::try_from(index)
            .ok()
            .and_then(|index| 1_u32.checked_shl(index))
            .is_some_and(|bit| self.0 & bit != 0)
    }

    /// Union with another set.
    #[must_use]
    pub(crate) const fn union(self, other: Self) -> Self {
        Self(self.0 | other.0)
    }

    fn fits(self, entries: usize) -> bool {
        entries >= Self::CAPACITY || self.0 >> entries == 0
    }
}

/// Evidence for one called position, in call order.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct CallEvidence {
    /// The called base on the read strand: A, C, G, T, or N when unresolved.
    pub(crate) base: char,
    pub(crate) profile: Option<EvidenceProfile>,
    pub(crate) mask: Option<CallMask>,
    pub(crate) vetoes: VetoSet,
}

/// Labels the core owns; a modality may not reuse them.
const CORE_REASONS: [&str; 4] = [
    "non_canonical_allele",
    "indel_length_exceeded",
    "outside_target_region",
    "read_end",
];

/// Whether the modality vouches for a read's calls up to its physical ends.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum ReadEnds {
    /// The read's ends bound its evidence: like an uninformative call, an end
    /// keeps the calls within `read_end_margin` of it from supporting a variant.
    Unvouched,
    /// The modality vouches for every call up to the read's ends, as for a
    /// reviewed consensus sequence; only uninformative calls bound the evidence.
    Vouched,
}

/// One read's evidence for the core.
#[derive(Debug, Clone)]
pub(crate) struct ReadEvidence {
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
    pub(crate) fn new(
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

    /// The same evidence with the modality vouching for the read's ends.
    pub(crate) fn with_vouched_ends(mut self) -> Self {
        self.ends = ReadEnds::Vouched;
        self
    }

    /// Whether the modality vouches for the read's calls up to its ends.
    pub(crate) fn ends(&self) -> ReadEnds {
        self.ends
    }

    /// Every call, in call order.
    pub(crate) fn calls(&self) -> &[CallEvidence] {
        &self.calls
    }

    /// The interval of calls the modality retains.
    pub(crate) fn informative(&self) -> Range<usize> {
        self.informative.clone()
    }

    /// The support-veto vocabulary, in reporting order.
    pub(crate) fn support_vetoes(&self) -> &[SupportVeto] {
        &self.support_vetoes
    }

    /// Clean evidence for `sequence`: one-hot profiles, no masks, no vetoes,
    /// every call informative.
    #[cfg(test)]
    pub(crate) fn clean(sequence: &str) -> Self {
        let calls = sequence
            .chars()
            .map(|base| {
                let mut weights = [0.0; 4];
                let profile = "ACGT".find(base).map(|channel| {
                    weights[channel] = 1.0;
                    EvidenceProfile { weights }
                });
                CallEvidence {
                    base,
                    profile,
                    mask: None,
                    vetoes: VetoSet::default(),
                }
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
    #[cfg(test)]
    pub(crate) fn with_call(mut self, index: usize, call: CallEvidence) -> Self {
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
        let call = CallEvidence {
            base: 'A',
            profile: None,
            mask: None,
            vetoes: VetoSet::default(),
        };
        assert!(ReadEvidence::new(vec![call; 2], 1..3, Vec::new()).is_err());
        let mut vetoed = call;
        vetoed.vetoes.insert(1);
        assert!(ReadEvidence::new(vec![vetoed], 0..1, vec![PEAK]).is_err());
        assert!(ReadEvidence::new(vec![call], 0..1, vec![PEAK, PEAK]).is_err());
        let core = SupportVeto {
            reason: EvidenceReason::new("read_end"),
            ..PEAK
        };
        assert!(ReadEvidence::new(vec![call], 0..1, vec![core]).is_err());
        let masked = CallEvidence {
            mask: Some(CallMask {
                alignment: MaskedAlignment::Unresolved,
                reason: EvidenceReason::new("Weak Signal"),
            }),
            ..call
        };
        assert!(ReadEvidence::new(vec![masked], 0..1, Vec::new()).is_err());
        vetoed = call;
        vetoed.vetoes.insert(0);
        assert!(ReadEvidence::new(vec![vetoed, call], 0..2, vec![PEAK]).is_ok());
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
