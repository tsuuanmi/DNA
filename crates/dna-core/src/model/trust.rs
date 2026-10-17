//! How far the core trusts each call of a read: the verdict its variant
//! calling and consensus already reach, exposed as diagnostic evidence for
//! evaluation and review. It is not a calling contract.

use dna_kernel::read_evidence::{CallMask, EvidenceReason};

/// What withholds trust from a call: the nearest uninformative call, or the
/// read's own end, within the read-end margin on one side.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Withholder {
    /// The read's own end, which its modality does not vouch for.
    ReadEnd,
    /// A call outside the read's informative interval.
    Untrimmed,
    /// A masked call that does not anchor the alignment, with its reason
    /// and onset as the modality gave them.
    Masked(CallMask),
}

/// The core's trust in one call of a read.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CallTrust {
    /// The call's mask, when it has one.
    pub mask: Option<CallMask>,
    /// The call reads no base (`N`).
    pub unresolved: bool,
    /// The support vetoes raised against the call, by label.
    pub vetoes: Vec<EvidenceReason>,
    /// The call is informative and outside the read-end margin of every
    /// uninformative call (`read_end`).
    pub trusted: bool,
    /// A bounded-run exemption lifts `read_end` from the call (ADR-0076).
    pub exempt: bool,
    /// What withholds trust, before then after the call in read order; empty
    /// when the call is trusted. An uninformative call is withheld by itself.
    pub withheld_by: Vec<Withholder>,
}
