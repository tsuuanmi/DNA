//! Signal-derived per-read callability: phase-state segments, a typed mask, and
//! the callable span.
//!
//! The callability core works on plain per-position records and never on the
//! chromatogram itself. A segment states how the measured signal behaves over a
//! run of call positions; the mask marks every position of a non-in-phase
//! segment; the callable span is the hull of the unmasked positions. These are
//! derived observations over immutable channels, loci, calls, and locus
//! evidence: a masked position keeps its original call and evidence.

use serde::Serialize;

/// Multi-channel evidence at one called position, in call order.
///
/// Channel order follows `Nucleotide::ALL`: A, C, G, T. `amplitudes` are the
/// non-negative corrected channel amplitudes at the position's event;
/// `coordinate` is the position's sample index, used only for spacing;
/// `primary` is the channel index of the primary call, `None` when the call is
/// unresolved.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct PositionEvidence {
    pub(crate) amplitudes: [f64; 4],
    pub(crate) coordinate: usize,
    pub(crate) primary: Option<usize>,
}

/// Per-position features derived from [`PositionEvidence`].
///
/// Every value is rounded to six decimals before it is stored or compared.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct PositionFeatures {
    /// Channel with the highest amplitude (ties favour the lower index), or
    /// `None` when no channel carries positive signal.
    pub(crate) primary_channel: Option<usize>,
    /// Channel with the second-highest amplitude, or `None` without a primary.
    pub(crate) secondary_channel: Option<usize>,
    /// Primary minus secondary normalized amplitude, in `[0, 1]`.
    pub(crate) dominance: f64,
    /// Secondary amplitude divided by primary amplitude, in `[0, 1]`.
    pub(crate) secondary_ratio: f64,
    /// Smallest call offset whose primary channel equals this position's
    /// secondary channel (the slippage "shadow"), searched in the fixed order
    /// `+1, -1, +2, -2, +3, -3`.
    pub(crate) shift_offset: Option<i8>,
    /// Larger relative deviation of the two neighbouring spacings from the
    /// local median spacing.
    pub(crate) spacing_deviation: f64,
    /// No positive signal, or a primary amplitude below the configured fraction
    /// of the read's median primary amplitude.
    pub(crate) weak: bool,
}

/// How the measured signal behaves over one run of call positions.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub(crate) enum PhaseState {
    /// One dominant ladder; the positions are callable.
    InPhase,
    /// Double peaks whose secondary channel follows a neighbouring primary
    /// channel: superimposed ladders offset by a few positions (slippage).
    Dephased,
    /// Double peaks without a consistent neighbour offset: mixed signal.
    Mixed,
    /// Absent or low primary signal.
    Weak,
    /// Irregular call spacing.
    Irregular,
}

impl PhaseState {
    /// Stable operational-log label.
    pub(crate) const fn label(self) -> &'static str {
        match self {
            Self::InPhase => "in_phase",
            Self::Dephased => "dephased",
            Self::Mixed => "mixed",
            Self::Weak => "weak",
            Self::Irregular => "irregular",
        }
    }
}

/// One maximal run of call positions sharing a phase state.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct PhaseSegment {
    pub(crate) call_start_0based: usize,
    pub(crate) call_end_0based_exclusive: usize,
    pub(crate) state: PhaseState,
    /// The segment starts in the window that follows a long repeat run.
    pub(crate) after_repeat: bool,
    /// Fraction of the segment's double peaks sharing the modal shift offset.
    pub(crate) coherence: Option<f64>,
    /// Most frequent shift offset among the segment's double peaks.
    pub(crate) modal_offset: Option<i8>,
}

/// Unit of a repeat run found in the read's own primary calls.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum RepeatUnit {
    /// One channel repeated: a homopolymer.
    Homopolymer(usize),
    /// Two alternating channels: a dinucleotide tandem repeat.
    Dinucleotide(usize, usize),
}

/// A run of at least the configured minimum length of one repeat unit.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct RepeatRun {
    pub(crate) call_start_0based: usize,
    pub(crate) call_end_0based_exclusive: usize,
    pub(crate) unit: RepeatUnit,
}

/// Complete signal-derived callability of one read.
///
/// `segments` partition `[0, call_count)` in order; `mask` has one entry per
/// call and is `None` exactly on in-phase positions; the callable span is the
/// interval from the first to the last unmasked position and is empty when
/// every position is masked.
#[derive(Debug, Clone)]
pub(crate) struct ReadCallability {
    pub(crate) repeats: Vec<RepeatRun>,
    pub(crate) segments: Vec<PhaseSegment>,
    pub(crate) mask: Vec<Option<PhaseState>>,
    pub(crate) callable_start_0based: usize,
    pub(crate) callable_end_0based_exclusive: usize,
}

impl ReadCallability {
    /// Number of unmasked calls.
    pub(crate) fn callable_count(&self) -> usize {
        self.mask.iter().filter(|state| state.is_none()).count()
    }

    /// Number of masked calls.
    pub(crate) fn masked_count(&self) -> usize {
        self.mask.len() - self.callable_count()
    }

    /// Number of segments in `state`.
    pub(crate) fn segment_count(&self, state: PhaseState) -> usize {
        self.segments
            .iter()
            .filter(|segment| segment.state == state)
            .count()
    }

    /// A read whose every call is in phase, for tests of consuming stages.
    #[cfg(test)]
    pub(crate) fn in_phase(positions: usize) -> Self {
        Self {
            repeats: Vec::new(),
            segments: vec![PhaseSegment {
                call_start_0based: 0,
                call_end_0based_exclusive: positions,
                state: PhaseState::InPhase,
                after_repeat: false,
                coherence: None,
                modal_offset: None,
            }],
            mask: vec![None; positions],
            callable_start_0based: 0,
            callable_end_0based_exclusive: positions,
        }
    }
}
