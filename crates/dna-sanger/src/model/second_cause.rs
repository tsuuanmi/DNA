//! Why a call shows a second peak (SRS-CALL-021): a label from the peak's
//! own shape, kept beside the core's records and never read by the core.

use dna_kernel::read_evidence::SecondBase;
use serde::Serialize;

/// The cause of a call's second peak, from its shape against the call's
/// neighbours, or of a call masked as mixed whose own second peak is small.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum SecondCause {
    /// The second peak barely stands above the call's lowest channel.
    Noise,
    /// The second channel's signal is the tail of a neighbouring call of that
    /// base: no local maximum in the call's window, or one shifted towards
    /// that neighbour.
    NeighbourTail,
    /// A second peak shifted off the call's centre, towards no neighbour of
    /// its base.
    OffCentre,
    /// A centred second peak much wider than the read's primary peaks there:
    /// a dye artefact.
    DyeBlob,
    /// A centred second peak of a neighbouring call's base: a shadow of
    /// slippage.
    PhaseShift,
    /// A centred second peak of another base, within the run level of its
    /// key in the run signal: a shadow the run's reads share.
    Context,
    /// A centred second peak of another base that nothing above explains.
    SecondTemplate,
    /// A call masked as mixed by its segment whose own second peak is small.
    MixedStretch,
}

impl SecondCause {
    /// Stable label.
    #[must_use]
    pub const fn label(self) -> &'static str {
        match self {
            Self::Noise => "noise",
            Self::NeighbourTail => "neighbour_tail",
            Self::OffCentre => "off_centre",
            Self::DyeBlob => "dye_blob",
            Self::PhaseShift => "phase_shift",
            Self::Context => "context",
            Self::SecondTemplate => "second_template",
            Self::MixedStretch => "mixed_stretch",
        }
    }
}

/// A call's second peak and its cause.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SecondPeak {
    /// The second base on the read strand (the strongest other channel), its
    /// share of the call's two strongest bases (SRS-CALL-015), and whether
    /// the share exceeds `signal_processing.second_floor`. A call masked as
    /// mixed without a second base has base `N` and share zero.
    pub second: SecondBase,
    /// Why the call shows it.
    pub cause: SecondCause,
}
