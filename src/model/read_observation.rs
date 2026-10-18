//! Complete scientific observation produced from one independently processed Sanger read.

use dna_call::model::called_read::CalledRead;
use dna_call::read_call::PlacementRejection;
use dna_sanger::model::attachment::{SangerAttachment, SangerRejection};
use dna_sanger::model::callability::ReadRejection;
use dna_sanger::model::signal::ShortRead;

/// A Sanger read after reference placement: the core's modality-neutral
/// products and the Sanger evidence that reports join to them.
#[derive(Debug, Clone)]
pub(crate) struct ReadObservation {
    pub(crate) called: CalledRead,
    pub(crate) sanger: SangerAttachment,
}

/// A Sanger read a sample operation set aside.
///
/// Public through [`sample_analysis`](crate::sample_analysis).
#[derive(Debug, Clone)]
#[non_exhaustive]
pub struct RejectedRead {
    /// Reviewer-facing name of the rejected trace.
    pub input_name: String,
    /// SHA-256 content identity of the rejected trace.
    pub input_sha256: String,
    /// Why the read was set aside, with the Sanger evidence its stages
    /// produced.
    pub cause: RejectionCause,
}

/// Why a read was set aside.
#[derive(Debug, Clone)]
#[non_exhaustive]
pub enum RejectionCause {
    /// Fewer calls than one signal window (SRS-SAMPLE-030); no stage ran, so
    /// there is no Sanger evidence.
    TooShort(ShortRead),
    /// Too few callable calls (SRS-SAMPLE-027).
    Callability {
        /// Integrity and callability of the read.
        sanger: SangerRejection,
        /// Callable and minimum call counts.
        rejection: ReadRejection,
    },
    /// The core could not place the read (SRS-SAMPLE-029).
    Placement {
        /// Integrity and callability of the read.
        sanger: SangerRejection,
        /// Why the placement was rejected.
        rejection: PlacementRejection,
    },
}

impl RejectionCause {
    /// Integrity and callability of the read, absent for a read too short to
    /// process.
    #[must_use]
    pub const fn sanger(&self) -> Option<&SangerRejection> {
        match self {
            Self::TooShort(_) => None,
            Self::Callability { sanger, .. } | Self::Placement { sanger, .. } => Some(sanger),
        }
    }
}
