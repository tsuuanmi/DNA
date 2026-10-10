//! Complete scientific observation produced from one independently processed Sanger read.

use dna_core::model::called_read::CalledRead;
use dna_core::read_call::PlacementRejection;
use dna_sanger::model::attachment::{SangerAttachment, SangerRejection};
use dna_sanger::model::callability::ReadRejection;

/// A Sanger read after reference placement: the core's modality-neutral
/// products and the Sanger evidence that reports join to them.
#[derive(Debug, Clone)]
pub(crate) struct ReadObservation {
    pub(crate) called: CalledRead,
    pub(crate) sanger: SangerAttachment,
}

/// A Sanger read a sample operation set aside, with its Sanger evidence.
#[derive(Debug, Clone)]
pub(crate) struct RejectedRead {
    pub(crate) sanger: SangerRejection,
    pub(crate) cause: RejectionCause,
}

/// Why a read was set aside.
#[derive(Debug, Clone, Copy)]
pub(crate) enum RejectionCause {
    /// Too few callable calls (SRS-SAMPLE-027).
    Callability(ReadRejection),
    /// The core could not place the read (SRS-SAMPLE-029).
    Placement(PlacementRejection),
}
