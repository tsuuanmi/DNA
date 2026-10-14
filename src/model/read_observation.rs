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
///
/// Public through [`sample_analysis`](crate::sample_analysis).
#[derive(Debug, Clone)]
#[non_exhaustive]
pub struct RejectedRead {
    /// Name, identity, integrity, and callability of the rejected trace.
    pub sanger: SangerRejection,
    /// Why the read was set aside.
    pub cause: RejectionCause,
}

/// Why a read was set aside.
#[derive(Debug, Clone, Copy)]
#[non_exhaustive]
pub enum RejectionCause {
    /// Too few callable calls (SRS-SAMPLE-027).
    Callability(ReadRejection),
    /// The core could not place the read (SRS-SAMPLE-029).
    Placement(PlacementRejection),
}
