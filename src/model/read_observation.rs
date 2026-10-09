//! Complete scientific observation produced from one independently processed Sanger read.

use dna_core::model::called_read::CalledRead;
use dna_sanger::model::attachment::SangerAttachment;

/// A Sanger read after reference placement: the core's modality-neutral
/// products and the Sanger evidence that reports join to them.
#[derive(Debug, Clone)]
pub(crate) struct ReadObservation {
    pub(crate) called: CalledRead,
    pub(crate) sanger: SangerAttachment,
}
