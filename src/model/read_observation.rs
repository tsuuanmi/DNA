//! Complete scientific observation produced from one independently processed Sanger read.

use crate::model::attachment::SangerAttachment;
use crate::model::called_read::CalledRead;

/// A Sanger read after reference placement: the core's modality-neutral
/// products and the Sanger evidence that reports join to them.
#[derive(Debug, Clone)]
pub(crate) struct ReadObservation {
    pub(crate) called: CalledRead,
    pub(crate) sanger: SangerAttachment,
}
