//! One read after the core caller: modality-neutral evidence, placement, and calls.

use crate::model::alignment::Alignment;
use crate::model::variant::VariantCallingResult;
use crate::read_evidence::ReadEvidence;

/// Modality-neutral products of one independently placed read.
///
/// This is the core's boundary between one-read calling and sample-level
/// aggregation (ADR-0069). The source filename is retained only as
/// reviewer-facing provenance and never constrains orientation, covered region,
/// or cross-read reconciliation.
#[derive(Debug, Clone)]
pub(crate) struct CalledRead {
    pub(crate) input_name: String,
    pub(crate) input_sha256: String,
    pub(crate) reference_sha256: String,
    pub(crate) configuration_sha256: String,
    pub(crate) evidence: ReadEvidence,
    pub(crate) alignment: Alignment,
    pub(crate) variants: VariantCallingResult,
}
