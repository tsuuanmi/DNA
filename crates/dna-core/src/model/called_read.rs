//! One read after the core caller: modality-neutral evidence, placement, and calls.

use crate::model::alignment::Alignment;
use crate::model::variant::VariantCallingResult;
use dna_kernel::read_evidence::ReadEvidence;

/// Modality-neutral products of one independently placed read.
///
/// This is the core's boundary between one-read calling and sample-level
/// aggregation (ADR-0069). The source filename is retained only as
/// reviewer-facing provenance and never constrains orientation, covered region,
/// or cross-read reconciliation.
#[derive(Debug, Clone)]
pub struct CalledRead {
    /// Reviewer-facing name of the read's source.
    pub input_name: String,
    /// SHA-256 content identity of the read's source.
    pub input_sha256: String,
    /// SHA-256 of the reference sequence the read was placed on.
    pub reference_sha256: String,
    /// SHA-256 of the configuration the read was called with.
    pub configuration_sha256: String,
    /// The modality evidence the core consumed.
    pub evidence: ReadEvidence,
    /// Selected placement of the read on the reference.
    pub alignment: Alignment,
    /// Variants called from the read, with eligibility.
    pub variants: VariantCallingResult,
}
