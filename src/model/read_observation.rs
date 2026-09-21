//! Complete scientific observation produced from one independently processed read.

use crate::model::alignment::Alignment;
use crate::model::basecalls::BaseCalls;
use crate::model::quality::QualityControlResult;
use crate::model::signal::DNAAnalysis;
use crate::model::variant::VariantCallingResult;

/// Immutable read-level products after evidence-driven reference placement.
///
/// This is the boundary between one-read processing and sample-level
/// reconciliation. Placement is already derived from alignment evidence. The
/// source filename is retained only as reviewer-facing provenance and never
/// constrains orientation, covered region, or cross-read reconciliation.
#[derive(Debug, Clone)]
pub(crate) struct ReadObservation {
    pub(crate) input_name: String,
    pub(crate) input_sha256: String,
    pub(crate) reference_sha256: String,
    pub(crate) configuration_sha256: String,
    pub(crate) calls: BaseCalls,
    pub(crate) signal: DNAAnalysis,
    pub(crate) quality: QualityControlResult,
    pub(crate) alignment: Alignment,
    pub(crate) variants: VariantCallingResult,
}
