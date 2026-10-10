//! Normalized and configured-filtered primary-sequence SNVs and small indels.

mod config;

pub use config::{RawVariantCallingConfig, VariantCallingConfig};

mod anchor;
pub(crate) mod eligibility;
mod extract;
mod filter;
mod mapping;

use crate::model::alignment::Alignment;
use crate::model::variant::VariantCallingResult;
use dna_kernel::error::Result;
use dna_kernel::model::reference::Reference;
use dna_kernel::read_evidence::ReadEvidence;

/// Extracts, anchors, and filters primary-sequence differences from one read's
/// alignment and modality evidence; `regions` are the target profile's
/// inclusive 1-based reportable regions.
pub(crate) fn call(
    alignment: &Alignment,
    reference: &Reference,
    evidence: &ReadEvidence,
    config: &VariantCallingConfig,
    regions: &[[usize; 2]],
) -> Result<VariantCallingResult> {
    let eligibility = eligibility::ReadEligibility::new(evidence, config);
    let extracted = extract::call(alignment, reference, &eligibility, config)?;
    filter::apply(
        extracted,
        evidence,
        &eligibility,
        alignment.orientation,
        regions,
    )
}
