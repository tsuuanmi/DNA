//! Normalized and configured-filtered primary-sequence SNVs and small indels.

pub mod components;
pub mod model;
pub mod read_call;

mod config;

pub use config::{RawVariantCallingConfig, VariantCallingConfig};
pub use eligibility::{Eligibility, read_trust};

mod anchor;
mod eligibility;
mod extract;
mod filter;
mod mapping;

use crate::model::variant::VariantCallingResult;
use dna_align::model::alignment::Alignment;
use dna_kernel::error::Result;
use dna_kernel::model::reference::Reference;
use dna_kernel::read_evidence::ReadEvidence;

/// Extracts, anchors, and filters primary-sequence differences from one read's
/// alignment and modality evidence; `regions` are the target profile's
/// inclusive 1-based reportable regions.
///
/// # Errors
///
/// Returns `VariantError` when a difference cannot be anchored or mapped to
/// its calls.
pub fn call(
    alignment: &Alignment,
    reference: &Reference,
    evidence: &ReadEvidence,
    config: &VariantCallingConfig,
    regions: &[[usize; 2]],
) -> Result<VariantCallingResult> {
    let eligibility = Eligibility::new(evidence, config);
    let extracted = extract::call(alignment, reference, &eligibility, config)?;
    filter::apply(
        extracted,
        evidence,
        &eligibility,
        alignment.orientation(),
        regions,
        reference,
    )
}
