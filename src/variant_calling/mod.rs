//! Normalized and configured-filtered primary-sequence SNVs and small indels.

mod anchor;
mod eligibility;
mod extract;
mod filter;
mod mapping;

use crate::config::VariantCallingConfig;
use crate::error::Result;
use crate::model::alignment::Alignment;
use crate::model::basecalls::BaseCalls;
use crate::model::callability::ReadCallability;
use crate::model::quality::QualityControlResult;
use crate::model::reference::Reference;
use crate::model::variant::VariantCallingResult;

/// Extracts, anchors, and filters primary-sequence differences; `regions` are
/// the target profile's inclusive 1-based reportable regions.
pub(crate) fn call(
    alignment: &Alignment,
    reference: &Reference,
    calls: &BaseCalls,
    quality: &QualityControlResult,
    callability: &ReadCallability,
    config: &VariantCallingConfig,
    regions: &[[usize; 2]],
) -> Result<VariantCallingResult> {
    let eligibility = eligibility::ReadEligibility::new(quality, callability, config);
    let extracted = extract::call(alignment, reference, &eligibility, config)?;
    filter::apply(extracted, calls, quality, &eligibility, config, regions)
}
