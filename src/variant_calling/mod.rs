//! Normalized and configured-filtered primary-sequence SNVs and small indels.

mod extract;
mod filter;
mod mapping;
mod normalize;

use crate::config::VariantCallingConfig;
use crate::error::Result;
use crate::model::alignment::Alignment;
use crate::model::basecalls::BaseCalls;
use crate::model::quality::QualityControlResult;
use crate::model::reference::Reference;
use crate::model::variant::VariantCallingResult;

/// Extracts, normalizes, and filters primary-sequence differences.
pub(crate) fn call(
    alignment: &Alignment,
    reference: &Reference,
    calls: &BaseCalls,
    quality: &QualityControlResult,
    config: &VariantCallingConfig,
) -> Result<VariantCallingResult> {
    let extracted = extract::call(alignment, reference, config)?;
    filter::apply(extracted, calls, quality, config)
}
