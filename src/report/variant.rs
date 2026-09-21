//! Projection of variant-associated calls into concise signal records.

use crate::error::{Error, Result};
use crate::model::alignment::Orientation;
use crate::model::basecalls::BaseCalls;
use crate::model::quality::QualityControlResult;
use crate::model::result::{PeakHeightsResult, VariantCallResult, VariantResult};
use crate::model::variant::{Variant, VariantCallMapping};

/// Projects normalized variants and joins their original calls to essential evidence.
pub(super) fn project(
    variants: Vec<Variant>,
    calls: &BaseCalls,
    quality: &QualityControlResult,
    orientation: Orientation,
) -> Result<Vec<VariantResult>> {
    variants
        .into_iter()
        .map(|variant| project_variant(variant, calls, quality, orientation))
        .collect()
}

fn project_variant(
    variant: Variant,
    calls: &BaseCalls,
    quality: &QualityControlResult,
    orientation: Orientation,
) -> Result<VariantResult> {
    let projected_calls = variant
        .calls
        .into_iter()
        .map(|mapping| project_call(mapping, calls, quality, orientation))
        .collect::<Result<Vec<_>>>()?;
    Ok(VariantResult {
        position: variant.position_1based,
        reference: variant.reference,
        alternate: variant.alternate,
        kind: variant.kind,
        calls: projected_calls,
    })
}

fn project_call(
    mapping: VariantCallMapping,
    calls: &BaseCalls,
    quality: &QualityControlResult,
    orientation: Orientation,
) -> Result<VariantCallResult> {
    let index = mapping.call_index_0based;
    let call = calls
        .calls
        .get(index)
        .ok_or_else(|| Error::Report(format!("variant references missing call index {index}")))?;
    let score = quality.per_call.get(index).ok_or_else(|| {
        Error::Report(format!("variant references missing quality index {index}"))
    })?;
    if call.index_0based != index || score.index_0based != index {
        return Err(Error::Report(format!(
            "variant call index {index} does not match call/quality records"
        )));
    }
    let primary = call.primary_peak_evidence.as_ref().ok_or_else(|| {
        Error::Report(format!(
            "variant call index {index} lacks primary-event peak evidence"
        ))
    })?;
    Ok(VariantCallResult {
        role: mapping.role,
        base: orientation.reference_base(call.primary),
        peaks: PeakHeightsResult::from(orientation.reference_peak_heights(primary.channel_heights)),
        quality: score.relative_quality_score,
    })
}
