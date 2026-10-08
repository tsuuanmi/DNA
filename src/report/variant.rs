//! Projection of variant-associated calls into concise signal records.

use crate::error::{ReportError, Result};
use crate::model::alignment::Orientation;
use crate::model::basecalls::BaseCalls;
use crate::model::quality::QualityControlResult;
use crate::model::reference_call;
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
    let evidence = reference_call::resolve(calls, quality, orientation, mapping.call_index_0based)
        .map_err(ReportError::CallEvidence)?;
    Ok(VariantCallResult {
        role: mapping.role,
        base: evidence.base,
        peaks: PeakHeightsResult::from(evidence.peak_heights),
        quality: evidence.quality,
    })
}
