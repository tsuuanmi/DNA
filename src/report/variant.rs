//! Projection of variant-associated calls into concise signal records.

use crate::error::{ReportError, Result};
use crate::model::alignment::Orientation;
use crate::model::basecalls::BaseCalls;
use crate::model::quality::QualityControlResult;
use crate::model::reference_call;
use crate::model::result::{PeakHeightsResult, VariantCallResult, VariantResult};
use crate::model::variant::Variant;

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
    let projected_calls =
        reference_call::resolve_public_calls(calls, quality, orientation, &variant.calls)
            .map_err(ReportError::CallEvidence)?
            .into_iter()
            .map(|call| VariantCallResult {
                role: call.mapping.role,
                base: call.evidence.base,
                peaks: PeakHeightsResult::from(call.evidence.peak_heights),
                quality: call.evidence.quality,
            })
            .collect();
    Ok(VariantResult {
        position: variant.position_1based,
        reference: variant.reference,
        alternate: variant.alternate,
        kind: variant.kind,
        calls: projected_calls,
    })
}
