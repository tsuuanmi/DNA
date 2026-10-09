//! Projection of variant-associated calls into concise signal records.

use crate::error::{ReportError, Result};
use crate::model::alignment::Orientation;
use crate::model::read_observation::SangerAttachment;
use crate::model::reference_call;
use crate::model::result::{PeakHeightsResult, VariantCallResult, VariantResult};
use crate::model::variant::Variant;
use crate::read_evidence::ReadEvidence;
use crate::report::sanger_call;

/// Projects normalized variants and joins their original calls to essential evidence.
pub(super) fn project(
    variants: Vec<Variant>,
    evidence: &ReadEvidence,
    sanger: &SangerAttachment,
    orientation: Orientation,
) -> Result<Vec<VariantResult>> {
    variants
        .into_iter()
        .map(|variant| project_variant(variant, evidence, sanger, orientation))
        .collect()
}

fn project_variant(
    variant: Variant,
    evidence: &ReadEvidence,
    sanger: &SangerAttachment,
    orientation: Orientation,
) -> Result<VariantResult> {
    let projected_calls =
        reference_call::resolve_public_calls(evidence, orientation, &variant.calls)
            .map_err(ReportError::CallEvidence)?
            .into_iter()
            .map(|call| {
                let joined = sanger_call::evidence(
                    &sanger.calls,
                    &sanger.quality,
                    orientation,
                    call.mapping.call_index_0based,
                )
                .map_err(ReportError::CallEvidence)?;
                Ok(VariantCallResult {
                    role: call.mapping.role,
                    base: call.base,
                    peaks: PeakHeightsResult::from(joined.peak_heights),
                    quality: joined.quality,
                })
            })
            .collect::<Result<_>>()?;
    Ok(VariantResult {
        position: variant.position_1based,
        reference: variant.reference,
        alternate: variant.alternate,
        kind: variant.kind,
        calls: projected_calls,
    })
}
