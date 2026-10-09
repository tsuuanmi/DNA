//! Projection of variant-associated calls into concise signal records.

use crate::model::result::{PeakHeightsResult, VariantCallResult, VariantResult};
use crate::report::sanger_call;
use dna_core::model::alignment::Orientation;
use dna_core::model::reference_call;
use dna_core::model::variant::Variant;
use dna_kernel::error::{ReportError, Result};
use dna_kernel::read_evidence::ReadEvidence;
use dna_sanger::model::attachment::SangerAttachment;

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
