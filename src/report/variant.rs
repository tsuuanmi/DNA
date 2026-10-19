//! Projection of variant-associated calls into concise signal records.

use crate::error::{CallEvidenceError, ReportError, Result};
use crate::model::result::{PeakHeightsResult, VariantCallResult, VariantResult};
use crate::model::sanger_read::SangerCallSignal;
use crate::report::sanger_call;
use dna_align::model::alignment::Orientation;
use dna_call::model::reference_call;
use dna_call::model::variant::Variant;
use dna_kernel::read_evidence::ReadEvidence;

/// Projects normalized variants and joins their original calls to essential evidence.
pub(super) fn project(
    variants: Vec<Variant>,
    evidence: &ReadEvidence,
    sanger: &[SangerCallSignal],
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
    sanger: &[SangerCallSignal],
    orientation: Orientation,
) -> Result<VariantResult> {
    let projected_calls =
        reference_call::resolve_public_calls(evidence, orientation, &variant.calls)
            .map_err(ReportError::CallEvidence)?
            .into_iter()
            .map(|call| {
                let index = call.mapping.call_index_0based;
                let signal = sanger
                    .get(index)
                    .copied()
                    .ok_or(CallEvidenceError::MissingCall { index })
                    .map_err(ReportError::CallEvidence)?;
                let joined = sanger_call::on_reference(signal, orientation, index)
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
