//! Resolution of variant call mappings to reference-strand called bases.

use crate::model::alignment::Orientation;
use crate::model::variant::{VariantCallMapping, VariantCallRole};
use dna_kernel::error::CallEvidenceError;
use dna_kernel::model::nucleotide::is_canonical;
use dna_kernel::read_evidence::ReadEvidence;

/// One variant call mapping with its reference-strand called base.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PublicCall {
    /// The variant's mapping to the source call.
    pub mapping: VariantCallMapping,
    /// Called base on the reference strand.
    pub base: char,
}

/// Resolves the public calls of one variant against the read's evidence.
///
/// A flanking call whose base is unresolved (`N`) has no called base to report,
/// so it is omitted rather than fabricated. Supporting calls are canonical by
/// construction and must resolve, and at least one call must remain. Modality
/// evidence for the remaining calls is joined by the report.
///
/// # Errors
///
/// Returns `CallEvidenceError` when a mapping names a missing call, a
/// supporting call is unresolved, or no call resolves.
pub fn resolve_public_calls(
    evidence: &ReadEvidence,
    orientation: Orientation,
    mappings: &[VariantCallMapping],
) -> Result<Vec<PublicCall>, CallEvidenceError> {
    let mut public = Vec::with_capacity(mappings.len());
    for &mapping in mappings {
        match resolve(evidence, orientation, mapping.call_index_0based) {
            Ok(base) => public.push(PublicCall { mapping, base }),
            Err(CallEvidenceError::UnresolvedCall { .. })
                if mapping.role == VariantCallRole::Flanking => {}
            Err(error) => return Err(error),
        }
    }
    if public.is_empty() {
        return Err(CallEvidenceError::NoResolvedCalls);
    }
    Ok(public)
}

/// Requires a canonical called base at `index` and projects it onto the
/// reference strand.
fn resolve(
    evidence: &ReadEvidence,
    orientation: Orientation,
    index: usize,
) -> Result<char, CallEvidenceError> {
    let call = evidence
        .calls()
        .get(index)
        .ok_or(CallEvidenceError::MissingCall { index })?;
    if !is_canonical(call.base) {
        return Err(CallEvidenceError::UnresolvedCall { index });
    }
    Ok(orientation.reference_base(call.base))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn mapping(role: VariantCallRole, index: usize) -> VariantCallMapping {
        VariantCallMapping {
            role,
            call_index_0based: index,
            reference_position_0based: Some(index),
        }
    }

    #[test]
    fn projects_reverse_reads_onto_the_reference_strand() {
        let evidence = ReadEvidence::clean("A");
        assert_eq!(resolve(&evidence, Orientation::Forward, 0), Ok('A'));
        assert_eq!(resolve(&evidence, Orientation::Reverse, 0), Ok('T'));
    }

    #[test]
    fn omits_unresolved_flanks_from_public_calls() {
        let evidence = ReadEvidence::clean("NGN");
        let mappings = [
            mapping(VariantCallRole::Flanking, 0),
            mapping(VariantCallRole::Flanking, 1),
            mapping(VariantCallRole::Flanking, 2),
        ];

        assert_eq!(
            resolve_public_calls(&evidence, Orientation::Forward, &mappings),
            Ok(vec![PublicCall {
                mapping: mappings[1],
                base: 'G',
            }])
        );
    }

    #[test]
    fn requires_supporting_calls_and_at_least_one_public_call() {
        let evidence = ReadEvidence::clean("NN");

        assert_eq!(
            resolve_public_calls(
                &evidence,
                Orientation::Forward,
                &[mapping(VariantCallRole::Supporting, 0)],
            ),
            Err(CallEvidenceError::UnresolvedCall { index: 0 })
        );
        assert_eq!(
            resolve_public_calls(
                &evidence,
                Orientation::Forward,
                &[
                    mapping(VariantCallRole::Flanking, 0),
                    mapping(VariantCallRole::Flanking, 1),
                ],
            ),
            Err(CallEvidenceError::NoResolvedCalls)
        );
        assert_eq!(
            resolve(&evidence, Orientation::Forward, 2),
            Err(CallEvidenceError::MissingCall { index: 2 })
        );
    }
}
