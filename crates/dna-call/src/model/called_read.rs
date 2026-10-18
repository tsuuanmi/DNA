//! One read after the core caller: modality-neutral evidence, placement, and calls.

use std::collections::BTreeSet;

use crate::model::variant::{ExcludedVariant, VariantCallingResult};
use dna_align::model::alignment::{Alignment, ReferenceSegment};
use dna_kernel::read_evidence::ReadEvidence;

/// Modality-neutral products of one independently placed read.
///
/// This is the core's boundary between one-read calling and sample-level
/// aggregation (ADR-0069). The source filename is retained only as
/// reviewer-facing provenance and never constrains orientation, covered region,
/// or cross-read reconciliation.
#[derive(Debug, Clone)]
pub struct CalledRead {
    /// Reviewer-facing name of the read's source.
    pub input_name: String,
    /// SHA-256 content identity of the read's source.
    pub input_sha256: String,
    /// SHA-256 of the reference sequence the read was placed on.
    pub reference_sha256: String,
    /// SHA-256 of the configuration the read was called with.
    pub configuration_sha256: String,
    /// The modality evidence the core consumed.
    pub evidence: ReadEvidence,
    /// Selected placement of the read on the reference.
    pub alignment: Alignment,
    /// Variants called from the read, with eligibility.
    pub variants: VariantCallingResult,
}

impl CalledRead {
    /// The read with the positions of its excluded variants removed from its
    /// callable reference segments: a substitution's position, a deletion's
    /// deleted positions, and both positions beside an insertion's junction.
    /// A position DNA excluded a call at must not read as the reference
    /// there (SRS-VAR-019); the call workflow applies this to its records.
    #[must_use]
    pub fn without_excluded_positions(mut self) -> Self {
        // Every observed variant excluded from reporting, and every candidate
        // excluded before it could be observed (an unresolved allele, an indel
        // no flank reads).
        let mut excluded = self
            .variants
            .excluded
            .iter()
            .flat_map(|candidate| candidate.undecided.iter().copied())
            .collect::<BTreeSet<_>>();
        for observed in &self.variants.observed {
            if !observed.exclusion_reasons.is_empty() {
                excluded.extend(ExcludedVariant::of(&observed.variant, Vec::new()).undecided);
            }
        }
        self.alignment.callable_segments = self
            .alignment
            .callable_segments
            .iter()
            .flat_map(|segment| {
                let mut pieces = Vec::new();
                let mut start = segment.start_0based;
                for position in excluded
                    .range(segment.start_0based..segment.end_0based_exclusive)
                    .copied()
                {
                    if position > start {
                        pieces.push(ReferenceSegment {
                            start_0based: start,
                            end_0based_exclusive: position,
                        });
                    }
                    start = position + 1;
                }
                if segment.end_0based_exclusive > start {
                    pieces.push(ReferenceSegment {
                        start_0based: start,
                        end_0based_exclusive: segment.end_0based_exclusive,
                    });
                }
                pieces
            })
            .collect();
        self
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::variant::{ObservedVariant, Variant, VariantExclusionReason, VariantKind};
    use dna_align::{RawAlignmentConfig, align_best};
    use dna_kernel::model::reference::{Reference, ReferenceTopology};

    fn excluded(
        position_1based: usize,
        reference: &str,
        alternate: &str,
        kind: VariantKind,
    ) -> ObservedVariant {
        ObservedVariant {
            variant: Variant {
                contig: "ref".into(),
                position_1based,
                reference: reference.into(),
                alternate: alternate.into(),
                kind,
                calls: Vec::new(),
            },
            exclusion_reasons: vec![VariantExclusionReason::ReadEnd],
        }
    }

    #[test]
    fn excluded_variants_leave_their_positions_out_of_the_callable_segments()
    -> Result<(), Box<dyn std::error::Error>> {
        let sequence = "TTGACGTCAGTACGATCGTACCTGAGTACGAGGACT";
        let reference = Reference {
            name: "ref".into(),
            sequence: sequence.into(),
            topology: ReferenceTopology::Linear,
            sequence_sha256: String::new(),
        };
        let config = RawAlignmentConfig {
            match_score: 3,
            mismatch_score: -5,
            ambiguous_score: 0,
            gap_open_score: -10,
            gap_extension_score: -4,
            minimum_callable_bases: 1,
            minimum_identity: 0.8,
        }
        .validate()?;
        let evidence = ReadEvidence::clean(sequence);
        let alignment = align_best(&evidence, &reference, &config)?;
        let read = CalledRead {
            input_name: "r".into(),
            input_sha256: "r".into(),
            reference_sha256: String::new(),
            configuration_sha256: String::new(),
            evidence,
            alignment,
            variants: VariantCallingResult {
                reported: Vec::new(),
                // An SNV at 5, a deletion of 11-12 anchored at 10, and an
                // insertion after 20 (1-based), all excluded.
                observed: vec![
                    excluded(5, "C", "T", VariantKind::Snv),
                    excluded(10, "GTA", "G", VariantKind::Del),
                    excluded(20, "A", "AT", VariantKind::Ins),
                ],
                excluded: Vec::new(),
            },
        };
        let segments = read
            .without_excluded_positions()
            .alignment
            .callable_segments
            .iter()
            .map(|segment| (segment.start_0based, segment.end_0based_exclusive))
            .collect::<Vec<_>>();
        assert_eq!(segments, [(0, 4), (5, 10), (12, 19), (21, 36)]);
        Ok(())
    }
}
