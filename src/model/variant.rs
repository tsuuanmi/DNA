//! Normalized primary-sequence differences and original-call mappings.

use serde::Serialize;

use crate::read_evidence::EvidenceReason;

/// Supported primary-difference type.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize)]
#[serde(rename_all = "UPPERCASE")]
pub(crate) enum VariantKind {
    /// Single-nucleotide substitution.
    Snv,
    /// Insertion relative to the reference.
    Ins,
    /// Deletion relative to the reference.
    Del,
}

impl VariantKind {
    /// Stable uppercase label used by reports and operational logs.
    pub(crate) const fn label(self) -> &'static str {
        match self {
            Self::Snv => "SNV",
            Self::Ins => "INS",
            Self::Del => "DEL",
        }
    }
}

/// How an original trace call relates to a reported difference.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize)]
#[serde(rename_all = "snake_case")]
pub(crate) enum VariantCallRole {
    /// A query call directly contributes an observed alternate base.
    Supporting,
    /// A reference-aligned query call bounds an indel.
    Flanking,
}

/// Mapping from a variant-associated call to the selected alignment.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) struct VariantCallMapping {
    pub(crate) role: VariantCallRole,
    pub(crate) call_index_0based: usize,
    /// Absent only for inserted query calls, which have no reference base.
    pub(crate) reference_position_0based: Option<usize>,
}

/// One normalized reportable primary-sequence difference.
#[derive(Debug, Clone)]
pub(crate) struct Variant {
    pub(crate) contig: String,
    pub(crate) position_1based: usize,
    pub(crate) reference: String,
    pub(crate) alternate: String,
    pub(crate) kind: VariantKind,
    pub(crate) calls: Vec<VariantCallMapping>,
}

/// Stable reason a primary-difference candidate was not reportable.
///
/// The core owns the first four reasons. Every other reason comes from the
/// read's modality evidence and is reported verbatim; the core never
/// interprets it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum VariantExclusionReason {
    /// At least one changed base was not canonical A/C/G/T.
    NonCanonicalAllele,
    /// An insertion or deletion exceeded the configured length cap.
    IndelLengthExceeded,
    /// The normalized anchor was outside every target-profile region.
    OutsideTargetRegion,
    /// A mapped call lies within the read-end margin of an uninformative call.
    ReadEnd,
    /// A support veto or mask reason supplied by the read's modality.
    Evidence(EvidenceReason),
}

impl VariantExclusionReason {
    /// Stable published label.
    pub(crate) const fn label(self) -> &'static str {
        match self {
            Self::NonCanonicalAllele => "non_canonical_allele",
            Self::IndelLengthExceeded => "indel_length_exceeded",
            Self::OutsideTargetRegion => "outside_target_region",
            Self::ReadEnd => "read_end",
            Self::Evidence(reason) => reason.label(),
        }
    }
}

impl Serialize for VariantExclusionReason {
    fn serialize<S: serde::Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        serializer.serialize_str(self.label())
    }
}

/// Concise diagnostic for one excluded candidate, intentionally without alleles.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct ExcludedVariant {
    pub(crate) contig: String,
    pub(crate) position_1based: Option<usize>,
    pub(crate) kind: VariantKind,
    pub(crate) reasons: Vec<VariantExclusionReason>,
}

/// One normalized canonical variant observed before configured eligibility filtering.
#[derive(Debug, Clone)]
pub(crate) struct ObservedVariant {
    pub(crate) variant: Variant,
    pub(crate) exclusion_reasons: Vec<VariantExclusionReason>,
}

impl ObservedVariant {
    /// Whether this normalized observation is reportable under the active config.
    pub(crate) fn eligible(&self) -> bool {
        self.exclusion_reasons.is_empty()
    }
}

/// Variant stage output.
#[derive(Debug, Clone)]
pub(crate) struct VariantCallingResult {
    pub(crate) reported: Vec<Variant>,
    pub(crate) observed: Vec<ObservedVariant>,
    pub(crate) excluded: Vec<ExcludedVariant>,
}

impl VariantCallingResult {
    /// Number of candidates excluded across extraction and configured filtering.
    pub(crate) fn excluded_count(&self) -> usize {
        self.excluded.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Pins every published exclusion-reason label in its serialized form.
    #[test]
    fn serializes_every_exclusion_reason_with_its_label() -> serde_json::Result<()> {
        let reasons = [
            (
                VariantExclusionReason::NonCanonicalAllele,
                "non_canonical_allele",
            ),
            (
                VariantExclusionReason::IndelLengthExceeded,
                "indel_length_exceeded",
            ),
            (
                VariantExclusionReason::OutsideTargetRegion,
                "outside_target_region",
            ),
            (VariantExclusionReason::ReadEnd, "read_end"),
            (
                VariantExclusionReason::Evidence(EvidenceReason::new("mixed_supporting_dna")),
                "mixed_supporting_dna",
            ),
        ];
        for (reason, label) in reasons {
            assert_eq!(reason.label(), label);
            assert_eq!(serde_json::to_string(&reason)?, format!("\"{label}\""));
        }
        Ok(())
    }
}
