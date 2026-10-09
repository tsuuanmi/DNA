//! Normalized primary-sequence differences and original-call mappings.

use dna_kernel::variant as public;
use serde::Serialize;

use dna_kernel::read_evidence::EvidenceReason;

/// Supported primary-difference type.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize)]
#[serde(rename_all = "UPPERCASE")]
pub enum VariantKind {
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
pub enum VariantCallRole {
    /// A query call directly contributes an observed alternate base.
    Supporting,
    /// A reference-aligned query call bounds an indel.
    Flanking,
}

/// Mapping from a variant-associated call to the selected alignment.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct VariantCallMapping {
    /// Whether the call supports or flanks the variant.
    pub role: VariantCallRole,
    /// Index of the source call in the read (0-based).
    pub call_index_0based: usize,
    /// Absent only for inserted query calls, which have no reference base.
    pub(crate) reference_position_0based: Option<usize>,
}

/// One normalized reportable primary-sequence difference.
#[derive(Debug, Clone)]
pub struct Variant {
    /// Reference record the variant is placed on.
    pub contig: String,
    /// 1-based position of the first reference allele base.
    pub position_1based: usize,
    /// Anchored reference allele.
    pub reference: String,
    /// Anchored alternate allele.
    pub alternate: String,
    /// Variant type.
    pub kind: VariantKind,
    /// Mappings to the read calls that support or flank the variant.
    pub calls: Vec<VariantCallMapping>,
}

/// Stable reason a primary-difference candidate was not reportable.
///
/// The core owns the first four reasons. Every other reason comes from the
/// read's modality evidence and is reported verbatim; the core never
/// interprets it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VariantExclusionReason {
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
pub struct ObservedVariant {
    /// The observed variant.
    pub variant: Variant,
    /// Why the variant is not reportable; empty when eligible.
    pub exclusion_reasons: Vec<VariantExclusionReason>,
}

impl ObservedVariant {
    /// Whether this normalized observation is reportable under the active config.
    #[must_use]
    pub fn eligible(&self) -> bool {
        self.exclusion_reasons.is_empty()
    }
}

/// Variant stage output.
#[derive(Debug, Clone)]
pub struct VariantCallingResult {
    /// Eligible variants, sorted and deduplicated.
    pub reported: Vec<Variant>,
    /// Every observed variant with its eligibility.
    pub observed: Vec<ObservedVariant>,
    pub(crate) excluded: Vec<ExcludedVariant>,
}

impl VariantCallingResult {
    /// Number of candidates excluded across extraction and configured filtering.
    #[must_use]
    pub fn excluded_count(&self) -> usize {
        self.excluded.len()
    }
}

impl From<&Variant> for public::Variant {
    /// Projects an internal called variant into the public variant boundary.
    fn from(variant: &Variant) -> Self {
        Self {
            contig: variant.contig.clone(),
            position_1based: variant.position_1based,
            reference: variant.reference.clone(),
            alternate: variant.alternate.clone(),
            kind: match variant.kind {
                VariantKind::Snv => public::VariantKind::Snv,
                VariantKind::Ins => public::VariantKind::Ins,
                VariantKind::Del => public::VariantKind::Del,
            },
        }
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
