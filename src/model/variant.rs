//! Normalized primary-sequence differences and original-call mappings.

use serde::Serialize;

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
    pub(crate) role: VariantCallRole,
    pub(crate) call_index_0based: usize,
    /// Absent only for inserted query calls, which have no reference base.
    pub(crate) reference_position_0based: Option<usize>,
}

/// One normalized reportable primary-sequence difference.
#[derive(Debug, Clone)]
pub struct Variant {
    pub(crate) contig: String,
    pub(crate) position_1based: usize,
    pub(crate) reference: String,
    pub(crate) alternate: String,
    pub(crate) kind: VariantKind,
    pub(crate) calls: Vec<VariantCallMapping>,
}

/// Stable reason a primary-difference candidate was not reportable.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum VariantExclusionReason {
    /// At least one changed base was not canonical A/C/G/T.
    NonCanonicalAllele,
    /// An insertion or deletion exceeded the configured length cap.
    IndelLengthExceeded,
    /// The normalized anchor was outside every configured region.
    OutsideConfiguredRegion,
    /// At least one supporting call was below the configured peak floor.
    PeakBelowMinimum,
    /// At least one supporting call did not strictly exceed the quality threshold.
    RelativeQualityNotAboveThreshold,
    /// An SNV supporting call retained more than one co-localized qualifying channel.
    #[serde(rename = "mixed_supporting_dna")]
    MixedSupportingDNA,
}

impl VariantExclusionReason {
    /// Stable operational-log label.
    pub(crate) const fn label(self) -> &'static str {
        match self {
            Self::NonCanonicalAllele => "non_canonical_allele",
            Self::IndelLengthExceeded => "indel_length_exceeded",
            Self::OutsideConfiguredRegion => "outside_configured_region",
            Self::PeakBelowMinimum => "peak_below_minimum",
            Self::RelativeQualityNotAboveThreshold => "relative_quality_not_above_threshold",
            Self::MixedSupportingDNA => "mixed_supporting_dna",
        }
    }
}

/// Concise diagnostic for one excluded candidate, intentionally without alleles.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ExcludedVariant {
    pub(crate) contig: String,
    pub(crate) position_1based: Option<usize>,
    pub(crate) kind: VariantKind,
    pub(crate) reasons: Vec<VariantExclusionReason>,
}

/// One normalized canonical variant observed before configured eligibility filtering.
#[derive(Debug, Clone)]
pub struct ObservedVariant {
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
pub struct VariantCallingResult {
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
