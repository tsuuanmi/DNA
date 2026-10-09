//! Serializable `dna.variants/v1` contract: the core's calls from
//! modality-neutral evidence, without modality evidence.

use serde::Serialize;

use crate::model::result::AlignmentResult;
use crate::model::sample_result::{SampleNotationResult, SampleProvenanceResult};
use crate::model::variant::{VariantExclusionReason, VariantKind};

/// Successful variants document.
#[derive(Debug, Serialize)]
pub(crate) struct VariantsResult {
    pub(crate) schema_version: &'static str,
    pub(crate) sample_id: String,
    pub(crate) provenance: SampleProvenanceResult,
    pub(crate) reads: Vec<ReadVariantsResult>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) notation: Option<SampleNotationResult>,
}

/// One read's placement and observed variants.
#[derive(Debug, Serialize)]
pub(crate) struct ReadVariantsResult {
    pub(crate) name: String,
    pub(crate) sha256: String,
    pub(crate) alignment: AlignmentResult,
    pub(crate) variants: Vec<ObservedVariantResult>,
}

/// One observed variant of a read with its eligibility.
#[derive(Debug, Serialize)]
pub(crate) struct ObservedVariantResult {
    pub(crate) position: usize,
    pub(crate) reference: String,
    pub(crate) alternate: String,
    pub(crate) kind: VariantKind,
    pub(crate) eligible: bool,
    pub(crate) exclusion_reasons: Vec<VariantExclusionReason>,
}
