//! Serializable `dna.notation/v1` contract: post-calling notation of a
//! variants document and its conformance findings.

use serde::Serialize;

use crate::model::result::{PluginResult, ProfileResult, ReferenceResult};
use crate::model::sample_result::SampleNotationResult;

/// Successful notation document.
#[derive(Debug, Serialize)]
pub(crate) struct NotationResult {
    pub(crate) schema_version: &'static str,
    pub(crate) sample_id: String,
    pub(crate) provenance: NotationProvenanceResult,
    pub(crate) notation: SampleNotationResult,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) conformance: Option<ConformanceResult>,
}

/// Identities of the inputs and plugins that produced the notation.
#[derive(Debug, Serialize)]
pub(crate) struct NotationProvenanceResult {
    pub(crate) reference: ReferenceResult,
    pub(crate) configuration_sha256: String,
    pub(crate) profile: ProfileResult,
    pub(crate) plugins: Vec<PluginResult>,
    pub(crate) source: SourceResult,
}

/// The variants document the notation was derived from.
#[derive(Debug, Serialize)]
pub(crate) struct SourceResult {
    pub(crate) schema_version: &'static str,
    pub(crate) sha256: String,
}

/// The profile's notation conventions and the calls that do not follow them.
#[derive(Debug, Serialize)]
pub(crate) struct ConformanceResult {
    pub(crate) rules: Vec<&'static str>,
    pub(crate) findings: Vec<FindingResult>,
}

/// One represented variant of one read that does not follow one rule.
#[derive(Debug, Serialize)]
pub(crate) struct FindingResult {
    pub(crate) rule: &'static str,
    pub(crate) read: String,
    /// The variant's rendered calls.
    pub(crate) calls: Vec<String>,
}
