//! Serializable `dna.basecalls/v2` reference-free result contract.

use serde::Serialize;

use crate::model::result::{DNAQualityResult, InputResult, IntervalResult};

/// Successful reference-free basecall document.
#[derive(Debug, Serialize)]
pub struct BasecallResult {
    pub(crate) schema_version: &'static str,
    pub(crate) provenance: BasecallProvenanceResult,
    pub(crate) read: BasecallReadResult,
    pub(crate) signal_quality: DNAQualityResult,
    pub(crate) warnings: BasecallWarningSummaryResult,
}

/// Deterministic input identities for a basecall operation.
#[derive(Debug, Serialize)]
pub struct BasecallProvenanceResult {
    pub(crate) input: InputResult,
    pub(crate) configuration_sha256: String,
}

/// Called sequences and the retained primary interval.
#[derive(Debug, Serialize)]
pub struct BasecallReadResult {
    pub(crate) call_count: usize,
    pub(crate) primary: String,
    pub(crate) ambiguity: String,
    pub(crate) retained: String,
    pub(crate) trim: IntervalResult,
}

/// Public non-fatal basecall counts.
#[derive(Debug, Serialize)]
pub struct BasecallWarningSummaryResult {
    pub(crate) unresolved_primary_calls: usize,
    pub(crate) multi_channel_unresolved_calls: usize,
    pub(crate) vendor_disagreements: usize,
    pub(crate) ploc_vendor_length_mismatches: usize,
    pub(crate) clipped_channel_samples: usize,
}
