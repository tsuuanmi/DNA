//! The run context contract `dna.run_context/v1` (PROP-0014) as written, and
//! the record of how an analysis used one.

use serde::Serialize;

use super::result::{PluginResult, ProfileResult, ReferenceResult};

/// A batch's run context.
#[derive(Debug, Serialize)]
pub(crate) struct RunContextResult {
    pub(crate) schema_version: &'static str,
    pub(crate) run_id: String,
    pub(crate) provenance: RunContextProvenanceResult,
    /// The content SHA-256 of every read that contributed, sorted.
    pub(crate) reads: Vec<String>,
    /// Reads set aside or not placed.
    pub(crate) rejected_reads: usize,
    /// Samples that contributed reads.
    pub(crate) samples: usize,
    pub(crate) length: LengthSectionResult,
}

/// What a run context was made with.
#[derive(Debug, Serialize)]
pub(crate) struct RunContextProvenanceResult {
    pub(crate) reference: ReferenceResult,
    pub(crate) configuration_sha256: String,
    pub(crate) profile: ProfileResult,
    pub(crate) plugins: Vec<PluginResult>,
    /// Listed input files left out: missing, not regular, empty, or
    /// unreadable.
    pub(crate) skipped: Vec<SkippedResult>,
}

/// A listed input file a run context left out.
#[derive(Debug, Serialize)]
pub(crate) struct SkippedResult {
    pub(crate) sample_id: String,
    pub(crate) name: String,
    pub(crate) reason: &'static str,
}

/// The length section: the batch's cross-strand pairs of the two-strand
/// rule's views, per run.
#[derive(Debug, Serialize)]
pub(crate) struct LengthSectionResult {
    pub(crate) method: u32,
    /// The profile's recorded anchored windows the pairs come from.
    pub(crate) windows: Vec<[usize; 2]>,
    /// The forward strand's τ the differences were normalised with.
    pub(crate) noise_forward: [f64; 2],
    /// The reverse strand's.
    pub(crate) noise_reverse: [f64; 2],
    /// The normalised differences' bin edges.
    pub(crate) edges: Vec<f64>,
    pub(crate) runs: Vec<LengthRunStatisticsResult>,
}

/// One run's additive statistics.
#[derive(Debug, Serialize)]
pub(crate) struct LengthRunStatisticsResult {
    pub(crate) run: usize,
    pub(crate) pairs: usize,
    pub(crate) sum_z2: f64,
    pub(crate) sum_z4: f64,
    pub(crate) histogram: Vec<usize>,
}

/// How an analysis used a run context, in its provenance.
#[derive(Debug, Serialize)]
pub(crate) struct ContextUseResult {
    pub(crate) run_id: String,
    /// The document's SHA-256.
    pub(crate) sha256: Option<String>,
    /// Whether the sample's own pairs were taken out.
    pub(crate) left_out: bool,
    /// τ's scale for run 1 and run 2; null where the batch's pairs were too
    /// few.
    pub(crate) length_scale: [Option<f64>; 2],
}

impl ContextUseResult {
    pub(crate) fn of(used: &crate::model::run_context::ContextUse) -> Self {
        Self {
            run_id: used.run_id.clone(),
            sha256: used.sha256.clone(),
            left_out: used.left_out,
            length_scale: used.scale,
        }
    }
}
