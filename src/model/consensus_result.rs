//! Serializable `dna.consensus/v1` result contract (PROP-0003).

use serde::Serialize;

use crate::model::result::IntervalResult;
use crate::model::sample_result::{ReadRejectionResult, SampleProvenanceResult};

/// Top-level sample consensus document.
#[derive(Debug, Serialize)]
pub(crate) struct ConsensusResult {
    pub(crate) schema_version: &'static str,
    pub(crate) sample_id: String,
    pub(crate) provenance: SampleProvenanceResult,
    pub(crate) reads: Vec<ConsensusReadResult>,
    pub(crate) rejected_reads: Vec<ConsensusRejectedReadResult>,
    pub(crate) segments: Vec<ConsensusSegmentResult>,
    pub(crate) sites: Vec<ConsensusSiteResult>,
    pub(crate) summary: ConsensusSummaryResult,
}

/// One admitted read, in the sample's read order.
#[derive(Debug, Serialize)]
pub(crate) struct ConsensusReadResult {
    pub(crate) name: String,
    pub(crate) sha256: String,
}

/// One read set aside before the consensus.
#[derive(Debug, Serialize)]
pub(crate) struct ConsensusRejectedReadResult {
    pub(crate) name: String,
    pub(crate) sha256: String,
    pub(crate) rejection: ReadRejectionResult,
}

/// One consensus segment; `name` is its consensus FASTA record identifier.
#[derive(Debug, Serialize)]
pub(crate) struct ConsensusSegmentResult {
    pub(crate) name: String,
    pub(crate) reference: IntervalResult,
    pub(crate) sequence: String,
}

/// One notable site: a reference interval and its decision.
#[derive(Debug, Serialize)]
pub(crate) struct ConsensusSiteResult {
    /// First reference position (1-based, inclusive).
    pub(crate) start: usize,
    /// Last reference position (1-based, inclusive).
    pub(crate) end: usize,
    /// Reference sequence over the interval.
    pub(crate) reference: String,
    /// Decided sequence over the interval, insertions at its edges included
    /// (empty when every position is deleted), or `None` when undecided.
    pub(crate) call: Option<String>,
    pub(crate) state: &'static str,
    pub(crate) supporting_reads: Vec<String>,
    pub(crate) opposing_reads: Vec<String>,
    pub(crate) uninformative_reads: Vec<String>,
    /// The call's runs and how each length is known, for a stretch decided
    /// run by run (ADR-0074); empty otherwise.
    pub(crate) runs: Vec<ConsensusRunResult>,
}

/// One run of a call: its base, its length, and how the length is known.
#[derive(Debug, Serialize)]
pub(crate) struct ConsensusRunResult {
    pub(crate) base: char,
    pub(crate) length: usize,
    pub(crate) length_evidence: &'static str,
}

/// Decision counts over the segments.
#[derive(Debug, Serialize)]
pub(crate) struct ConsensusSummaryResult {
    pub(crate) called_positions: usize,
    pub(crate) contested_sites: usize,
    pub(crate) unresolved_positions: usize,
    pub(crate) phase_loss_runs: usize,
    pub(crate) reference_frame_runs: usize,
}
