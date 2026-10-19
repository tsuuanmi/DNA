//! Serializable `dna.signal/v1` run signal contract (PROP-0009).

use serde::Serialize;

use crate::model::result::PluginResult;

/// A run's second-channel signal.
#[derive(Debug, Serialize)]
pub(crate) struct RunSignalResult {
    pub(crate) schema_version: &'static str,
    pub(crate) run_id: String,
    pub(crate) provenance: RunSignalProvenanceResult,
    /// Content SHA-256 of every trace of the run, sorted.
    pub(crate) traces: Vec<String>,
    /// Reads that contributed calls.
    pub(crate) reads: usize,
    /// Reads rejected before their evidence (too short, too few callable
    /// calls).
    pub(crate) rejected_reads: usize,
    /// Calls that contributed.
    pub(crate) calls: u64,
    /// Second-channel shares by key, sorted by key.
    pub(crate) keys: Vec<RunSignalKeyResult>,
}

/// Inputs and selection a run signal was made with.
#[derive(Debug, Serialize)]
pub(crate) struct RunSignalProvenanceResult {
    pub(crate) configuration_sha256: String,
    pub(crate) plugins: Vec<PluginResult>,
    /// Version of the call selection.
    pub(crate) selection: u32,
    /// Calls of context on each side of a key's call.
    pub(crate) flank: usize,
    /// Steps of a share from 0 to 1.
    pub(crate) resolution: u16,
}

/// One key: the read-strand context, the call's base, its second base, and
/// the counts of calls by second-channel share in `resolution` steps.
#[derive(Debug, Serialize)]
pub(crate) struct RunSignalKeyResult {
    pub(crate) context: String,
    pub(crate) base: char,
    pub(crate) second: char,
    /// `[share_permille, calls]` pairs, by share.
    pub(crate) shares: Vec<[u64; 2]>,
}
