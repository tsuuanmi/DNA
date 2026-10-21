//! A run context in memory (`dna.run_context/v1`, PROP-0014), as the
//! command builds it and the input adapter reads it back.

use std::collections::BTreeSet;

use dna_consensus::noise::LengthStatistics;

/// The length section's method version.
pub const LENGTH_METHOD: u32 = 1;

/// A run context: a batch's own statistics.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct RunContext {
    /// The caller's batch identifier.
    pub run_id: String,
    /// The content SHA-256 of every read that contributed, sorted.
    pub reads: BTreeSet<String>,
    /// Reads the batch's samples set aside or the core could not place.
    pub rejected_reads: usize,
    /// The SHA-256 of the reference the reads were placed against.
    pub reference_sha256: String,
    /// The profile's recorded anchored windows the length pairs come from,
    /// each `[first, last]` (1-based, inclusive).
    pub windows: Vec<[usize; 2]>,
    /// The length section.
    pub length: LengthStatistics,
    /// The SHA-256 of the document it was read from, when it was.
    pub sha256: Option<String>,
}

/// How one analysis used a run context, for its provenance.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct ContextUse {
    /// The context's batch identifier.
    pub run_id: String,
    /// The SHA-256 of its document, when it was read from one.
    pub sha256: Option<String>,
    /// Whether the sample's own pairs were taken out.
    pub left_out: bool,
    /// τ's scale for run 1 and run 2; `None` where the batch's pairs, the
    /// sample's own out, were too few.
    pub scale: [Option<f64>; 2],
}
