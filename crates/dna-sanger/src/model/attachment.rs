//! Sanger evidence kept beside the core's records, which reports join to them
//! by read identity and call index (ADR-0069).

use crate::model::basecalls::BaseCalls;
use crate::model::callability::{ReadCallability, ReadRejection};
use crate::model::locus_evidence::LocusEvidence;
use crate::model::quality::QualityControlResult;
use crate::model::signal::{SangerIntegrity, SignalAnalysis};

/// Sanger products of one read that the core never reads; reports join them
/// to core records by read identity and call index (ADR-0069).
#[derive(Debug, Clone)]
pub struct SangerAttachment {
    /// Signal-derived base calls.
    pub calls: BaseCalls,
    /// Locus evidence, integrity, and noisy regions.
    pub signal: SignalAnalysis,
    /// Phase-state segments, mask, and callable span.
    pub callability: ReadCallability,
    /// Relative quality and the trim interval.
    pub quality: QualityControlResult,
}

impl SangerAttachment {
    /// Locus evidence of call `index`, when its record carries that index.
    #[must_use]
    pub fn locus(&self, index: usize) -> Option<&LocusEvidence> {
        self.signal
            .loci
            .get(index)
            .filter(|locus| locus.call_index_0based == index)
    }

    /// Whether call `index` lies in a merged candidate-noisy region.
    #[must_use]
    pub fn in_noisy_region(&self, index: usize) -> bool {
        self.signal.noisy_regions.iter().any(|region| {
            region.call_start_0based <= index && index < region.call_end_0based_exclusive
        })
    }
}

/// Sanger evidence of a read with too few callable calls, which never reaches
/// the core.
#[derive(Debug, Clone)]
pub struct SangerRejection {
    /// Reviewer-facing name of the rejected trace.
    pub input_name: String,
    /// SHA-256 content identity of the rejected trace.
    pub input_sha256: String,
    /// Trace-integrity evidence of the rejected read.
    pub integrity: SangerIntegrity,
    /// Callability of the rejected read.
    pub callability: ReadCallability,
    /// Why the read was rejected.
    pub rejection: ReadRejection,
}
