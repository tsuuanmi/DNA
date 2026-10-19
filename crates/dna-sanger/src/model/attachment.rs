//! Sanger evidence kept beside the core's records, which reports join to them
//! by read identity and call index (ADR-0069).

use crate::model::basecalls::BaseCalls;
use crate::model::callability::ReadCallability;
use crate::model::locus_evidence::LocusEvidence;
use crate::model::quality::QualityControlResult;
use crate::model::second_cause::SecondPeak;
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
    /// Each call's second peak and its cause (SRS-CALL-021), by call index.
    pub second: Vec<Option<SecondPeak>>,
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

/// Sanger evidence of a rejected read that reached callability: one with too
/// few callable calls, which never reaches the core, or one the core could
/// not place. A read too short for signal processing has none.
#[derive(Debug, Clone)]
pub struct SangerRejection {
    /// Trace-integrity evidence of the rejected read.
    pub integrity: SangerIntegrity,
    /// Callability of the rejected read.
    pub callability: ReadCallability,
}
