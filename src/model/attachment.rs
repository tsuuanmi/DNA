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
pub(crate) struct SangerAttachment {
    pub(crate) calls: BaseCalls,
    pub(crate) signal: SignalAnalysis,
    pub(crate) callability: ReadCallability,
    pub(crate) quality: QualityControlResult,
}

impl SangerAttachment {
    /// Locus evidence of call `index`, when its record carries that index.
    pub(crate) fn locus(&self, index: usize) -> Option<&LocusEvidence> {
        self.signal
            .loci
            .get(index)
            .filter(|locus| locus.call_index_0based == index)
    }

    /// Whether call `index` lies in a merged candidate-noisy region.
    pub(crate) fn in_noisy_region(&self, index: usize) -> bool {
        self.signal.noisy_regions.iter().any(|region| {
            region.call_start_0based <= index && index < region.call_end_0based_exclusive
        })
    }
}

/// Sanger evidence of a read with too few callable calls, which never reaches
/// the core.
#[derive(Debug, Clone)]
pub(crate) struct SangerRejection {
    pub(crate) input_name: String,
    pub(crate) input_sha256: String,
    pub(crate) integrity: SangerIntegrity,
    pub(crate) callability: ReadCallability,
    pub(crate) rejection: ReadRejection,
}
