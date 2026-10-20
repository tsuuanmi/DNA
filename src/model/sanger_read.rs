//! The Sanger evidence of one read that the analysis document reports beside
//! the core's products, as `dna.read_evidence/v1` carries it (PROP-0011).

use crate::model::basecall_result::BasecallWarningSummaryResult;
use crate::model::result::{CallabilityResult, IntervalResult, SignalQualityResult};

/// One read's Sanger evidence as the documents report it: projected from the
/// Sanger products of a trace, or read back from its read evidence document.
#[derive(Debug)]
pub(crate) struct SangerRead {
    /// The relative-quality trim interval.
    pub(crate) trim: IntervalResult,
    /// The read callability view.
    pub(crate) callability: CallabilityResult,
    /// Trace integrity and merged noisy regions.
    pub(crate) signal_quality: SignalQualityResult,
    /// The modality's warning counts.
    pub(crate) warnings: BasecallWarningSummaryResult,
    /// Each call's signal, in call order.
    pub(crate) calls: Vec<SangerCallSignal>,
}

/// The signal of one call a variant record joins.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct SangerCallSignal {
    /// A/C/G/T channel heights at the primary event, on the read strand;
    /// `None` when the call has no primary-event peak evidence.
    pub(crate) peaks: Option<[i32; 4]>,
    /// The call's relative quality.
    pub(crate) quality: u8,
}
