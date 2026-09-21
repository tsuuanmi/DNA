//! Typed assembly of the reference-free basecall result.

use crate::config::Config;
use crate::error::{Error, Result};
use crate::model::basecall_result::{
    BasecallProvenanceResult, BasecallReadResult, BasecallResult, BasecallWarningSummaryResult,
};
use crate::model::basecalls::BaseCalls;
use crate::model::quality::QualityControlResult;
use crate::model::result::{InputResult, IntervalResult};
use crate::model::signal::DNAAnalysis;
use crate::model::trace::Chromatogram;
use crate::report::signal;

/// Inputs consumed to build one immutable basecall document.
pub(crate) struct CompletedBasecall {
    pub(crate) config: Config,
    pub(crate) trace: Chromatogram,
    pub(crate) calls: BaseCalls,
    pub(crate) signal: DNAAnalysis,
    pub(crate) quality: QualityControlResult,
}

/// Builds `dna.basecalls/v2` without filesystem side effects.
pub(crate) fn build(completed: CompletedBasecall) -> Result<BasecallResult> {
    let CompletedBasecall {
        config,
        trace,
        calls,
        signal: signal_analysis,
        quality,
    } = completed;
    let call_count = calls.len();
    let ambiguity = calls
        .calls
        .iter()
        .map(|call| call.ambiguity)
        .collect::<String>();
    if calls.primary_sequence.len() != call_count || ambiguity.len() != call_count {
        return Err(Error::Report(
            "basecall sequence lengths do not match call count".into(),
        ));
    }
    if quality.trim_start_0based > quality.trim_end_0based_exclusive
        || quality.trim_end_0based_exclusive > call_count
        || quality.retained_sequence
            != calls.primary_sequence[quality.trim_start_0based..quality.trim_end_0based_exclusive]
    {
        return Err(Error::Report(
            "basecall retained sequence does not match trim bounds".into(),
        ));
    }
    let unresolved_primary_calls = calls
        .calls
        .iter()
        .filter(|call| call.primary == 'N')
        .count();
    let multi_channel_unresolved_calls = calls
        .calls
        .iter()
        .filter(|call| call.ambiguity == 'N' && call.qualifying_channels.len() > 2)
        .count();
    let vendor_disagreements = calls
        .calls
        .iter()
        .filter(|call| call.vendor_agrees == Some(false))
        .count();
    let ploc_vendor_length_mismatches = signal_analysis.integrity.vendor_length_mismatch_count();
    let clipped_channel_samples = signal_analysis.integrity.clipped_channel_samples;
    let signal_quality = signal::project(signal_analysis);

    Ok(BasecallResult {
        schema_version: "dna.basecalls/v2",
        provenance: BasecallProvenanceResult {
            input: InputResult {
                sha256: trace.source_sha256,
            },
            configuration_sha256: config.source_sha256,
        },
        read: BasecallReadResult {
            call_count,
            primary: calls.primary_sequence,
            ambiguity,
            retained: quality.retained_sequence,
            trim: IntervalResult {
                start: quality.trim_start_0based,
                end: quality.trim_end_0based_exclusive,
            },
        },
        signal_quality,
        warnings: BasecallWarningSummaryResult {
            unresolved_primary_calls,
            multi_channel_unresolved_calls,
            vendor_disagreements,
            ploc_vendor_length_mismatches,
            clipped_channel_samples,
        },
    })
}
