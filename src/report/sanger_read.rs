//! Projection of one read's Sanger products into the Sanger evidence the
//! analysis and read evidence documents report (PROP-0011).

use crate::error::{ReportError, Result};
use crate::model::basecall_result::BasecallWarningSummaryResult;
use crate::model::result::IntervalResult;
use crate::model::sanger_read::SangerRead;
use crate::report::{callability, sanger_call, signal};
use dna_sanger::model::basecalls::BaseCalls;
use dna_sanger::model::callability::ReadCallability;
use dna_sanger::model::quality::QualityControlResult;
use dna_sanger::model::signal::{SangerIntegrity, SignalAnalysis};

/// Projects the trim, callability, signal quality, warning counts, and every
/// call's signal of one read.
pub(crate) fn project(
    calls: &BaseCalls,
    signal: SignalAnalysis,
    read_callability: &ReadCallability,
    quality: &QualityControlResult,
) -> Result<SangerRead> {
    if quality.per_call.len() != calls.len() {
        return Err(
            ReportError::Inconsistent("Sanger calls and quality differ in call count").into(),
        );
    }
    let warnings = warnings(calls, &signal.integrity);
    let call_signals = (0..calls.len())
        .map(|index| sanger_call::signal(calls, quality, index).map_err(ReportError::CallEvidence))
        .collect::<std::result::Result<Vec<_>, _>>()?;
    Ok(SangerRead {
        trim: IntervalResult {
            start: quality.trim_start_0based,
            end: quality.trim_end_0based_exclusive,
        },
        callability: callability::project(read_callability),
        signal_quality: signal::project(signal),
        warnings,
        calls: call_signals,
    })
}

/// The modality's warning counts of one read.
pub(super) fn warnings(
    calls: &BaseCalls,
    integrity: &SangerIntegrity,
) -> BasecallWarningSummaryResult {
    BasecallWarningSummaryResult {
        unresolved_primary_calls: calls
            .calls
            .iter()
            .filter(|call| call.primary == 'N')
            .count(),
        multi_channel_unresolved_calls: calls
            .calls
            .iter()
            .filter(|call| call.ambiguity == 'N' && call.qualifying_channels.len() > 2)
            .count(),
        vendor_disagreements: calls
            .calls
            .iter()
            .filter(|call| call.vendor_agrees == Some(false))
            .count(),
        ploc_vendor_length_mismatches: integrity.vendor_length_mismatch_count(),
        clipped_channel_samples: integrity.clipped_channel_samples,
    }
}
