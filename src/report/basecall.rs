//! Typed assembly of the reference-free basecall result.

use crate::config::Config;
use crate::error::{ReportError, Result};
use crate::model::basecall_result::{BasecallProvenanceResult, BasecallReadResult, BasecallResult};
use crate::model::result::{InputResult, IntervalResult};
use crate::model::sanger_read;
use crate::report::json::project_plugins;
use dna_kernel::plugin::PluginDescriptor;
use dna_sanger::model::basecalls::BaseCalls;
use dna_sanger::model::callability::ReadCallability;
use dna_sanger::model::quality::QualityControlResult;
use dna_sanger::model::sanger::Chromatogram;
use dna_sanger::model::signal::SignalAnalysis;

/// Inputs consumed to build one immutable basecall document.
pub(crate) struct CompletedBasecall {
    pub(crate) config: Config,
    pub(crate) trace: Chromatogram,
    pub(crate) calls: BaseCalls,
    pub(crate) signal: SignalAnalysis,
    pub(crate) callability: ReadCallability,
    pub(crate) quality: QualityControlResult,
    /// Plugins of the workflow, in execution order.
    pub(crate) plugins: &'static [&'static PluginDescriptor],
}

/// Builds `dna.basecalls/v3` without filesystem side effects.
pub(crate) fn build(completed: CompletedBasecall) -> Result<BasecallResult> {
    let CompletedBasecall {
        config,
        trace,
        calls,
        signal: signal_analysis,
        callability: read_callability,
        quality,
        plugins,
    } = completed;
    let call_count = calls.len();
    let ambiguity = calls
        .calls
        .iter()
        .map(|call| call.ambiguity)
        .collect::<String>();
    if calls.primary_sequence.len() != call_count || ambiguity.len() != call_count {
        return Err(
            ReportError::Inconsistent("basecall sequence lengths do not match call count").into(),
        );
    }
    if quality.trim_start_0based > quality.trim_end_0based_exclusive
        || quality.trim_end_0based_exclusive > call_count
        || quality.retained_sequence
            != calls.primary_sequence[quality.trim_start_0based..quality.trim_end_0based_exclusive]
    {
        return Err(ReportError::Inconsistent(
            "basecall retained sequence does not match trim bounds",
        )
        .into());
    }
    let warnings = sanger_read::warnings(&calls, &signal_analysis.integrity);
    let signal_quality = sanger_read::signal_quality(&signal_analysis);

    Ok(BasecallResult {
        schema_version: "dna.basecalls/v3",
        provenance: BasecallProvenanceResult {
            input: InputResult {
                sha256: trace.source_sha256,
            },
            configuration_sha256: config.source_sha256,
            plugins: project_plugins(plugins),
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
            callability: sanger_read::callability(&read_callability),
        },
        signal_quality,
        warnings,
    })
}
