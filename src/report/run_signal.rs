//! `dna.signal/v1` projection of a run signal (PROP-0009).

use dna_kernel::plugin::PluginDescriptor;
use dna_sanger::run_signal::RunSignal;

use crate::model::run_signal_result::{
    RunSignalKeyResult, RunSignalProvenanceResult, RunSignalResult,
};
use crate::report::json::project_plugins;

/// A built run signal and the inputs it was made with.
pub(crate) struct CompletedRunSignal {
    pub(crate) signal: RunSignal,
    pub(crate) rejected_reads: usize,
    pub(crate) configuration_sha256: String,
    pub(crate) plugins: &'static [&'static PluginDescriptor],
}

/// Builds `dna.signal/v1` without filesystem side effects.
pub(crate) fn build(completed: CompletedRunSignal) -> RunSignalResult {
    let CompletedRunSignal {
        signal,
        rejected_reads,
        configuration_sha256,
        plugins,
    } = completed;
    RunSignalResult {
        schema_version: "dna.signal/v1",
        run_id: signal.run_id,
        provenance: RunSignalProvenanceResult {
            configuration_sha256,
            plugins: project_plugins(plugins),
            selection: signal.selection,
            flank: signal.flank,
            resolution: signal.resolution,
        },
        traces: signal.traces.into_iter().collect(),
        reads: signal.reads,
        rejected_reads,
        calls: signal.calls,
        keys: signal
            .keys
            .into_iter()
            .map(|(key, shares)| RunSignalKeyResult {
                context: key.context,
                base: key.base,
                second: key.second,
                shares: shares
                    .into_iter()
                    .map(|(share, calls)| [u64::from(share), calls])
                    .collect(),
            })
            .collect(),
    }
}
