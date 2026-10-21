//! `dna.run_context/v1` construction (PROP-0014).

use dna_kernel::model::reference::Reference;
use dna_kernel::plugin::PluginDescriptor;
use dna_kernel::profile::ProfileIdentity;

use crate::input::run_context::SCHEMA_VERSION;
use crate::model::result::ReferenceResult;
use crate::model::run_context::{LENGTH_METHOD, RunContext};
use crate::model::run_context_result::{
    LengthRunStatisticsResult, LengthSectionResult, RunContextProvenanceResult, RunContextResult,
    SkippedResult,
};

use super::json::{project_plugins, project_profile};

/// A completed run context with what it was made with.
pub(crate) struct CompletedRunContext {
    pub(crate) context: RunContext,
    /// Samples that contributed reads.
    pub(crate) samples: usize,
    /// Listed input files left out, in list order.
    pub(crate) skipped: Vec<Skipped>,
    pub(crate) reference: Reference,
    pub(crate) profile: ProfileIdentity,
    pub(crate) configuration_sha256: String,
    /// Plugins of the workflow, in execution order.
    pub(crate) plugins: &'static [&'static PluginDescriptor],
}

/// A listed input file the command left out, and why.
pub(crate) struct Skipped {
    pub(crate) sample_id: String,
    /// The file name.
    pub(crate) name: String,
    /// `missing`, `not_a_file`, `empty`, or `unreadable`.
    pub(crate) reason: &'static str,
}

/// Builds `dna.run_context/v1` without filesystem side effects.
pub(crate) fn build(completed: CompletedRunContext) -> RunContextResult {
    let CompletedRunContext {
        context,
        samples,
        skipped,
        reference,
        profile,
        configuration_sha256,
        plugins,
    } = completed;
    let length = context.length;
    RunContextResult {
        schema_version: SCHEMA_VERSION,
        run_id: context.run_id,
        provenance: RunContextProvenanceResult {
            reference: ReferenceResult {
                name: reference.name,
                topology: reference.topology,
                sha256: reference.sequence_sha256,
            },
            configuration_sha256,
            profile: project_profile(profile),
            plugins: project_plugins(plugins),
            skipped: skipped
                .into_iter()
                .map(|skipped| SkippedResult {
                    sample_id: skipped.sample_id,
                    name: skipped.name,
                    reason: skipped.reason,
                })
                .collect(),
        },
        reads: context.reads.into_iter().collect(),
        rejected_reads: context.rejected_reads,
        samples,
        length: LengthSectionResult {
            method: LENGTH_METHOD,
            windows: context.windows,
            noise_forward: length.forward,
            noise_reverse: length.reverse,
            edges: length.edges,
            runs: length
                .runs
                .into_iter()
                .enumerate()
                .map(|(index, run)| LengthRunStatisticsResult {
                    run: index + 1,
                    pairs: run.pairs,
                    sum_z2: run.sum_z2,
                    sum_z4: run.sum_z4,
                    histogram: run.histogram,
                })
                .collect(),
        },
    }
}
