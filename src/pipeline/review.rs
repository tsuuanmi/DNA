//! Sample review orchestration and publication (PROP-0004).

use std::collections::BTreeMap;
use std::path::Path;
use std::time::Instant;

use crate::error::Result;
use crate::input::signal as input_signal;
use crate::model::review_result::{
    ComparedResult, ComponentResult, ItemKind, LabelResult, MixtureResult,
};
use crate::operation_log::OperationLog;
use crate::pipeline::path;
use crate::pipeline::plugins as plugin;
use crate::pipeline::represent;
use crate::report::{self, CompletedReview, Reviewed};
use crate::sample_analysis::{self, ConsensusCall, Outcome, ReviewCheck, SampleAnalysis};
use dna_call::components::ComponentCall;
use dna_kernel::variant::Variant;
use dna_sanger::run_signal::RunSignal;

use super::sample::validate_sample_id;
use super::{Operation, SampleInput};

/// Runs one sample review operation with one sample-level append-only log.
pub(crate) fn run(
    sample_id: &str,
    input: SampleInput<'_>,
    reference: &Path,
    signal: Option<&Path>,
    config_path: &Path,
) -> Result<()> {
    validate_sample_id(sample_id)?;
    Operation::begin(sample_id, "review_failed", || {
        tracing::info!(
            event = "review_started",
            version = env!("CARGO_PKG_VERSION"),
            sample_id = ?sample_id,
            traces = input.len(),
            evidence = matches!(input, SampleInput::Evidence(_)),
            reference_path = ?reference.display().to_string(),
        );
    })?
    .run(|log, started| {
        review(
            sample_id,
            input,
            reference,
            signal,
            config_path,
            log,
            started,
        )
    })
}

fn review(
    sample_id: &str,
    input: SampleInput<'_>,
    reference: &Path,
    signal: Option<&Path>,
    config_path: &Path,
    log: &OperationLog,
    started: Instant,
) -> Result<()> {
    let stage = tracing::info_span!("input_loading").entered();
    let super::LoadedSample {
        loaded: inputs,
        reads,
    } = super::load_sample(input, reference, config_path)?;
    let signal = signal.map(input_signal::load).transpose()?;
    if let Some((signal, sha256)) = &signal {
        tracing::info!(
            event = "review_signal_loaded",
            run_id = ?signal.run_id,
            sha256 = %sha256,
            reads = signal.reads,
            keys = signal.keys.len(),
        );
    }
    let output = path::review_output(sample_id)?;
    tracing::info!(
        event = "review_inputs_loaded",
        traces = reads.len(),
        reference_sha256 = %inputs.reference.sequence_sha256,
        config_sha256 = %inputs.config.source_sha256,
        profile_id = ?inputs.profile.identity().id,
        output_path = ?output.display().to_string(),
    );
    let reference = inputs.reference.clone();
    let profile = inputs.profile.clone();
    let configuration_sha256 = inputs.config.source_sha256.clone();

    drop(stage);
    let outcome = sample_analysis::analyze(reads.reads(), inputs, log, super::read_workers())?;

    let stage = tracing::info_span!("review").entered();
    let stage_started = Instant::now();
    let mut read_mixtures = BTreeMap::new();
    let outcome = match outcome {
        Outcome::Rejected { rejections, .. } => Reviewed::Rejected(rejections),
        Outcome::Analyzed(analysis) => {
            read_mixtures = mixtures(&analysis)?;
            let call = analysis.call(sample_id);
            if let Err(error) = &call {
                tracing::warn!(event = "review_call_failed", error = ?error.to_string());
            }
            Reviewed::Analyzed { analysis, call }
        }
    };
    let plugins = plugin::review(profile.notation.is_some(), profile.empop);
    let result = report::build_review(CompletedReview {
        sample_id: sample_id.to_owned(),
        reference,
        profile,
        configuration_sha256,
        outcome,
        plugins,
        signal,
        mixtures: read_mixtures,
    })?;
    let judged = result
        .items
        .iter()
        .flat_map(|item| &item.labels)
        .filter_map(|label| match label {
            LabelResult::PossibleHeteroplasmy { reads, .. } => Some(reads),
            _ => None,
        })
        .flatten()
        .map(|read| read.calls.is_some())
        .collect::<Vec<_>>();
    tracing::info!(
        event = "review_completed",
        elapsed_ms = stage_started.elapsed().as_millis(),
        status = ?result.status,
        reads = result.reads.len(),
        items = result.items.len(),
        checks = result
            .items
            .iter()
            .filter(|item| item.kind == ItemKind::Check)
            .count(),
        clear = result.summary.clear,
        heteroplasmy_reads_run_level = judged.iter().filter(|&&level| level).count(),
        heteroplasmy_reads_floor = judged.iter().filter(|&&level| !level).count(),
    );

    drop(stage);
    let _stage = tracing::info_span!("reporting").entered();
    let bytes = report::serialize(&result)?;
    tracing::info!(
        event = "review_ready_for_publication",
        total_elapsed_ms = started.elapsed().as_millis(),
        schema = result.schema_version,
        output_path = ?output.display().to_string(),
        bytes = bytes.len(),
    );
    log.sync()?;
    report::publish(&output, &bytes)
}

/// Each admitted read's phase mixtures (ADR-0082), by content SHA-256: every
/// component placed and called as a read of its own, its eligible edits
/// written in the profile notation over the reference interval it covers,
/// and the reference positions of the calls the components do not explain.
fn mixtures(analysis: &SampleAnalysis) -> Result<BTreeMap<String, Vec<MixtureResult>>> {
    let mut found = BTreeMap::new();
    for read in analysis.components() {
        let positions = analysis
            .called
            .iter()
            .find(|called| called.input_sha256 == read.sha256)
            .map(|called| {
                called
                    .alignment
                    .columns()
                    .iter()
                    .filter_map(|column| {
                        Some((
                            column.original_call_index_0based()?,
                            column.reference_index_0based()? + 1,
                        ))
                    })
                    .collect::<BTreeMap<_, _>>()
            })
            .unwrap_or_default();
        let mixtures = read
            .mixtures
            .iter()
            .map(|mixture| {
                // Every component the modality asserts, by offset and
                // fraction; an explained mixture's placed ones also compared
                // with the reference.
                let components = mixture
                    .mixture
                    .components
                    .iter()
                    .map(|component| {
                        let compared = mixture
                            .components
                            .iter()
                            .find(|placed| placed.offset == component.offset)
                            .filter(|_| mixture.mixture.explained)
                            .map(|placed| {
                                Ok::<_, crate::error::Error>(ComparedResult {
                                    first: placed.span[0],
                                    last: placed.span[1],
                                    sequence: placed.sequence.clone(),
                                    calls: notation_texts(analysis, placed)?,
                                })
                            })
                            .transpose()?;
                        Ok(ComponentResult {
                            offset: component.offset,
                            fraction: round(component.fraction),
                            compared,
                        })
                    })
                    .collect::<Result<Vec<_>>>()?;
                Ok(MixtureResult {
                    explained: mixture.mixture.explained,
                    components,
                    unexplained: mixture
                        .mixture
                        .unexplained
                        .iter()
                        .filter_map(|call| positions.get(call).copied())
                        .collect::<std::collections::BTreeSet<_>>()
                        .into_iter()
                        .collect(),
                })
            })
            .collect::<Result<Vec<_>>>()?;
        found.insert(read.sha256, mixtures);
    }
    Ok(found)
}

/// A fraction to three decimals, as every DNA document writes it.
fn round(value: f64) -> f64 {
    (value * 1000.0).round() / 1000.0
}

/// A component's eligible edits in the profile notation, as the read
/// representation writes a read's (ADR-0060); none under a profile without
/// notation.
fn notation_texts(analysis: &SampleAnalysis, component: &ComponentCall) -> Result<Vec<String>> {
    let loaded = &analysis.loaded;
    let variants = component
        .variants
        .iter()
        .filter(|observed| observed.exclusion_reasons.is_empty())
        .map(|observed| Variant::from(&observed.variant))
        .collect();
    // A span whose first position follows its last crosses the origin of a
    // circular reference.
    let [first, last] = component.span;
    let mut covered = Vec::with_capacity(2);
    if first <= last {
        covered.push(first - 1..last);
    } else {
        covered.extend([first - 1..loaded.reference.sequence.len(), 0..last]);
    }
    let Some(notation) = represent::represent_variants(
        std::iter::once((String::new(), variants, covered)),
        &loaded.reference,
        &loaded.profile,
    )?
    else {
        return Ok(Vec::new());
    };
    Ok(report::notation_calls(&loaded.reference, &notation)?
        .into_keys()
        .map(|call| call.text)
        .collect())
}

impl SampleAnalysis {
    /// The items of this analysis's review list that a reviewer must look at
    /// (`check` items), given the call of its consensus and the run signal
    /// it is judged against, if any, by the rules `dna review` applies,
    /// without building or publishing the document.
    ///
    /// # Errors
    ///
    /// Returns [`ReportError`](crate::error::ReportError) when the analysis's
    /// records disagree, as `dna review` would.
    pub fn review_checks(
        &self,
        call: &ConsensusCall,
        signal: Option<&RunSignal>,
    ) -> Result<Vec<ReviewCheck>> {
        Ok(report::review_checks(self, call, signal)?
            .into_iter()
            .map(|(first, last, labels)| ReviewCheck {
                first,
                last,
                labels,
            })
            .collect())
    }
}
