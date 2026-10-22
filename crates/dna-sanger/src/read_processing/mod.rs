//! Shared reference-free read processing for reusable scientific capabilities.

use dna_kernel::plugin::{Contract, PluginDescriptor, PluginFamily};
mod config;
pub(crate) mod evidence;
mod stage_log;

pub use config::{RawSangerConfig, RawSangerEvidenceConfig, SangerConfig, SangerEvidenceConfig};
pub use evidence::{read_evidence, reason};

use std::time::Instant;

use crate::basecalling;
use crate::callability;
use crate::error::{CallabilityError, Result};
use crate::model::attachment::SangerRejection;
use crate::model::basecalls::BaseCalls;
use crate::model::callability::{PhaseState, ReadCallability, ReadRejection};
use crate::model::quality::QualityControlResult;
use crate::model::sanger::Chromatogram;
use crate::model::second_cause::SecondPeak;
use crate::model::signal::{ShortRead, SignalAnalysis};
use crate::quality_control;
use crate::second_cause;
use crate::signal_processing;
use dna_kernel::read_evidence::{PhaseMixture, TractLengths};
use dna_trace::decompose::TraceCall;

/// Warning counts shared by command-level operational summaries.
pub struct ReadWarnings {
    /// Calls without a resolved primary base.
    pub unresolved_primary_calls: usize,
    /// Calls with more than two qualifying channels.
    pub multi_channel_unresolved_calls: usize,
    /// Calls whose primary base differs from the vendor base.
    pub vendor_disagreements: usize,
    /// Vendor series whose length differs from the locus count.
    pub locus_vendor_length_mismatches: usize,
    /// Channel samples at the signed 16-bit limits.
    pub clipped_channel_samples: usize,
}

impl ReadWarnings {
    /// The warning counts of a read's calls and signal.
    fn of(calls: &BaseCalls, signal: &SignalAnalysis) -> Self {
        Self {
            unresolved_primary_calls: calls
                .calls
                .iter()
                .filter(|call| call.primary == 'N')
                .count(),
            multi_channel_unresolved_calls: calls
                .calls
                .iter()
                .filter(|call| call.qualifying_channels.len() > 2 && call.ambiguity == 'N')
                .count(),
            vendor_disagreements: calls
                .calls
                .iter()
                .filter(|call| call.vendor_agrees == Some(false))
                .count(),
            locus_vendor_length_mismatches: signal.integrity.vendor_length_mismatch_count(),
            clipped_channel_samples: signal.integrity.clipped_channel_samples,
        }
    }
}

/// Scientific read products shared by reference-free and reference-guided paths.
pub struct ProcessedRead {
    /// Signal-derived base calls.
    pub calls: BaseCalls,
    /// Signal evidence of the read.
    pub signal: SignalAnalysis,
    /// Callability of the read.
    pub callability: ReadCallability,
    /// Relative quality and trim interval.
    pub quality: QualityControlResult,
    /// Each call's second peak and its cause (SRS-CALL-021), by call index.
    pub second: Vec<Option<SecondPeak>>,
    /// The read's phase mixtures after repeat runs (SRS-CALL-023).
    pub mixtures: Vec<PhaseMixture>,
    /// The read's tracts and their length decodes (PROP-0012 increment 2).
    pub tracts: Vec<TractLengths>,
    /// Warning counts for operational summaries.
    pub warnings: ReadWarnings,
}

/// A read after the stages that decide whether it can be analyzed at all.
pub struct PreparedRead {
    pub(crate) calls: BaseCalls,
    /// Signal evidence of the read.
    pub signal: SignalAnalysis,
    /// Callability of the read.
    pub callability: ReadCallability,
    pub(crate) warnings: ReadWarnings,
}

impl PreparedRead {
    /// Why the read cannot be analyzed, or `None` when it has enough callable
    /// calls.
    #[must_use]
    pub fn rejection(&self, config: &SangerConfig) -> Option<ReadRejection> {
        let callable_calls = self.callability.callable_count();
        let minimum_callable_calls = config.callability.minimum_callable_calls;
        (callable_calls < minimum_callable_calls).then_some(ReadRejection {
            callable_calls,
            minimum_callable_calls,
        })
    }
}

/// Runs the scientific stages that require no reference, emitting one
/// `tracing` stage span and completion event per stage; a read shorter than
/// one signal window, or with too few callable calls, fails typed.
///
/// # Errors
///
/// Returns the failure of the first stage that fails, including
/// `SignalError::TooFewCalls` and `CallabilityError::TooFewCallableCalls`.
pub fn process(trace: &Chromatogram, config: &SangerConfig) -> Result<ProcessedRead> {
    let prepared = prepare(trace, config)?;
    finish(trace, prepared, config)
}

/// What the Sanger stages make of one trace in a multi-read operation: its
/// products, or why it is set aside, with the evidence its stages produced.
pub enum Admission {
    /// The read can be analyzed.
    Processed(Box<ProcessedRead>),
    /// The read has fewer calls than one signal window (SRS-SIG-007); no
    /// stage ran on it.
    TooShort(ShortRead),
    /// Too few of the read's calls are callable.
    TooFewCallable {
        /// Integrity and callability of the read.
        sanger: Box<SangerRejection>,
        /// Callable and minimum call counts.
        rejection: ReadRejection,
    },
}

/// Runs the reference-free stages of one trace as a multi-read operation
/// does: a read too short for one signal window, or with too few callable
/// calls, is an [`Admission`] that sets it aside rather than an error.
///
/// # Errors
///
/// Returns the failure of a stage other than those two rejections.
pub fn admit(trace: &Chromatogram, config: &SangerConfig) -> Result<Admission> {
    let prepared = match prepare(trace, config) {
        Ok(prepared) => prepared,
        Err(error) => {
            return ShortRead::of(&error).map(Admission::TooShort).ok_or(error);
        }
    };
    if let Some(rejection) = prepared.rejection(config) {
        return Ok(Admission::TooFewCallable {
            sanger: Box::new(SangerRejection {
                integrity: prepared.signal.integrity,
                callability: prepared.callability,
            }),
            rejection,
        });
    }
    Ok(Admission::Processed(Box::new(finish(
        trace, prepared, config,
    )?)))
}

/// Runs basecalling, signal processing, and callability; a read shorter than
/// one signal window fails typed before basecalling, and
/// [`ShortRead::of`](crate::model::signal::ShortRead::of) reads that failure.
///
/// # Errors
///
/// Returns `SignalError::TooFewCalls` for a read shorter than one signal
/// window, or the failure of basecalling, signal processing, or callability.
pub fn prepare(trace: &Chromatogram, config: &SangerConfig) -> Result<PreparedRead> {
    if let Err(error) =
        signal_processing::check_length(trace.call_count(), &config.signal_processing)
    {
        // The rule belongs to signal processing, which the failure names.
        let _stage = tracing::info_span!("signal_processing").entered();
        return Err(error);
    }
    let stage = tracing::info_span!("basecalling").entered();
    let stage_started = Instant::now();
    let calls = basecalling::call(trace, &config.basecalling)?;
    stage_log::basecalling(&calls, config, stage_started);

    drop(stage);
    let stage = tracing::info_span!("signal_processing").entered();
    let stage_started = Instant::now();
    let signal = signal_processing::analyze(trace, &calls, &config.signal_processing)?;
    stage_log::signal_processing(&signal, config, stage_started);

    drop(stage);
    let stage = tracing::info_span!("callability").entered();
    let stage_started = Instant::now();
    let callability = callability::analyze(trace, &calls, &signal, &config.callability)?;
    stage_log::callability(&calls, &callability, config, stage_started);

    drop(stage);
    let warnings = ReadWarnings::of(&calls, &signal);
    Ok(PreparedRead {
        calls,
        signal,
        callability,
        warnings,
    })
}

/// Runs quality control on a prepared read; a read with too few callable
/// calls fails typed.
///
/// # Errors
///
/// Returns `CallabilityError::TooFewCallableCalls` for a read with too few
/// callable calls, or the quality-control failure.
pub fn finish(
    trace: &Chromatogram,
    prepared: PreparedRead,
    config: &SangerConfig,
) -> Result<ProcessedRead> {
    if let Some(rejection) = prepared.rejection(config) {
        return Err(CallabilityError::TooFewCallableCalls {
            callable: rejection.callable_calls,
            minimum: rejection.minimum_callable_calls,
        }
        .into());
    }
    let PreparedRead {
        calls,
        signal,
        callability,
        warnings,
    } = prepared;
    let _stage = tracing::info_span!("quality_control").entered();
    let stage_started = Instant::now();
    let quality = quality_control::analyze(trace, &calls, &callability, &config.quality_control)?;
    stage_log::quality_control(&calls, &quality, stage_started);

    let second = second_cause::classify(&calls, &signal, &callability, &config.signal_processing);
    let mixtures = phase_mixtures(&calls, &signal, &callability, &quality, config);
    let tracts = crate::tract::tracts(trace, &calls, &callability, &config.trace);
    Ok(ProcessedRead {
        calls,
        signal,
        callability,
        quality,
        second,
        mixtures,
        tracts,
        warnings,
    })
}

/// Sanger ABIF modality: basecalling, signal processing, callability, quality
/// control, and the evidence adapter.
pub const PLUGIN: PluginDescriptor = PluginDescriptor {
    id: "sanger",
    family: PluginFamily::Modality,
    version: 23,
    provides: &[Contract::ReadEvidence],
    requires: &[],
    config_sections: &[
        "basecalling",
        "signal_processing",
        "callability",
        "quality_control",
        "sanger_evidence",
        "trace",
    ],
};

/// The read's phase mixtures (SRS-CALL-023), decomposed by `dna-trace`
/// (ADR-0084) from the read's calls, their normalised locus profiles, the
/// callability mask, its dephased segments, and the trim start.
fn phase_mixtures(
    calls: &BaseCalls,
    signal: &SignalAnalysis,
    callability: &ReadCallability,
    quality: &QualityControlResult,
    config: &SangerConfig,
) -> Vec<PhaseMixture> {
    let trace = calls
        .calls
        .iter()
        .enumerate()
        .map(|(index, call)| TraceCall {
            primary: call.primary,
            single: call.qualifying_channels.len() == 1,
            profile: signal
                .loci
                .get(index)
                .filter(|locus| locus.call_index_0based == index)
                .and_then(|locus| locus.profile)
                .map(|profile| profile.weights),
            in_phase: callability.mask.get(index).is_some_and(Option::is_none),
        })
        .collect::<Vec<_>>();
    let dephased = callability
        .segments
        .iter()
        .filter(|segment| segment.state == PhaseState::Dephased)
        .map(|segment| segment.call_start_0based..segment.call_end_0based_exclusive)
        .collect::<Vec<_>>();
    dna_trace::decompose::mixtures(
        &trace,
        &dephased,
        quality.trim_start_0based,
        config.signal_processing.decompose_unexplained_quantile,
        config.signal_processing.decompose_minimum_calls,
        config.signal_processing.decompose_onset_window,
    )
    .unwrap_or_default()
}
