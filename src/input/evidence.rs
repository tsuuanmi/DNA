//! Read evidence documents (`dna.read_evidence/v1`) as input (PROP-0011):
//! one read's evidence, read back into the kernel's `ReadEvidence` through
//! its own constructors, with the Sanger evidence the documents report
//! beside the core's products; or a read the Sanger plugin set aside, with
//! its cause.

use std::fs;
use std::path::{Path, PathBuf};

use serde::Deserialize;

use crate::config::Config;
use crate::error::{
    CallabilityError, KernelError, ReadEvidenceError, Result, SampleError, SangerError, SignalError,
};
use crate::model::basecall_result::BasecallWarningSummaryResult;
use crate::model::result::{
    CallabilityResult, CallabilitySegmentResult, IntervalResult, NoisyRegionResult,
    SignalQualityResult, TraceIntegrityResult,
};
use crate::model::sample_result::ReadRejectionResult;
use crate::model::sanger_read::{RejectedSanger, Rejection, SangerCallSignal, SangerRead};
use dna_kernel::model::reference::Reference;
use dna_kernel::plugin::PluginDescriptor;
use dna_kernel::profile::Profile;
use dna_kernel::read_evidence::{
    Bound, CallEvidence, EvidenceProfile, MixtureComponent, PhaseMixture, Placement, PostRepeatRun,
    ReadEvidence, Reasons, SecondBase, UnreliableSpan, Verdict,
};
use dna_kernel::read_evidence::{
    CleanedLengths, LeakageShare, LengthMarginal, Margin, StructureFit, TractLengths,
};
use dna_sanger::model::callability::PhaseState;
use dna_sanger::model::second_cause::{SecondCause, SecondPeak};
use dna_sanger::read_processing;

use super::{load_config, load_profile, load_reference, require_regular_file};

/// Largest accepted read evidence document.
const MAX_EVIDENCE_BYTES: u64 = 16 * 1024 * 1024;

/// The contract this reader accepts.
const SCHEMA_VERSION: &str = "dna.read_evidence/v1";

/// Scientific inputs for analyzing one recorded read.
pub(crate) struct EvidenceInputs {
    pub(crate) config: Config,
    pub(crate) profile: Profile,
    pub(crate) reference: Reference,
    pub(crate) read: RecordedRead,
}

/// Scientific inputs for one sample from its reads' recorded evidence.
pub(crate) struct SampleInputs {
    pub(crate) config: Config,
    pub(crate) profile: Profile,
    pub(crate) reference: Reference,
    /// The reads, in document order.
    pub(crate) reads: Vec<Recorded>,
}

/// One read as its read evidence document records it.
pub(crate) enum Recorded {
    /// A read the Sanger plugin hands to the core.
    Read(Box<RecordedRead>),
    /// A read the Sanger plugin set aside, with its cause.
    Rejected(Box<Rejection>),
}

/// One read the Sanger plugin hands to the core, as its document records it.
pub(crate) struct RecordedRead {
    /// Reviewer-facing name of the input: the file name of the trace the
    /// evidence was written from.
    pub(crate) input_name: String,
    /// SHA-256 of the trace the evidence was written from.
    pub(crate) input_sha256: String,
    /// The read's evidence as the core receives it.
    pub(crate) evidence: ReadEvidence,
    /// The read's Sanger evidence as the documents report it.
    pub(crate) sanger: SangerRead,
}

/// What a document must have been written under: the loaded
/// configuration, its rejection minimums, and this build's evidence plugins.
struct Expected<'a> {
    configuration_sha256: &'a str,
    minimum_calls: usize,
    minimum_callable_calls: usize,
    plugins: &'a [&'a PluginDescriptor],
}

impl<'a> Expected<'a> {
    fn of(config: &'a Config, plugins: &'a [&'a PluginDescriptor]) -> Self {
        Self {
            configuration_sha256: &config.source_sha256,
            minimum_calls: config.sanger.signal_processing.window_size_bases,
            minimum_callable_calls: config.sanger.minimum_callable_calls(),
            plugins,
        }
    }
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct RawDocument {
    schema_version: String,
    provenance: RawProvenance,
    read: Option<RawRead>,
    signal_quality: Option<RawSignalQuality>,
    warnings: Option<RawWarnings>,
    calls: Option<Vec<RawCall>>,
    rejection: Option<RawRejection>,
    integrity: Option<RawIntegrity>,
    callability: Option<RawCallability>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct RawProvenance {
    input: RawInput,
    configuration_sha256: String,
    plugins: Vec<RawPlugin>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct RawInput {
    name: String,
    sha256: String,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct RawPlugin {
    id: String,
    family: String,
    version: u32,
}

#[derive(Clone, Copy, Deserialize)]
#[serde(tag = "reason", rename_all = "snake_case", deny_unknown_fields)]
enum RawRejection {
    TooShort {
        calls: usize,
        minimum_calls: usize,
    },
    CallableCallsBelowMinimum {
        callable_calls: usize,
        minimum_callable_calls: usize,
    },
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct RawRead {
    call_count: usize,
    trim: RawInterval,
    informative: RawInterval,
    ends: RawEnds,
    reasons: Vec<String>,
    callability: RawCallability,
    mixtures: Vec<RawMixture>,
    unreliable: Vec<RawUnreliable>,
    tracts: Vec<RawTract>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct RawTract {
    onset: usize,
    calls: RawInterval,
    extent: usize,
    run: char,
    mark: Option<char>,
    candidates: Vec<RawStructureFit>,
    residual: f64,
    sigma: f64,
    autocorrelation: f64,
    first: RawMarginal,
    total: RawMarginal,
    downstream: String,
    cleaned: Option<RawCleaned>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct RawStructureFit {
    mark: Option<char>,
    residual: f64,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct RawMarginal {
    lengths: Vec<usize>,
    shares: Vec<f64>,
    errors: Vec<f64>,
    dominant: Option<usize>,
    next: Option<usize>,
    z: Option<f64>,
    z_adjusted: Option<f64>,
    margin: RawMargin,
}

#[derive(Clone, Copy, Deserialize)]
#[serde(rename_all = "snake_case")]
enum RawMargin {
    Clear,
    NearEqual,
    Undefined,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct RawCleaned {
    leakage: Option<Vec<RawLeakage>>,
    residual: f64,
    first: RawMarginal,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct RawLeakage {
    offset: i32,
    share: f64,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct RawUnreliable {
    calls: RawInterval,
    reason: String,
}

#[derive(Clone, Copy, Deserialize)]
#[serde(deny_unknown_fields)]
struct RawInterval {
    start: usize,
    end: usize,
}

#[derive(Clone, Copy, Deserialize)]
#[serde(rename_all = "snake_case")]
enum RawEnds {
    Unvouched,
    Vouched,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct RawCallability {
    callable_span: RawInterval,
    segments: Vec<RawSegment>,
    masked_calls: usize,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct RawSegment {
    calls: RawInterval,
    state: RawPhaseState,
    after_repeat: bool,
    shadow_offsets: Option<Vec<i8>>,
}

#[derive(Clone, Copy, Deserialize)]
#[serde(rename_all = "snake_case")]
enum RawPhaseState {
    InPhase,
    Dephased,
    Mixed,
    Weak,
    Irregular,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct RawMixture {
    calls: RawInterval,
    run_calls: usize,
    explained: bool,
    residual: f64,
    components: Vec<RawComponent>,
    unexplained: Vec<usize>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct RawComponent {
    offset: i32,
    fraction: f64,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct RawSignalQuality {
    integrity: RawIntegrity,
    noisy_regions: Vec<RawNoisyRegion>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct RawIntegrity {
    ploc_count: usize,
    duplicate_ploc_positions: usize,
    vendor_primary_count: Option<usize>,
    vendor_quality_count: Option<usize>,
    minimum_ploc_spacing: Option<usize>,
    median_ploc_spacing: Option<f64>,
    maximum_ploc_spacing: Option<usize>,
    clipped_channel_samples: usize,
    maximum_to_median_event_signal_ratio: Option<f64>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct RawNoisyRegion {
    calls: RawInterval,
    samples: RawInterval,
    minimum_primary_snr: f64,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct RawWarnings {
    unresolved_primary_calls: usize,
    multi_channel_unresolved_calls: usize,
    vendor_disagreements: usize,
    ploc_vendor_length_mismatches: usize,
    clipped_channel_samples: usize,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct RawCall {
    index: usize,
    base: char,
    profile: Option<[f64; 4]>,
    placement: RawPlacement,
    identity: RawVerdict,
    indel: RawVerdict,
    bound: RawBound,
    post_repeat: RawPostRepeat,
    second: Option<RawSecond>,
    second_unexplained: bool,
    reference_only: bool,
    signal: RawSignal,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct RawPlacement {
    kind: RawPlacementKind,
    reason: Option<String>,
}

#[derive(Clone, Copy, Deserialize)]
#[serde(rename_all = "snake_case")]
enum RawPlacementKind {
    Scored,
    Anchoring,
    Unresolved,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct RawVerdict {
    edge: bool,
    withheld: Vec<String>,
}

#[derive(Clone, Copy, Deserialize)]
#[serde(rename_all = "snake_case")]
enum RawBound {
    InPhase,
    Anchored,
    Withheld,
}

#[derive(Clone, Copy, Deserialize)]
#[serde(deny_unknown_fields)]
struct RawPostRepeat {
    call: bool,
    bound: bool,
}

#[derive(Clone, Copy, Deserialize)]
#[serde(deny_unknown_fields)]
struct RawSecond {
    base: char,
    share: f64,
    stands_out: bool,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct RawSignal {
    scan: usize,
    channels: RawPeaks,
    peaks: Option<RawPeaks>,
    relative_quality: u8,
    cause: Option<RawCause>,
}

#[derive(Clone, Copy, Deserialize)]
#[serde(deny_unknown_fields)]
struct RawPeaks {
    #[serde(rename = "A")]
    a: i32,
    #[serde(rename = "C")]
    c: i32,
    #[serde(rename = "G")]
    g: i32,
    #[serde(rename = "T")]
    t: i32,
}

impl RawPeaks {
    const fn channels(self) -> [i32; 4] {
        [self.a, self.c, self.g, self.t]
    }
}

/// The second-peak causes the modality gives a call from its shape.
#[derive(Clone, Copy, Deserialize)]
#[serde(rename_all = "snake_case")]
enum RawCause {
    Noise,
    NeighbourTail,
    OffCentre,
    DyeBlob,
    PhaseShift,
    SecondTemplate,
    MixedStretch,
}

/// Validates paths, loads configuration, profile, and reference, and reads
/// the read evidence document `document`, which `plugins` must have written
/// under the loaded configuration. A document of a read the Sanger plugin set
/// aside fails as the analysis of its trace does.
pub(crate) fn load_analysis(
    document: &Path,
    reference_path: &Path,
    config_path: &Path,
    plugins: &[&PluginDescriptor],
) -> Result<EvidenceInputs> {
    require_regular_file(document, "read evidence")?;
    require_regular_file(reference_path, "reference")?;
    let config = load_config(config_path)?;
    let profile = load_profile(&config)?;
    let read = match read_document(document, &Expected::of(&config, plugins))? {
        Recorded::Read(read) => *read,
        Recorded::Rejected(rejection) => return Err(set_aside(&rejection.cause)),
    };
    let reference = load_reference(reference_path, &profile)?;
    Ok(EvidenceInputs {
        config,
        profile,
        reference,
        read,
    })
}

/// Validates paths, loads configuration, profile, and reference, and reads
/// one sample's read evidence `documents`, each of which `plugins` must have
/// written under the loaded configuration.
pub(crate) fn load_sample(
    documents: &[PathBuf],
    reference_path: &Path,
    config_path: &Path,
    plugins: &[&PluginDescriptor],
) -> Result<SampleInputs> {
    if documents.is_empty() {
        return Err(SampleError::NoTraces.into());
    }
    for document in documents {
        require_regular_file(document, "read evidence")?;
    }
    require_regular_file(reference_path, "reference")?;
    let config = load_config(config_path)?;
    let profile = load_profile(&config)?;
    let reads = {
        let expected = Expected::of(&config, plugins);
        documents
            .iter()
            .map(|document| read_document(document, &expected))
            .collect::<Result<Vec<_>>>()?
    };
    let reference = load_reference(reference_path, &profile)?;
    Ok(SampleInputs {
        config,
        profile,
        reference,
        reads,
    })
}

/// Loads configuration, profile, and reference, and reads one sample's read
/// evidence `documents` from memory, each named for its errors, each of
/// which `plugins` must have written under the loaded configuration.
pub(crate) fn load_sample_documents(
    documents: &[(String, Vec<u8>)],
    reference_path: &Path,
    config_path: &Path,
    plugins: &[&PluginDescriptor],
) -> Result<SampleInputs> {
    if documents.is_empty() {
        return Err(SampleError::NoTraces.into());
    }
    require_regular_file(reference_path, "reference")?;
    let config = load_config(config_path)?;
    let profile = load_profile(&config)?;
    let reads = {
        let expected = Expected::of(&config, plugins);
        documents
            .iter()
            .map(|(name, bytes)| {
                let length = bytes.len() as u64;
                if length == 0 || length > MAX_EVIDENCE_BYTES {
                    return Err(ReadEvidenceError::FileSize {
                        bytes: length,
                        maximum: MAX_EVIDENCE_BYTES,
                    }
                    .into());
                }
                parse(bytes, Path::new(name), &expected)
            })
            .collect::<Result<Vec<_>>>()?
    };
    let reference = load_reference(reference_path, &profile)?;
    Ok(SampleInputs {
        config,
        profile,
        reference,
        reads,
    })
}

/// Reads one bounded document from disk into the read it records.
fn read_document(document: &Path, expected: &Expected<'_>) -> Result<Recorded> {
    let read_error = |source| KernelError::Read {
        kind: "read evidence",
        path: document.to_path_buf(),
        source,
    };
    let length = fs::metadata(document).map_err(read_error)?.len();
    if length == 0 || length > MAX_EVIDENCE_BYTES {
        return Err(ReadEvidenceError::FileSize {
            bytes: length,
            maximum: MAX_EVIDENCE_BYTES,
        }
        .into());
    }
    let bytes = fs::read(document).map_err(read_error)?;
    parse(&bytes, document, expected)
}

/// The failure the analysis of a set-aside read's trace returns.
fn set_aside(rejection: &ReadRejectionResult) -> crate::error::Error {
    match *rejection {
        ReadRejectionResult::TooShort {
            calls,
            minimum_calls,
        } => SangerError::from(SignalError::TooFewCalls {
            calls,
            window: minimum_calls,
        })
        .into(),
        ReadRejectionResult::CallableCallsBelowMinimum {
            callable_calls,
            minimum_callable_calls,
        } => SangerError::from(CallabilityError::TooFewCallableCalls {
            callable: callable_calls,
            minimum: minimum_callable_calls,
        })
        .into(),
        _ => inconsistent("a recorded rejection is not the Sanger plugin's"),
    }
}

/// Reads one document into the read it records.
fn parse(bytes: &[u8], path: &Path, expected: &Expected<'_>) -> Result<Recorded> {
    let raw: RawDocument =
        serde_json::from_slice(bytes).map_err(|error| ReadEvidenceError::Parse {
            path: path.to_path_buf(),
            source: Box::new(error),
        })?;
    if raw.schema_version != SCHEMA_VERSION {
        return Err(ReadEvidenceError::UnsupportedSchema {
            found: raw.schema_version,
        }
        .into());
    }
    check_provenance(&raw.provenance, expected)?;
    let RawDocument {
        provenance,
        read,
        signal_quality,
        warnings,
        calls,
        rejection,
        integrity,
        callability,
        ..
    } = raw;
    let RawInput { name, sha256 } = provenance.input;
    match (read, signal_quality, warnings, calls, rejection) {
        (Some(read), Some(signal_quality), Some(warnings), Some(calls), None)
            if integrity.is_none() && callability.is_none() =>
        {
            let (evidence, sanger) = recorded_read(&read, &signal_quality, &warnings, &calls)?;
            Ok(Recorded::Read(Box::new(RecordedRead {
                input_name: name,
                input_sha256: sha256,
                evidence,
                sanger,
            })))
        }
        (None, None, None, None, Some(rejection)) => {
            let (cause, sanger) = recorded_rejection(rejection, integrity, callability, expected)?;
            Ok(Recorded::Rejected(Box::new(Rejection {
                input_name: name,
                input_sha256: sha256,
                cause,
                sanger,
            })))
        }
        _ => Err(inconsistent(
            "a document records neither a read with its calls nor a rejection alone",
        )),
    }
}

/// Requires a trace name and identity, the loaded configuration, and this
/// build's evidence plugins.
fn check_provenance(provenance: &RawProvenance, expected: &Expected<'_>) -> Result<()> {
    let name = &provenance.input.name;
    if name.is_empty() || Path::new(name).file_name() != Some(name.as_ref()) {
        return Err(inconsistent("provenance.input.name is not a file name"));
    }
    let sha256 = &provenance.input.sha256;
    if sha256.len() != 64
        || !sha256
            .bytes()
            .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte))
    {
        return Err(inconsistent("provenance.input.sha256 is not a SHA-256"));
    }
    if provenance.configuration_sha256 != expected.configuration_sha256 {
        return Err(ReadEvidenceError::Configuration {
            found: provenance.configuration_sha256.clone(),
            expected: expected.configuration_sha256.to_owned(),
        }
        .into());
    }
    let plugins = expected.plugins;
    let same = provenance.plugins.len() == plugins.len()
        && provenance
            .plugins
            .iter()
            .zip(plugins)
            .all(|(found, plugin)| {
                found.id == plugin.id
                    && found.family == plugin.family.label()
                    && found.version == plugin.version
            });
    if !same {
        let list = |items: Vec<String>| items.join(",");
        return Err(ReadEvidenceError::Plugins {
            found: list(
                provenance
                    .plugins
                    .iter()
                    .map(|plugin| format!("{}@{}", plugin.id, plugin.version))
                    .collect(),
            ),
            expected: list(
                plugins
                    .iter()
                    .map(|plugin| format!("{}@{}", plugin.id, plugin.version))
                    .collect(),
            ),
        }
        .into());
    }
    Ok(())
}

/// A rejection the Sanger plugin records, checked against the loaded
/// configuration's minimums and the evidence its stages produced.
fn recorded_rejection(
    rejection: RawRejection,
    integrity: Option<RawIntegrity>,
    callability: Option<RawCallability>,
    expected: &Expected<'_>,
) -> Result<(ReadRejectionResult, Option<RejectedSanger>)> {
    match rejection {
        RawRejection::TooShort {
            calls,
            minimum_calls,
        } => {
            if minimum_calls != expected.minimum_calls || calls >= minimum_calls {
                return Err(inconsistent(
                    "a too-short read does not fall short of the configured signal window",
                ));
            }
            if integrity.is_some() || callability.is_some() {
                return Err(inconsistent(
                    "a too-short read records integrity or callability",
                ));
            }
            Ok((
                ReadRejectionResult::TooShort {
                    calls,
                    minimum_calls,
                },
                None,
            ))
        }
        RawRejection::CallableCallsBelowMinimum {
            callable_calls,
            minimum_callable_calls,
        } => {
            if minimum_callable_calls != expected.minimum_callable_calls
                || callable_calls >= minimum_callable_calls
            {
                return Err(inconsistent(
                    "a read rejected for its callable calls reaches the configured minimum",
                ));
            }
            let (Some(integrity), Some(raw)) = (integrity, callability) else {
                return Err(inconsistent(
                    "a read rejected for its callable calls lacks its integrity or callability",
                ));
            };
            let count = raw.segments.last().map_or(0, |segment| segment.calls.end);
            let callability = callability_result(&raw, count)?;
            if callable_calls + callability.masked_calls != count {
                return Err(inconsistent(
                    "callable and masked calls do not add up to the read's calls",
                ));
            }
            Ok((
                ReadRejectionResult::CallableCallsBelowMinimum {
                    callable_calls,
                    minimum_callable_calls,
                },
                Some(RejectedSanger {
                    integrity: integrity_result(&integrity),
                    callability,
                }),
            ))
        }
    }
}

/// The kernel evidence and the Sanger evidence of a read the Sanger plugin
/// hands to the core.
fn recorded_read(
    read: &RawRead,
    signal_quality: &RawSignalQuality,
    warnings: &RawWarnings,
    calls: &[RawCall],
) -> Result<(ReadEvidence, SangerRead)> {
    let count = read.call_count;
    if calls.len() != count {
        return Err(inconsistent(
            "read.call_count differs from the number of calls",
        ));
    }
    for (position, call) in calls.iter().enumerate() {
        if call.index != position {
            return Err(ReadEvidenceError::CallOrder {
                position,
                index: call.index,
            }
            .into());
        }
    }
    let sanger = sanger_read(read, signal_quality, warnings, calls)?;
    let reason_order = read
        .reasons
        .iter()
        .map(|label| reason(label))
        .collect::<Result<Vec<_>>>()?;
    let call_evidence = calls
        .iter()
        .map(call_evidence)
        .collect::<Result<Vec<_>>>()?;
    let informative = interval(read.informative);
    let evidence = ReadEvidence::new(
        call_evidence,
        informative.start..informative.end,
        reason_order,
    )?;
    let evidence = match read.ends {
        RawEnds::Unvouched => evidence,
        RawEnds::Vouched => evidence.with_vouched_ends(),
    };
    let evidence = evidence.with_mixtures(read.mixtures.iter().map(mixture).collect())?;
    let unreliable = read
        .unreliable
        .iter()
        .map(|span| {
            Ok(UnreliableSpan {
                calls: span.calls.start..span.calls.end,
                reason: reason(&span.reason)?,
            })
        })
        .collect::<Result<Vec<_>>>()?;
    let evidence = evidence.with_unreliable(unreliable)?;
    let evidence = evidence.with_tracts(read.tracts.iter().map(tract).collect())?;
    Ok((evidence, sanger))
}

/// The Sanger evidence the documents report, checked against the calls it
/// summarizes.
fn sanger_read(
    read: &RawRead,
    signal_quality: &RawSignalQuality,
    warnings: &RawWarnings,
    calls: &[RawCall],
) -> Result<SangerRead> {
    let count = read.call_count;
    let within = |interval: RawInterval| interval.start <= interval.end && interval.end <= count;
    if !within(read.trim)
        || !signal_quality
            .noisy_regions
            .iter()
            .all(|region| within(region.calls) && region.samples.start <= region.samples.end)
    {
        return Err(inconsistent("an interval of the read exceeds its calls"));
    }
    let callability = callability_result(&read.callability, count)?;
    let integrity = &signal_quality.integrity;
    let mismatches = [
        integrity.vendor_primary_count,
        integrity.vendor_quality_count,
    ]
    .into_iter()
    .flatten()
    .filter(|&vendor| vendor != integrity.ploc_count)
    .count();
    let unresolved = calls.iter().filter(|call| call.base == 'N').count();
    if warnings.unresolved_primary_calls != unresolved
        || warnings.ploc_vendor_length_mismatches != mismatches
        || warnings.clipped_channel_samples != integrity.clipped_channel_samples
        || warnings.multi_channel_unresolved_calls > count
        || warnings.vendor_disagreements > count
    {
        return Err(inconsistent(
            "warnings disagree with the calls or the trace integrity",
        ));
    }
    Ok(SangerRead {
        trim: interval(read.trim),
        callability,
        signal_quality: SignalQualityResult {
            integrity: integrity_result(integrity),
            noisy_regions: signal_quality
                .noisy_regions
                .iter()
                .map(|region| NoisyRegionResult {
                    calls: interval(region.calls),
                    samples: interval(region.samples),
                    minimum_primary_snr: region.minimum_primary_snr,
                })
                .collect(),
        },
        warnings: BasecallWarningSummaryResult {
            unresolved_primary_calls: warnings.unresolved_primary_calls,
            multi_channel_unresolved_calls: warnings.multi_channel_unresolved_calls,
            vendor_disagreements: warnings.vendor_disagreements,
            ploc_vendor_length_mismatches: warnings.ploc_vendor_length_mismatches,
            clipped_channel_samples: warnings.clipped_channel_samples,
        },
        calls: calls
            .iter()
            .map(|call| {
                Ok(SangerCallSignal {
                    scan: call.signal.scan,
                    channels: call.signal.channels.channels(),
                    peaks: call.signal.peaks.map(RawPeaks::channels),
                    quality: call.signal.relative_quality,
                    second: second_peak(call.second, call.signal.cause)?,
                })
            })
            .collect::<Result<_>>()?,
    })
}

/// The callability view of a read of `count` calls, whose segments must
/// partition its calls in order.
fn callability_result(raw: &RawCallability, count: usize) -> Result<CallabilityResult> {
    let span = raw.callable_span;
    let mut next = 0;
    for segment in &raw.segments {
        if segment.calls.start != next || segment.calls.end <= segment.calls.start {
            return Err(inconsistent(
                "the callability segments do not partition the read's calls",
            ));
        }
        next = segment.calls.end;
    }
    if next != count || span.start > span.end || span.end > count || raw.masked_calls > count {
        return Err(inconsistent(
            "the callability does not cover exactly the read's calls",
        ));
    }
    Ok(CallabilityResult {
        callable_span: interval(span),
        segments: raw
            .segments
            .iter()
            .map(|segment| CallabilitySegmentResult {
                calls: interval(segment.calls),
                state: match segment.state {
                    RawPhaseState::InPhase => PhaseState::InPhase,
                    RawPhaseState::Dephased => PhaseState::Dephased,
                    RawPhaseState::Mixed => PhaseState::Mixed,
                    RawPhaseState::Weak => PhaseState::Weak,
                    RawPhaseState::Irregular => PhaseState::Irregular,
                },
                after_repeat: segment.after_repeat,
                shadow_offsets: segment.shadow_offsets.clone(),
            })
            .collect(),
        masked_calls: raw.masked_calls,
    })
}

fn integrity_result(integrity: &RawIntegrity) -> TraceIntegrityResult {
    TraceIntegrityResult {
        ploc_count: integrity.ploc_count,
        duplicate_ploc_positions: integrity.duplicate_ploc_positions,
        vendor_primary_count: integrity.vendor_primary_count,
        vendor_quality_count: integrity.vendor_quality_count,
        minimum_ploc_spacing: integrity.minimum_ploc_spacing,
        median_ploc_spacing: integrity.median_ploc_spacing,
        maximum_ploc_spacing: integrity.maximum_ploc_spacing,
        clipped_channel_samples: integrity.clipped_channel_samples,
        maximum_to_median_event_signal_ratio: integrity.maximum_to_median_event_signal_ratio,
    }
}

/// A call's second peak (SRS-CALL-021) from its second base and the cause
/// the modality gives it: a peak whose share exceeds the floor stands out
/// with its shape's cause; a call masked as mixed whose own second peak is
/// small is a `mixed_stretch`, with its second base if it has one.
fn second_peak(second: Option<RawSecond>, cause: Option<RawCause>) -> Result<Option<SecondPeak>> {
    let Some(cause) = cause else {
        if second.is_some_and(|second| second.stands_out) {
            return Err(inconsistent(
                "a second base stands out without the cause of its peak",
            ));
        }
        return Ok(None);
    };
    let cause = match cause {
        RawCause::Noise => SecondCause::Noise,
        RawCause::NeighbourTail => SecondCause::NeighbourTail,
        RawCause::OffCentre => SecondCause::OffCentre,
        RawCause::DyeBlob => SecondCause::DyeBlob,
        RawCause::PhaseShift => SecondCause::PhaseShift,
        RawCause::SecondTemplate => SecondCause::SecondTemplate,
        RawCause::MixedStretch => SecondCause::MixedStretch,
    };
    let stands_out = cause != SecondCause::MixedStretch;
    match second {
        Some(second) if second.stands_out == stands_out => Ok(Some(SecondPeak {
            base: second.base,
            share: second.share,
            stands_out,
            cause,
        })),
        None if !stands_out => Ok(Some(SecondPeak {
            base: 'N',
            share: 0.0,
            stands_out,
            cause,
        })),
        _ => Err(inconsistent(
            "the cause of a second peak disagrees with the call's second base",
        )),
    }
}

/// One call's kernel evidence, through the kernel's own types.
fn call_evidence(call: &RawCall) -> Result<CallEvidence> {
    if !matches!(call.base, 'A' | 'C' | 'G' | 'T' | 'N') {
        return Err(inconsistent("a call's base is not A, C, G, T, or N"));
    }
    let profile = call
        .profile
        .map(|weights| {
            weights
                .iter()
                .all(|weight| weight.is_finite() && (0.0..=1.0).contains(weight))
                .then_some(EvidenceProfile { weights })
                .ok_or_else(|| inconsistent("a profile weight is not in [0, 1]"))
        })
        .transpose()?;
    let placement = match (call.placement.kind, call.placement.reason.as_deref()) {
        (RawPlacementKind::Scored, None) => Placement::Scored,
        (RawPlacementKind::Anchoring, Some(label)) => Placement::Anchoring(reason(label)?),
        (RawPlacementKind::Unresolved, Some(label)) => Placement::Unresolved(reason(label)?),
        _ => {
            return Err(inconsistent(
                "a masked placement has no reason, or a scored one has one",
            ));
        }
    };
    let second = call
        .second
        .map(|second| {
            (matches!(second.base, 'A' | 'C' | 'G' | 'T')
                && second.share.is_finite()
                && (0.0..=0.5).contains(&second.share))
            .then_some(SecondBase {
                base: second.base,
                share: second.share,
                stands_out: second.stands_out,
            })
            .ok_or_else(|| inconsistent("a second base or its share is out of range"))
        })
        .transpose()?;
    Ok(CallEvidence {
        base: call.base,
        profile,
        placement,
        identity: verdict(&call.identity)?,
        indel: verdict(&call.indel)?,
        bound: match call.bound {
            RawBound::InPhase => Bound::InPhase,
            RawBound::Anchored => Bound::Anchored,
            RawBound::Withheld => Bound::Withheld,
        },
        post_repeat: PostRepeatRun {
            call: call.post_repeat.call,
            bound: call.post_repeat.bound,
        },
        second,
        second_unexplained: call.second_unexplained,
        reference_only: call.reference_only,
    })
}

/// A verdict with its reasons in the order raised; a repeated reason is not
/// the contract's.
fn verdict(raw: &RawVerdict) -> Result<Verdict> {
    let mut withheld = Reasons::default();
    for (position, label) in raw.withheld.iter().enumerate() {
        if raw.withheld[..position].contains(label) {
            return Err(inconsistent("a verdict lists a reason twice"));
        }
        withheld.insert(reason(label)?)?;
    }
    Ok(Verdict {
        edge: raw.edge,
        withheld,
    })
}

/// The Sanger reason a label names.
fn reason(label: &str) -> Result<dna_kernel::read_evidence::EvidenceReason> {
    read_processing::reason(label).ok_or_else(|| {
        ReadEvidenceError::UnknownReason {
            label: label.to_owned(),
        }
        .into()
    })
}

fn mixture(raw: &RawMixture) -> PhaseMixture {
    PhaseMixture {
        call_start: raw.calls.start,
        call_end: raw.calls.end,
        run_calls: raw.run_calls,
        components: raw
            .components
            .iter()
            .map(|component| MixtureComponent {
                offset: component.offset,
                fraction: component.fraction,
            })
            .collect(),
        residual: raw.residual,
        unexplained: raw.unexplained.clone(),
        explained: raw.explained,
    }
}

fn tract(raw: &RawTract) -> TractLengths {
    TractLengths {
        onset: raw.onset,
        call_start: raw.calls.start,
        extent: raw.extent,
        window: raw.calls.end.saturating_sub(raw.calls.start),
        run: raw.run,
        mark: raw.mark,
        candidates: raw
            .candidates
            .iter()
            .map(|candidate| StructureFit {
                mark: candidate.mark,
                residual: candidate.residual,
            })
            .collect(),
        residual: raw.residual,
        sigma: raw.sigma,
        autocorrelation: raw.autocorrelation,
        first: marginal(&raw.first),
        total: marginal(&raw.total),
        downstream: raw.downstream.clone(),
        cleaned: raw.cleaned.as_ref().map(|cleaned| CleanedLengths {
            leakage: cleaned.leakage.as_ref().map(|shares| {
                shares
                    .iter()
                    .map(|leak| LeakageShare {
                        offset: leak.offset,
                        share: leak.share,
                    })
                    .collect()
            }),
            residual: cleaned.residual,
            first: marginal(&cleaned.first),
        }),
    }
}

fn marginal(raw: &RawMarginal) -> LengthMarginal {
    LengthMarginal {
        lengths: raw.lengths.clone(),
        shares: raw.shares.clone(),
        errors: raw.errors.clone(),
        dominant: raw.dominant,
        next: raw.next,
        z: raw.z,
        z_adjusted: raw.z_adjusted,
        margin: match raw.margin {
            RawMargin::Clear => Margin::Clear,
            RawMargin::NearEqual => Margin::NearEqual,
            RawMargin::Undefined => Margin::Undefined,
        },
    }
}

const fn interval(raw: RawInterval) -> IntervalResult {
    IntervalResult {
        start: raw.start,
        end: raw.end,
    }
}

fn inconsistent(reason: &'static str) -> crate::error::Error {
    ReadEvidenceError::Inconsistent(reason).into()
}

#[cfg(test)]
mod tests {
    use dna_kernel::read_evidence::{EvidenceReason, ReadEnds};
    use dna_sanger::read_processing::PLUGIN as SANGER;

    use super::*;
    use crate::error::{Error, SangerError};

    const EXAMPLE: &str =
        include_str!("../../docs/reference/examples/read-evidence-v1.example.json");
    const CONFIGURATION: &str = "07e5ef7cf9306eb6c41f21454f49a3d91e404e0b355d34a9a14defa5d4c2a675";

    const EXPECTED: Expected<'static> = Expected {
        configuration_sha256: CONFIGURATION,
        minimum_calls: 10,
        minimum_callable_calls: 20,
        plugins: &[&SANGER],
    };

    fn recorded(document: &str) -> Result<Recorded> {
        parse(
            document.as_bytes(),
            Path::new("trace.evidence.json"),
            &EXPECTED,
        )
    }

    fn read(document: &str) -> Result<RecordedRead> {
        match recorded(document)? {
            Recorded::Read(read) => Ok(*read),
            Recorded::Rejected(_) => Err(inconsistent("the document records a rejection")),
        }
    }

    fn rejection(document: &str) -> Result<Rejection> {
        match recorded(document)? {
            Recorded::Rejected(rejection) => Ok(*rejection),
            Recorded::Read(_) => Err(inconsistent("the document records a read")),
        }
    }

    fn edited(edit: impl FnOnce(&mut serde_json::Value)) -> String {
        let mut value: serde_json::Value = serde_json::from_str(EXAMPLE).unwrap_or_default();
        edit(&mut value);
        value.to_string()
    }

    fn failure(document: &str) -> Option<ReadEvidenceError> {
        match recorded(document) {
            Err(Error::Sanger(SangerError::ReadEvidence(failure))) => Some(failure),
            _ => None,
        }
    }

    // SRS-API-018
    #[test]
    fn reads_the_example_back_into_the_kernel_evidence() -> Result<()> {
        let example = read(EXAMPLE)?;
        let evidence = &example.evidence;
        // The read is named after its trace, whatever the document's name.
        assert_eq!(example.input_name, "trace.ab1");
        assert_eq!(evidence.calls().len(), 6);
        assert_eq!(evidence.informative(), 0..5);
        assert_eq!(evidence.ends(), ReadEnds::Unvouched);
        assert_eq!(
            evidence
                .reason_order()
                .iter()
                .map(|reason| reason.label())
                .collect::<Vec<_>>(),
            [
                "peak_below_minimum",
                "relative_quality_not_above_threshold",
                "mixed_supporting_dna"
            ]
        );
        // Profile weights come back at full precision.
        assert_eq!(
            evidence.calls()[0].profile.map(|profile| profile.weights),
            Some([98.0 / 101.0, 1.0 / 101.0, 1.0 / 101.0, 1.0 / 101.0])
        );
        assert_eq!(
            evidence.calls()[4]
                .placement
                .reason()
                .map(EvidenceReason::label),
            Some("mixed_signal")
        );
        assert_eq!(evidence.mixtures().len(), 1);
        assert_eq!(example.sanger.calls[0].peaks, Some([812, 6, 10, 4]));
        assert_eq!(example.sanger.calls[5].peaks, None);
        assert_eq!(example.sanger.trim.end, 5);
        // The review's Sanger evidence: scan, window peaks, and second peaks.
        assert_eq!(example.sanger.calls[0].scan, 12);
        assert_eq!(example.sanger.calls[0].channels, [812, 6, 10, 4]);
        // Shares come back at full precision.
        let second = example.sanger.calls[4]
            .second
            .ok_or_else(|| inconsistent("no peak"))?;
        assert_eq!(
            (second.base, second.share, second.stands_out, second.cause),
            ('A', 168.0 / 808.0, true, SecondCause::NeighbourTail)
        );
        assert_eq!(example.sanger.calls[0].second, None);
        // A call masked as mixed whose own second peak is small.
        let mixed = read(&edited(|value| {
            value["calls"][4]["second"]["stands_out"] = false.into();
            value["calls"][4]["signal"]["cause"] = "mixed_stretch".into();
        }))?;
        let mixed = mixed.sanger.calls[4]
            .second
            .ok_or_else(|| inconsistent("no peak"))?;
        assert_eq!(
            (mixed.base, mixed.stands_out, mixed.cause),
            ('A', false, SecondCause::MixedStretch)
        );
        Ok(())
    }

    /// The example as the document of a read rejected for too few callable
    /// calls, with the example's integrity and callability.
    fn rejected(edit: impl FnOnce(&mut serde_json::Value)) -> String {
        edited(|value| {
            let integrity = value["signal_quality"]["integrity"].take();
            let callability = value["read"]["callability"].take();
            let mut document = serde_json::json!({
                "schema_version": value["schema_version"].take(),
                "provenance": value["provenance"].take(),
                "rejection": {
                    "reason": "callable_calls_below_minimum",
                    "callable_calls": 4,
                    "minimum_callable_calls": 20
                },
                "integrity": integrity,
                "callability": callability,
            });
            edit(&mut document);
            *value = document;
        })
    }

    // SRS-API-017 SRS-API-019
    #[test]
    fn reads_a_rejected_read_back_with_its_cause() -> Result<()> {
        let read = rejection(&rejected(|_| {}))?;
        assert_eq!(read.input_name, "trace.ab1");
        assert!(matches!(
            read.cause,
            ReadRejectionResult::CallableCallsBelowMinimum {
                callable_calls: 4,
                minimum_callable_calls: 20
            }
        ));
        let sanger = read
            .sanger
            .ok_or_else(|| inconsistent("no Sanger evidence"))?;
        assert_eq!(sanger.callability.masked_calls, 2);
        assert_eq!(sanger.integrity.ploc_count, 6);
        let short = rejection(&rejected(|value| {
            value["rejection"] =
                serde_json::json!({"reason": "too_short", "calls": 6, "minimum_calls": 10});
            value.as_object_mut().map(|object| {
                object.remove("integrity");
                object.remove("callability")
            });
        }))?;
        assert!(matches!(
            short.cause,
            ReadRejectionResult::TooShort {
                calls: 6,
                minimum_calls: 10
            }
        ));
        assert!(short.sanger.is_none());
        // The analysis of a rejected read fails as that of its trace does.
        assert!(matches!(
            set_aside(&short.cause),
            Error::Sanger(SangerError::Signal(SignalError::TooFewCalls {
                calls: 6,
                window: 10
            }))
        ));
        Ok(())
    }

    // SRS-API-019
    #[test]
    fn refuses_a_rejection_the_configuration_or_the_read_contradicts() {
        type Edit = fn(&mut serde_json::Value);
        let cases: [(&str, Edit); 5] = [
            ("another minimum", |value| {
                value["rejection"]["minimum_callable_calls"] = 21.into();
            }),
            ("enough callable calls", |value| {
                value["rejection"]["callable_calls"] = 20.into();
            }),
            ("calls that do not add up", |value| {
                value["rejection"]["callable_calls"] = 3.into();
            }),
            ("no callability", |value| {
                value
                    .as_object_mut()
                    .map(|object| object.remove("callability"));
            }),
            ("calls beside the rejection", |value| {
                value["calls"] = serde_json::json!([]);
            }),
        ];
        for (case, edit) in cases {
            assert!(
                matches!(
                    failure(&rejected(edit)),
                    Some(ReadEvidenceError::Inconsistent(_))
                ),
                "{case}"
            );
        }
        let unknown = rejected(|value| value["rejection"]["reason"] = "ambiguous_placement".into());
        assert!(matches!(
            failure(&unknown),
            Some(ReadEvidenceError::Parse { .. })
        ));
    }

    // SRS-API-018
    #[test]
    fn refuses_another_contract_configuration_or_plugin() {
        let other = edited(|value| value["schema_version"] = "dna.read_evidence/v0".into());
        assert!(matches!(
            failure(&other),
            Some(ReadEvidenceError::UnsupportedSchema { .. })
        ));
        let configuration =
            edited(|value| value["provenance"]["configuration_sha256"] = "0".repeat(64).into());
        assert!(matches!(
            failure(&configuration),
            Some(ReadEvidenceError::Configuration { .. })
        ));
        let plugin = edited(|value| value["provenance"]["plugins"][0]["version"] = 1.into());
        assert!(matches!(
            failure(&plugin),
            Some(ReadEvidenceError::Plugins { .. })
        ));
        let unknown = edited(|value| value["calls"][0]["placement"]["extra"] = true.into());
        assert!(matches!(
            failure(&unknown),
            Some(ReadEvidenceError::Parse { .. })
        ));
    }

    // SRS-API-018
    #[test]
    fn refuses_unknown_reasons_and_calls_out_of_order() {
        let reason = edited(|value| value["read"]["reasons"][0] = "made_up_reason".into());
        assert!(matches!(
            failure(&reason),
            Some(ReadEvidenceError::UnknownReason { label }) if label == "made_up_reason"
        ));
        // A core label is no Sanger reason.
        let mask = edited(|value| value["calls"][5]["placement"]["reason"] = "read_end".into());
        assert!(matches!(
            failure(&mask),
            Some(ReadEvidenceError::UnknownReason { .. })
        ));
        let order = edited(|value| {
            value["calls"][0]["index"] = 1.into();
            value["calls"][1]["index"] = 0.into();
        });
        assert!(matches!(
            failure(&order),
            Some(ReadEvidenceError::CallOrder {
                position: 0,
                index: 1
            })
        ));
    }

    // SRS-API-018
    #[test]
    fn refuses_values_the_contract_or_the_kernel_rejects() {
        type Edit = fn(&mut serde_json::Value);
        let cases: [(&str, Edit); 9] = [
            ("scored call with a reason", |value| {
                value["calls"][0]["placement"]["reason"] = "weak_signal".into();
            }),
            ("weight above one", |value| {
                value["calls"][0]["profile"][0] = 1.5.into();
            }),
            ("warnings that miscount", |value| {
                value["warnings"]["unresolved_primary_calls"] = 0.into();
            }),
            ("trim past the calls", |value| {
                value["read"]["trim"]["end"] = 7.into();
            }),
            ("repeated verdict reason", |value| {
                value["calls"][3]["identity"]["withheld"] =
                    serde_json::json!(["peak_below_minimum", "peak_below_minimum"]);
            }),
            ("call count", |value| {
                value["read"]["call_count"] = 5.into();
            }),
            ("trace path for a name", |value| {
                value["provenance"]["input"]["name"] = "runs/trace.ab1".into();
            }),
            ("standing second base without a cause", |value| {
                value["calls"][4]["signal"]
                    .as_object_mut()
                    .map(|signal| signal.remove("cause"));
            }),
            ("segments with a gap", |value| {
                value["read"]["callability"]["segments"][1]["calls"]["start"] = 5.into();
            }),
        ];
        for (case, edit) in cases {
            assert!(
                matches!(
                    failure(&edited(edit)),
                    Some(ReadEvidenceError::Inconsistent(_))
                ),
                "{case}"
            );
        }
        // The kernel's own checks: a verdict reason the order does not list,
        // a mixture past the calls.
        let unlisted = edited(|value| {
            value["read"]["reasons"] = serde_json::json!(["peak_below_minimum"]);
        });
        assert!(matches!(read(&unlisted), Err(Error::Kernel(_))));
        let mixture = edited(|value| value["read"]["mixtures"][0]["calls"]["end"] = 9.into());
        assert!(matches!(read(&mixture), Err(Error::Kernel(_))));
    }
}
