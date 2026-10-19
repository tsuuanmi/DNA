//! Read evidence documents (`dna.read_evidence/v1`) as input to `analyze`
//! (PROP-0011): one read's evidence, read back into the kernel's
//! `ReadEvidence` through its own constructors, with the Sanger evidence the
//! analysis document reports beside the core's products.

use std::fs;
use std::path::Path;

use serde::Deserialize;

use crate::config::Config;
use crate::error::{KernelError, ReadEvidenceError, Result};
use crate::model::basecall_result::BasecallWarningSummaryResult;
use crate::model::result::{
    CallabilityResult, CallabilitySegmentResult, IntervalResult, NoisyRegionResult,
    SignalQualityResult, TraceIntegrityResult,
};
use crate::model::sanger_read::{SangerCallSignal, SangerRead};
use dna_kernel::model::reference::Reference;
use dna_kernel::plugin::PluginDescriptor;
use dna_kernel::profile::Profile;
use dna_kernel::read_evidence::{
    Bound, CallEvidence, EvidenceProfile, MixtureComponent, PhaseMixture, Placement, PostRepeatRun,
    ReadEvidence, Reasons, SecondBase, Verdict,
};
use dna_sanger::model::callability::PhaseState;
use dna_sanger::read_processing;

use super::{load_config, load_profile, load_reference, require_regular_file};

/// Largest accepted read evidence document.
const MAX_EVIDENCE_BYTES: u64 = 64 * 1024 * 1024;

/// The contract this reader accepts.
const SCHEMA_VERSION: &str = "dna.read_evidence/v1";

/// Scientific inputs for analyzing one recorded read.
pub(crate) struct EvidenceInputs {
    pub(crate) config: Config,
    pub(crate) profile: Profile,
    pub(crate) reference: Reference,
    pub(crate) read: RecordedRead,
}

/// One read as its read evidence document records it.
pub(crate) struct RecordedRead {
    /// Reviewer-facing name of the input: the document's file name.
    pub(crate) input_name: String,
    /// SHA-256 of the trace the evidence was written from.
    pub(crate) input_sha256: String,
    /// The read's evidence as the core receives it.
    pub(crate) evidence: ReadEvidence,
    /// The read's Sanger evidence as the analysis document reports it.
    pub(crate) sanger: SangerRead,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct RawDocument {
    schema_version: String,
    provenance: RawProvenance,
    read: RawRead,
    signal_quality: RawSignalQuality,
    warnings: RawWarnings,
    calls: Vec<RawCall>,
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
    sha256: String,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct RawPlugin {
    id: String,
    family: String,
    version: u32,
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
    #[serde(rename = "scan")]
    _scan: usize,
    #[serde(rename = "height")]
    _height: u32,
    peaks: Option<RawPeaks>,
    relative_quality: u8,
    #[serde(rename = "cause")]
    _cause: Option<RawCause>,
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

/// The second-peak causes the contract names; read only to check the shape.
#[derive(Deserialize)]
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
/// under the loaded configuration.
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
    let read = parse(&bytes, document, &config.source_sha256, plugins)?;
    let reference = load_reference(reference_path, &profile)?;
    Ok(EvidenceInputs {
        config,
        profile,
        reference,
        read,
    })
}

/// Reads one document into the read it records.
fn parse(
    bytes: &[u8],
    path: &Path,
    configuration_sha256: &str,
    plugins: &[&PluginDescriptor],
) -> Result<RecordedRead> {
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
    check_provenance(&raw.provenance, configuration_sha256, plugins)?;
    let count = raw.read.call_count;
    if raw.calls.len() != count {
        return Err(inconsistent(
            "read.call_count differs from the number of calls",
        ));
    }
    for (position, call) in raw.calls.iter().enumerate() {
        if call.index != position {
            return Err(ReadEvidenceError::CallOrder {
                position,
                index: call.index,
            }
            .into());
        }
    }
    let sanger = sanger_read(&raw, count)?;
    let reason_order = raw
        .read
        .reasons
        .iter()
        .map(|label| reason(label))
        .collect::<Result<Vec<_>>>()?;
    let calls = raw
        .calls
        .iter()
        .map(call_evidence)
        .collect::<Result<Vec<_>>>()?;
    let informative = interval(raw.read.informative);
    let evidence = ReadEvidence::new(calls, informative.start..informative.end, reason_order)?;
    let evidence = match raw.read.ends {
        RawEnds::Unvouched => evidence,
        RawEnds::Vouched => evidence.with_vouched_ends(),
    };
    let evidence = evidence.with_mixtures(raw.read.mixtures.iter().map(mixture).collect())?;
    Ok(RecordedRead {
        input_name: path
            .file_name()
            .map(|name| name.to_string_lossy().into_owned())
            .unwrap_or_default(),
        input_sha256: raw.provenance.input.sha256,
        evidence,
        sanger,
    })
}

/// Requires a trace identity, the loaded configuration, and this build's
/// evidence plugins.
fn check_provenance(
    provenance: &RawProvenance,
    configuration_sha256: &str,
    plugins: &[&PluginDescriptor],
) -> Result<()> {
    let sha256 = &provenance.input.sha256;
    if sha256.len() != 64
        || !sha256
            .bytes()
            .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte))
    {
        return Err(inconsistent("provenance.input.sha256 is not a SHA-256"));
    }
    if provenance.configuration_sha256 != configuration_sha256 {
        return Err(ReadEvidenceError::Configuration {
            found: provenance.configuration_sha256.clone(),
            expected: configuration_sha256.to_owned(),
        }
        .into());
    }
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

/// The Sanger evidence the analysis document reports, checked against the
/// calls it summarizes.
fn sanger_read(raw: &RawDocument, count: usize) -> Result<SangerRead> {
    let within = |interval: RawInterval| interval.start <= interval.end && interval.end <= count;
    let callability = &raw.read.callability;
    if !within(raw.read.trim)
        || !within(callability.callable_span)
        || !callability
            .segments
            .iter()
            .all(|segment| within(segment.calls))
        || callability.masked_calls > count
        || !raw
            .signal_quality
            .noisy_regions
            .iter()
            .all(|region| within(region.calls) && region.samples.start <= region.samples.end)
    {
        return Err(inconsistent("an interval of the read exceeds its calls"));
    }
    let integrity = &raw.signal_quality.integrity;
    let warnings = &raw.warnings;
    let mismatches = [
        integrity.vendor_primary_count,
        integrity.vendor_quality_count,
    ]
    .into_iter()
    .flatten()
    .filter(|&vendor| vendor != integrity.ploc_count)
    .count();
    let unresolved = raw.calls.iter().filter(|call| call.base == 'N').count();
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
        trim: interval(raw.read.trim),
        callability: CallabilityResult {
            callable_span: interval(callability.callable_span),
            segments: callability
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
            masked_calls: callability.masked_calls,
        },
        signal_quality: SignalQualityResult {
            integrity: TraceIntegrityResult {
                ploc_count: integrity.ploc_count,
                duplicate_ploc_positions: integrity.duplicate_ploc_positions,
                vendor_primary_count: integrity.vendor_primary_count,
                vendor_quality_count: integrity.vendor_quality_count,
                minimum_ploc_spacing: integrity.minimum_ploc_spacing,
                median_ploc_spacing: integrity.median_ploc_spacing,
                maximum_ploc_spacing: integrity.maximum_ploc_spacing,
                clipped_channel_samples: integrity.clipped_channel_samples,
                maximum_to_median_event_signal_ratio: integrity
                    .maximum_to_median_event_signal_ratio,
            },
            noisy_regions: raw
                .signal_quality
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
        calls: raw
            .calls
            .iter()
            .map(|call| SangerCallSignal {
                peaks: call
                    .signal
                    .peaks
                    .map(|peaks| [peaks.a, peaks.c, peaks.g, peaks.t]),
                quality: call.signal.relative_quality,
            })
            .collect(),
    })
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

    fn read(document: &str) -> Result<RecordedRead> {
        parse(
            document.as_bytes(),
            Path::new("trace.evidence.json"),
            CONFIGURATION,
            &[&SANGER],
        )
    }

    fn edited(edit: impl FnOnce(&mut serde_json::Value)) -> String {
        let mut value: serde_json::Value = serde_json::from_str(EXAMPLE).unwrap_or_default();
        edit(&mut value);
        value.to_string()
    }

    fn failure(document: &str) -> Option<ReadEvidenceError> {
        match read(document) {
            Err(Error::Sanger(SangerError::ReadEvidence(failure))) => Some(failure),
            _ => None,
        }
    }

    // SRS-API-018
    #[test]
    fn reads_the_example_back_into_the_kernel_evidence() -> Result<()> {
        let read = read(EXAMPLE)?;
        let evidence = &read.evidence;
        assert_eq!(read.input_name, "trace.evidence.json");
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
        assert_eq!(read.sanger.calls[0].peaks, Some([812, 6, 10, 4]));
        assert_eq!(read.sanger.calls[5].peaks, None);
        assert_eq!(read.sanger.trim.end, 5);
        Ok(())
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
        let cases: [(&str, Edit); 6] = [
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
