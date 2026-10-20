//! Typed assembly of `dna.read_evidence/v1`: one read's evidence exactly as
//! the modality hands it to the core, with the Sanger signal each verdict
//! rests on, or why the modality set the read aside (PROP-0011).

use crate::config::Config;
use crate::error::{ReportError, Result};
use crate::model::evidence_result::{
    BoundResult, CallEvidenceResult, CallSignalResult, Ends, EvidenceDocument, EvidenceInputResult,
    EvidenceProvenanceResult, EvidenceReadResult, EvidenceResult, MixtureComponentResult,
    MixtureEvidenceResult, PostRepeatResult, RejectedEvidenceResult, SecondBaseResult,
    UnreliableEvidenceResult, VerdictResult,
};
use crate::model::result::{IntervalResult, PeakHeightsResult, PlacementKind, PlacementResult};
use crate::model::sanger_read::{Rejection, SangerRead};
use crate::report::json::project_plugins;
use dna_kernel::plugin::PluginDescriptor;
use dna_kernel::read_evidence::{
    Bound, EvidenceReason, Placement, ReadEnds, ReadEvidence, Verdict,
};

/// The schema every read evidence document names.
const SCHEMA_VERSION: &str = "dna.read_evidence/v1";

/// Inputs consumed to build one evidence document.
pub(crate) struct CompletedEvidence {
    pub(crate) config: Config,
    /// The trace's file name and content SHA-256.
    pub(crate) input: (String, String),
    pub(crate) outcome: Evidenced,
    /// Plugins of the workflow, in execution order.
    pub(crate) plugins: &'static [&'static PluginDescriptor],
}

/// What the modality made of the read.
pub(crate) enum Evidenced {
    /// The read's evidence as the core receives it, and its Sanger evidence
    /// as the documents report it.
    Read {
        sanger: Box<SangerRead>,
        evidence: Box<ReadEvidence>,
    },
    /// The modality set the read aside.
    Rejected(Box<Rejection>),
}

/// Builds `dna.read_evidence/v1` without filesystem side effects.
pub(crate) fn build(completed: CompletedEvidence) -> Result<EvidenceDocument> {
    let CompletedEvidence {
        config,
        input: (name, sha256),
        outcome,
        plugins,
    } = completed;
    let provenance = EvidenceProvenanceResult {
        input: EvidenceInputResult { name, sha256 },
        configuration_sha256: config.source_sha256,
        plugins: project_plugins(plugins),
    };
    let (sanger, evidence) = match outcome {
        Evidenced::Read { sanger, evidence } => (*sanger, *evidence),
        Evidenced::Rejected(rejection) => {
            let Rejection { cause, sanger, .. } = *rejection;
            let (integrity, callability) = sanger
                .map(|sanger| (sanger.integrity, sanger.callability))
                .unzip();
            return Ok(EvidenceDocument::Rejected(RejectedEvidenceResult {
                schema_version: SCHEMA_VERSION,
                provenance,
                rejection: cause,
                integrity,
                callability,
            }));
        }
    };
    let count = evidence.calls().len();
    if sanger.calls.len() != count {
        return Err(
            ReportError::Inconsistent("evidence and basecalls differ in call count").into(),
        );
    }
    let informative = evidence.informative();
    let results = evidence
        .calls()
        .iter()
        .enumerate()
        .map(|(index, call)| {
            let signal = sanger.calls[index];
            CallEvidenceResult {
                index,
                base: call.base,
                // At full precision: alignment scores quantize the weights
                // finer than three decimals, so a rounded profile could
                // place or call the read differently when read back.
                profile: call.profile.map(|profile| profile.weights),
                placement: placement(call.placement),
                identity: verdict(call.identity),
                indel: verdict(call.indel),
                bound: match call.bound {
                    Bound::InPhase => BoundResult::InPhase,
                    Bound::Anchored => BoundResult::Anchored,
                    Bound::Withheld => BoundResult::Withheld,
                },
                post_repeat: PostRepeatResult {
                    call: call.post_repeat.call,
                    bound: call.post_repeat.bound,
                },
                // At full precision, as every other value read back: the
                // review judges each share against its threshold.
                second: call.second.map(|second| SecondBaseResult {
                    base: second.base,
                    share: second.share,
                    stands_out: second.stands_out,
                }),
                second_unexplained: call.second_unexplained,
                reference_only: call.reference_only,
                signal: CallSignalResult {
                    scan: signal.scan,
                    channels: PeakHeightsResult::from(signal.channels),
                    peaks: signal.peaks.map(PeakHeightsResult::from),
                    relative_quality: signal.quality,
                    cause: signal.second.map(|peak| peak.cause.label()),
                },
            }
        })
        .collect();
    Ok(EvidenceDocument::Read(EvidenceResult {
        schema_version: SCHEMA_VERSION,
        provenance,
        read: EvidenceReadResult {
            call_count: count,
            trim: sanger.trim,
            informative: IntervalResult {
                start: informative.start,
                end: informative.end,
            },
            ends: match evidence.ends() {
                ReadEnds::Unvouched => Ends::Unvouched,
                ReadEnds::Vouched => Ends::Vouched,
            },
            reasons: evidence
                .reason_order()
                .iter()
                .copied()
                .map(EvidenceReason::label)
                .collect(),
            callability: sanger.callability,
            mixtures: evidence
                .mixtures()
                .iter()
                .map(|mixture| MixtureEvidenceResult {
                    calls: IntervalResult {
                        start: mixture.call_start,
                        end: mixture.call_end,
                    },
                    run_calls: mixture.run_calls,
                    explained: mixture.explained,
                    residual: mixture.residual,
                    components: mixture
                        .components
                        .iter()
                        .map(|component| MixtureComponentResult {
                            offset: component.offset,
                            fraction: component.fraction,
                        })
                        .collect(),
                    unexplained: mixture.unexplained.clone(),
                })
                .collect(),
            unreliable: evidence
                .unreliable()
                .iter()
                .map(|span| UnreliableEvidenceResult {
                    calls: IntervalResult {
                        start: span.calls.start,
                        end: span.calls.end,
                    },
                    reason: span.reason.label(),
                })
                .collect(),
        },
        signal_quality: sanger.signal_quality,
        warnings: sanger.warnings,
        calls: results,
    }))
}

/// The kernel placement as its kind and the modality's reason.
pub(super) fn placement(placement: Placement) -> PlacementResult {
    let kind = match placement {
        Placement::Scored => PlacementKind::Scored,
        Placement::Anchoring(_) => PlacementKind::Anchoring,
        Placement::Unresolved(_) => PlacementKind::Unresolved,
    };
    PlacementResult {
        kind,
        reason: placement.reason().map(EvidenceReason::label),
    }
}

/// A verdict's edge flag and reasons, in the order the modality raised them.
fn verdict(verdict: Verdict) -> VerdictResult {
    VerdictResult {
        edge: verdict.edge,
        withheld: verdict.withheld.iter().map(EvidenceReason::label).collect(),
    }
}
