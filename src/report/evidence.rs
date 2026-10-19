//! Typed assembly of `dna.read_evidence/v1`: one read's evidence exactly as
//! the modality hands it to the core, with the Sanger signal each verdict
//! rests on (PROP-0011).

use crate::config::Config;
use crate::error::{ReportError, Result};
use crate::model::evidence_result::{
    BoundResult, CallEvidenceResult, CallSignalResult, Ends, EvidenceProvenanceResult,
    EvidenceReadResult, EvidenceResult, MixtureComponentResult, MixtureEvidenceResult,
    PostRepeatResult, SecondBaseResult, VerdictResult,
};
use crate::model::result::{InputResult, IntervalResult, PlacementKind, PlacementResult};
use crate::report::callability;
use crate::report::json::project_plugins;
use dna_kernel::plugin::PluginDescriptor;
use dna_kernel::read_evidence::{
    Bound, EvidenceReason, Placement, ReadEnds, ReadEvidence, Verdict,
};
use dna_sanger::model::basecalls::BaseCalls;
use dna_sanger::model::callability::ReadCallability;
use dna_sanger::model::quality::QualityControlResult;
use dna_sanger::model::sanger::Chromatogram;
use dna_sanger::model::second_cause::SecondPeak;

/// Inputs consumed to build one evidence document.
pub(crate) struct CompletedEvidence {
    pub(crate) config: Config,
    pub(crate) trace: Chromatogram,
    pub(crate) calls: BaseCalls,
    pub(crate) callability: ReadCallability,
    pub(crate) quality: QualityControlResult,
    pub(crate) second: Vec<Option<SecondPeak>>,
    pub(crate) evidence: ReadEvidence,
    /// Plugins of the workflow, in execution order.
    pub(crate) plugins: &'static [&'static PluginDescriptor],
}

/// Builds `dna.read_evidence/v1` without filesystem side effects.
pub(crate) fn build(completed: CompletedEvidence) -> Result<EvidenceResult> {
    let CompletedEvidence {
        config,
        trace,
        calls,
        callability: read_callability,
        quality,
        second,
        evidence,
        plugins,
    } = completed;
    let count = evidence.calls().len();
    if calls.calls.len() != count || quality.per_call.len() != count {
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
            let basecall = &calls.calls[index];
            let height = basecall
                .peaks
                .iter()
                .map(|peak| peak.height)
                .max()
                .unwrap_or(0);
            CallEvidenceResult {
                index,
                base: call.base,
                profile: call.profile.map(|profile| profile.weights.map(round)),
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
                second: call.second.map(|second| SecondBaseResult {
                    base: second.base,
                    share: round(second.share),
                    stands_out: second.stands_out,
                }),
                second_unexplained: call.second_unexplained,
                reference_only: call.reference_only,
                signal: CallSignalResult {
                    scan: basecall.locus_position_0based,
                    height: u32::try_from(height).unwrap_or(0),
                    relative_quality: quality.per_call[index].relative_quality_score,
                    cause: second
                        .get(index)
                        .copied()
                        .flatten()
                        .map(|peak| peak.cause.label()),
                },
            }
        })
        .collect();
    Ok(EvidenceResult {
        schema_version: "dna.read_evidence/v1",
        provenance: EvidenceProvenanceResult {
            input: InputResult {
                sha256: trace.source_sha256,
            },
            configuration_sha256: config.source_sha256,
            plugins: project_plugins(plugins),
        },
        read: EvidenceReadResult {
            call_count: count,
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
            callability: callability::project(&read_callability),
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
                    residual: round(mixture.residual),
                    components: mixture
                        .components
                        .iter()
                        .map(|component| MixtureComponentResult {
                            offset: component.offset,
                            fraction: round(component.fraction),
                        })
                        .collect(),
                    unexplained: mixture.unexplained.clone(),
                })
                .collect(),
        },
        calls: results,
    })
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

/// A share or weight to three decimals, as every DNA document writes them.
fn round(value: f64) -> f64 {
    (value * 1000.0).round() / 1000.0
}
