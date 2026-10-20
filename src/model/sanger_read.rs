//! The Sanger evidence of one read that the documents report beside the
//! core's products, as `dna.read_evidence/v1` carries it (PROP-0011), and its
//! one projection from the Sanger plugin's products.

use crate::error::{CallEvidenceError, ReportError, Result};
use crate::model::basecall_result::BasecallWarningSummaryResult;
use crate::model::read_observation::{RejectedRead, RejectionCause};
use crate::model::result::{
    CallabilityResult, CallabilitySegmentResult, IntervalResult, NoisyRegionResult,
    SignalQualityResult, TraceIntegrityResult,
};
use crate::model::sample_result::ReadRejectionResult;
use dna_call::read_call::PlacementRejection;
use dna_sanger::model::attachment::SangerAttachment;
use dna_sanger::model::basecalls::BaseCalls;
use dna_sanger::model::callability::{PhaseState, ReadCallability};
use dna_sanger::model::quality::QualityControlResult;
use dna_sanger::model::second_cause::SecondPeak;
use dna_sanger::model::signal::{SangerIntegrity, SignalAnalysis};
use dna_sanger::read_processing::ReadWarnings;

/// One read's Sanger evidence as the documents report it: projected from the
/// Sanger products of a trace, or read back from its read evidence document.
#[derive(Debug, Clone)]
pub(crate) struct SangerRead {
    /// The relative-quality trim interval.
    pub(crate) trim: IntervalResult,
    /// The read callability view.
    pub(crate) callability: CallabilityResult,
    /// Trace integrity and merged noisy regions.
    pub(crate) signal_quality: SignalQualityResult,
    /// The modality's warning counts.
    pub(crate) warnings: BasecallWarningSummaryResult,
    /// Each call's signal, in call order.
    pub(crate) calls: Vec<SangerCallSignal>,
}

impl SangerRead {
    /// Projects the trim, callability, signal quality, warning counts, and
    /// every call's signal of one read.
    pub(crate) fn project(
        calls: &BaseCalls,
        signal: &SignalAnalysis,
        read_callability: &ReadCallability,
        quality: &QualityControlResult,
        second: &[Option<SecondPeak>],
    ) -> Result<Self> {
        if quality.per_call.len() != calls.len() {
            return Err(
                ReportError::Inconsistent("Sanger calls and quality differ in call count").into(),
            );
        }
        let call_signals = (0..calls.len())
            .map(|index| {
                SangerCallSignal::of(calls, quality, second, index)
                    .map_err(ReportError::CallEvidence)
            })
            .collect::<std::result::Result<Vec<_>, _>>()?;
        Ok(Self {
            trim: IntervalResult {
                start: quality.trim_start_0based,
                end: quality.trim_end_0based_exclusive,
            },
            callability: callability(read_callability),
            signal_quality: signal_quality(signal),
            warnings: warnings(calls, &signal.integrity),
            calls: call_signals,
        })
    }

    /// Projects the Sanger products a sample operation keeps of an admitted
    /// read.
    pub(crate) fn of(sanger: &SangerAttachment) -> Result<Self> {
        Self::project(
            &sanger.calls,
            &sanger.signal,
            &sanger.callability,
            &sanger.quality,
            &sanger.second,
        )
    }

    /// Whether call `index` lies in a merged candidate-noisy region.
    pub(crate) fn in_noisy_region(&self, index: usize) -> bool {
        self.signal_quality
            .noisy_regions
            .iter()
            .any(|region| region.calls.start <= index && index < region.calls.end)
    }

    /// The read's warning counts for the operational warning summary.
    pub(crate) const fn read_warnings(&self) -> ReadWarnings {
        let warnings = &self.warnings;
        ReadWarnings {
            unresolved_primary_calls: warnings.unresolved_primary_calls,
            multi_channel_unresolved_calls: warnings.multi_channel_unresolved_calls,
            vendor_disagreements: warnings.vendor_disagreements,
            locus_vendor_length_mismatches: warnings.ploc_vendor_length_mismatches,
            clipped_channel_samples: warnings.clipped_channel_samples,
        }
    }

    /// Each call's selected A/C/G/T peak heights in its window, in call
    /// order, as the run signal counts them.
    pub(crate) fn heights(&self) -> Vec<[i32; 4]> {
        self.calls.iter().map(|call| call.channels).collect()
    }
}

/// The signal of one call the documents join to the core's records.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct SangerCallSignal {
    /// The trace sample of the call's locus.
    pub(crate) scan: usize,
    /// The selected A/C/G/T peak heights in the call's window, on the read
    /// strand.
    pub(crate) channels: [i32; 4],
    /// A/C/G/T channel heights at the primary event, on the read strand;
    /// `None` when the call has no primary-event peak evidence.
    pub(crate) peaks: Option<[i32; 4]>,
    /// The call's relative quality.
    pub(crate) quality: u8,
    /// The call's second peak and its cause (SRS-CALL-021), when it has one.
    pub(crate) second: Option<SecondPeak>,
}

impl SangerCallSignal {
    /// The signal of call `index` from its call, quality, and second-peak
    /// records, on the read strand.
    pub(crate) fn of(
        calls: &BaseCalls,
        quality: &QualityControlResult,
        second: &[Option<SecondPeak>],
        index: usize,
    ) -> std::result::Result<Self, CallEvidenceError> {
        let call = calls
            .calls
            .get(index)
            .ok_or(CallEvidenceError::MissingCall { index })?;
        let score = quality
            .per_call
            .get(index)
            .ok_or(CallEvidenceError::MissingQuality { index })?;
        if call.index_0based != index || score.index_0based != index {
            return Err(CallEvidenceError::IndexMismatch { index });
        }
        Ok(Self {
            scan: call.locus_position_0based,
            channels: call.peaks.map(|peak| peak.height),
            peaks: call
                .primary_peak_evidence
                .as_ref()
                .map(|primary| primary.channel_heights),
            quality: score.relative_quality_score,
            second: second.get(index).copied().flatten(),
        })
    }
}

/// A read a sample operation set aside, as the sample documents report it:
/// projected from the rejection of a trace, or read back from a read
/// evidence document.
#[derive(Debug, Clone)]
pub(crate) struct Rejection {
    /// Reviewer-facing name of the input: the trace's file name.
    pub(crate) input_name: String,
    /// SHA-256 content identity of the trace.
    pub(crate) input_sha256: String,
    /// Why the read was set aside, with the measured value and its minimum.
    pub(crate) cause: ReadRejectionResult,
    /// Integrity and callability of the read, absent for a read too short
    /// to process.
    pub(crate) sanger: Option<RejectedSanger>,
}

impl Rejection {
    /// Projects a read a sample operation set aside: its cause and, when it
    /// reached callability, its integrity and callability.
    pub(crate) fn of(read: &RejectedRead) -> Self {
        let cause = match read.cause {
            RejectionCause::TooShort(short) => ReadRejectionResult::TooShort {
                calls: short.calls,
                minimum_calls: short.minimum_calls,
            },
            RejectionCause::Callability { rejection, .. } => {
                ReadRejectionResult::CallableCallsBelowMinimum {
                    callable_calls: rejection.callable_calls,
                    minimum_callable_calls: rejection.minimum_callable_calls,
                }
            }
            RejectionCause::Placement { rejection, .. } => placement(rejection),
        };
        Self {
            input_name: read.input_name.clone(),
            input_sha256: read.input_sha256.clone(),
            cause,
            sanger: read.cause.sanger().map(|sanger| RejectedSanger {
                integrity: integrity(&sanger.integrity),
                callability: callability(&sanger.callability),
            }),
        }
    }
}

/// Integrity and callability of a rejected read that reached callability.
#[derive(Debug, Clone)]
pub(crate) struct RejectedSanger {
    pub(crate) integrity: TraceIntegrityResult,
    pub(crate) callability: CallabilityResult,
}

/// Why the core could not place a read.
pub(crate) const fn placement(rejection: PlacementRejection) -> ReadRejectionResult {
    match rejection {
        PlacementRejection::CallableColumnsBelowMinimum {
            callable_columns,
            minimum,
        } => ReadRejectionResult::CallableColumnsBelowMinimum {
            callable_columns,
            minimum_callable_columns: minimum,
        },
        PlacementRejection::IdentityBelowMinimum { identity, minimum } => {
            ReadRejectionResult::IdentityBelowMinimum {
                identity,
                minimum_identity: minimum,
            }
        }
        PlacementRejection::AmbiguousPlacement => ReadRejectionResult::AmbiguousPlacement,
        PlacementRejection::OrientationTie => ReadRejectionResult::OrientationTie,
    }
}

/// Projects the callable span, ordered phase segments with the shadow offsets
/// of dephased segments, and masked-call count; per-position features, shadow
/// shares, and the mask itself stay internal.
pub(crate) fn callability(callability: &ReadCallability) -> CallabilityResult {
    CallabilityResult {
        callable_span: IntervalResult {
            start: callability.callable_start_0based,
            end: callability.callable_end_0based_exclusive,
        },
        segments: callability
            .segments
            .iter()
            .map(|segment| CallabilitySegmentResult {
                calls: IntervalResult {
                    start: segment.call_start_0based,
                    end: segment.call_end_0based_exclusive,
                },
                state: segment.state,
                after_repeat: segment.after_repeat,
                shadow_offsets: segment
                    .shadow
                    .filter(|_| segment.state == PhaseState::Dephased)
                    .map(|shadow| shadow.offsets().collect()),
            })
            .collect(),
        masked_calls: callability.masked_count(),
    }
}

/// Projects merged noisy regions while omitting internal rolling windows.
pub(crate) fn signal_quality(signal: &SignalAnalysis) -> SignalQualityResult {
    SignalQualityResult {
        integrity: integrity(&signal.integrity),
        noisy_regions: signal
            .noisy_regions
            .iter()
            .map(|region| NoisyRegionResult {
                calls: IntervalResult {
                    start: region.call_start_0based,
                    end: region.call_end_0based_exclusive,
                },
                samples: IntervalResult {
                    start: region.sample_start_0based,
                    end: region.sample_end_0based_exclusive,
                },
                minimum_primary_snr: region.minimum_primary_snr,
            })
            .collect(),
    }
}

/// Projects immutable trace-integrity evidence.
pub(crate) fn integrity(integrity: &SangerIntegrity) -> TraceIntegrityResult {
    TraceIntegrityResult {
        ploc_count: integrity.locus_count,
        duplicate_ploc_positions: integrity.duplicate_locus_count,
        vendor_primary_count: integrity.vendor_primary_count,
        vendor_quality_count: integrity.vendor_quality_count,
        minimum_ploc_spacing: integrity.minimum_locus_spacing,
        median_ploc_spacing: integrity.median_locus_spacing,
        maximum_ploc_spacing: integrity.maximum_locus_spacing,
        clipped_channel_samples: integrity.clipped_channel_samples,
        maximum_to_median_event_signal_ratio: integrity.maximum_to_median_event_signal_ratio,
    }
}

/// The modality's warning counts of one read.
pub(crate) fn warnings(
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
