//! Shared projection of signal-derived read callability.

use crate::model::callability::{PhaseState, ReadCallability};
use crate::model::result::{CallabilityResult, CallabilitySegmentResult, IntervalResult};

/// Projects the callable span, ordered phase segments with the shadow offsets
/// of dephased segments, and masked-call count; per-position features, shadow
/// shares, and the mask itself stay internal.
pub(super) fn project(callability: &ReadCallability) -> CallabilityResult {
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
