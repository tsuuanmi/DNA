//! The per-position mask and callable span, with their invariants.

use crate::error::{CallabilityError, Result};
use crate::model::callability::{PhaseSegment, PhaseState, ReadCallability, RepeatRun};

/// Assembles the read's callability and checks that the segments partition the
/// read in order.
pub(super) fn build(
    repeats: Vec<RepeatRun>,
    segments: Vec<PhaseSegment>,
    positions: usize,
) -> Result<ReadCallability> {
    let mut expected_start = 0;
    for segment in &segments {
        if segment.call_start_0based != expected_start
            || segment.call_end_0based_exclusive <= segment.call_start_0based
        {
            return Err(
                CallabilityError::Inconsistent("segments do not partition the read").into(),
            );
        }
        expected_start = segment.call_end_0based_exclusive;
    }
    if expected_start != positions {
        return Err(CallabilityError::Inconsistent("segments do not cover every call").into());
    }
    let mut mask = vec![None; positions];
    for segment in &segments {
        if segment.state != PhaseState::InPhase {
            mask[segment.call_start_0based..segment.call_end_0based_exclusive]
                .fill(Some(segment.state));
        }
    }
    let callable_start_0based = mask.iter().position(Option::is_none).unwrap_or(positions);
    let callable_end_0based_exclusive = mask
        .iter()
        .rposition(Option::is_none)
        .map_or(callable_start_0based, |index| index + 1);
    Ok(ReadCallability {
        repeats,
        segments,
        mask,
        callable_start_0based,
        callable_end_0based_exclusive,
    })
}

#[cfg(test)]
mod tests {
    use crate::error::Error;

    use super::*;

    fn segment(start: usize, end: usize, state: PhaseState) -> PhaseSegment {
        PhaseSegment {
            call_start_0based: start,
            call_end_0based_exclusive: end,
            state,
            after_repeat: false,
            coherence: None,
            modal_offset: None,
        }
    }

    #[test]
    fn masks_every_position_outside_in_phase_segments() -> Result<()> {
        let callability = build(
            Vec::new(),
            vec![
                segment(0, 2, PhaseState::Weak),
                segment(2, 6, PhaseState::InPhase),
                segment(6, 7, PhaseState::Dephased),
                segment(7, 9, PhaseState::InPhase),
                segment(9, 10, PhaseState::Mixed),
            ],
            10,
        )?;
        assert_eq!(
            callability.mask,
            [
                Some(PhaseState::Weak),
                Some(PhaseState::Weak),
                None,
                None,
                None,
                None,
                Some(PhaseState::Dephased),
                None,
                None,
                Some(PhaseState::Mixed),
            ]
        );
        assert_eq!(
            (
                callability.callable_start_0based,
                callability.callable_end_0based_exclusive
            ),
            (2, 9)
        );
        assert_eq!(callability.callable_count(), 6);
        assert_eq!(callability.masked_count(), 4);
        assert_eq!(callability.segment_count(PhaseState::InPhase), 2);
        Ok(())
    }

    #[test]
    fn reports_an_empty_callable_span_when_everything_is_masked() -> Result<()> {
        let callability = build(Vec::new(), vec![segment(0, 3, PhaseState::Mixed)], 3)?;
        assert_eq!(
            (
                callability.callable_start_0based,
                callability.callable_end_0based_exclusive
            ),
            (3, 3)
        );
        assert_eq!(callability.callable_count(), 0);
        Ok(())
    }

    #[test]
    fn rejects_segments_that_do_not_partition_the_read() {
        for segments in [
            vec![
                segment(0, 2, PhaseState::InPhase),
                segment(3, 4, PhaseState::Weak),
            ],
            vec![segment(0, 2, PhaseState::InPhase)],
            vec![segment(0, 0, PhaseState::InPhase)],
        ] {
            assert!(matches!(
                build(Vec::new(), segments, 4),
                Err(Error::Callability(CallabilityError::Inconsistent(_)))
            ));
        }
    }
}
