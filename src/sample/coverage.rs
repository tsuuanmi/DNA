//! Run-length encoded sample coverage topology from independently placed reads.

use std::collections::BTreeMap;

use crate::error::{Error, Result};
use crate::model::alignment::Orientation;
use crate::model::sample_evidence::{SampleCoverageEvidence, SampleReadEvidence};

#[derive(Debug, Default, Clone, Copy)]
struct CoverageEvent {
    forward_starts: usize,
    forward_ends: usize,
    reverse_starts: usize,
    reverse_ends: usize,
}

/// Summarizes reference-coordinate coverage depth without interpreting nucleotide agreement.
pub(super) fn summarize(reads: &[SampleReadEvidence]) -> Result<Vec<SampleCoverageEvidence>> {
    let mut events: BTreeMap<usize, CoverageEvent> = BTreeMap::new();

    for read in reads {
        validate_segments(read)?;
        for segment in &read.alignment.reference_segments {
            match read.alignment.orientation {
                Orientation::Forward => {
                    let start = events.entry(segment.start_0based).or_default();
                    start.forward_starts =
                        start.forward_starts.checked_add(1).ok_or_else(|| {
                            Error::Sample("forward coverage start count overflow".into())
                        })?;
                    let end = events.entry(segment.end_0based_exclusive).or_default();
                    end.forward_ends = end.forward_ends.checked_add(1).ok_or_else(|| {
                        Error::Sample("forward coverage end count overflow".into())
                    })?;
                }
                Orientation::Reverse => {
                    let start = events.entry(segment.start_0based).or_default();
                    start.reverse_starts =
                        start.reverse_starts.checked_add(1).ok_or_else(|| {
                            Error::Sample("reverse coverage start count overflow".into())
                        })?;
                    let end = events.entry(segment.end_0based_exclusive).or_default();
                    end.reverse_ends = end.reverse_ends.checked_add(1).ok_or_else(|| {
                        Error::Sample("reverse coverage end count overflow".into())
                    })?;
                }
            }
        }
    }

    let mut coverage = Vec::new();
    let mut previous = None;
    let mut forward_depth = 0usize;
    let mut reverse_depth = 0usize;

    for (position, event) in events {
        if let Some(start) = previous
            && start < position
            && (forward_depth > 0 || reverse_depth > 0)
        {
            push_segment(&mut coverage, start, position, forward_depth, reverse_depth)?;
        }

        forward_depth = forward_depth
            .checked_sub(event.forward_ends)
            .ok_or_else(|| Error::Sample("forward coverage end exceeds active depth".into()))?;
        reverse_depth = reverse_depth
            .checked_sub(event.reverse_ends)
            .ok_or_else(|| Error::Sample("reverse coverage end exceeds active depth".into()))?;
        forward_depth = forward_depth
            .checked_add(event.forward_starts)
            .ok_or_else(|| Error::Sample("forward coverage depth overflow".into()))?;
        reverse_depth = reverse_depth
            .checked_add(event.reverse_starts)
            .ok_or_else(|| Error::Sample("reverse coverage depth overflow".into()))?;
        previous = Some(position);
    }

    if forward_depth != 0 || reverse_depth != 0 {
        return Err(Error::Sample(
            "coverage sweep ended with active reference segments".into(),
        ));
    }

    Ok(coverage)
}

fn validate_segments(read: &SampleReadEvidence) -> Result<()> {
    if read.alignment.reference_segments.is_empty() {
        return Err(Error::Sample(format!(
            "read {} has no mapped reference coverage",
            read.input_sha256
        )));
    }
    let mut segments = read.alignment.reference_segments.clone();
    segments.sort_by_key(|segment| (segment.start_0based, segment.end_0based_exclusive));
    let mut previous_end = None;
    for segment in segments {
        if segment.start_0based >= segment.end_0based_exclusive {
            return Err(Error::Sample(format!(
                "read {} contains empty or reversed reference coverage",
                read.input_sha256
            )));
        }
        if previous_end.is_some_and(|end| segment.start_0based < end) {
            return Err(Error::Sample(format!(
                "read {} contains overlapping reference segments",
                read.input_sha256
            )));
        }
        previous_end = Some(segment.end_0based_exclusive);
    }
    Ok(())
}

fn push_segment(
    output: &mut Vec<SampleCoverageEvidence>,
    start_0based: usize,
    end_0based_exclusive: usize,
    forward_depth: usize,
    reverse_depth: usize,
) -> Result<()> {
    let read_depth = forward_depth
        .checked_add(reverse_depth)
        .ok_or_else(|| Error::Sample("coverage read depth overflow".into()))?;
    if read_depth == 0 || start_0based >= end_0based_exclusive {
        return Err(Error::Sample("invalid non-empty coverage segment".into()));
    }

    if let Some(previous) = output.last_mut()
        && previous.end_0based_exclusive == start_0based
        && previous.read_depth == read_depth
        && previous.forward_depth == forward_depth
        && previous.reverse_depth == reverse_depth
    {
        previous.end_0based_exclusive = end_0based_exclusive;
        return Ok(());
    }

    output.push(SampleCoverageEvidence {
        start_0based,
        end_0based_exclusive,
        read_depth,
        forward_depth,
        reverse_depth,
    });
    Ok(())
}

#[cfg(test)]
mod tests {
    use crate::model::alignment::{Orientation, ReferenceSegment};
    use crate::model::sample_evidence::{SampleReadAlignmentEvidence, SampleReadEvidence};
    use crate::model::signal::TraceIntegrity;

    use super::*;

    fn read(id: &str, orientation: Orientation, segments: &[(usize, usize)]) -> SampleReadEvidence {
        SampleReadEvidence {
            input_name: format!("{id}.ab1"),
            input_sha256: id.into(),
            integrity: TraceIntegrity {
                ploc_count: 20,
                vendor_primary_count: None,
                vendor_quality_count: None,
                minimum_ploc_spacing: Some(4),
                median_ploc_spacing: Some(4.0),
                maximum_ploc_spacing: Some(4),
                clipped_channel_samples: 0,
                maximum_to_median_event_signal_ratio: Some(1.0),
            },
            alignment: SampleReadAlignmentEvidence {
                orientation,
                callable_bases: 20,
                identity: 1.0,
                gap_opens: 0,
                unresolved_bases: 0,
                reference_segments: segments
                    .iter()
                    .map(|&(start_0based, end_0based_exclusive)| ReferenceSegment {
                        start_0based,
                        end_0based_exclusive,
                    })
                    .collect(),
                wraps_origin: segments.len() > 1,
            },
        }
    }

    #[test]
    fn run_length_encodes_total_and_orientation_depth() -> Result<()> {
        let reads = vec![
            read("a", Orientation::Forward, &[(0, 10)]),
            read("b", Orientation::Reverse, &[(5, 15)]),
            read("c", Orientation::Forward, &[(15, 20)]),
        ];

        assert_eq!(
            summarize(&reads)?,
            vec![
                SampleCoverageEvidence {
                    start_0based: 0,
                    end_0based_exclusive: 5,
                    read_depth: 1,
                    forward_depth: 1,
                    reverse_depth: 0,
                },
                SampleCoverageEvidence {
                    start_0based: 5,
                    end_0based_exclusive: 10,
                    read_depth: 2,
                    forward_depth: 1,
                    reverse_depth: 1,
                },
                SampleCoverageEvidence {
                    start_0based: 10,
                    end_0based_exclusive: 15,
                    read_depth: 1,
                    forward_depth: 0,
                    reverse_depth: 1,
                },
                SampleCoverageEvidence {
                    start_0based: 15,
                    end_0based_exclusive: 20,
                    read_depth: 1,
                    forward_depth: 1,
                    reverse_depth: 0,
                },
            ]
        );
        Ok(())
    }

    #[test]
    fn merges_adjacent_runs_with_identical_depth() -> Result<()> {
        let reads = vec![
            read("a", Orientation::Forward, &[(0, 5)]),
            read("b", Orientation::Forward, &[(5, 10)]),
        ];

        assert_eq!(
            summarize(&reads)?,
            vec![SampleCoverageEvidence {
                start_0based: 0,
                end_0based_exclusive: 10,
                read_depth: 1,
                forward_depth: 1,
                reverse_depth: 0,
            }]
        );
        Ok(())
    }

    #[test]
    fn keeps_origin_wrapping_segments_as_linearized_coverage() -> Result<()> {
        let reads = vec![read("a", Orientation::Reverse, &[(90, 100), (0, 10)])];

        assert_eq!(
            summarize(&reads)?,
            vec![
                SampleCoverageEvidence {
                    start_0based: 0,
                    end_0based_exclusive: 10,
                    read_depth: 1,
                    forward_depth: 0,
                    reverse_depth: 1,
                },
                SampleCoverageEvidence {
                    start_0based: 90,
                    end_0based_exclusive: 100,
                    read_depth: 1,
                    forward_depth: 0,
                    reverse_depth: 1,
                },
            ]
        );
        Ok(())
    }

    #[test]
    fn rejects_missing_reference_segments() {
        let reads = vec![read("a", Orientation::Forward, &[])];
        assert!(summarize(&reads).is_err());
    }

    #[test]
    fn rejects_overlapping_segments_within_one_read() {
        let reads = vec![read("a", Orientation::Forward, &[(0, 10), (5, 15)])];
        assert!(summarize(&reads).is_err());
    }
}
