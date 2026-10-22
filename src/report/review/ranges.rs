//! Reference ranges of the review: 1-based inclusive `[first, last]`
//! pairs, merged, clipped to the reportable regions, and tested for a
//! position.

use dna_align::model::alignment::ReferenceSegment;

/// A 0-based half-open reference segment as a pair.
pub(super) const fn span(segment: &ReferenceSegment) -> (usize, usize) {
    (segment.start_0based, segment.end_0based_exclusive)
}

/// 0-based half-open segments as 1-based inclusive ranges, empty ones left
/// out.
pub(super) fn segments(segments: impl IntoIterator<Item = (usize, usize)>) -> Vec<[usize; 2]> {
    segments
        .into_iter()
        .filter(|(start, end)| start < end)
        .map(|(start, end)| [start + 1, end])
        .collect()
}

/// 1-based inclusive ranges sorted and merged where they overlap or touch.
pub(super) fn merge(mut ranges: Vec<[usize; 2]>) -> Vec<[usize; 2]> {
    ranges.sort_unstable();
    let mut merged: Vec<[usize; 2]> = Vec::with_capacity(ranges.len());
    for [start, end] in ranges {
        match merged.last_mut() {
            Some(last) if start <= last[1] + 1 => last[1] = last[1].max(end),
            _ => merged.push([start, end]),
        }
    }
    merged
}

/// The profile's regions within a reference of `length` bases, merged.
pub(super) fn clip(regions: &[[usize; 2]], length: usize) -> Vec<[usize; 2]> {
    merge(
        regions
            .iter()
            .map(|&[start, end]| [start.max(1), end.min(length)])
            .filter(|[start, end]| start <= end)
            .collect(),
    )
}

/// Whether merged `ranges` hold `position`.
pub(super) fn inside(ranges: &[[usize; 2]], position: usize) -> bool {
    ranges
        .iter()
        .any(|&[start, end]| start <= position && position <= end)
}

/// The first and last positions of `first..=last` inside merged `regions`.
pub(super) fn within(regions: &[[usize; 2]], first: usize, last: usize) -> Option<[usize; 2]> {
    let clipped = regions
        .iter()
        .map(|&[start, end]| [start.max(first), end.min(last)])
        .filter(|[start, end]| start <= end)
        .collect::<Vec<_>>();
    Some([clipped.first()?[0], clipped.last()?[1]])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ranges_merge_clip_and_contain() {
        assert_eq!(segments([(0, 10), (5, 5), (12, 20)]), [[1, 10], [13, 20]]);
        assert_eq!(
            merge(vec![[13, 20], [1, 10], [11, 11], [30, 31]]),
            [[1, 11], [13, 20], [30, 31]]
        );
        assert_eq!(
            clip(&[[16024, 16365], [73, 340], [1, 50_000]], 400),
            [[1, 400]]
        );
        assert_eq!(
            clip(&[[16024, 16365], [73, 340], [438, 576]], 16569),
            [[73, 340], [438, 576], [16024, 16365]]
        );
        assert_eq!(clip(&[[500, 600]], 400), Vec::<[usize; 2]>::new());
        let regions = [[73, 340], [438, 576]];
        assert!(inside(&regions, 73) && inside(&regions, 576));
        assert!(!inside(&regions, 341) && !inside(&regions, 437));
        assert_eq!(within(&regions, 300, 500), Some([300, 500]));
        assert_eq!(within(&regions, 60, 80), Some([73, 80]));
        assert_eq!(within(&regions, 341, 437), None);
    }
}
