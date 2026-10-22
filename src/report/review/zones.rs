//! The review's zones: the reportable positions as maximal stretches of one
//! kind a reviewer opens whole, from the profile's repeat windows, the callable
//! ranges, and the reads' masked stretches, with the check items each touches.

use super::ranges::inside;
use crate::model::review_result::{ItemKind, ItemResult, ReviewReadResult, ZoneKind, ZoneResult};

/// The reportable positions within `regions` as zones a reviewer opens
/// whole. Each position takes the first kind that holds, in [`ZoneKind`]
/// order: inside a repeat window of the profile (1-based inclusive, with its
/// name), outside the `callable` ranges, inside a read's masked stretch past a
/// repeat, inside another masked stretch, else clean. A zone is a maximal
/// stretch of one kind (and window) within one region; it counts the `check`
/// items touching it. Every boundary comes from the profile and the sample's
/// reads, none from a width.
pub(super) fn zones(
    regions: &[[usize; 2]],
    windows: &[(usize, usize, &str)],
    callable: &[[usize; 2]],
    reads: &[ReviewReadResult],
    items: &[ItemResult],
) -> Vec<ZoneResult> {
    let masked = reads
        .iter()
        .flat_map(|read| read.masked.iter().flatten())
        .collect::<Vec<_>>();
    let masks = |position: usize, after_repeat: bool| {
        masked.iter().any(|span| {
            (!after_repeat || span.after_repeat) && span.first <= position && position <= span.last
        })
    };
    let mut zones: Vec<ZoneResult> = Vec::new();
    for &[start, end] in regions {
        for position in start..=end {
            let window = windows
                .iter()
                .find(|(first, last, _)| *first <= position && position <= *last);
            let (kind, name) = match window {
                Some(&(_, _, name)) => (ZoneKind::Repeat, Some(name)),
                None if !inside(callable, position) => (ZoneKind::NotCallable, None),
                None if masks(position, true) => (ZoneKind::Dephased, None),
                None if masks(position, false) => (ZoneKind::Masked, None),
                None => (ZoneKind::Clean, None),
            };
            match zones.last_mut() {
                Some(zone)
                    if position != start && zone.kind == kind && zone.name.as_deref() == name =>
                {
                    zone.last = position;
                }
                _ => zones.push(ZoneResult {
                    kind,
                    first: position,
                    last: position,
                    name: name.map(str::to_owned),
                    checks: 0,
                }),
            }
        }
    }
    for zone in &mut zones {
        zone.checks = items
            .iter()
            .filter(|item| {
                item.kind == ItemKind::Check
                    && item.first.position <= zone.last
                    && zone.first <= item.last.position
            })
            .count();
    }
    zones
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::result::PhaseStateResult as PhaseState;
    use crate::model::review_result::{LabelResult, MaskedSpanResult, ReviewStatus, ZoneCounts};
    use crate::report::review::labels::summarize;
    use crate::report::review::test_support::{F, item};

    // SRS-API-014
    #[test]
    fn zones_follow_the_profile_windows_and_what_the_reads_show() {
        let span = |first, last, state, after_repeat| MaskedSpanResult {
            first,
            last,
            state,
            after_repeat,
        };
        let reads = [ReviewReadResult {
            name: "a".to_owned(),
            sha256: String::new(),
            strand: Some(F),
            admitted: true,
            rejection: None,
            masked: Some(vec![
                span(6, 6, PhaseState::Dephased, true),
                span(7, 7, PhaseState::Mixed, false),
            ]),
            unreliable: Some(Vec::new()),
            mixtures: Some(Vec::new()),
        }];
        let check = item(ItemKind::Check, vec![LabelResult::Undecided]);
        let variant = item(ItemKind::Variant, Vec::new());
        let zones = zones(
            &[[1, 10], [11, 13]],
            &[(4, 5, "W")],
            &[[1, 8], [11, 13]],
            &reads,
            &[check, variant],
        );
        let zone = |kind, first, last, name: Option<&str>, checks| ZoneResult {
            kind,
            first,
            last,
            name: name.map(str::to_owned),
            checks,
        };
        assert_eq!(
            zones,
            [
                // Only the check item at position 1 counts.
                zone(ZoneKind::Clean, 1, 3, None, 1),
                zone(ZoneKind::Repeat, 4, 5, Some("W"), 0),
                zone(ZoneKind::Dephased, 6, 6, None, 0),
                zone(ZoneKind::Masked, 7, 7, None, 0),
                zone(ZoneKind::Clean, 8, 8, None, 0),
                zone(ZoneKind::NotCallable, 9, 10, None, 0),
                // A zone never crosses into the next region.
                zone(ZoneKind::Clean, 11, 13, None, 0),
            ]
        );
        let summary = summarize(&[], ReviewStatus::Called, &[[1, 13]], &zones, &[], 0);
        assert_eq!(summary.zones, ZoneCounts { total: 7, check: 1 });
    }
}
