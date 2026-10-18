//! Review workload: how much of each sample DNA leaves to a reviewer, and
//! which differences from the truth no reviewer would see.
//!
//! The review items are `dna review`'s own (`SampleAnalysis::review_checks`).
//! The HV regions split into windows of a fixed size aligned to each region's
//! first position, the last window of a region possibly shorter; a reviewer
//! looks at every window a check item touches, so the positions to review are
//! those windows' positions. A `haplotype_lod` difference block that overlaps
//! no check item is unreviewed, an error no reviewer would see. Forced windows
//! are always reviewed; which windows to force is a policy decision, learned
//! on one run and judged on a hold-out or later batches.

use std::collections::BTreeSet;

use crate::compare::HV_REGIONS;
use crate::haplotype::Block;
use crate::zones::{Instance, Kind, Zones};

/// One window of an HV region, 1-based and inclusive.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct Window {
    pub(crate) region: &'static str,
    pub(crate) start: i64,
    pub(crate) end: i64,
}

impl Window {
    /// Whether the window holds any position of `first..=last`.
    pub(crate) const fn overlaps(&self, first: i64, last: i64) -> bool {
        first <= self.end && self.start <= last
    }

    pub(crate) const fn length(&self) -> i64 {
        self.end - self.start + 1
    }
}

/// The windows of `size` positions over every HV region, in region order.
pub(crate) fn windows(size: i64) -> Vec<Window> {
    let size = size.max(1);
    HV_REGIONS
        .iter()
        .flat_map(|&(region, first, last)| {
            (0..)
                .map(move |index| first + index * size)
                .take_while(move |start| *start <= last)
                .map(move |start| Window {
                    region,
                    start,
                    end: (start + size - 1).min(last),
                })
        })
        .collect()
}

/// One sample's review workload.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct Workload {
    /// Check items in the review list.
    pub(crate) checks: usize,
    /// Positions of the windows a check item touches.
    pub(crate) positions: i64,
    /// The same with every forced window reviewed.
    pub(crate) positions_forced: i64,
    /// Difference blocks no check item overlaps.
    pub(crate) unreviewed: Vec<String>,
    /// The same, also leaving out blocks a forced window overlaps.
    pub(crate) unreviewed_forced: Vec<String>,
    /// Per window, whether a check item touches it.
    pub(crate) checked: Vec<bool>,
    /// Per window, whether an unreviewed block overlaps it.
    pub(crate) escaped: Vec<bool>,
    /// Region of each unreviewed block, by its first anchor.
    pub(crate) regions: Vec<&'static str>,
    /// The zone instances a check item touches.
    pub(crate) zones: BTreeSet<Instance>,
    /// Check items that touch more than one zone instance.
    pub(crate) spanning: usize,
    /// Zone kind of each unreviewed block, by its first anchor.
    pub(crate) zone_kinds: Vec<Kind>,
}

/// The workload of one sample from its check items (`first..=last`,
/// 1-based) and its `haplotype_lod` difference blocks.
pub(crate) fn of(
    checks: &[(i64, i64)],
    blocks: &[Block],
    windows: &[Window],
    forced: &[(i64, i64)],
    zones: Option<&Zones>,
) -> Workload {
    let touched = |window: &Window| {
        checks
            .iter()
            .any(|&(first, last)| window.overlaps(first, last))
    };
    let is_forced = |window: &Window| {
        forced
            .iter()
            .any(|&(first, last)| window.overlaps(first, last))
    };
    let checked = windows.iter().map(touched).collect::<Vec<_>>();
    let positions = windows
        .iter()
        .zip(&checked)
        .filter(|(_, checked)| **checked)
        .map(|(window, _)| window.length())
        .sum();
    let positions_forced = windows
        .iter()
        .zip(&checked)
        .filter(|(window, checked)| **checked || is_forced(window))
        .map(|(window, _)| window.length())
        .sum();
    let reviewed = |block: &Block| {
        checks
            .iter()
            .any(|&(first, last)| block.overlaps(first, last))
    };
    let unreviewed = blocks
        .iter()
        .filter(|block| !reviewed(block))
        .collect::<Vec<_>>();
    let unreviewed_forced = unreviewed
        .iter()
        .filter(|block| {
            !forced
                .iter()
                .any(|&(first, last)| block.overlaps(first, last))
        })
        .map(|block| block.text())
        .collect();
    let escaped = windows
        .iter()
        .map(|window| {
            unreviewed
                .iter()
                .any(|block| window.overlaps(block.start, block.end))
        })
        .collect();
    let regions = unreviewed
        .iter()
        .filter_map(|block| {
            HV_REGIONS
                .iter()
                .find(|(_, first, last)| (*first..=*last).contains(&block.start))
                .map(|(region, _, _)| *region)
        })
        .collect();
    let touched_zones = checks
        .iter()
        .map(|&(first, last)| {
            zones
                .map(|zones| zones.touched(first, last))
                .unwrap_or_default()
        })
        .collect::<Vec<_>>();
    let zone_kinds = unreviewed
        .iter()
        .filter_map(|block| zones?.kind(block.start))
        .collect();
    Workload {
        zones: touched_zones.iter().flatten().copied().collect(),
        spanning: touched_zones.iter().filter(|zones| zones.len() > 1).count(),
        zone_kinds,
        checks: checks.len(),
        positions,
        positions_forced,
        unreviewed: unreviewed.iter().map(|block| block.text()).collect(),
        unreviewed_forced,
        checked,
        escaped,
        regions,
    }
}

/// Parses a forced-window list: one `start-end` per line, 1-based and
/// inclusive; blank lines and `#` comments are ignored.
pub(crate) fn parse_forced(text: &str) -> Result<Vec<(i64, i64)>, String> {
    text.lines()
        .enumerate()
        .map(|(index, line)| (index + 1, line.trim()))
        .filter(|(_, line)| !line.is_empty() && !line.starts_with('#'))
        .map(|(number, line)| {
            let window = line.split_once('-').and_then(|(start, end)| {
                Some((
                    start.trim().parse::<i64>().ok()?,
                    end.trim().parse::<i64>().ok()?,
                ))
            });
            match window {
                Some((start, end)) if 1 <= start && start <= end => Ok((start, end)),
                _ => Err(format!("line {number}: expected `start-end`, got `{line}`")),
            }
        })
        .collect()
}

/// The upper bound of the Wilson 95 % interval of `events` in `trials`.
pub(crate) fn wilson_upper(events: usize, trials: usize) -> f64 {
    if trials == 0 {
        return 1.0;
    }
    let z = crate::estimate::Z_95;
    let n = trials as f64;
    let p = events as f64 / n;
    let centre = p + z * z / (2.0 * n);
    let spread = z * (p * (1.0 - p) / n + z * z / (4.0 * n * n)).sqrt();
    ((centre + spread) / (1.0 + z * z / n)).min(1.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn block(start: i64, end: i64) -> Block {
        Block::fixture(start, end, "A", "C")
    }

    #[test]
    fn windows_align_to_each_region_start_and_end_short() {
        let windows = windows(10);
        assert_eq!(
            (windows[0].region, windows[0].start, windows[0].end),
            ("HV1", 16024, 16033)
        );
        let last_hv1 = windows.iter().rfind(|window| window.region == "HV1");
        assert_eq!(
            last_hv1.map(|window| (window.start, window.end)),
            Some((16364, 16365))
        );
        let first_hv2 = windows.iter().find(|window| window.region == "HV2");
        assert_eq!(
            first_hv2.map(|window| (window.start, window.end)),
            Some((73, 82))
        );
        assert_eq!(
            windows.iter().map(Window::length).sum::<i64>(),
            (16365 - 16024 + 1) + (340 - 73 + 1) + (576 - 438 + 1)
        );
    }

    #[test]
    fn a_block_beside_a_check_item_is_unreviewed_and_marks_its_windows() {
        let windows = windows(10);
        let zones = Zones::new(
            &[(16024, 16365)],
            &[(16180, 16193)],
            crate::zones::Widths {
                after: 20,
                edge: 40,
                overrun: 2,
            },
        );
        let workload = of(
            &[(16030, 16031)],
            &[block(16031, 16031), block(16040, 16041)],
            &windows,
            &[],
            Some(&zones),
        );
        assert_eq!(workload.checks, 1);
        assert_eq!(workload.zones, BTreeSet::from([(0, Kind::Start)]));
        assert_eq!(workload.spanning, 0);
        assert_eq!(workload.zone_kinds, [Kind::Start]);
        assert_eq!(workload.positions, 10);
        assert_eq!(workload.unreviewed, ["16040-16041:A/C"]);
        assert_eq!(workload.regions, ["HV1"]);
        assert!(workload.checked[0] && !workload.escaped[0]);
        assert!(workload.escaped[1] && !workload.checked[1]);
    }

    #[test]
    fn forced_windows_count_as_reviewed() {
        let windows = windows(10);
        let workload = of(
            &[],
            &[block(16185, 16186)],
            &windows,
            &[(16184, 16193)],
            None,
        );
        assert!(workload.zones.is_empty() && workload.zone_kinds.is_empty());
        assert_eq!(workload.positions, 0);
        assert_eq!(workload.positions_forced, 10);
        assert_eq!(workload.unreviewed.len(), 1);
        assert!(workload.unreviewed_forced.is_empty());
        assert_eq!(
            parse_forced("# hotspots\n16184-16193\n\n303-340\n"),
            Ok(vec![(16184, 16193), (303, 340)])
        );
        assert!(parse_forced("16193-16184\n").is_err());
    }

    #[test]
    fn the_wilson_bound_holds_its_events() {
        assert_eq!(wilson_upper(0, 0), 1.0);
        let bound = wilson_upper(5, 100);
        assert!(bound > 0.05 && bound < 0.12, "{bound}");
    }
}
