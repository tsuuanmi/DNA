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
//!
//! A sample, a region of the zones, or a zone that no check item touches
//! passes without review; [`passes`] counts them over a set of samples, with
//! the passed ones that carry an unreviewed block, het or not.

use std::collections::BTreeSet;

use serde::Serialize;

use crate::compare::{self, HV_REGIONS};
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
    /// Each unreviewed block's zone region and whether it is het.
    pub(crate) escapes: Vec<Escape>,
}

/// One unreviewed difference block, as the pass measures count it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct Escape {
    /// The zone region that holds the block's first anchor, when one does.
    pub(crate) region: Option<usize>,
    /// Whether the truth writes a heteroplasmy in the block ([`het`]).
    pub(crate) het: bool,
}

/// Whether the truth writes a heteroplasmy in `block`: an IUPAC code of two
/// or three bases (`R`, `Y`, `K`, `M`, `S`, `W`, `B`, `D`, `H`, `V`); `N`,
/// any base, is not one.
pub(crate) fn het(block: &Block) -> bool {
    block.truth.char_indices().any(|(index, symbol)| {
        matches!(
            compare::bases(&block.truth[index..index + symbol.len_utf8()]).len(),
            2 | 3
        )
    })
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
        escapes: unreviewed
            .iter()
            .map(|block| Escape {
                region: zones.and_then(|zones| zones.region(block.start)),
                het: het(block),
            })
            .collect(),
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

/// Passed samples or regions that carry an unreviewed block: in all, with a
/// block that is not het, and with het blocks only.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize)]
pub(crate) struct Escaped {
    pub(crate) total: usize,
    pub(crate) non_het: usize,
    pub(crate) het_only: usize,
}

impl Escaped {
    /// Counts one passed sample or region with its unreviewed `escapes`.
    fn add<'e>(&mut self, escapes: impl IntoIterator<Item = &'e Escape>) {
        let (mut any, mut non_het) = (false, false);
        for escape in escapes {
            any = true;
            non_het |= !escape.het;
        }
        self.total += usize::from(any);
        self.non_het += usize::from(non_het);
        self.het_only += usize::from(any && !non_het);
    }
}

/// What passes without review over a set of samples: the samples, the
/// regions of the zones, and the zones that no check item touches.
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub(crate) struct Passes {
    /// The samples counted.
    pub(crate) samples: usize,
    /// Samples without a check item (auto-pass).
    pub(crate) auto_pass: AutoPass,
    /// One region per sample and region of the zones.
    pub(crate) regions: Regions,
    /// One zone per sample and zone instance.
    pub(crate) zones: Checked,
}

/// The samples without a check item, and those of them with an unreviewed
/// block.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
pub(crate) struct AutoPass {
    pub(crate) samples: usize,
    pub(crate) unreviewed: Escaped,
}

/// Units in all, and those a check item touches.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
pub(crate) struct Checked {
    pub(crate) total: usize,
    pub(crate) checked: usize,
}

/// Regions in all, those with a check item in one of their zones, and the
/// others (passed) that hold an unreviewed block's first anchor; overall and
/// per region.
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub(crate) struct Regions {
    pub(crate) total: usize,
    pub(crate) checked: usize,
    pub(crate) unreviewed: Escaped,
    pub(crate) by_region: Vec<Region>,
}

/// One region of the zones over the samples counted.
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub(crate) struct Region {
    pub(crate) name: String,
    pub(crate) total: usize,
    pub(crate) checked: usize,
    pub(crate) unreviewed: Escaped,
}

/// What passes without review over `workloads`, in the regions and zone
/// instances of `zones`; without zones there are no regions or zones.
pub(crate) fn passes<'w>(
    workloads: impl IntoIterator<Item = &'w Workload>,
    zones: Option<&Zones>,
) -> Passes {
    let mut by_region = zones
        .into_iter()
        .flat_map(Zones::names)
        .map(|region| Region {
            name: region.to_owned(),
            total: 0,
            checked: 0,
            unreviewed: Escaped::default(),
        })
        .collect::<Vec<_>>();
    let mut samples = 0;
    let mut auto_pass = AutoPass {
        samples: 0,
        unreviewed: Escaped::default(),
    };
    let mut zones_checked = 0;
    for workload in workloads {
        samples += 1;
        zones_checked += workload.zones.len();
        if workload.checks == 0 {
            auto_pass.samples += 1;
            auto_pass.unreviewed.add(&workload.escapes);
        }
        for (index, region) in by_region.iter_mut().enumerate() {
            region.total += 1;
            if workload.zones.iter().any(|&(found, _)| found == index) {
                region.checked += 1;
            } else {
                region.unreviewed.add(
                    workload
                        .escapes
                        .iter()
                        .filter(|escape| escape.region == Some(index)),
                );
            }
        }
    }
    let unreviewed = Escaped {
        total: by_region.iter().map(|region| region.unreviewed.total).sum(),
        non_het: by_region
            .iter()
            .map(|region| region.unreviewed.non_het)
            .sum(),
        het_only: by_region
            .iter()
            .map(|region| region.unreviewed.het_only)
            .sum(),
    };
    Passes {
        samples,
        auto_pass,
        regions: Regions {
            total: by_region.iter().map(|region| region.total).sum(),
            checked: by_region.iter().map(|region| region.checked).sum(),
            unreviewed,
            by_region,
        },
        zones: Checked {
            total: samples * zones.map_or(0, |zones| zones.instances().len()),
            checked: zones_checked,
        },
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
    fn a_block_is_het_when_the_truth_writes_a_two_or_three_base_code() {
        let het = |truth: &str| het(&Block::fixture(1, 1, "A", truth));
        assert!(het("ACY") && het("R") && het("B") && het("V"));
        assert!(!het("ACGT") && !het("N") && !het(""));
    }

    #[test]
    fn passes_count_samples_regions_and_zones_without_a_check_item() {
        let windows = windows(10);
        // HV3 holds no repeat window here, so it has three zone kinds.
        let zones = Zones::new(
            &[(16024, 16365), (73, 340), (438, 576)],
            &[(16184, 16193), (303, 315)],
            crate::zones::Widths {
                after: 20,
                edge: 40,
                overrun: 2,
            },
        );
        let of = |checks: &[(i64, i64)], blocks: &[Block]| {
            of(checks, blocks, &windows, &[], Some(&zones))
        };
        let workloads = [
            // Auto-pass with a non-het block in HV1 and a het one in HV2.
            of(
                &[],
                &[
                    Block::fixture(16100, 16100, "A", "C"),
                    Block::fixture(100, 100, "A", "Y"),
                ],
            ),
            // Auto-pass with a het block only.
            of(&[], &[Block::fixture(200, 200, "C", "R")]),
            // Auto-pass without a difference.
            of(&[], &[]),
            // HV1 checked, its block reviewed; HV2 and HV3 pass with a het
            // and an `N` block, and a block outside every region counts for
            // none.
            of(
                &[(16030, 16030)],
                &[
                    Block::fixture(16030, 16030, "A", "C"),
                    Block::fixture(150, 150, "A", "Y"),
                    Block::fixture(500, 500, "A", "N"),
                    Block::fixture(600, 600, "A", "C"),
                ],
            ),
            // HV2 checked: its unreviewed block is in a checked region.
            of(&[(310, 310)], &[Block::fixture(120, 120, "A", "C")]),
        ];
        let passes = passes(&workloads, Some(&zones));
        let escaped = |total, non_het, het_only| Escaped {
            total,
            non_het,
            het_only,
        };
        assert_eq!(passes.samples, 5);
        assert_eq!(
            passes.auto_pass,
            AutoPass {
                samples: 3,
                unreviewed: escaped(2, 1, 1)
            }
        );
        let region = |region: &str, checked, unreviewed| Region {
            name: region.to_owned(),
            total: 5,
            checked,
            unreviewed,
        };
        assert_eq!(
            passes.regions,
            Regions {
                total: 15,
                checked: 2,
                unreviewed: escaped(5, 2, 3),
                by_region: vec![
                    region("HV1", 1, escaped(1, 1, 0)),
                    region("HV2", 1, escaped(3, 0, 3)),
                    region("HV3", 0, escaped(1, 1, 0)),
                ],
            }
        );
        assert_eq!(
            passes.zones,
            Checked {
                total: 5 * 13,
                checked: 2
            }
        );
        let bare = super::passes(&workloads, None);
        assert_eq!((bare.regions.total, bare.zones.total), (0, 0));
        assert_eq!(bare.auto_pass, passes.auto_pass);
    }

    #[test]
    fn the_wilson_bound_holds_its_events() {
        assert_eq!(wilson_upper(0, 0), 1.0);
        let bound = wilson_upper(5, 100);
        assert!(bound > 0.05 && bound < 0.12, "{bound}");
    }
}
