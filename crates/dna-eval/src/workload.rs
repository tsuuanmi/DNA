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
//! the passed ones that carry an unreviewed block by the [`Class`] of its
//! blocks.

use std::collections::BTreeSet;

use serde::Serialize;

use crate::compare::{self, HV_REGIONS};
use crate::haplotype::Block;
use crate::recheck;
use crate::zones::{Instance, Zones};

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
    /// Each unreviewed block's zone instance and class.
    pub(crate) escapes: Vec<Escape>,
}

/// One unreviewed difference block, as the zone and pass measures count it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct Escape {
    /// The block's first and last anchors, 1-based.
    pub(crate) first: i64,
    pub(crate) last: i64,
    /// The zone instance that holds the block's first anchor, when one does.
    pub(crate) zone: Option<Instance>,
    /// What the block is, as the pass measures split it.
    pub(crate) class: Class,
}

/// What an unreviewed difference block is, as the pass measures split it,
/// from the least to the most severe. A block is `Recheck` when the re-check
/// list names it, else `TruthN` when the truth writes `N` in it or every
/// difference is an insertion past the truth's interval (the truth records
/// nothing there, [`Block::past_interval`]), else `Het`
/// when the truth writes a heteroplasmy in it ([`het`]), else `NonHet`: a
/// DNA error that neither a heteroplasmy, a truth `N`, nor the re-check list
/// explains. A passed unit counts under its most severe block.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) enum Class {
    Recheck,
    Het,
    TruthN,
    NonHet,
}

impl Class {
    /// The class's label in `unreviewed.tsv`.
    pub(crate) const fn label(self) -> &'static str {
        match self {
            Self::Recheck => "recheck",
            Self::Het => "het",
            Self::TruthN => "truth_n",
            Self::NonHet => "non_het",
        }
    }

    /// The class of `block`, listed when it overlaps one of `listed`, its
    /// sample's spans on the re-check list.
    pub(crate) fn of(block: &Block, listed: &[(i64, i64)]) -> Self {
        if recheck::listed(block, listed) {
            Self::Recheck
        } else if block.truth.contains('N') || block.past_interval() {
            Self::TruthN
        } else if het(block) {
            Self::Het
        } else {
            Self::NonHet
        }
    }
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
/// 1-based), its `haplotype_lod` difference blocks, and its spans on the
/// re-check list.
pub(crate) fn of(
    checks: &[(i64, i64)],
    blocks: &[Block],
    windows: &[Window],
    forced: &[(i64, i64)],
    zones: Option<&Zones>,
    listed: &[(i64, i64)],
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
    Workload {
        zones: touched_zones.iter().flatten().copied().collect(),
        spanning: touched_zones.iter().filter(|zones| zones.len() > 1).count(),
        escapes: unreviewed
            .iter()
            .map(|block| Escape {
                first: block.start,
                last: block.end,
                zone: zones.and_then(|zones| zones.instance(block.start)),
                class: Class::of(block, listed),
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

/// Passed samples, regions, or zones that carry an unreviewed block: in all,
/// and by the most severe [`Class`] of their blocks: with a `NonHet` block,
/// else with a truth `N`, else het only, else on the re-check list only.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize)]
pub(crate) struct Escaped {
    pub(crate) total: usize,
    pub(crate) non_het: usize,
    pub(crate) het_only: usize,
    pub(crate) truth_n_only: usize,
    pub(crate) recheck_only: usize,
}

impl Escaped {
    /// Counts one passed unit with its unreviewed `escapes`.
    fn add<'e>(&mut self, escapes: impl IntoIterator<Item = &'e Escape>) {
        let Some(class) = escapes.into_iter().map(|escape| escape.class).max() else {
            return;
        };
        self.total += 1;
        *match class {
            Class::NonHet => &mut self.non_het,
            Class::TruthN => &mut self.truth_n_only,
            Class::Het => &mut self.het_only,
            Class::Recheck => &mut self.recheck_only,
        } += 1;
    }

    /// The sum of `counts`.
    pub(crate) fn sum<'e>(counts: impl IntoIterator<Item = &'e Self>) -> Self {
        counts.into_iter().fold(Self::default(), |sum, count| Self {
            total: sum.total + count.total,
            non_het: sum.non_het + count.non_het,
            het_only: sum.het_only + count.het_only,
            truth_n_only: sum.truth_n_only + count.truth_n_only,
            recheck_only: sum.recheck_only + count.recheck_only,
        })
    }
}

/// What passes without review over a set of samples: the samples, the
/// regions of the zones, and the zone instances that no check item touches.
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub(crate) struct Passes {
    /// The samples counted.
    pub(crate) samples: usize,
    /// Samples without a check item (auto-pass).
    pub(crate) auto_pass: AutoPass,
    /// One region per sample and region of the zones.
    pub(crate) regions: Regions,
    /// One zone per sample and zone instance.
    pub(crate) zones: ZonePasses,
}

/// The samples without a check item, and those of them with an unreviewed
/// block.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
pub(crate) struct AutoPass {
    pub(crate) samples: usize,
    pub(crate) unreviewed: Escaped,
}

/// Regions in all, those with a check item in one of their zones, and the
/// others (passed) that hold an unreviewed block's first anchor; overall and
/// per region.
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub(crate) struct Regions {
    pub(crate) total: usize,
    pub(crate) checked: usize,
    pub(crate) unreviewed: Escaped,
    pub(crate) by_region: Vec<Unit>,
}

/// Zone instances in all, those a check item touches, and the others
/// (passed) that hold an unreviewed block's first anchor; overall and per
/// instance, in the order of [`Zones::instances`].
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub(crate) struct ZonePasses {
    pub(crate) total: usize,
    pub(crate) checked: usize,
    pub(crate) unreviewed: Escaped,
    pub(crate) by_zone: Vec<Unit>,
}

/// One region or zone instance over the samples counted.
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub(crate) struct Unit {
    pub(crate) name: String,
    pub(crate) total: usize,
    pub(crate) checked: usize,
    pub(crate) unreviewed: Escaped,
}

impl Unit {
    fn new(name: String) -> Self {
        Self {
            name,
            total: 0,
            checked: 0,
            unreviewed: Escaped::default(),
        }
    }

    /// Counts one sample: checked when `checked`, else passed with the
    /// unreviewed `escapes` that start in the unit.
    fn count<'e>(&mut self, checked: bool, escapes: impl IntoIterator<Item = &'e Escape>) {
        self.total += 1;
        if checked {
            self.checked += 1;
        } else {
            self.unreviewed.add(escapes);
        }
    }
}

/// What passes without review over `workloads`, in the regions and zone
/// instances of `zones`; without zones there are no regions or zones.
pub(crate) fn passes<'w>(
    workloads: impl IntoIterator<Item = &'w Workload>,
    zones: Option<&Zones>,
) -> Passes {
    let instances = zones.map(Zones::instances).unwrap_or_default();
    let mut by_region = zones
        .into_iter()
        .flat_map(Zones::names)
        .map(|name| Unit::new(name.to_owned()))
        .collect::<Vec<_>>();
    let mut by_zone = zones
        .into_iter()
        .flat_map(|zones| instances.iter().map(|&instance| zones.label(instance)))
        .map(Unit::new)
        .collect::<Vec<_>>();
    let mut samples = 0;
    let mut auto_pass = AutoPass {
        samples: 0,
        unreviewed: Escaped::default(),
    };
    for workload in workloads {
        samples += 1;
        if workload.checks == 0 {
            auto_pass.samples += 1;
            auto_pass.unreviewed.add(&workload.escapes);
        }
        for (index, region) in by_region.iter_mut().enumerate() {
            region.count(
                workload.zones.iter().any(|&(found, _)| found == index),
                workload
                    .escapes
                    .iter()
                    .filter(|escape| escape.zone.is_some_and(|(found, _)| found == index)),
            );
        }
        for (&instance, zone) in instances.iter().zip(&mut by_zone) {
            zone.count(
                workload.zones.contains(&instance),
                workload
                    .escapes
                    .iter()
                    .filter(|escape| escape.zone == Some(instance)),
            );
        }
    }
    Passes {
        samples,
        auto_pass,
        regions: Regions {
            total: by_region.iter().map(|region| region.total).sum(),
            checked: by_region.iter().map(|region| region.checked).sum(),
            unreviewed: Escaped::sum(by_region.iter().map(|region| &region.unreviewed)),
            by_region,
        },
        zones: ZonePasses {
            total: by_zone.iter().map(|zone| zone.total).sum(),
            checked: by_zone.iter().map(|zone| zone.checked).sum(),
            unreviewed: Escaped::sum(by_zone.iter().map(|zone| &zone.unreviewed)),
            by_zone,
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
    use crate::zones::Kind;

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
            &[],
        );
        assert_eq!(workload.checks, 1);
        assert_eq!(workload.zones, BTreeSet::from([(0, Kind::Start)]));
        assert_eq!(workload.spanning, 0);
        assert_eq!(
            workload.escapes,
            [Escape {
                first: 16040,
                last: 16041,
                zone: Some((0, Kind::Start)),
                class: Class::NonHet
            }]
        );
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
            &[],
        );
        assert!(workload.zones.is_empty());
        assert_eq!(workload.escapes[0].zone, None);
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
            of(checks, blocks, &windows, &[], Some(&zones), &[])
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
            // and a truth `N` block, and a block outside every region counts
            // for none.
            of(
                &[(16030, 16030)],
                &[
                    Block::fixture(16030, 16030, "A", "C"),
                    Block::fixture(150, 150, "A", "Y"),
                    Block::fixture(500, 500, "A", "N"),
                    Block::fixture(600, 600, "A", "C"),
                ],
            ),
            // HV2 checked at its window: an unreviewed block in the window
            // is in a checked zone, one in the interior passes its zone but
            // not its region.
            of(
                &[(310, 310)],
                &[
                    Block::fixture(120, 120, "A", "C"),
                    Block::fixture(305, 305, "A", "C"),
                ],
            ),
        ];
        let passes = passes(&workloads, Some(&zones));
        let escaped = |non_het, het_only, truth_n_only| Escaped {
            total: non_het + het_only + truth_n_only,
            non_het,
            het_only,
            truth_n_only,
            recheck_only: 0,
        };
        assert_eq!(passes.samples, 5);
        assert_eq!(
            passes.auto_pass,
            AutoPass {
                samples: 3,
                unreviewed: escaped(1, 1, 0)
            }
        );
        let unit = |name: &str, checked, unreviewed| Unit {
            name: name.to_owned(),
            total: 5,
            checked,
            unreviewed,
        };
        assert_eq!(
            passes.regions,
            Regions {
                total: 15,
                checked: 2,
                unreviewed: escaped(1, 3, 1),
                by_region: vec![
                    unit("HV1", 1, escaped(1, 0, 0)),
                    unit("HV2", 1, escaped(0, 3, 0)),
                    unit("HV3", 0, escaped(0, 0, 1)),
                ],
            }
        );
        let none = escaped(0, 0, 0);
        assert_eq!(
            passes.zones,
            ZonePasses {
                total: 5 * 13,
                checked: 2,
                unreviewed: escaped(2, 3, 1),
                by_zone: vec![
                    unit("window@HV1", 0, none),
                    unit("after@HV1", 0, none),
                    unit("start@HV1", 1, none),
                    unit("end@HV1", 0, none),
                    unit("interior@HV1", 0, escaped(1, 0, 0)),
                    unit("window@HV2", 1, none),
                    unit("after@HV2", 0, none),
                    unit("start@HV2", 0, escaped(0, 1, 0)),
                    unit("end@HV2", 0, none),
                    unit("interior@HV2", 0, escaped(1, 2, 0)),
                    unit("start@HV3", 0, none),
                    unit("end@HV3", 0, none),
                    unit("interior@HV3", 0, escaped(0, 0, 1)),
                ],
            }
        );
        let bare = super::passes(&workloads, None);
        assert_eq!((bare.regions.total, bare.zones.total), (0, 0));
        assert_eq!(bare.auto_pass, passes.auto_pass);
    }

    #[test]
    fn a_block_takes_one_class_and_a_unit_its_most_severe() {
        use Class::{Het, NonHet, Recheck, TruthN};
        let class = |truth: &str, listed: &[(i64, i64)]| {
            Class::of(&Block::fixture(100, 101, "AC", truth), listed)
        };
        // The re-check list comes first, then a truth `N`, then het.
        assert_eq!(class("AC", &[(101, 105)]), Recheck);
        assert_eq!(class("NY", &[(90, 100)]), Recheck);
        assert_eq!(class("NY", &[(102, 105)]), TruthN);
        assert_eq!(class("AY", &[]), Het);
        assert_eq!(class("GG", &[]), NonHet);
        // An insertion past the truth's interval: the truth records nothing
        // there. A block that also differs at the last base, holds more than
        // it, or stays inside the interval is DNA's.
        let open = |start, dna: &str, truth: &str| {
            Class::of(&Block::open_fixture(start, 573, dna, truth), &[])
        };
        assert_eq!(open(573, "CCC", "C"), TruthN);
        assert_eq!(open(573, "C", "CC"), TruthN);
        assert_eq!(open(573, "ACC", "C"), NonHet);
        assert_eq!(open(572, "CCCC", "CC"), NonHet);
        assert_eq!(
            Class::of(&Block::fixture(573, 573, "CCC", "C"), &[]),
            NonHet
        );
        // A passed unit counts once, under its most severe block.
        let escape = |class| Escape {
            first: 1,
            last: 1,
            zone: None,
            class,
        };
        let count = |classes: &[Class]| {
            let mut escaped = Escaped::default();
            escaped.add(&classes.iter().copied().map(escape).collect::<Vec<_>>());
            escaped
        };
        let only = |non_het, het_only, truth_n_only, recheck_only| Escaped {
            total: 1,
            non_het,
            het_only,
            truth_n_only,
            recheck_only,
        };
        assert_eq!(count(&[Recheck, Het, TruthN, NonHet]), only(1, 0, 0, 0));
        assert_eq!(count(&[Recheck, Het, TruthN]), only(0, 0, 1, 0));
        assert_eq!(count(&[Recheck, Het]), only(0, 1, 0, 0));
        assert_eq!(count(&[Recheck, Recheck]), only(0, 0, 0, 1));
        assert_eq!(count(&[]), Escaped::default());
        // A listed block passes as re-check only; unlisted it is non-het.
        let windows = windows(10);
        let block = [Block::fixture(573, 573, "CCCCCCCA", "CCCCCCA")];
        let listed = of(&[], &block, &windows, &[], None, &[(573, 573)]);
        let unlisted = of(&[], &block, &windows, &[], None, &[(16057, 16057)]);
        assert_eq!(listed.escapes[0].class, Recheck);
        assert_eq!(unlisted.escapes[0].class, NonHet);
        assert_eq!(listed.unreviewed, unlisted.unreviewed);
        let auto_pass = |workload: &Workload| passes([workload], None).auto_pass.unreviewed;
        assert_eq!(auto_pass(&listed), only(0, 0, 0, 1));
        assert_eq!(auto_pass(&unlisted), only(1, 0, 0, 0));
    }

    #[test]
    fn the_wilson_bound_holds_its_events() {
        assert_eq!(wilson_upper(0, 0), 1.0);
        let bound = wilson_upper(5, 100);
        assert!(bound > 0.05 && bound < 0.12, "{bound}");
    }
}
