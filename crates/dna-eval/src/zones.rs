//! Review zones: the functional parts of the reportable regions a reviewer
//! opens as a whole, so the workload counts zones rather than positions.
//!
//! Every position of a region has one zone kind, by precedence: a repeat
//! window of the profile's nomenclature windows, then the positions within
//! `after` of a window on either side (the harness does not see a read's
//! direction), then the first and last `edge` positions of the region (the
//! amplicon's start and end), then the interior. A zone instance is a kind in
//! one region, whether or not its positions are contiguous. A check item
//! counts for every instance it touches, except that an item overlapping a
//! repeat window counts its positions within `overrun` of that window for
//! the window: the consensus often decides a repeat as one stretch a
//! position or two past it.

use std::collections::BTreeSet;

use crate::compare::HV_REGIONS;

/// The widths that place the zones, in positions; recorded in `run.json`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct Widths {
    /// Positions on each side of a repeat window that follow it.
    pub(crate) after: i64,
    /// Positions at each end of a region that are its amplicon's start or end.
    pub(crate) edge: i64,
    /// Positions past a repeat window that an item overlapping the window
    /// still counts for the window.
    pub(crate) overrun: i64,
}

/// What part of a region a position lies in, in precedence order.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub(crate) enum Kind {
    Window,
    After,
    Start,
    End,
    Interior,
}

impl Kind {
    pub(crate) const ALL: [Self; 5] = [
        Self::Window,
        Self::After,
        Self::Start,
        Self::End,
        Self::Interior,
    ];

    pub(crate) const fn label(self) -> &'static str {
        match self {
            Self::Window => "window",
            Self::After => "after",
            Self::Start => "start",
            Self::End => "end",
            Self::Interior => "interior",
        }
    }
}

/// One zone instance: a kind in the region of that index.
pub(crate) type Instance = (usize, Kind);

/// One reportable region with the kind of each of its positions.
#[derive(Debug, Clone, PartialEq, Eq)]
struct Region {
    name: String,
    start: i64,
    /// The kind of `start + index`.
    kinds: Vec<Kind>,
    /// The repeat windows inside it, 1-based and inclusive.
    windows: Vec<(i64, i64)>,
}

/// The zones of a profile's reportable regions.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct Zones {
    regions: Vec<Region>,
    overrun: i64,
}

impl Zones {
    /// The zones of `regions` and repeat `windows` (both 1-based and
    /// inclusive), placed with `widths`. A region takes the name of the HV
    /// region that holds its first position, else its bounds.
    pub(crate) fn new(regions: &[(i64, i64)], windows: &[(i64, i64)], widths: Widths) -> Self {
        let regions = regions
            .iter()
            .map(|&(start, end)| {
                let inside = windows
                    .iter()
                    .copied()
                    .filter(|&(first, last)| first <= end && start <= last)
                    .collect::<Vec<_>>();
                let kind = |position: i64| {
                    if inside
                        .iter()
                        .any(|&(first, last)| (first..=last).contains(&position))
                    {
                        Kind::Window
                    } else if inside.iter().any(|&(first, last)| {
                        (first - widths.after..first).contains(&position)
                            || (last + 1..=last + widths.after).contains(&position)
                    }) {
                        Kind::After
                    } else if position < start + widths.edge {
                        Kind::Start
                    } else if position > end - widths.edge {
                        Kind::End
                    } else {
                        Kind::Interior
                    }
                };
                let name = HV_REGIONS
                    .iter()
                    .find(|(_, first, last)| (*first..=*last).contains(&start))
                    .map_or_else(
                        || format!("{start}-{end}"),
                        |(name, _, _)| (*name).to_owned(),
                    );
                Region {
                    name,
                    start,
                    kinds: (start..=end).map(kind).collect(),
                    windows: inside,
                }
            })
            .collect();
        Self {
            regions,
            overrun: widths.overrun,
        }
    }

    /// The instances the item `first..=last` (1-based) touches.
    pub(crate) fn touched(&self, first: i64, last: i64) -> BTreeSet<Instance> {
        let mut instances = BTreeSet::new();
        for (index, region) in self.regions.iter().enumerate() {
            let end = region.start + i64::try_from(region.kinds.len()).unwrap_or(0) - 1;
            let (from, to) = (first.max(region.start), last.min(end));
            let held = region
                .windows
                .iter()
                .filter(|&&(start, stop)| first <= stop && start <= last)
                .map(|&(start, stop)| (start - self.overrun, stop + self.overrun))
                .collect::<Vec<_>>();
            for position in from..=to {
                let kind = if held
                    .iter()
                    .any(|&(start, stop)| (start..=stop).contains(&position))
                {
                    Kind::Window
                } else {
                    region.kinds[usize::try_from(position - region.start).unwrap_or(0)]
                };
                instances.insert((index, kind));
            }
        }
        instances
    }

    /// The instance that holds `position`, when it lies in a region.
    pub(crate) fn instance(&self, position: i64) -> Option<Instance> {
        self.regions.iter().enumerate().find_map(|(index, region)| {
            let offset = usize::try_from(position - region.start).ok()?;
            Some((index, *region.kinds.get(offset)?))
        })
    }

    /// The regions' names, in region order.
    pub(crate) fn names(&self) -> impl Iterator<Item = &str> {
        self.regions.iter().map(|region| region.name.as_str())
    }

    /// `kind@region`, as `scores.tsv` writes an instance.
    pub(crate) fn label(&self, (region, kind): Instance) -> String {
        format!(
            "{}@{}",
            kind.label(),
            self.regions.get(region).map_or("?", |region| &region.name)
        )
    }

    /// Every instance of the zones, in region then kind order.
    pub(crate) fn instances(&self) -> Vec<Instance> {
        self.regions
            .iter()
            .enumerate()
            .flat_map(|(index, region)| {
                Kind::ALL
                    .into_iter()
                    .filter(|kind| region.kinds.contains(kind))
                    .map(move |kind| (index, kind))
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const WIDTHS: Widths = Widths {
        after: 20,
        edge: 40,
        overrun: 2,
    };

    fn hv2() -> Zones {
        Zones::new(&[(73, 340)], &[(303, 315)], WIDTHS)
    }

    #[test]
    fn positions_take_the_kind_of_highest_precedence() {
        let zones = hv2();
        let kind = |position| zones.instance(position).map(|(_, kind)| kind);
        assert_eq!(kind(73), Some(Kind::Start));
        assert_eq!(kind(112), Some(Kind::Start));
        assert_eq!(kind(113), Some(Kind::Interior));
        assert_eq!(kind(282), Some(Kind::Interior));
        assert_eq!(kind(283), Some(Kind::After));
        assert_eq!(kind(303), Some(Kind::Window));
        assert_eq!(kind(315), Some(Kind::Window));
        // After the window outranks the region's end.
        assert_eq!(kind(335), Some(Kind::After));
        assert_eq!(kind(336), Some(Kind::End));
        assert_eq!(kind(340), Some(Kind::End));
        assert_eq!(kind(341), None);
        assert_eq!(zones.label((0, Kind::After)), "after@HV2");
        assert_eq!(zones.instance(72), None);
        assert_eq!(zones.instance(340), Some((0, Kind::End)));
        assert_eq!(zones.names().collect::<Vec<_>>(), ["HV2"]);
        assert_eq!(
            zones.instances(),
            [
                (0, Kind::Window),
                (0, Kind::After),
                (0, Kind::Start),
                (0, Kind::End),
                (0, Kind::Interior)
            ]
        );
    }

    #[test]
    fn a_window_item_counts_its_overrun_for_the_window() {
        let zones = hv2();
        // 303-316 runs one past the window: the window alone.
        assert_eq!(zones.touched(303, 316), BTreeSet::from([(0, Kind::Window)]));
        assert_eq!(zones.touched(301, 317), BTreeSet::from([(0, Kind::Window)]));
        // Three past it reaches the after zone.
        assert_eq!(
            zones.touched(303, 318),
            BTreeSet::from([(0, Kind::Window), (0, Kind::After)])
        );
        // An item that does not overlap the window keeps its own kind.
        assert_eq!(zones.touched(316, 316), BTreeSet::from([(0, Kind::After)]));
        // Every instance an item touches counts, and positions outside the
        // regions none.
        assert_eq!(
            zones.touched(60, 120),
            BTreeSet::from([(0, Kind::Start), (0, Kind::Interior)])
        );
        assert!(zones.touched(1, 72).is_empty());
    }
}
