//! The target geometry the harness compares and counts in, read once from
//! the configuration's target profile; nothing else in the harness holds a
//! reference coordinate.
//!
//! - **Regions:** the profile's reportable regions (`variant_calling.regions`),
//!   in the profile's order, named by the laboratory's convention
//!   ([`conventions::region_name`]). They bound the comparison, the callable
//!   yield, the review windows, and the review zones.
//! - **Windows:** the repeat windows a discordant position is classified by
//!   (`class_regions`, `regions`, the written forms of the discrepancy
//!   causes): the profile's nomenclature windows, then its length tracts
//!   that lie outside every nomenclature window, each named by its profile
//!   name in lower case with every other character than a letter or digit
//!   written `_` (`HVS-II` is `hvs_ii`). A position in none is `other`.

use std::fmt;

use dna_kernel::profile::Profile;

use crate::conventions;

/// A named stretch of the reference, 1-based and inclusive.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct Span {
    pub(crate) name: String,
    pub(crate) first: i64,
    pub(crate) last: i64,
}

impl Span {
    pub(crate) fn contains(&self, position: i64) -> bool {
        (self.first..=self.last).contains(&position)
    }
}

/// The window that classifies a position: one of the geometry's windows, by
/// its name, or none of them. Windows order by name, `other` last.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub(crate) enum Window {
    Named(String),
    Other,
}

impl Window {
    pub(crate) fn label(&self) -> &str {
        match self {
            Self::Named(name) => name,
            Self::Other => "other",
        }
    }
}

impl fmt::Display for Window {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(self.label())
    }
}

/// The regions and windows of one target profile.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct Geometry {
    regions: Vec<Span>,
    /// The nomenclature windows' bounds: the repeat windows of the zones.
    repeats: Vec<(i64, i64)>,
    windows: Vec<Span>,
}

/// A profile coordinate as the harness's signed positions.
fn position(value: usize) -> i64 {
    i64::try_from(value).unwrap_or(i64::MAX)
}

/// A profile name as a window label: lower case, every other character than
/// an ASCII letter or digit written `_`.
fn label(name: &str) -> String {
    name.chars()
        .map(|symbol| {
            if symbol.is_ascii_alphanumeric() {
                symbol.to_ascii_lowercase()
            } else {
                '_'
            }
        })
        .collect()
}

impl Geometry {
    /// The geometry of `profile`.
    pub(crate) fn of(profile: &Profile) -> Self {
        let regions = profile
            .regions
            .iter()
            .enumerate()
            .map(|(index, &[first, last])| Span {
                name: conventions::region_name(index),
                first: position(first),
                last: position(last),
            })
            .collect();
        let mut windows = profile
            .windows
            .iter()
            .map(|window| Span {
                name: label(&window.name),
                first: position(window.start + 1),
                last: position(window.end()),
            })
            .collect::<Vec<_>>();
        let tracts = profile
            .length_tracts
            .iter()
            .map(|tract| Span {
                name: label(&tract.name),
                first: position(tract.window[0]),
                last: position(tract.window[1]),
            })
            .filter(|tract| {
                !windows
                    .iter()
                    .any(|window| tract.first <= window.last && window.first <= tract.last)
            })
            .collect::<Vec<_>>();
        let repeats = windows
            .iter()
            .map(|window| (window.first, window.last))
            .collect();
        windows.extend(tracts);
        Self {
            regions,
            repeats,
            windows,
        }
    }

    /// The reportable regions, in the profile's order.
    pub(crate) fn regions(&self) -> &[Span] {
        &self.regions
    }

    /// The nomenclature windows' bounds, in the profile's order: the repeat
    /// windows the review zones are placed around.
    pub(crate) fn repeats(&self) -> &[(i64, i64)] {
        &self.repeats
    }

    /// The regions' bounds, in the profile's order.
    pub(crate) fn bounds(&self) -> Vec<(i64, i64)> {
        self.regions
            .iter()
            .map(|region| (region.first, region.last))
            .collect()
    }

    /// The repeat windows' bounds, in the profile's order.
    pub(crate) fn window_bounds(&self) -> impl Iterator<Item = (Window, i64, i64)> + '_ {
        self.windows.iter().map(|window| {
            (
                Window::Named(window.name.clone()),
                window.first,
                window.last,
            )
        })
    }

    /// The region that holds the 1-based reference `position`.
    pub(crate) fn region_at(&self, position: i64) -> Option<&Span> {
        self.regions.iter().find(|region| region.contains(position))
    }

    /// The window that holds the 1-based reference `position`.
    pub(crate) fn window_at(&self, position: i64) -> Window {
        self.windows
            .iter()
            .find(|window| window.contains(position))
            .map_or(Window::Other, |window| Window::Named(window.name.clone()))
    }

    /// Every window, `other` included, in label order: the columns of the
    /// summary's window tables.
    pub(crate) fn all_windows(&self) -> Vec<Window> {
        let mut all = self
            .windows
            .iter()
            .map(|window| Window::Named(window.name.clone()))
            .collect::<Vec<_>>();
        all.sort();
        all.push(Window::Other);
        all
    }

    /// The window a label names, `other` included.
    pub(crate) fn window_named(&self, label: &str) -> Option<Window> {
        self.all_windows()
            .into_iter()
            .find(|window| window.label() == label)
    }
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;

    /// The geometry of the repository's rCRS profile, loaded once.
    pub(crate) fn rcrs() -> &'static Geometry {
        static GEOMETRY: std::sync::OnceLock<Geometry> = std::sync::OnceLock::new();
        GEOMETRY.get_or_init(|| {
            let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
                .join("../../config/profiles/human-mtdna-rcrs.toml");
            match Profile::load(&path) {
                Ok(profile) => Geometry::of(&profile),
                Err(error) => panic!("{}: {error}", path.display()),
            }
        })
    }

    /// The regions of the rCRS profile that `bounds` name, in that order.
    pub(crate) fn spans(bounds: &[(i64, i64)]) -> Vec<Span> {
        bounds
            .iter()
            .map(|&(first, last)| {
                rcrs()
                    .regions()
                    .iter()
                    .find(|region| (region.first, region.last) == (first, last))
                    .cloned()
                    .unwrap_or_else(|| panic!("no rCRS region {first}-{last}"))
            })
            .collect()
    }

    #[test]
    fn reads_regions_and_windows_from_the_profile() {
        let geometry = rcrs();
        assert_eq!(
            geometry
                .regions()
                .iter()
                .map(|region| (region.name.as_str(), region.first, region.last))
                .collect::<Vec<_>>(),
            [("HV1", 16024, 16365), ("HV2", 73, 340), ("HV3", 438, 576)]
        );
        assert_eq!(
            geometry
                .window_bounds()
                .map(|(window, first, last)| (window.label().to_owned(), first, last))
                .collect::<Vec<_>>(),
            [
                ("hvs_ii".to_owned(), 303, 315),
                ("hvs_iii".to_owned(), 512, 524),
                ("hvs_i".to_owned(), 16180, 16193),
                ("hvs3_c_tract".to_owned(), 568, 573),
            ]
        );
        assert_eq!(geometry.window_at(309).label(), "hvs_ii");
        assert_eq!(geometry.window_at(16195), Window::Other);
        assert_eq!(
            geometry.region_at(16193).map(|region| region.name.as_str()),
            Some("HV1")
        );
        assert!(geometry.region_at(600).is_none());
        assert_eq!(
            geometry
                .all_windows()
                .iter()
                .map(Window::label)
                .collect::<Vec<_>>(),
            ["hvs3_c_tract", "hvs_i", "hvs_ii", "hvs_iii", "other"]
        );
    }
}
