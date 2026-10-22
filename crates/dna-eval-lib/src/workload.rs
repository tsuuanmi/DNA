//! Review workload: how much of each sample DNA leaves to a reviewer, and
//! which differences from the truth no reviewer would see.
//!
//! The review items are `dna review`'s own (`SampleAnalysis::review_list`).
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
//! blocks. The counting itself is [`PassSample`]'s, which `dna-eval diff`
//! also counts a run's files with, so a run and its comparison count the
//! same passes.

use std::collections::{BTreeMap, BTreeSet};

use serde::{Deserialize, Serialize};

use crate::compare;
use crate::geometry::Span;
use crate::haplotype::Block;
use crate::zones::{Instance, Zones};

/// One window of an HV region, 1-based and inclusive.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct Window {
    pub(crate) region: String,
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

/// The windows of `size` positions over every region, in region order.
pub(crate) fn windows(size: i64, regions: &[Span]) -> Vec<Window> {
    let size = size.max(1);
    regions
        .iter()
        .flat_map(|region| {
            (0..)
                .map(move |index| region.first + index * size)
                .take_while(move |start| *start <= region.last)
                .map(move |start| Window {
                    region: region.name.clone(),
                    start,
                    end: (start + size - 1).min(region.last),
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
    /// Difference blocks with a differing position no check item reviews.
    pub(crate) unreviewed: Vec<String>,
    /// The same, also leaving out blocks a forced window overlaps.
    pub(crate) unreviewed_forced: Vec<String>,
    /// Per window, whether a check item touches it.
    pub(crate) checked: Vec<bool>,
    /// Per window, whether an unreviewed block overlaps it.
    pub(crate) escaped: Vec<bool>,
    /// Region of each unreviewed block, by its first anchor.
    pub(crate) regions: Vec<String>,
    /// The zone instances a check item touches.
    pub(crate) zones: BTreeSet<Instance>,
    /// Check items that touch more than one zone instance.
    pub(crate) spanning: usize,
    /// Each unreviewed block's zones with their classes, and each
    /// unreviewed omission.
    pub(crate) escapes: Vec<Escape>,
}

/// One unreviewed difference: a block's unreviewed positions in one zone,
/// or a truth allele the report omits, as the zone and pass measures count
/// it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct Escape {
    /// The block's first and last anchors, or the omitted allele's anchor,
    /// 1-based.
    pub(crate) first: i64,
    pub(crate) last: i64,
    /// The block's text, or the omitted allele's.
    pub(crate) text: String,
    /// The zone instance that holds these positions, when one does.
    pub(crate) zone: Option<Instance>,
    /// What it is, as the pass measures split it.
    pub(crate) class: Class,
    /// What it would be without the re-check list.
    pub(crate) unlisted: Class,
}

/// What an unreviewed difference is, as the pass measures split it, from the
/// least to the most severe; a passed unit counts under its most severe one.
///
/// A difference block is decided per differing position (a rule position of
/// the `haplotype_lod` comparison inside it), and takes the most severe class
/// of its unreviewed positions. A position the re-check list names is
/// `Recheck`. Otherwise the whole block is `TruthN` when every difference is
/// an insertion past the truth's interval ([`Block::past_interval`]), and
/// `LengthNote` when it differs only in the length of a profile length tract,
/// by one unit; else a position is `TruthN` where the truth writes `N`, `Het`
/// where it writes a heteroplasmy, and `NonHet` anywhere else: a DNA error
/// that neither explains.
///
/// A truth allele inside the truth's analysed interval that DNA's callable
/// intervals leave out is omitted from the report: `OmittedNotAnalysed` in a
/// stretch DNA did not analyse (a poor read start or end), `OmittedTract` in
/// a profile length tract under a `tract_length` note whose reading is the
/// truth's length, `NoteMisread` there when the reading differs, and
/// `OmittedOther` anywhere else.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) enum Class {
    LengthNote,
    OmittedTract,
    NoteMisread,
    Recheck,
    Het,
    TruthN,
    OmittedNotAnalysed,
    OmittedOther,
    NonHet,
}

impl Class {
    /// Every class, from the least to the most severe.
    pub(crate) const ALL: [Self; 9] = [
        Self::LengthNote,
        Self::OmittedTract,
        Self::NoteMisread,
        Self::Recheck,
        Self::Het,
        Self::TruthN,
        Self::OmittedNotAnalysed,
        Self::OmittedOther,
        Self::NonHet,
    ];

    /// The class a label of `unreviewed.tsv` names.
    pub(crate) fn parse(label: &str) -> Option<Self> {
        Self::ALL.into_iter().find(|class| class.label() == label)
    }

    /// The class's label in `unreviewed.tsv`.
    pub(crate) const fn label(self) -> &'static str {
        match self {
            Self::LengthNote => "length_note",
            Self::OmittedTract => "omitted_tract",
            Self::NoteMisread => "note_misread",
            Self::Recheck => "recheck",
            Self::Het => "het",
            Self::TruthN => "truth_n",
            Self::OmittedNotAnalysed => "omitted_not_analysed",
            Self::OmittedOther => "omitted_other",
            Self::NonHet => "non_het",
        }
    }

    /// The class of a truth symbol at a differing position: `N` written
    /// there, a heteroplasmy (an IUPAC code of two or three bases), or
    /// anything else.
    pub(crate) fn of_truth(symbols: &str) -> Self {
        if symbols.contains('N') {
            Self::TruthN
        } else if symbols.char_indices().any(|(index, symbol)| {
            matches!(
                compare::bases(&symbols[index..index + symbol.len_utf8()]).len(),
                2 | 3
            )
        }) {
            Self::Het
        } else {
            Self::NonHet
        }
    }
}

/// One differing position of a block: its anchor, the reference stretch a
/// check item must touch to review it (the anchor, or the reference run an
/// indel there lies in), and the class of the truth's symbol there.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct Position {
    pub(crate) anchor: i64,
    pub(crate) reach: (i64, i64),
    pub(crate) truth: Class,
}

/// A truth allele DNA's report omits: its text, anchor, the stretch a check
/// item must touch to review it, and its omission class.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct Omission {
    pub(crate) text: String,
    pub(crate) anchor: i64,
    pub(crate) reach: (i64, i64),
    pub(crate) class: Class,
}

/// What a sample's workload is measured on besides its check items: its
/// `haplotype_lod` difference blocks with each one's differing positions,
/// and the truth alleles its report omits.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub(crate) struct Differences<'b> {
    pub(crate) blocks: Vec<(&'b Block, Vec<Position>)>,
    pub(crate) omissions: Vec<Omission>,
}

/// A homopolymer tract of the profile with length heteroplasmy, 1-based
/// and inclusive: a length difference there is noted, never gated.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct Tract {
    /// The profile's name of the tract.
    pub(crate) name: String,
    pub(crate) first: i64,
    pub(crate) last: i64,
    pub(crate) base: char,
}

impl Tract {
    /// Whether `block` differs only in the length of a run of the tract:
    /// both sides read the same sequence once the tract's base is taken out
    /// (the structure equal: no `T310C`, no `16189C` one side lacks), and
    /// the first run whose length differs overlaps the tract, starting no
    /// later than its last position (a `16183C` run reaches before it). A block
    /// may reach past the tract, to where its two sides' lengths meet again
    /// or its interval ends; a run that differs after the tract's own
    /// (`315.1C`) is no tract length.
    pub(crate) fn length_only(&self, block: &Block) -> bool {
        let rest = |text: &str| {
            text.chars()
                .filter(|&symbol| symbol != self.base)
                .collect::<String>()
        };
        if block.dna.len() == block.truth.len() || rest(&block.dna) != rest(&block.truth) {
            return false;
        }
        // The first symbol where the sides differ lies in a run of the
        // tract's base on the longer side; that run's start, counted from
        // the block's first anchor, is where the differing run begins.
        let longer = if block.dna.len() > block.truth.len() {
            block.dna.as_bytes()
        } else {
            block.truth.as_bytes()
        };
        let shorter = if block.dna.len() > block.truth.len() {
            block.truth.as_bytes()
        } else {
            block.dna.as_bytes()
        };
        let common = longer
            .iter()
            .zip(shorter)
            .take_while(|(left, right)| left == right)
            .count();
        let base = u8::try_from(self.base).unwrap_or(0);
        let Some(mut start) = (common < longer.len()).then_some(common) else {
            return false;
        };
        if longer[start] != base {
            return false;
        }
        while start > 0 && longer[start - 1] == base {
            start -= 1;
        }
        let mut end = common;
        while end + 1 < longer.len() && longer[end + 1] == base {
            end += 1;
        }
        let at = |index: usize| block.start + i64::try_from(index).unwrap_or(i64::MAX);
        // The longer side's run lies at least one base past the shorter's, so
        // its first and last bases bound the run on the reference.
        at(start) <= self.last && self.first <= at(end)
    }

    /// How many more of the tract's base DNA reads in `block` than the
    /// truth.
    pub(crate) fn excess(&self, block: &Block) -> i64 {
        let count = |text: &str| {
            i64::try_from(text.chars().filter(|&symbol| symbol == self.base).count())
                .unwrap_or(i64::MAX)
        };
        count(&block.dna) - count(&block.truth)
    }
}

/// The workload of one sample from its check items (`first..=last`,
/// 1-based), its differences, its spans on the re-check list, and the
/// profile's length `tracts`. A position or an omission is reviewed when a
/// check item touches its reach; a block is unreviewed when any of its
/// positions is, and counts once per zone that holds an unreviewed position.
pub(crate) fn of(
    checks: &[(i64, i64)],
    differences: &Differences<'_>,
    windows: &[Window],
    forced: &[(i64, i64)],
    zones: Option<&Zones>,
    (listed, tracts): (&[(i64, i64)], &[Tract]),
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
    let reviewed =
        |(first, last): (i64, i64)| checks.iter().any(|&(from, to)| from <= last && first <= to);
    let on_list = |anchor: i64| {
        listed
            .iter()
            .any(|&(first, last)| (first..=last).contains(&anchor))
    };
    let zone_of = |anchor: i64| zones.and_then(|zones| zones.instance(anchor));
    let mut differences_out = Vec::new();
    let mut unreviewed = Vec::new();
    for (block, block_positions) in &differences.blocks {
        let open = block_positions
            .iter()
            .filter(|position| !reviewed(position.reach))
            .collect::<Vec<_>>();
        if open.is_empty() {
            continue;
        }
        unreviewed.push(*block);
        let whole = if block.past_interval() {
            Some(Class::TruthN)
        } else if tracts
            .iter()
            .any(|tract| tract.length_only(block) && tract.excess(block).abs() == 1)
        {
            Some(Class::LengthNote)
        } else {
            None
        };
        let mut by_zone = BTreeMap::<Option<Instance>, (Class, Class)>::new();
        for position in open {
            let unlisted = whole.unwrap_or(position.truth);
            let class = if on_list(position.anchor) {
                Class::Recheck
            } else {
                unlisted
            };
            let entry = by_zone
                .entry(zone_of(position.anchor))
                .or_insert((class, unlisted));
            *entry = (entry.0.max(class), entry.1.max(unlisted));
        }
        for (zone, (class, unlisted)) in by_zone {
            differences_out.push(Escape {
                first: block.start,
                last: block.end,
                text: block.text(),
                zone,
                class,
                unlisted,
            });
        }
    }
    for omission in &differences.omissions {
        if reviewed(omission.reach) {
            continue;
        }
        differences_out.push(Escape {
            first: omission.anchor,
            last: omission.anchor,
            text: omission.text.clone(),
            zone: zone_of(omission.anchor),
            class: if on_list(omission.anchor) {
                Class::Recheck
            } else {
                omission.class
            },
            unlisted: omission.class,
        });
    }
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
            windows
                .iter()
                .find(|window| window.overlaps(block.start, block.start))
                .map(|window| window.region.clone())
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
        escapes: differences_out,
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

/// Passed samples, regions, or zones that carry an unreviewed difference:
/// in all, and by the most severe [`Class`] of their differences: with a
/// `NonHet` block, else with an omission elsewhere, else with an omission
/// in a stretch DNA did not analyse, else with a truth `N`, else het only,
/// else on the re-check list only, else a misread tract note, else omitted
/// tract lengths only, else noted tract lengths only.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
pub(crate) struct Escaped {
    pub(crate) total: usize,
    pub(crate) non_het: usize,
    /// The omission classes are 0 in a run that predates them.
    #[serde(default)]
    pub(crate) omitted_other: usize,
    #[serde(default)]
    pub(crate) omitted_not_analysed: usize,
    pub(crate) het_only: usize,
    pub(crate) truth_n_only: usize,
    pub(crate) recheck_only: usize,
    #[serde(default)]
    pub(crate) note_misread_only: usize,
    #[serde(default)]
    pub(crate) omitted_tract_only: usize,
    pub(crate) length_note_only: usize,
}

impl Escaped {
    /// Counts one passed unit with the classes of its unreviewed blocks.
    fn add(&mut self, classes: impl IntoIterator<Item = Class>) {
        let Some(class) = classes.into_iter().max() else {
            return;
        };
        self.total += 1;
        *match class {
            Class::NonHet => &mut self.non_het,
            Class::OmittedOther => &mut self.omitted_other,
            Class::OmittedNotAnalysed => &mut self.omitted_not_analysed,
            Class::TruthN => &mut self.truth_n_only,
            Class::Het => &mut self.het_only,
            Class::Recheck => &mut self.recheck_only,
            Class::NoteMisread => &mut self.note_misread_only,
            Class::OmittedTract => &mut self.omitted_tract_only,
            Class::LengthNote => &mut self.length_note_only,
        } += 1;
    }

    /// The sum of `counts`.
    pub(crate) fn sum<'e>(counts: impl IntoIterator<Item = &'e Self>) -> Self {
        counts.into_iter().fold(Self::default(), |sum, count| Self {
            total: sum.total + count.total,
            non_het: sum.non_het + count.non_het,
            omitted_other: sum.omitted_other + count.omitted_other,
            omitted_not_analysed: sum.omitted_not_analysed + count.omitted_not_analysed,
            het_only: sum.het_only + count.het_only,
            truth_n_only: sum.truth_n_only + count.truth_n_only,
            recheck_only: sum.recheck_only + count.recheck_only,
            note_misread_only: sum.note_misread_only + count.note_misread_only,
            omitted_tract_only: sum.omitted_tract_only + count.omitted_tract_only,
            length_note_only: sum.length_note_only + count.length_note_only,
        })
    }
}

/// What passes without review over a set of samples: the samples, the
/// regions of the zones, and the zone instances that no check item touches.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
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
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub(crate) struct AutoPass {
    pub(crate) samples: usize,
    pub(crate) unreviewed: Escaped,
}

/// Regions in all, those with a check item in one of their zones, and the
/// others (passed) that hold an unreviewed block's first anchor; overall and
/// per region.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub(crate) struct Regions {
    pub(crate) total: usize,
    pub(crate) checked: usize,
    pub(crate) unreviewed: Escaped,
    pub(crate) by_region: Vec<Unit>,
}

/// Zone instances in all, those a check item touches, and the others
/// (passed) that hold an unreviewed block's first anchor; overall and per
/// instance, in the order of [`Zones::instances`].
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub(crate) struct ZonePasses {
    pub(crate) total: usize,
    pub(crate) checked: usize,
    pub(crate) unreviewed: Escaped,
    pub(crate) by_zone: Vec<Unit>,
}

/// One region or zone instance over the samples counted.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub(crate) struct Unit {
    pub(crate) name: String,
    pub(crate) total: usize,
    pub(crate) checked: usize,
    pub(crate) unreviewed: Escaped,
}

impl Unit {
    pub(crate) fn new(name: String) -> Self {
        Self {
            name,
            total: 0,
            checked: 0,
            unreviewed: Escaped::default(),
        }
    }

    /// Counts `sample` in the unit `key`: checked when a check item of the
    /// sample touches it, else passed with the sample's unreviewed blocks
    /// that start in it.
    pub(crate) fn add<K: Ord>(&mut self, sample: &PassSample<K>, key: &K) {
        self.total += 1;
        if sample.zones.contains(key) {
            self.checked += 1;
        } else {
            self.unreviewed.add(
                sample
                    .escapes
                    .iter()
                    .filter(|(zone, _)| zone.as_ref() == Some(key))
                    .map(|&(_, class)| class),
            );
        }
    }
}

impl AutoPass {
    /// Counts `sample`: an auto-pass sample when no check item is in its
    /// review list, with all its unreviewed blocks.
    pub(crate) fn add<K>(&mut self, sample: &PassSample<K>) {
        if sample.clear {
            self.samples += 1;
            self.unreviewed
                .add(sample.escapes.iter().map(|&(_, class)| class));
        }
    }
}

/// One sample as the pass measures count it, its units keyed by `K`: whether
/// its review list has no check item, the units its check items touch, and
/// each unreviewed block's unit (the one that holds its first anchor, if
/// any) and class.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct PassSample<K> {
    pub(crate) clear: bool,
    pub(crate) zones: BTreeSet<K>,
    pub(crate) escapes: Vec<(Option<K>, Class)>,
}

impl Workload {
    /// The workload as the pass measures count it, its units mapped by
    /// `key`.
    fn pass_sample<K: Ord>(&self, key: impl Fn(Instance) -> K) -> PassSample<K> {
        PassSample {
            clear: self.checks == 0,
            zones: self.zones.iter().map(|&instance| key(instance)).collect(),
            escapes: self
                .escapes
                .iter()
                .map(|escape| (escape.zone.map(&key), escape.class))
                .collect(),
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
        let sample = workload.pass_sample(|instance| instance);
        auto_pass.add(&sample);
        let by_index = workload.pass_sample(|(index, _)| index);
        for (index, region) in by_region.iter_mut().enumerate() {
            region.add(&by_index, &index);
        }
        for (instance, zone) in instances.iter().zip(&mut by_zone) {
            zone.add(&sample, instance);
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

    /// The blocks with one differing position each, at the block's first
    /// anchor, of the class of the block's truth symbols.
    fn differences(blocks: &[Block]) -> Differences<'_> {
        Differences {
            blocks: blocks
                .iter()
                .map(|block| {
                    (
                        block,
                        vec![Position {
                            anchor: block.start,
                            reach: (block.start, block.start),
                            truth: Class::of_truth(&block.truth),
                        }],
                    )
                })
                .collect(),
            omissions: Vec::new(),
        }
    }

    #[test]
    fn windows_align_to_each_region_start_and_end_short() {
        let windows = windows(10, crate::geometry::tests::rcrs().regions());
        assert_eq!(
            (windows[0].region.as_str(), windows[0].start, windows[0].end),
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
        let windows = windows(10, crate::geometry::tests::rcrs().regions());
        let zones = Zones::new(
            &crate::geometry::tests::spans(&[(16024, 16365)]),
            &[(16180, 16193)],
            crate::zones::Widths {
                after: 20,
                edge: 40,
                overrun: 2,
            },
        );
        let blocks = [block(16031, 16031), block(16040, 16041)];
        let workload = of(
            &[(16030, 16031)],
            &differences(&blocks),
            &windows,
            &[],
            Some(&zones),
            (&[], &[]),
        );
        assert_eq!(workload.checks, 1);
        assert_eq!(workload.zones, BTreeSet::from([(0, Kind::Start)]));
        assert_eq!(workload.spanning, 0);
        assert_eq!(
            workload.escapes,
            [Escape {
                first: 16040,
                last: 16041,
                text: "16040-16041:A/C".to_owned(),
                zone: Some((0, Kind::Start)),
                class: Class::NonHet,
                unlisted: Class::NonHet,
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
        let windows = windows(10, crate::geometry::tests::rcrs().regions());
        let blocks = [block(16185, 16186)];
        let workload = of(
            &[],
            &differences(&blocks),
            &windows,
            &[(16184, 16193)],
            None,
            (&[], &[]),
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
    fn a_position_is_het_where_the_truth_writes_a_two_or_three_base_code() {
        let het = |truth: &str| Class::of_truth(truth) == Class::Het;
        assert!(het("ACY") && het("R") && het("B") && het("V"));
        assert!(!het("ACGT") && !het("N") && !het(""));
        assert_eq!(Class::of_truth("NY"), Class::TruthN);
    }

    /// A far check item no longer reviews an open block: each differing
    /// position needs its own check, and the block counts in the zone of
    /// every unreviewed position (correctness audit, F2).
    #[test]
    fn a_block_is_reviewed_position_by_position() {
        let windows = windows(10, crate::geometry::tests::rcrs().regions());
        let zones = Zones::new(
            &crate::geometry::tests::spans(&[(438, 576)]),
            &[(512, 524)],
            crate::zones::Widths {
                after: 20,
                edge: 40,
                overrun: 2,
            },
        );
        let open = Block::open_fixture(523, 576, "AC", "ACACAC");
        let position = |anchor, truth| Position {
            anchor,
            reach: (anchor, anchor),
            truth,
        };
        let differences = Differences {
            blocks: vec![(
                &open,
                vec![position(523, Class::NonHet), position(560, Class::Het)],
            )],
            omissions: Vec::new(),
        };
        // A check over 554-560 reviews 560 only: 523 stays unreviewed, in
        // the HVS-III window, and its class is its own position's.
        let far = of(
            &[(554, 560)],
            &differences,
            &windows,
            &[],
            Some(&zones),
            (&[], &[]),
        );
        assert_eq!(far.unreviewed.len(), 1);
        assert_eq!(
            far.escapes
                .iter()
                .map(|escape| (escape.zone, escape.class))
                .collect::<Vec<_>>(),
            [(Some((0, Kind::Window)), Class::NonHet)]
        );
        // Both positions checked: reviewed.
        let both = of(
            &[(520, 525), (554, 560)],
            &differences,
            &windows,
            &[],
            Some(&zones),
            (&[], &[]),
        );
        assert!(both.unreviewed.is_empty() && both.escapes.is_empty());
        // Unchecked, the block counts in both zones, each with its class.
        let none = of(&[], &differences, &windows, &[], Some(&zones), (&[], &[]));
        assert_eq!(
            none.escapes
                .iter()
                .map(|escape| (escape.zone, escape.class))
                .collect::<Vec<_>>(),
            [
                (Some((0, Kind::Window)), Class::NonHet),
                (Some((0, Kind::End)), Class::Het)
            ]
        );
        // An omitted truth allele is an escape unless a check touches its
        // reach.
        let omitted = Differences {
            blocks: Vec::new(),
            omissions: vec![Omission {
                text: "73G".to_owned(),
                anchor: 73,
                reach: (73, 73),
                class: Class::OmittedNotAnalysed,
            }],
        };
        let hv2 = windows.clone();
        let flagged = of(&[(70, 80)], &omitted, &hv2, &[], None, (&[], &[]));
        assert!(flagged.escapes.is_empty());
        let silent = of(&[], &omitted, &hv2, &[(73, 73)], None, (&[(73, 73)], &[]));
        assert_eq!(
            (silent.escapes[0].class, silent.escapes[0].unlisted),
            (Class::Recheck, Class::OmittedNotAnalysed)
        );
        assert!(silent.unreviewed.is_empty());
    }

    #[test]
    fn passes_count_samples_regions_and_zones_without_a_check_item() {
        let windows = windows(10, crate::geometry::tests::rcrs().regions());
        // HV3 holds no repeat window here, so it has three zone kinds.
        let zones = Zones::new(
            &crate::geometry::tests::spans(&[(16024, 16365), (73, 340), (438, 576)]),
            &[(16184, 16193), (303, 315)],
            crate::zones::Widths {
                after: 20,
                edge: 40,
                overrun: 2,
            },
        );
        let of = |checks: &[(i64, i64)], blocks: &[Block]| {
            of(
                checks,
                &differences(blocks),
                &windows,
                &[],
                Some(&zones),
                (&[], &[]),
            )
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
            ..Escaped::default()
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
        let windows = windows(10, crate::geometry::tests::rcrs().regions());
        let class_with = |block: Block, listed: &[(i64, i64)], tracts: &[Tract]| {
            let blocks = [block];
            of(
                &[],
                &differences(&blocks),
                &windows,
                &[],
                None,
                (listed, tracts),
            )
            .escapes[0]
                .class
        };
        let class = |truth: &str, listed: &[(i64, i64)]| {
            class_with(Block::fixture(100, 101, "AC", truth), listed, &[])
        };
        // A listed position is re-check, then a truth `N`, then het.
        assert_eq!(class("AC", &[(100, 105)]), Recheck);
        assert_eq!(class("NY", &[(90, 100)]), Recheck);
        assert_eq!(class("NY", &[(102, 105)]), TruthN);
        assert_eq!(class("AY", &[]), Het);
        assert_eq!(class("GG", &[]), NonHet);
        // An insertion past the truth's interval: the truth records nothing
        // there. A block that also differs at the last base, holds more than
        // it, or stays inside the interval is DNA's.
        let open = |start, dna: &str, truth: &str| {
            class_with(Block::open_fixture(start, 573, dna, truth), &[], &[])
        };
        assert_eq!(open(573, "CCC", "C"), TruthN);
        assert_eq!(open(573, "C", "CC"), TruthN);
        assert_eq!(open(573, "ACC", "C"), NonHet);
        assert_eq!(open(572, "CCCC", "CC"), NonHet);
        assert_eq!(
            class_with(Block::fixture(573, 573, "CCC", "C"), &[], &[]),
            NonHet
        );
        // A length difference of one unit inside a length tract, the
        // structure equal, is noted; two units, a T read as C, or a length
        // outside the tract are not.
        let tracts = [Tract {
            name: "hvs2".to_owned(),
            first: 303,
            last: 309,
            base: 'C',
        }];
        let tract = |start, end, dna: &str, truth: &str| {
            class_with(Block::fixture(start, end, dna, truth), &[], &tracts)
        };
        assert_eq!(tract(309, 309, "CC", "C"), Class::LengthNote);
        assert_eq!(tract(309, 309, "CCC", "C"), NonHet);
        assert_eq!(tract(309, 310, "CCC", "CT"), NonHet);
        assert_eq!(tract(315, 315, "CC", "C"), NonHet);
        // A block that reaches past the tract: run 1 differs, or run 2.
        assert_eq!(
            tract(309, 316, "CCCCTCCCCCCG", "CCCTCCCCCCG"),
            Class::LengthNote
        );
        assert_eq!(tract(309, 316, "CCCTCCCCCCCG", "CCCTCCCCCCG"), NonHet);
        // Starting before the tract, the run that differs starts in it.
        assert_eq!(
            tract(301, 309, "AACCCCCCCCC", "AACCCCCCCC"),
            Class::LengthNote
        );
        assert_eq!(tract(309, 309, "CC", "Y"), Het);
        assert_eq!(tracts[0].excess(&Block::fixture(309, 309, "CCC", "C")), 2);
        // A passed unit counts once, under its most severe difference.
        let count = |classes: &[Class]| {
            let mut escaped = Escaped::default();
            escaped.add(classes.iter().copied());
            escaped
        };
        let only = |non_het, het_only, truth_n_only, recheck_only| Escaped {
            total: 1,
            non_het,
            het_only,
            truth_n_only,
            recheck_only,
            ..Escaped::default()
        };
        assert_eq!(count(&[Recheck, Het, TruthN, NonHet]), only(1, 0, 0, 0));
        assert_eq!(count(&[Recheck, Het, TruthN]), only(0, 0, 1, 0));
        assert_eq!(count(&[Recheck, Het]), only(0, 1, 0, 0));
        assert_eq!(count(&[Recheck, Recheck]), only(0, 0, 0, 1));
        assert_eq!(count(&[]), Escaped::default());
        // An omission elsewhere outranks every class but non-het; an omitted
        // tract length and a noted one are the least severe.
        assert_eq!(
            count(&[TruthN, Class::OmittedOther]),
            Escaped {
                omitted_other: 1,
                ..only(0, 0, 0, 0)
            }
        );
        assert_eq!(
            count(&[Class::OmittedNotAnalysed, Het]),
            Escaped {
                omitted_not_analysed: 1,
                ..only(0, 0, 0, 0)
            }
        );
        assert_eq!(count(&[Class::LengthNote, Recheck]), only(0, 0, 0, 1));
        assert_eq!(
            count(&[Class::LengthNote, Class::OmittedTract]),
            Escaped {
                omitted_tract_only: 1,
                ..only(0, 0, 0, 0)
            }
        );
        assert_eq!(
            count(&[Class::LengthNote]),
            Escaped {
                length_note_only: 1,
                ..only(0, 0, 0, 0)
            }
        );
        // A listed position passes as re-check only; unlisted it is non-het,
        // and the escape keeps both classes.
        let block = [Block::fixture(573, 573, "CCCCCCCA", "CCCCCCA")];
        let listed = of(
            &[],
            &differences(&block),
            &windows,
            &[],
            None,
            (&[(573, 573)], &[]),
        );
        let unlisted = of(
            &[],
            &differences(&block),
            &windows,
            &[],
            None,
            (&[(16057, 16057)], &[]),
        );
        assert_eq!(
            (listed.escapes[0].class, listed.escapes[0].unlisted),
            (Recheck, NonHet)
        );
        assert_eq!(unlisted.escapes[0].class, NonHet);
        assert_eq!(listed.unreviewed, unlisted.unreviewed);
        let auto_pass = |workload: &Workload| passes([workload], None).auto_pass.unreviewed;
        assert_eq!(auto_pass(&listed), only(0, 0, 0, 1));
        assert_eq!(auto_pass(&unlisted), only(1, 0, 0, 0));
        for class in Class::ALL {
            assert_eq!(Class::parse(class.label()), Some(class));
        }
    }

    #[test]
    fn the_wilson_bound_holds_its_events() {
        assert_eq!(wilson_upper(0, 0), 1.0);
        let bound = wilson_upper(5, 100);
        assert!(bound > 0.05 && bound < 0.12, "{bound}");
    }
}
