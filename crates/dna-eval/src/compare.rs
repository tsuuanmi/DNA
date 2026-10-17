//! Concordance of DNA's calls with the reviewer truth of one sample.
//!
//! The rules adapt the pairwise comparison of the lab's `mtdna_raw` pipeline
//! (`src/core/comparison.py`), with DNA as tool A (the automated caller) and
//! the reviewer truth as tool B:
//!
//! - an allele is `(position, reference, alternate)` with `-` for an inserted
//!   or deleted base; positions are normalised (`16193.0` → `16193`), and
//!   alleles whose reference equals the alternate are dropped;
//! - `N/A` when either side has no allele;
//! - insertions and deletions at 309, 455, 463, and 573 (and their `.n`
//!   insertion positions) are left out ("special" positions);
//! - `N` when a position both sides call holds different allele sets;
//!   `comparison.py` also never compares 309.1, 309.2, and 573.1, but the
//!   insertions there are special positions already, so that rule adds
//!   nothing and is not repeated here;
//! - `N` when DNA alone calls an allele inside HV1, HV2, or HV3;
//! - `N` when the truth alone calls an allele inside HV1, HV2, or HV3 and
//!   inside DNA's callable intervals (or DNA has none);
//! - `Y` otherwise.
//!
//! Interval and region tests use an allele's anchor base, so `16193.1` lies in
//! `[16024, 16193]`.
//!
//! Two secondary measures change one rule each ([`Measure`]):
//!
//! - **strict:** a truth-only allele in an HV region also counts inside the
//!   truth's analysed intervals, so a truth allele where DNA has no callable
//!   base is missed rather than ignored;
//! - **IUPAC-compatible:** DNA's base matches a truth IUPAC code that admits
//!   it, like `iupac_bases_compatible` in `mtdna_raw`: two alleles at one
//!   position with the same reference match when their called bases share a
//!   possible base (`C` and `Y`), and where DNA calls nothing inside its
//!   callable intervals its base is the reference, which matches a truth IUPAC
//!   allele whose code admits the reference (`16189Y` over a reference `T`).
//!
//! Two haplotype measures apply the primary or the IUPAC-compatible rules
//! only to the alleles that do not describe the same sequence on both sides:
//! inside the comparable intervals (the HV regions within DNA's callable
//! intervals, or the HV regions when DNA reports none), each side's alleles,
//! special positions included, are applied to the reference, and alleles
//! where the two sequences agree leave the comparison
//! ([`haplotype`](crate::haplotype)). `16188.1C` and `16193.1C` in one C run,
//! or `248DEL` and `249DEL` over `AA`, are then concordant. Special-position
//! indel lengths cannot decide a difference, and a symbol matches a
//! compatible one only where a truth allele wrote it, as under the
//! IUPAC-compatible rules. The rules still decide the rest, special positions
//! left out, so a sample concordant under the base rules stays concordant.

use std::collections::{BTreeMap, BTreeSet};
use std::fmt;
use std::path::PathBuf;

use crate::haplotype::{Block, Edit, Haplotypes};

/// Positions whose insertions and deletions are left out of the comparison.
const SPECIAL_POSITIONS: [&str; 4] = ["455", "463", "573", "309"];
/// Hypervariable regions, closed at both ends.
const HV_REGIONS: [(&str, i64, i64); 3] =
    [("HV1", 16024, 16365), ("HV2", 73, 340), ("HV3", 438, 576)];
/// Single-base codes for bases A, C, G, T; any other base is an IUPAC code.
const DEFINITE: [&str; 4] = ["A", "C", "G", "T"];

/// One allele state at one position.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) struct Allele {
    /// Normalised position: `263`, `309.1`.
    pub(crate) position: String,
    /// Reference base, or `-` for an inserted base.
    pub(crate) reference: String,
    /// Called base, or `-` for a deleted base.
    pub(crate) alternate: String,
}

impl Allele {
    pub(crate) fn new(position: &str, reference: &str, alternate: &str) -> Self {
        Self {
            position: normalize_position(position),
            reference: reference.to_owned(),
            alternate: alternate.to_owned(),
        }
    }

    /// Insertion or deletion at a special position.
    fn is_special(&self) -> bool {
        let special = SPECIAL_POSITIONS.iter().any(|special| {
            self.position == *special
                || self
                    .position
                    .strip_prefix(special)
                    .is_some_and(|rest| rest.starts_with('.'))
        });
        special
            && (self.reference == "-"
                || self.alternate == "-"
                || self.reference.len() != self.alternate.len())
    }

    /// The reviewer notation of the allele: `263G`, `309.1C`, `249DEL`.
    pub(crate) fn text(&self) -> String {
        if self.alternate == "-" && self.reference != "-" {
            format!("{}DEL", self.position)
        } else {
            format!("{}{}", self.position, self.alternate)
        }
    }

    pub(crate) fn is_iupac(&self) -> bool {
        self.alternate != "-" && !DEFINITE.contains(&self.alternate.as_str())
    }

    /// Sort key: position order, then the alleles.
    fn order(&self) -> ((i64, i64), &str, &str) {
        (
            position_order(&self.position),
            &self.reference,
            &self.alternate,
        )
    }
}

/// Reference order of a position: anchor base, then insertion ordinal (`-1`
/// for the base itself).
fn position_order(position: &str) -> (i64, i64) {
    let base = anchor(position).unwrap_or(i64::MAX);
    let ordinal = position
        .split_once('.')
        .and_then(|(_, ordinal)| ordinal.parse::<i64>().ok())
        .unwrap_or(-1);
    (base, ordinal)
}

/// Positions in reference order.
pub(crate) fn sorted_positions(positions: &BTreeSet<String>) -> Vec<&str> {
    let mut sorted = positions.iter().map(String::as_str).collect::<Vec<_>>();
    sorted.sort_by_key(|position| position_order(position));
    sorted
}

/// The canonical text of a position: a whole number without separators or a
/// decimal part, otherwise the trimmed text without thousands separators.
pub(crate) fn normalize_position(text: &str) -> String {
    let cleaned = text.trim().replace(',', "");
    if let Ok(value) = cleaned.parse::<i64>() {
        return value.to_string();
    }
    match cleaned.parse::<f64>() {
        Ok(value) if value.is_finite() && value.fract() == 0.0 => format!("{value:.0}"),
        _ => cleaned,
    }
}

/// The anchor base of a position: `309` for `309.1`.
pub(crate) fn anchor(position: &str) -> Option<i64> {
    position.split('.').next()?.parse::<i64>().ok()
}

fn hv_region(position: &str) -> Option<&'static str> {
    let base = anchor(position)?;
    HV_REGIONS
        .iter()
        .find(|(_, start, end)| (*start..=*end).contains(&base))
        .map(|(name, _, _)| *name)
}

/// A set of comparison rules, declared in [`Measure::ALL`] order.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) enum Measure {
    /// The lab's rules.
    Primary,
    /// The lab's rules, with truth-only alleles also counted inside the
    /// truth's analysed intervals.
    Strict,
    /// The lab's rules, with IUPAC-compatible called bases matching.
    IupacCompatible,
    /// The lab's rules over the alleles that do not describe the same
    /// sequence.
    Haplotype,
    /// The IUPAC-compatible rules over the alleles that do not describe the
    /// same sequence, IUPAC codes that share a base counting as the same.
    HaplotypeIupacCompatible,
    /// As `haplotype`, but a truth heteroplasmy below the limit of detection
    /// compares IUPAC-compatibly ([`crate::lod`]).
    HaplotypeLod,
}

impl Measure {
    /// The number of measures.
    pub(crate) const COUNT: usize = 6;

    pub(crate) const ALL: [Self; Self::COUNT] = [
        Self::Primary,
        Self::Strict,
        Self::IupacCompatible,
        Self::Haplotype,
        Self::HaplotypeIupacCompatible,
        Self::HaplotypeLod,
    ];

    pub(crate) const fn label(self) -> &'static str {
        match self {
            Self::Primary => "primary",
            Self::Strict => "strict",
            Self::IupacCompatible => "iupac_compatible",
            Self::Haplotype => "haplotype",
            Self::HaplotypeIupacCompatible => "haplotype_iupac_compatible",
            Self::HaplotypeLod => "haplotype_lod",
        }
    }

    /// The name in Markdown tables.
    pub(crate) const fn title(self) -> &'static str {
        match self {
            Self::Primary => "Primary",
            Self::Strict => "Strict",
            Self::IupacCompatible => "IUPAC-compatible",
            Self::Haplotype => "Haplotype",
            Self::HaplotypeIupacCompatible => "Haplotype IUPAC-compatible",
            Self::HaplotypeLod => "Haplotype LoD",
        }
    }

    /// The measure's place in [`Measure::ALL`]: its declaration order.
    pub(crate) const fn index(self) -> usize {
        self as usize
    }

    /// The scores column that holds the measure's status.
    pub(crate) const fn column(self) -> &'static str {
        match self {
            Self::Primary => "status",
            Self::Strict => "strict_status",
            Self::IupacCompatible => "iupac_compatible_status",
            Self::Haplotype => "haplotype_status",
            Self::HaplotypeIupacCompatible => "haplotype_iupac_compatible_status",
            Self::HaplotypeLod => "haplotype_lod_status",
        }
    }

    /// Whether called bases match when their IUPAC codes share a base.
    const fn compatible(self) -> bool {
        matches!(self, Self::IupacCompatible | Self::HaplotypeIupacCompatible)
    }

    /// Whether alleles that describe the same sequence leave the comparison.
    pub(crate) const fn haplotype(self) -> bool {
        matches!(
            self,
            Self::Haplotype | Self::HaplotypeIupacCompatible | Self::HaplotypeLod
        )
    }

    /// Whether two alleles match under this measure, the truth IUPAC alleles
    /// at `below` (below the limit of detection) matching compatibly under
    /// `haplotype_lod`.
    fn matches(self, left: &Allele, right: &Allele, below: &BTreeSet<String>) -> bool {
        left.position == right.position
            && left.reference == right.reference
            && (left.alternate == right.alternate
                || ((self.lenient(left, below) || self.lenient(right, below))
                    && compatible(&left.alternate, &right.alternate)))
    }

    /// Whether `allele` compares IUPAC-compatibly under this measure.
    fn lenient(self, allele: &Allele, below: &BTreeSet<String>) -> bool {
        self.compatible()
            || (self == Self::HaplotypeLod && allele.is_iupac() && below.contains(&allele.position))
    }
}

/// The bases an IUPAC code admits; empty for anything else.
pub(crate) fn bases(code: &str) -> &'static [char] {
    match code {
        "A" => &['A'],
        "C" => &['C'],
        "G" => &['G'],
        "T" => &['T'],
        "R" => &['A', 'G'],
        "Y" => &['C', 'T'],
        "S" => &['C', 'G'],
        "W" => &['A', 'T'],
        "K" => &['G', 'T'],
        "M" => &['A', 'C'],
        "B" => &['C', 'G', 'T'],
        "D" => &['A', 'G', 'T'],
        "H" => &['A', 'C', 'T'],
        "V" => &['A', 'C', 'G'],
        "N" => &['A', 'C', 'G', 'T'],
        _ => &[],
    }
}

/// Whether two called bases share a possible base.
fn compatible(left: &str, right: &str) -> bool {
    let right = bases(right);
    bases(left).iter().any(|base| right.contains(base))
}

/// One side of a comparison: its alleles and intervals, closed and 1-based.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Side<'a> {
    pub(crate) alleles: &'a [Allele],
    pub(crate) intervals: &'a [(i64, i64)],
}

fn inside(position: &str, intervals: &[(i64, i64)]) -> bool {
    anchor(position).is_some_and(|base| {
        intervals
            .iter()
            .any(|(start, end)| (*start..=*end).contains(&base))
    })
}

/// The concordance of one sample.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) enum Concordance {
    /// `Y`.
    Concordant,
    /// `N`.
    Discordant,
    /// `N/A`: one side has no allele.
    NotApplicable,
}

/// A rule that makes a sample discordant.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) enum Reason {
    /// A position both sides call holds different alleles.
    SharedPosition,
    /// DNA alone calls an allele in an HV region.
    DnaOnly,
    /// The truth alone calls an allele in an HV region inside DNA's callable
    /// intervals.
    TruthOnly,
    /// Strict measure: the truth alone calls an allele in an HV region inside
    /// its analysed intervals but outside DNA's callable intervals.
    Uncovered,
}

impl Reason {
    pub(crate) const ALL: [Self; 4] = [
        Self::SharedPosition,
        Self::DnaOnly,
        Self::TruthOnly,
        Self::Uncovered,
    ];

    pub(crate) const fn label(self) -> &'static str {
        match self {
            Self::SharedPosition => "shared_position",
            Self::DnaOnly => "dna_only",
            Self::TruthOnly => "truth_only",
            Self::Uncovered => "uncovered",
        }
    }
}

/// The comparison of one sample.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct Comparison {
    pub(crate) concordance: Concordance,
    /// Distinct DNA alleles, special positions included.
    pub(crate) dna_alleles: usize,
    /// Distinct truth alleles, special positions included.
    pub(crate) truth_alleles: usize,
    /// DNA alleles absent from the truth, special positions left out.
    pub(crate) dna_unique: Vec<Allele>,
    /// Truth alleles absent from DNA, special positions left out.
    pub(crate) truth_unique: Vec<Allele>,
    /// Every rule that makes the sample discordant; empty unless `N`.
    pub(crate) reasons: BTreeSet<Reason>,
    /// The positions those rules name.
    pub(crate) positions: BTreeSet<String>,
    /// Whether every such position holds a truth IUPAC allele.
    pub(crate) iupac_only: bool,
    /// Truth-only alleles in an HV region inside the truth's analysed
    /// intervals but outside DNA's callable intervals: DNA's coverage gaps.
    pub(crate) uncovered: Vec<Allele>,
    /// Haplotype measures: the blocks where the two sequences differ that
    /// hold a discordant position, in reference order.
    pub(crate) blocks: Vec<Block>,
}

/// The intervals a haplotype measure compares as sequences: the HV regions
/// within DNA's callable intervals (merged where they touch), or the HV
/// regions when DNA reports none.
fn comparable(intervals: &[(i64, i64)]) -> Vec<(i64, i64)> {
    let regions = HV_REGIONS.map(|(_, start, end)| (start, end));
    if intervals.is_empty() {
        let mut regions = regions.to_vec();
        regions.sort_unstable();
        return regions;
    }
    let mut sorted = intervals.to_vec();
    sorted.sort_unstable();
    let mut merged = Vec::<(i64, i64)>::new();
    for (start, end) in sorted {
        match merged.last_mut() {
            Some(last) if start <= last.1.saturating_add(1) => last.1 = last.1.max(end),
            _ => merged.push((start, end)),
        }
    }
    let mut comparable = regions
        .iter()
        .flat_map(|&(region_start, region_end)| {
            merged.iter().filter_map(move |&(start, end)| {
                let (start, end) = (start.max(region_start), end.min(region_end));
                (start <= end).then_some((start, end))
            })
        })
        .collect::<Vec<_>>();
    comparable.sort_unstable();
    comparable
}

/// The allele as an edit of the reference, or `None` without an anchor base.
/// An insertion (`-` reference) needs an insertion ordinal and an inserted
/// base; any other allele replaces as many bases as its reference has. Its
/// symbols are `None` when it cannot be applied: an insertion without an
/// ordinal or a base, an ordinal on a reference base, or a reference or
/// alternate that is empty or not ASCII.
fn edit(allele: &Allele) -> Option<Edit<'_>> {
    let base = anchor(&allele.position)?;
    let ordinal = allele
        .position
        .split_once('.')
        .map(|(_, ordinal)| ordinal.parse::<u32>().ok());
    let readable = |text: &str| !text.is_empty() && text.is_ascii();
    let (ordinal, symbols) = match (ordinal, allele.reference == "-") {
        (Some(Some(ordinal)), true) if allele.alternate != "-" => {
            (Some(ordinal), Some(allele.alternate.as_str()))
        }
        (None, false) if allele.alternate == "-" => (None, Some("")),
        (None, false) => (None, Some(allele.alternate.as_str())),
        _ => (None, None),
    };
    let usable = readable(&allele.reference) && readable(&allele.alternate);
    Some(Edit {
        base,
        ordinal,
        span: i64::try_from(allele.reference.len()).unwrap_or(i64::MAX),
        symbols: symbols.filter(|_| usable),
        special: allele.is_special(),
        lenient: true,
    })
}

/// The reference bases an allele touches, from its anchor base.
fn bases_of(allele: &Allele) -> Option<(i64, i64)> {
    let edit = edit(allele)?;
    Some((edit.base, edit.last()))
}

/// Both sides' sequences over the comparable intervals (the HV regions
/// within DNA's callable intervals, or the HV regions when DNA reports none)
/// on `reference`, built once for the haplotype measures. With `below`, only
/// the truth IUPAC alleles at those positions write symbols that match
/// compatibly (`haplotype_lod`); without it, every truth allele does.
pub(crate) fn haplotypes(
    dna: Side<'_>,
    truth: Side<'_>,
    reference: &str,
    below: Option<&BTreeSet<String>>,
) -> Haplotypes {
    fn edits<'a>(alleles: &'a [Allele], lenient: &dyn Fn(&Allele) -> bool) -> Vec<Edit<'a>> {
        alleles
            .iter()
            .filter(|allele| allele.reference != allele.alternate)
            .filter_map(|allele| {
                edit(allele).map(|edit| Edit {
                    lenient: lenient(allele),
                    ..edit
                })
            })
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect()
    }
    let truth_lenient = |allele: &Allele| {
        below.is_none_or(|below| allele.is_iupac() && below.contains(&allele.position))
    };
    Haplotypes::new(
        reference,
        &comparable(dna.intervals),
        &edits(dna.alleles, &|_| true),
        &edits(truth.alleles, &truth_lenient),
    )
}

/// Compares DNA's side (its callable intervals) with the truth's side (its
/// analysed intervals) under `measure`; the haplotype measures read both
/// sides' sequences from `haplotypes`, and `haplotype_lod` compares the truth
/// IUPAC alleles at `below` compatibly.
pub(crate) fn compare(
    dna: Side<'_>,
    truth: Side<'_>,
    measure: Measure,
    haplotypes: &Haplotypes,
    below: &BTreeSet<String>,
) -> Comparison {
    let dna_all = distinct(dna.alleles);
    let truth_all = distinct(truth.alleles);
    // Special positions take part in the sequences, but not in the rules. A
    // truth symbol written by a truth allele also matches a compatible DNA
    // symbol, as two alleles at one position do under the IUPAC-compatible
    // rules; DNA's codes never match the truth's reference.
    let equivalence = measure.haplotype().then(|| {
        haplotypes.equivalence(
            (measure.compatible() || measure == Measure::HaplotypeLod)
                .then_some(compatible as fn(&str, &str) -> bool),
        )
    });
    let kept = |allele: &&Allele| {
        !allele.is_special()
            && !equivalence.as_ref().is_some_and(|equivalence| {
                bases_of(allele).is_some_and(|(first, last)| equivalence.explains(first, last))
            })
    };
    let dna_kept = dna_all.iter().filter(kept).collect::<Vec<_>>();
    let truth_kept = truth_all.iter().filter(kept).collect::<Vec<_>>();
    let unique = |left: &[&Allele], right: &[&Allele]| {
        left.iter()
            .filter(|allele| {
                !right
                    .iter()
                    .any(|other| measure.matches(allele, other, below))
            })
            .map(|allele| (*allele).clone())
            .collect::<Vec<_>>()
    };
    let dna_unique = unique(&dna_kept, &truth_kept);
    let dna_by_position = by_position(&dna_kept);
    let mut truth_unique = unique(&truth_kept, &dna_kept);
    if measure.compatible() || measure == Measure::HaplotypeLod {
        // Where DNA calls nothing but reads the position, its base is the
        // reference; that matches an IUPAC code that admits it.
        truth_unique.retain(|allele| {
            !(allele.is_iupac()
                && measure.lenient(allele, below)
                && !dna_by_position.contains_key(&allele.position)
                && inside(&allele.position, dna.intervals)
                && compatible(&allele.reference, &allele.alternate))
        });
    }

    let mut reasons = BTreeSet::new();
    let mut positions = BTreeSet::new();
    let truth_by_position = by_position(&truth_kept);
    for (position, alleles) in &dna_by_position {
        let Some(others) = truth_by_position.get(position) else {
            continue;
        };
        let covered = |left: &BTreeSet<&Allele>, right: &BTreeSet<&Allele>| {
            left.iter().all(|allele| {
                right
                    .iter()
                    .any(|other| measure.matches(allele, other, below))
            })
        };
        if !(covered(alleles, others) && covered(others, alleles)) {
            reasons.insert(Reason::SharedPosition);
            positions.insert(position.clone());
        }
    }
    for allele in dna_unique
        .iter()
        .filter(|allele| hv_region(&allele.position).is_some())
    {
        reasons.insert(Reason::DnaOnly);
        positions.insert(allele.position.clone());
    }
    let mut uncovered = Vec::new();
    for allele in truth_unique
        .iter()
        .filter(|allele| hv_region(&allele.position).is_some())
    {
        if dna.intervals.is_empty() || inside(&allele.position, dna.intervals) {
            reasons.insert(Reason::TruthOnly);
            positions.insert(allele.position.clone());
        } else if inside(&allele.position, truth.intervals) {
            uncovered.push(allele.clone());
            if measure == Measure::Strict {
                reasons.insert(Reason::Uncovered);
                positions.insert(allele.position.clone());
            }
        }
    }

    let concordance = if dna_all.is_empty() || truth_all.is_empty() {
        reasons.clear();
        positions.clear();
        Concordance::NotApplicable
    } else if reasons.is_empty() {
        Concordance::Concordant
    } else {
        Concordance::Discordant
    };
    let iupac_only = !positions.is_empty()
        && positions.iter().all(|position| {
            truth_kept
                .iter()
                .any(|allele| allele.position == *position && allele.is_iupac())
        });
    let blocks = equivalence
        .map(|equivalence| {
            equivalence
                .blocks()
                .filter(|block| {
                    positions.iter().any(|position| {
                        anchor(position).is_some_and(|base| block.overlaps(base, base))
                    })
                })
                .cloned()
                .collect()
        })
        .unwrap_or_default();
    Comparison {
        concordance,
        dna_alleles: dna_all.len(),
        truth_alleles: truth_all.len(),
        dna_unique,
        truth_unique,
        reasons,
        positions,
        iupac_only,
        uncovered,
        blocks,
    }
}

/// Distinct alleles that differ from the reference, in position order.
fn distinct(alleles: &[Allele]) -> Vec<Allele> {
    let mut kept = alleles
        .iter()
        .filter(|allele| allele.reference != allele.alternate)
        .cloned()
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect::<Vec<_>>();
    kept.sort_by(|left, right| left.order().cmp(&right.order()));
    kept
}

fn by_position<'a>(alleles: &[&'a Allele]) -> BTreeMap<String, BTreeSet<&'a Allele>> {
    let mut grouped = BTreeMap::<String, BTreeSet<&Allele>>::new();
    for allele in alleles {
        grouped
            .entry(allele.position.clone())
            .or_default()
            .insert(allele);
    }
    grouped
}

/// Alleles in the reviewer notation, space-separated, in position order.
pub(crate) fn texts(alleles: &[Allele]) -> String {
    let mut sorted = alleles.iter().collect::<Vec<_>>();
    sorted.sort_by(|left, right| left.order().cmp(&right.order()));
    sorted
        .iter()
        .map(|allele| allele.text())
        .collect::<Vec<_>>()
        .join(" ")
}

/// The window a discordant position falls in, for the breakdown.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) enum Region {
    /// HVS-III AC repeat, 513–526.
    AcRepeat,
    /// HVS-I poly-C tract, 16180–16195.
    Hv1PolyC,
    /// HVS-II poly-C tract, 302–316.
    Hv2PolyC,
    /// HVS-III poly-C tract, 560–576.
    Hv3PolyC,
    /// Anywhere else.
    Other,
}

impl Region {
    /// Every region, in label order.
    pub(crate) const ALL: [Self; 5] = [
        Self::AcRepeat,
        Self::Hv1PolyC,
        Self::Hv2PolyC,
        Self::Hv3PolyC,
        Self::Other,
    ];

    pub(crate) const fn label(self) -> &'static str {
        match self {
            Self::AcRepeat => "ac_repeat",
            Self::Hv1PolyC => "hv1_polyc",
            Self::Hv2PolyC => "hv2_polyc",
            Self::Hv3PolyC => "hv3_polyc",
            Self::Other => "other",
        }
    }

    /// The repeat windows with their inclusive 1-based bounds.
    pub(crate) const WINDOWS: [(Self, i64, i64); 4] = [
        (Self::Hv2PolyC, 302, 316),
        (Self::Hv1PolyC, 16180, 16195),
        (Self::Hv3PolyC, 560, 576),
        (Self::AcRepeat, 513, 526),
    ];

    /// The window of a position's anchor base.
    pub(crate) fn of(position: &str) -> Self {
        anchor(position).map_or(Self::Other, Self::at)
    }

    /// The window of a 1-based reference position.
    pub(crate) fn at(position: i64) -> Self {
        Self::WINDOWS
            .iter()
            .find(|(_, first, last)| (*first..=*last).contains(&position))
            .map_or(Self::Other, |(window, _, _)| *window)
    }
}

/// Why a sample has no comparison.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) enum Failure {
    /// A trace could not be read for duplicate removal.
    Trace { path: PathBuf, message: String },
    /// DNA returned an error.
    Dna(String),
    /// DNA panicked; this is a DNA defect.
    Panic,
    /// The target profile declares no notation, so there is nothing to compare.
    Notation,
    /// The truth has no record of the sample.
    TruthMissing,
    /// DNA analysed the sample on another reference sequence than the one
    /// the harness compares on.
    Reference { analysed: String, compared: String },
}

impl Failure {
    /// Stable label of the failure kind.
    pub(crate) const fn kind(&self) -> &'static str {
        match self {
            Self::Trace { .. } => "trace",
            Self::Dna(_) => "dna",
            Self::Panic => "panic",
            Self::Notation => "notation",
            Self::TruthMissing => "truth_missing",
            Self::Reference { .. } => "reference",
        }
    }
}

impl fmt::Display for Failure {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Trace { path, message } => {
                write!(formatter, "cannot read trace {}: {message}", path.display())
            }
            Self::Dna(message) => formatter.write_str(message),
            Self::Panic => formatter.write_str("DNA panicked"),
            Self::Notation => formatter.write_str("the target profile declares no notation"),
            Self::TruthMissing => formatter.write_str("the truth has no record of the sample"),
            Self::Reference { analysed, compared } => write!(
                formatter,
                "DNA analysed reference sequence {analysed}, the harness compares on {compared}"
            ),
        }
    }
}

impl std::error::Error for Failure {}

#[cfg(test)]
pub(crate) mod tests {
    use std::path::Path;
    use std::sync::OnceLock;

    use dna_kernel::model::reference::ReferenceTopology;

    use super::*;

    /// The rCRS sequence of the repository's reference.
    pub(crate) fn rcrs() -> &'static str {
        static SEQUENCE: OnceLock<String> = OnceLock::new();
        SEQUENCE.get_or_init(|| {
            let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../references/rCRS.fasta");
            match dna_kernel::reference::load(&path, ReferenceTopology::Linear) {
                Ok(reference) => reference.sequence,
                Err(error) => panic!("cannot load {}: {error}", path.display()),
            }
        })
    }

    /// Alleles from reviewer tokens: `263G`, `309.1C` (insertion), `249DEL`,
    /// with the reference base given after a colon where it matters.
    pub(crate) fn alleles(tokens: &[&str]) -> Vec<Allele> {
        tokens
            .iter()
            .map(|token| {
                let (token, reference) = token.split_once(':').unwrap_or((token, "N"));
                if let Some(position) = token.strip_suffix("DEL") {
                    Allele::new(position, reference, "-")
                } else {
                    let split = token.len() - 1;
                    let (position, alternate) = token.split_at(split);
                    if position.contains('.') {
                        Allele::new(position, "-", alternate)
                    } else {
                        Allele::new(position, reference, alternate)
                    }
                }
            })
            .collect()
    }

    const FULL: [(i64, i64); 3] = [(16024, 16365), (73, 340), (438, 576)];

    /// The primary measure, with the truth analysed over the full HV regions.
    fn primary(dna: &[Allele], intervals: &[(i64, i64)], truth: &[Allele]) -> Comparison {
        measure(dna, intervals, truth, Measure::Primary)
    }

    fn measure(
        dna: &[Allele],
        intervals: &[(i64, i64)],
        truth: &[Allele],
        measure: Measure,
    ) -> Comparison {
        compared(
            Side {
                alleles: dna,
                intervals,
            },
            Side {
                alleles: truth,
                intervals: &FULL,
            },
            measure,
        )
    }

    /// `compare` with both sides' sequences built on rCRS, every truth IUPAC
    /// allele at or above the limit of detection.
    fn compared(dna: Side<'_>, truth: Side<'_>, measure: Measure) -> Comparison {
        let below = BTreeSet::new();
        let lod = (measure == Measure::HaplotypeLod).then_some(&below);
        compare(
            dna,
            truth,
            measure,
            &haplotypes(dna, truth, rcrs(), lod),
            &below,
        )
    }

    #[test]
    fn identical_alleles_are_concordant() {
        let calls = alleles(&["263G", "309.1C", "315.1C", "16189C"]);
        let comparison = primary(&calls, &FULL, &calls);
        assert_eq!(comparison.concordance, Concordance::Concordant);
        assert!(comparison.dna_unique.is_empty() && comparison.reasons.is_empty());
    }

    #[test]
    fn an_empty_side_is_not_applicable() {
        let calls = alleles(&["263G"]);
        assert_eq!(
            primary(&[], &FULL, &calls).concordance,
            Concordance::NotApplicable
        );
        assert_eq!(
            primary(&calls, &FULL, &[]).concordance,
            Concordance::NotApplicable
        );
        // An allele equal to the reference is not an allele.
        let reference_state = [Allele::new("263", "A", "A")];
        assert_eq!(
            primary(&reference_state, &FULL, &calls).concordance,
            Concordance::NotApplicable
        );
    }

    #[test]
    fn special_indels_are_left_out_but_special_snps_are_not() {
        // 309.2C, 455.1T, 573DEL are special; 309T is a substitution.
        let dna = alleles(&["263G", "309.1C", "309.2C", "455.1T", "573DEL"]);
        let truth = alleles(&["263G", "309.1C"]);
        assert_eq!(
            primary(&dna, &FULL, &truth).concordance,
            Concordance::Concordant
        );
        let truth = alleles(&["263G", "309T"]);
        let comparison = primary(&alleles(&["263G"]), &FULL, &truth);
        assert_eq!(comparison.concordance, Concordance::Discordant);
        assert_eq!(comparison.reasons, BTreeSet::from([Reason::TruthOnly]));
        // The positions `comparison.py` never compares (309.1, 309.2, 573.1)
        // hold insertions, which are special, so differing alleles there stay
        // concordant without a rule of their own.
        let dna = alleles(&["263G", "309.1T", "309.2T", "573.1T"]);
        let truth = alleles(&["263G", "309.1C", "309.2C", "573.1C"]);
        let comparison = primary(&dna, &FULL, &truth);
        assert_eq!(comparison.concordance, Concordance::Concordant);
        assert!(comparison.positions.is_empty());
    }

    #[test]
    fn a_shared_position_with_other_alleles_is_discordant() {
        let comparison = primary(
            &alleles(&["263G", "16189C"]),
            &FULL,
            &alleles(&["263G", "16189Y"]),
        );
        assert_eq!(comparison.concordance, Concordance::Discordant);
        assert!(comparison.reasons.contains(&Reason::SharedPosition));
        assert_eq!(comparison.positions, BTreeSet::from(["16189".to_owned()]));
        assert!(comparison.iupac_only);
        assert_eq!(comparison.dna_unique, alleles(&["16189C"]));
        assert_eq!(comparison.truth_unique, alleles(&["16189Y"]));
    }

    #[test]
    fn dna_only_alleles_count_only_inside_hv_regions() {
        let truth = alleles(&["263G"]);
        // 750G lies outside every HV region.
        assert_eq!(
            primary(&alleles(&["263G", "750G"]), &FULL, &truth).concordance,
            Concordance::Concordant
        );
        // 16390A lies outside HV1 (16024-16365); 16300A inside.
        assert_eq!(
            primary(&alleles(&["263G", "16390A"]), &FULL, &truth).concordance,
            Concordance::Concordant
        );
        let comparison = primary(&alleles(&["263G", "16300A"]), &FULL, &truth);
        assert_eq!(comparison.reasons, BTreeSet::from([Reason::DnaOnly]));
        assert!(!comparison.iupac_only);
    }

    #[test]
    fn truth_only_alleles_count_inside_dna_intervals() {
        let dna = alleles(&["263G"]);
        let truth = alleles(&["263G", "16300A"]);
        // Outside DNA's callable intervals: ignored.
        assert_eq!(
            primary(&dna, &[(73, 340)], &truth).concordance,
            Concordance::Concordant
        );
        // Inside, or DNA reports no interval: discordant.
        assert_eq!(
            primary(&dna, &[(73, 340), (16024, 16365)], &truth).concordance,
            Concordance::Discordant
        );
        assert_eq!(
            primary(&dna, &[], &truth).concordance,
            Concordance::Discordant
        );
    }

    #[test]
    fn insertions_belong_to_their_anchor_base() {
        let dna = alleles(&["263G"]);
        let truth = alleles(&["263G", "16193.1C"]);
        assert_eq!(
            primary(&dna, &[(73, 340), (16024, 16193)], &truth).concordance,
            Concordance::Discordant
        );
        assert_eq!(
            primary(&dna, &[(73, 340), (16024, 16192)], &truth).concordance,
            Concordance::Concordant
        );
    }

    #[test]
    fn references_take_part_in_the_allele() {
        let comparison = primary(
            &alleles(&["263G:A", "249DEL:A"]),
            &FULL,
            &alleles(&["263G:C", "249DEL:A"]),
        );
        assert_eq!(comparison.concordance, Concordance::Discordant);
        assert_eq!(comparison.positions, BTreeSet::from(["263".to_owned()]));
    }

    #[test]
    fn normalises_positions_and_renders_notation() {
        assert_eq!(normalize_position(" 16,193 "), "16193");
        assert_eq!(normalize_position("16193.0"), "16193");
        assert_eq!(normalize_position("309.1"), "309.1");
        assert_eq!(normalize_position("309.10"), "309.10");
        assert_eq!(normalize_position("x1"), "x1");
        let calls = alleles(&["16189C", "249DEL", "309.2C", "309.1C", "73G"]);
        assert_eq!(texts(&calls), "73G 249DEL 309.1C 309.2C 16189C");
        assert_eq!(Region::of("309.1"), Region::Hv2PolyC);
        assert_eq!(Region::of("x"), Region::Other);
    }

    #[test]
    fn the_strict_measure_counts_truth_alleles_dna_does_not_cover() {
        let dna = alleles(&["263G"]);
        let truth = alleles(&["263G", "16300A"]);
        let intervals = [(73, 340)];
        let primary = measure(&dna, &intervals, &truth, Measure::Primary);
        assert_eq!(primary.concordance, Concordance::Concordant);
        assert_eq!(primary.uncovered, alleles(&["16300A"]));
        let strict = measure(&dna, &intervals, &truth, Measure::Strict);
        assert_eq!(strict.concordance, Concordance::Discordant);
        assert_eq!(strict.reasons, BTreeSet::from([Reason::Uncovered]));
        // Outside the truth's analysed intervals it is not missed either.
        let narrow = compared(
            Side {
                alleles: &dna,
                intervals: &intervals,
            },
            Side {
                alleles: &truth,
                intervals: &[(73, 340), (16024, 16200)],
            },
            Measure::Strict,
        );
        assert_eq!(narrow.concordance, Concordance::Concordant);
        assert!(narrow.uncovered.is_empty());
    }

    #[test]
    fn the_iupac_compatible_measure_matches_a_base_the_code_admits() {
        let truth = alleles(&["263G", "16189Y", "16311R"]);
        // C is one of Y's bases; A is one of R's.
        let dna = alleles(&["263G", "16189C", "16311A"]);
        assert_eq!(
            measure(&dna, &FULL, &truth, Measure::Primary).concordance,
            Concordance::Discordant
        );
        assert_eq!(
            measure(&dna, &FULL, &truth, Measure::IupacCompatible).concordance,
            Concordance::Concordant
        );
        // T is not one of R's bases, and a deletion matches no base.
        let dna = alleles(&["263G", "16189C", "16311T"]);
        assert_eq!(
            measure(&dna, &FULL, &truth, Measure::IupacCompatible).concordance,
            Concordance::Discordant
        );
        // No DNA allele where DNA reads the reference: the reference A is one
        // of R's bases, C is not; outside DNA's intervals nothing is known.
        let dna = alleles(&["263G", "16189C"]);
        let truth = alleles(&["263G", "16189Y", "16311R:A"]);
        assert_eq!(
            measure(&dna, &FULL, &truth, Measure::IupacCompatible).concordance,
            Concordance::Concordant
        );
        let truth = alleles(&["263G", "16189Y", "16311R:C"]);
        assert_eq!(
            measure(&dna, &FULL, &truth, Measure::IupacCompatible).concordance,
            Concordance::Discordant
        );
        assert_eq!(
            measure(
                &dna,
                &[],
                &alleles(&["263G", "16311R:A"]),
                Measure::IupacCompatible
            )
            .concordance,
            Concordance::Discordant
        );
        assert!(!compatible("-", "N"));
        assert!(compatible("N", "Y"));
        assert!(!compatible("X", "A"));
    }

    /// Concordance under primary, haplotype, and haplotype IUPAC-compatible.
    fn concordances(dna: &[&str], intervals: &[(i64, i64)], truth: &[&str]) -> [Concordance; 3] {
        let (dna, truth) = (alleles(dna), alleles(truth));
        [
            Measure::Primary,
            Measure::Haplotype,
            Measure::HaplotypeIupacCompatible,
        ]
        .map(|of| measure(&dna, intervals, &truth, of).concordance)
    }

    const Y: Concordance = Concordance::Concordant;
    const N: Concordance = Concordance::Discordant;

    #[test]
    fn the_haplotype_measures_accept_other_placements_of_one_sequence() {
        for (dna, truth) in [
            // HVS-I poly-C: an insertion anywhere in the C run 16190-16193.
            (&["263G", "16190.1C"][..], &["263G", "16193.1C"][..]),
            // HVS-II: a C after 310T lengthens the run 311-315.
            (&["263G", "310.1C"], &["263G", "315.1C"]),
            // A deletion in AA (248-249).
            (&["249DEL", "263G"], &["248DEL", "263G"]),
            // The AC repeat 514-524 loses one AC unit.
            (&["263G", "514DEL", "515DEL"], &["263G", "523DEL", "524DEL"]),
            // An insertion and a deletion against substitutions (the HVS-I
            // nomenclature form 16183C 16184A 16189C).
            (
                &["263G", "16182.1C", "16189C", "16193DEL"],
                &["263G", "16183C", "16184A", "16189C"],
            ),
            // The HVS-III C run 568-573: 573.1C is special, 570.1C is not.
            (&["263G", "570.1C"], &["263G", "573.1C"]),
            // A T after the special 455 with a C lost, against 456T.
            (&["263G", "455.1T", "459DEL"], &["263G", "456T"]),
        ] {
            assert_eq!(
                concordances(dna, &FULL, truth),
                [N, Y, Y],
                "{dna:?} {truth:?}"
            );
        }
    }

    #[test]
    fn the_haplotype_measures_keep_real_differences() {
        // One more C before the T and one fewer after it: equal lengths, a
        // different sequence.
        let dna = alleles(&["263G", "16188.1C", "16188.2T"]);
        let truth = alleles(&["263G", "16190T", "16193.1C", "16193.2C"]);
        let comparison = measure(&dna, &FULL, &truth, Measure::Haplotype);
        assert_eq!(comparison.concordance, N);
        assert_eq!(
            comparison
                .blocks
                .iter()
                .map(Block::text)
                .collect::<Vec<_>>(),
            ["16188-16193:CCTTCCCC/CTTCCCCC"]
        );
        assert_eq!(comparison.dna_unique, alleles(&["16188.1C", "16188.2T"]));
        // An insertion against a substitution that lengthens the run is
        // another sequence.
        assert_eq!(
            concordances(&["263G", "16194C"], &FULL, &["263G", "16193.1C"]),
            [N, N, N]
        );
        // A truth IUPAC code where DNA calls a base it admits, in another
        // placement: only the IUPAC-compatible haplotype measure accepts it.
        assert_eq!(
            concordances(
                &["263G", "16190.1T"],
                &FULL,
                &["263G", "16191Y", "16193.1C"]
            ),
            [N, N, Y]
        );
    }

    #[test]
    fn the_haplotype_measures_keep_the_base_rules_elsewhere() {
        // Concordant samples stay concordant, special-position lengths
        // included; an empty side stays not applicable.
        assert_eq!(
            concordances(
                &["263G", "309.1C", "315.1C"],
                &FULL,
                &["263G", "309.1C", "309.2C", "315.1C"]
            ),
            [Y, Y, Y]
        );
        assert_eq!(
            concordances(&[], &FULL, &["263G"]),
            [Concordance::NotApplicable; 3]
        );
        // A DNA-only allele outside DNA's callable intervals is not judged
        // as a sequence, and counts as on the primary measure.
        let comparison = measure(
            &alleles(&["263G", "16300A"]),
            &[(73, 340)],
            &alleles(&["263G"]),
            Measure::Haplotype,
        );
        assert_eq!(comparison.reasons, BTreeSet::from([Reason::DnaOnly]));
        assert!(comparison.blocks.is_empty());
        // A different reference written for the same base is the same
        // sequence.
        assert_eq!(concordances(&["263G:A"], &FULL, &["263G:C"]), [N, Y, Y]);
    }

    #[test]
    fn a_callable_interval_edge_leaves_a_split_length_change_open() {
        // DNA reads up to 16191 and places the extra C at 16190; the truth's
        // 16193.1C lies outside DNA's intervals and is ignored, so the
        // difference cannot close inside the comparable interval.
        let intervals = [(73, 340), (16024, 16191)];
        let comparison = measure(
            &alleles(&["263G", "16190.1C"]),
            &intervals,
            &alleles(&["263G", "16193.1C"]),
            Measure::Haplotype,
        );
        assert_eq!(comparison.concordance, N);
        assert_eq!(
            comparison
                .blocks
                .iter()
                .map(Block::text)
                .collect::<Vec<_>>(),
            ["16190-16191+:CCC/CC"]
        );
        assert_eq!(
            concordances(&["263G", "16190.1C"], &FULL, &["263G", "16193.1C"]),
            [N, Y, Y]
        );
    }

    #[test]
    fn measures_are_declared_in_their_order() {
        for (index, measure) in Measure::ALL.into_iter().enumerate() {
            assert_eq!(measure.index(), index);
        }
    }

    #[test]
    fn dna_codes_never_match_the_truth_reference() {
        // DNA's Y or N where the truth writes nothing: discordant on the
        // IUPAC-compatible rules, so on the haplotype measures too.
        for dna in ["16304Y:T", "16304N:T"] {
            let dna = alleles(&["263G", dna]);
            let truth = alleles(&["263G"]);
            for of in Measure::ALL {
                assert_eq!(measure(&dna, &FULL, &truth, of).concordance, N, "{of:?}");
            }
        }
        // Both write the position: codes that share a base match on both
        // IUPAC-compatible measures, as two alleles do.
        assert_eq!(
            concordances(&["263G", "16304Y:T"], &FULL, &["263G", "16304C:T"]),
            [N, N, Y]
        );
    }

    #[test]
    fn several_reference_bases_are_one_edit() {
        let with = |mut alleles: Vec<Allele>, allele: Allele| {
            alleles.push(allele);
            alleles
        };
        let deletion = with(alleles(&["263G"]), Allele::new("514", "CA", "-"));
        let haplotype = |dna: &[&str], truth: &[Allele]| {
            measure(&alleles(dna), &FULL, truth, Measure::Haplotype).concordance
        };
        // Two bases deleted at once are the two one-base deletions, here or
        // one AC unit later, but not one of them.
        assert_eq!(haplotype(&["263G", "514DEL", "515DEL"], &deletion), Y);
        assert_eq!(haplotype(&["263G", "523DEL", "524DEL"], &deletion), Y);
        assert_eq!(haplotype(&["263G", "514DEL"], &deletion), N);
        // Two substituted bases are the two substitutions.
        let mnv = with(alleles(&["263G"]), Allele::new("16300", "AC", "GT"));
        assert_eq!(haplotype(&["263G", "16300G", "16301T"], &mnv), Y);
        assert_eq!(haplotype(&["263G", "16300G"], &mnv), N);
    }

    #[test]
    fn a_special_length_difference_does_not_hide_a_placement() {
        // The 303-309 run differs by a special 309.2C, and the C after 310T
        // is written at another place. The special length is left out of
        // the rules, so only the placement remains: one sequence.
        assert_eq!(
            concordances(
                &["263G", "309.1C", "315.1C"],
                &FULL,
                &["263G", "309.1C", "309.2C", "310.1C"]
            ),
            [N, Y, Y]
        );
        // A real difference after it stays one.
        assert_eq!(
            concordances(
                &["263G", "309.1C", "315.1C", "320T"],
                &FULL,
                &["263G", "309.1C", "309.2C", "310.1C"]
            ),
            [N, N, N]
        );
    }

    #[test]
    fn comparable_intervals_are_hv_regions_within_callable_intervals() {
        // The tests' reference loaded.
        assert_eq!(rcrs().len(), 16_569);
        assert_eq!(comparable(&[]), [(73, 340), (438, 576), (16024, 16365)]);
        assert_eq!(
            comparable(&[(16300, 16400), (1, 100), (101, 200), (500, 600)]),
            [(73, 200), (500, 576), (16300, 16365)]
        );
    }

    /// Per-base alleles of anchored variants, as DNA's notation writes them.
    fn per_base(variants: &[dna_kernel::variant::Variant]) -> Vec<Allele> {
        use dna_kernel::variant::VariantKind;
        let mut alleles = Vec::new();
        for variant in variants {
            let position = variant.position_1based;
            match variant.kind {
                VariantKind::Ins => {
                    for (ordinal, base) in variant.alternate.chars().skip(1).enumerate() {
                        let position = format!("{position}.{}", ordinal + 1);
                        alleles.push(Allele::new(&position, "-", &base.to_string()));
                    }
                }
                VariantKind::Del => {
                    for (offset, base) in variant.reference.chars().enumerate().skip(1) {
                        let position = (position + offset).to_string();
                        alleles.push(Allele::new(&position, &base.to_string(), "-"));
                    }
                }
                _ => alleles.push(Allele::new(
                    &position.to_string(),
                    &variant.reference,
                    &variant.alternate,
                )),
            }
        }
        alleles
    }

    /// DNA's post-calling plugins as the oracle: every single-base edit in
    /// and around the HVS-II and HVS-III repeat windows is normalized and
    /// represented by `dna-normalize`, which proves the represented variants
    /// rebuild the source haplotype. The haplotype comparison must find the
    /// source and represented alleles the same sequence, and any two edits
    /// the same exactly when `dna-normalize`'s rebuilt sequences are equal.
    #[test]
    fn agrees_with_dna_normalize_on_sequence_equivalence() -> Result<(), Box<dyn std::error::Error>>
    {
        use dna_kernel::profile::Profile;
        use dna_kernel::variant::{CalledVariantSet, Variant, VariantKind};
        use dna_normalize::variant_nomenclature::{apply_with, from_normalization};
        use dna_normalize::variant_normalization::{NormalizationPolicy, normalize_with};

        let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
        let reference = dna_kernel::reference::load(
            &root.join("references/rCRS.fasta"),
            ReferenceTopology::Linear,
        )?;
        let profile = Profile::load(&root.join("config/profiles/human-mtdna-rcrs.toml"))?;
        let sequence = reference.sequence.as_str();
        let base = |position: usize| sequence[position - 1..position].to_owned();
        let mut edits = Vec::new();
        for (start, end) in [(300, 318), (510, 527)] {
            for position in start..=end {
                for inserted in ["A", "C", "G", "T"] {
                    edits.push(Variant {
                        contig: reference.name.clone(),
                        position_1based: position,
                        reference: base(position),
                        alternate: format!("{}{inserted}", base(position)),
                        kind: VariantKind::Ins,
                    });
                    if base(position) != inserted {
                        edits.push(Variant {
                            contig: reference.name.clone(),
                            position_1based: position,
                            reference: base(position),
                            alternate: inserted.to_owned(),
                            kind: VariantKind::Snv,
                        });
                    }
                }
                edits.push(Variant {
                    contig: reference.name.clone(),
                    position_1based: position - 1,
                    reference: format!("{}{}", base(position - 1), base(position)),
                    alternate: base(position - 1),
                    kind: VariantKind::Del,
                });
            }
        }
        // Intervals around the edits; the sequences agree outside them.
        let intervals = [(295, 325), (505, 532)];
        let built = |left: &[Allele], right: &[Allele]| {
            let side = |alleles| Side {
                alleles,
                intervals: &intervals,
            };
            haplotypes(side(left), side(right), sequence, None).equivalence(None)
        };
        let mut rebuilt = Vec::new();
        let mut moved = 0;
        for edit in &edits {
            let called = CalledVariantSet {
                reference: reference.identity(),
                variants: vec![edit.clone()],
            };
            let normalized =
                normalize_with(&reference, &called, NormalizationPolicy::RightAligned)?;
            let represented = apply_with(&reference, &profile, from_normalization(&normalized))?;
            let source = per_base(&represented.source_variants);
            let target = per_base(&represented.represented_variants);
            moved += usize::from(source != target);
            let equivalence = built(&source, &target);
            assert_eq!(equivalence.blocks().count(), 0, "{edit:?}");
            for allele in source.iter().chain(&target) {
                let (first, last) = bases_of(allele).ok_or("an anchor base")?;
                assert!(equivalence.explains(first, last), "{allele:?}");
            }
            rebuilt.push((source, represented.alternate_sequence));
        }
        assert!(moved > 50, "{moved} edits change their representation");
        for (left, left_sequence) in &rebuilt {
            for (right, right_sequence) in &rebuilt {
                assert_eq!(
                    built(left, right).blocks().count() == 0,
                    left_sequence == right_sequence,
                    "{left:?} {right:?}"
                );
            }
        }
        Ok(())
    }
}
