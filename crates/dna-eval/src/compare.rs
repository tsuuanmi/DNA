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
//! - `N` when a position both sides call holds different allele sets
//!   (except 309.1, 309.2, 573.1);
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

use std::collections::{BTreeMap, BTreeSet};
use std::fmt;
use std::path::PathBuf;

/// Positions whose insertions and deletions are left out of the comparison.
const SPECIAL_POSITIONS: [&str; 4] = ["455", "463", "573", "309"];
/// Shared positions never compared.
const ALWAYS_CONCORDANT: [&str; 3] = ["573.1", "309.1", "309.2"];
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

    fn is_iupac(&self) -> bool {
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
fn anchor(position: &str) -> Option<i64> {
    position.split('.').next()?.parse::<i64>().ok()
}

fn hv_region(position: &str) -> Option<&'static str> {
    let base = anchor(position)?;
    HV_REGIONS
        .iter()
        .find(|(_, start, end)| (*start..=*end).contains(&base))
        .map(|(name, _, _)| *name)
}

/// A set of comparison rules.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) enum Measure {
    /// The lab's rules.
    Primary,
    /// The lab's rules, with truth-only alleles also counted inside the
    /// truth's analysed intervals.
    Strict,
    /// The lab's rules, with IUPAC-compatible called bases matching.
    IupacCompatible,
}

impl Measure {
    pub(crate) const ALL: [Self; 3] = [Self::Primary, Self::Strict, Self::IupacCompatible];

    pub(crate) const fn label(self) -> &'static str {
        match self {
            Self::Primary => "primary",
            Self::Strict => "strict",
            Self::IupacCompatible => "iupac_compatible",
        }
    }

    /// The measure's place in [`Measure::ALL`].
    pub(crate) const fn index(self) -> usize {
        match self {
            Self::Primary => 0,
            Self::Strict => 1,
            Self::IupacCompatible => 2,
        }
    }

    /// The scores column that holds the measure's status.
    pub(crate) const fn column(self) -> &'static str {
        match self {
            Self::Primary => "status",
            Self::Strict => "strict_status",
            Self::IupacCompatible => "iupac_compatible_status",
        }
    }

    /// Whether two alleles match under this measure.
    fn matches(self, left: &Allele, right: &Allele) -> bool {
        left.position == right.position
            && left.reference == right.reference
            && (left.alternate == right.alternate
                || (self == Self::IupacCompatible && compatible(&left.alternate, &right.alternate)))
    }
}

/// The bases an IUPAC code admits; empty for anything else.
fn bases(code: &str) -> &'static [char] {
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
}

/// Compares DNA's side (its callable intervals) with the truth's side (its
/// analysed intervals) under `measure`.
pub(crate) fn compare(dna: Side<'_>, truth: Side<'_>, measure: Measure) -> Comparison {
    let dna_all = distinct(dna.alleles);
    let truth_all = distinct(truth.alleles);
    let dna_kept = dna_all
        .iter()
        .filter(|allele| !allele.is_special())
        .collect::<Vec<_>>();
    let truth_kept = truth_all
        .iter()
        .filter(|allele| !allele.is_special())
        .collect::<Vec<_>>();
    let unique = |left: &[&Allele], right: &[&Allele]| {
        left.iter()
            .filter(|allele| !right.iter().any(|other| measure.matches(allele, other)))
            .map(|allele| (*allele).clone())
            .collect::<Vec<_>>()
    };
    let dna_unique = unique(&dna_kept, &truth_kept);
    let dna_by_position = by_position(&dna_kept);
    let mut truth_unique = unique(&truth_kept, &dna_kept);
    if measure == Measure::IupacCompatible {
        // Where DNA calls nothing but reads the position, its base is the
        // reference; that matches an IUPAC code that admits it.
        truth_unique.retain(|allele| {
            !(allele.is_iupac()
                && !dna_by_position.contains_key(&allele.position)
                && inside(&allele.position, dna.intervals)
                && compatible(&allele.reference, &allele.alternate))
        });
    }

    let mut reasons = BTreeSet::new();
    let mut positions = BTreeSet::new();
    let truth_by_position = by_position(&truth_kept);
    for (position, alleles) in &dna_by_position {
        if ALWAYS_CONCORDANT.contains(&position.as_str()) {
            continue;
        }
        let Some(others) = truth_by_position.get(position) else {
            continue;
        };
        let covered = |left: &BTreeSet<&Allele>, right: &BTreeSet<&Allele>| {
            left.iter()
                .all(|allele| right.iter().any(|other| measure.matches(allele, other)))
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

    /// The window of a position's anchor base.
    pub(crate) fn of(position: &str) -> Self {
        match anchor(position) {
            Some(302..=316) => Self::Hv2PolyC,
            Some(16180..=16195) => Self::Hv1PolyC,
            Some(560..=576) => Self::Hv3PolyC,
            Some(513..=526) => Self::AcRepeat,
            _ => Self::Other,
        }
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
        }
    }
}

impl std::error::Error for Failure {}

#[cfg(test)]
mod tests {
    use super::*;

    /// Alleles from reviewer tokens: `263G`, `309.1C` (insertion), `249DEL`,
    /// with the reference base given after a colon where it matters.
    fn alleles(tokens: &[&str]) -> Vec<Allele> {
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
        compare(
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
        let narrow = compare(
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
}
