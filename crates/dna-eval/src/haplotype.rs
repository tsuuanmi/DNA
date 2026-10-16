//! Haplotype equivalence of two edit sets over comparable intervals.
//!
//! Each side's edits are applied to the reference inside an interval: a
//! substitution replaces its bases, a deletion removes them, and insertions
//! (`309.1`, `309.2`, …) follow their anchor base in ordinal order. Symbols are
//! kept as written, so an IUPAC allele stays an IUPAC code in the sequence.
//!
//! The interval is walked base by base. Wherever both local sequences have
//! the same length after a base, they can be cut there: the two sequences
//! are equal exactly when every piece between consecutive cuts is equal.
//! Pieces that differ are the [`Block`]s; edits anchored in the other pieces
//! are written differently but describe the same sequence. A piece that
//! reaches the end of its interval with unequal lengths is `open`: the
//! difference may continue past the interval, where one side has no callable
//! base.
//!
//! Special edits (indels the comparison rules leave out) must not decide a
//! difference by themselves. Each interval is therefore built three ways:
//! each side with its own special edits, both with DNA's, and both with the
//! truth's. One way decides the whole interval: the first that leaves no
//! block over any non-special edit, or else the first with the fewest blocks
//! over one; an edit is explained when that way finds no block over it.
//!
//! Symbols are compared exactly or, when the caller gives a `compatible`
//! test, a truth symbol written by a lenient truth allele also matches a DNA
//! symbol it is compatible with. An interval where every way has two edits at one
//! base or insertion ordinal, or an edit that cannot be applied or leaves the
//! interval, is not decided.

use std::collections::BTreeMap;

/// One allele as an edit of the reference.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) struct Edit<'a> {
    /// Anchor base, 1-based.
    pub(crate) base: i64,
    /// Insertion ordinal after the anchor (`309.2` → 2); `None` for an edit
    /// of reference bases.
    pub(crate) ordinal: Option<u32>,
    /// Reference bases the edit replaces from `base` on; unused for
    /// insertions.
    pub(crate) span: i64,
    /// The ASCII symbols written there, empty for deleted bases; `None` when
    /// the allele cannot be applied.
    pub(crate) symbols: Option<&'a str>,
    /// Whether the comparison rules leave the edit out.
    pub(crate) special: bool,
    /// Whether a truth edit's symbols may match a compatible DNA symbol under
    /// an IUPAC-compatible symbol test.
    pub(crate) lenient: bool,
}

impl Edit<'_> {
    /// The last reference base the edit touches.
    pub(crate) const fn last(&self) -> i64 {
        if self.ordinal.is_some() || self.span < 1 {
            self.base
        } else {
            self.base + self.span - 1
        }
    }
}

/// Shown symbols of a block sequence; longer ones are shortened.
const SHOWN: usize = 24;
/// Symbols kept at each end of a shortened sequence.
const KEPT: usize = 10;

/// A stretch of the reference over which the two haplotypes differ.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct Block {
    /// First anchor base, 1-based.
    pub(crate) start: i64,
    /// Last anchor base, 1-based; insertions after it belong to the block.
    pub(crate) end: i64,
    /// DNA's local sequence.
    pub(crate) dna: String,
    /// The truth's local sequence.
    pub(crate) truth: String,
    /// Per truth symbol, whether a lenient truth allele wrote it.
    truth_written: Vec<bool>,
    /// The block reaches the end of its interval with unequal lengths.
    pub(crate) open: bool,
}

impl Block {
    /// `start-end:dna/truth`, `+` after `end` when open, `-` for an empty
    /// sequence, and a sequence over 24 symbols as its first and last 10 and
    /// its length: `CCCTCCCCCC..TTCTGGCCAC(40)`.
    pub(crate) fn text(&self) -> String {
        let sequence = |text: &str| match text.len() {
            0 => "-".to_owned(),
            length if length <= SHOWN => text.to_owned(),
            length => format!("{}..{}({length})", &text[..KEPT], &text[length - KEPT..]),
        };
        format!(
            "{}-{}{}:{}/{}",
            self.start,
            self.end,
            if self.open { "+" } else { "" },
            sequence(&self.dna),
            sequence(&self.truth)
        )
    }

    /// Whether the block holds any base of `first..=last`.
    pub(crate) fn overlaps(&self, first: i64, last: i64) -> bool {
        first <= self.end && self.start <= last
    }

    /// Whether the block differs under `compatible`.
    fn differs(&self, compatible: Option<fn(&str, &str) -> bool>) -> bool {
        self.open
            || self.dna.len() != self.truth.len()
            || (0..self.dna.len()).any(|index| {
                let (dna, truth) = (&self.dna[index..=index], &self.truth[index..=index]);
                dna != truth
                    && !(self.truth_written[index]
                        && compatible.is_some_and(|compatible| compatible(dna, truth)))
            })
    }
}

/// One comparable interval: the bases of its non-special edits on either
/// side, and the pieces whose symbols differ, per way it could be built
/// (each side's own special edits first).
#[derive(Debug, Clone, PartialEq, Eq)]
struct Interval {
    start: i64,
    end: i64,
    edits: Vec<(i64, i64)>,
    ways: Vec<Vec<Block>>,
}

/// The edits whose `special` flag is `keep`.
fn special<'a>(edits: &[Edit<'a>], keep: bool) -> Vec<Edit<'a>> {
    edits
        .iter()
        .filter(|edit| edit.special == keep)
        .copied()
        .collect()
}

/// Both sides' sequences over the comparable intervals, built once for every
/// measure.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub(crate) struct Haplotypes {
    intervals: Vec<Interval>,
}

impl Haplotypes {
    /// Builds both sides over each of `intervals` (closed, 1-based, disjoint,
    /// in order) on `reference`.
    pub(crate) fn new(
        reference: &str,
        intervals: &[(i64, i64)],
        dna: &[Edit<'_>],
        truth: &[Edit<'_>],
    ) -> Self {
        let intervals = intervals
            .iter()
            .filter_map(|&(start, end)| {
                let within = |edit: &&Edit<'_>| (start..=end).contains(&edit.base);
                let dna = dna.iter().filter(within).copied().collect::<Vec<_>>();
                let truth = truth.iter().filter(within).copied().collect::<Vec<_>>();
                let (dna_special, truth_special) = (special(&dna, true), special(&truth, true));
                let mut ways = vec![(dna.clone(), Vec::new(), truth.clone(), Vec::new())];
                if dna_special != truth_special {
                    ways.push((
                        dna.clone(),
                        Vec::new(),
                        special(&truth, false),
                        dna_special.clone(),
                    ));
                    ways.push((
                        special(&dna, false),
                        truth_special,
                        truth.clone(),
                        Vec::new(),
                    ));
                }
                let edits = special(&dna, false)
                    .iter()
                    .chain(&special(&truth, false))
                    .map(|edit| (edit.base, edit.last()))
                    .collect();
                let ways = ways
                    .into_iter()
                    .filter_map(|(dna, dna_borrowed, truth, truth_borrowed)| {
                        pieces(
                            reference,
                            (start, end),
                            [(&dna, &dna_borrowed), (&truth, &truth_borrowed)],
                        )
                    })
                    .collect::<Vec<_>>();
                (!ways.is_empty()).then_some(Interval {
                    start,
                    end,
                    edits,
                    ways,
                })
            })
            .collect();
        Self { intervals }
    }

    /// The differing blocks under `compatible` (exact symbols when `None`) of
    /// the way that decides each interval.
    pub(crate) fn equivalence(&self, compatible: Option<fn(&str, &str) -> bool>) -> Equivalence {
        Equivalence {
            intervals: self
                .intervals
                .iter()
                .map(|interval| {
                    let ways = interval.ways.iter().map(|blocks| {
                        blocks
                            .iter()
                            .filter(|block| block.differs(compatible))
                            .cloned()
                            .collect::<Vec<_>>()
                    });
                    // The first way with the fewest blocks over a non-special
                    // edit; `min_by_key` keeps the first of equal ones.
                    let blocks = ways
                        .min_by_key(|blocks| {
                            blocks
                                .iter()
                                .filter(|block| {
                                    interval
                                        .edits
                                        .iter()
                                        .any(|&(first, last)| block.overlaps(first, last))
                                })
                                .count()
                        })
                        .unwrap_or_default();
                    Decided {
                        start: interval.start,
                        end: interval.end,
                        blocks,
                    }
                })
                .collect(),
        }
    }
}

/// The differing blocks of the built intervals under one symbol test.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub(crate) struct Equivalence {
    intervals: Vec<Decided>,
}

/// One built interval and the differing blocks of the way that decides it.
#[derive(Debug, Clone, PartialEq, Eq)]
struct Decided {
    start: i64,
    end: i64,
    blocks: Vec<Block>,
}

impl Equivalence {
    /// Whether an edit over `first..=last` lies in a built interval whose
    /// deciding way finds no differing block over it: its side writes the
    /// same sequence as the other there.
    pub(crate) fn explains(&self, first: i64, last: i64) -> bool {
        self.intervals.iter().any(|interval| {
            interval.start <= first
                && last <= interval.end
                && !interval
                    .blocks
                    .iter()
                    .any(|block| block.overlaps(first, last))
        })
    }

    /// The differing blocks of each interval's deciding way, in reference
    /// order.
    pub(crate) fn blocks(&self) -> impl Iterator<Item = &Block> {
        self.intervals.iter().flat_map(|interval| &interval.blocks)
    }
}

/// What one side writes at one anchor base.
#[derive(Debug, Default)]
struct Column<'a> {
    /// The replacing symbols, empty for a deletion, and whether the side
    /// wrote them itself; `None` keeps the reference.
    base: Option<(&'a str, bool)>,
    /// Inserted symbols by ordinal, and whether the side wrote them itself.
    inserted: BTreeMap<u32, (&'a str, bool)>,
}

/// One side's own and borrowed edits by anchor base, or `None` when they
/// cannot be applied inside `start..=end`.
fn columns<'a>(
    (start, end): (i64, i64),
    own: &[Edit<'a>],
    borrowed: &[Edit<'a>],
) -> Option<BTreeMap<i64, Column<'a>>> {
    let mut columns = BTreeMap::<i64, Column<'a>>::new();
    let edits = own
        .iter()
        .map(|edit| (edit, edit.lenient))
        .chain(borrowed.iter().map(|edit| (edit, false)));
    for (edit, written) in edits {
        let symbols = edit.symbols?;
        if edit.last() > end || edit.base < start {
            return None;
        }
        let taken = match edit.ordinal {
            Some(ordinal) => columns
                .entry(edit.base)
                .or_default()
                .inserted
                .insert(ordinal, (symbols, written))
                .is_some(),
            None if edit.span < 1 => return None,
            // A substitution of several bases with as many symbols replaces
            // them one by one; anything else writes its symbols at the first
            // base and deletes the rest.
            None => (0..edit.span).any(|offset| {
                let replacing = if usize::try_from(edit.span).ok() == Some(symbols.len()) {
                    let index = usize::try_from(offset).unwrap_or(usize::MAX);
                    symbols.get(index..=index).unwrap_or("")
                } else if offset == 0 {
                    symbols
                } else {
                    ""
                };
                columns
                    .entry(edit.base + offset)
                    .or_default()
                    .base
                    .replace((replacing, written))
                    .is_some()
            }),
        };
        if taken {
            return None;
        }
    }
    Some(columns)
}

/// The pieces of one interval whose symbols differ, or `None` when a side
/// cannot be applied or the interval leaves the reference. Each side is its
/// own edits and the edits it borrows from the other side.
fn pieces(
    reference: &str,
    (start, end): (i64, i64),
    [(dna, dna_borrowed), (truth, truth_borrowed)]: [(&[Edit<'_>], &[Edit<'_>]); 2],
) -> Option<Vec<Block>> {
    let dna = columns((start, end), dna, dna_borrowed)?;
    let truth = columns((start, end), truth, truth_borrowed)?;
    let mut found = Vec::new();
    let mut first = start;
    let (mut local_dna, mut local_truth) = (Local::default(), Local::default());
    for base in start..=end {
        let index = usize::try_from(base.checked_sub(1)?).ok()?;
        let reference_base = reference.get(index..=index)?;
        for (columns, local) in [(&dna, &mut local_dna), (&truth, &mut local_truth)] {
            let column = columns.get(&base);
            let (symbols, written) = column
                .and_then(|column| column.base)
                .unwrap_or((reference_base, false));
            local.push(symbols, written);
            for (symbols, written) in column.iter().flat_map(|column| column.inserted.values()) {
                local.push(symbols, *written);
            }
        }
        if local_dna.text.len() == local_truth.text.len() {
            if local_dna.text != local_truth.text {
                found.push(Block::of(first, base, &local_dna, &local_truth, false));
            }
            local_dna = Local::default();
            local_truth = Local::default();
            first = base + 1;
        }
    }
    if first <= end {
        found.push(Block::of(first, end, &local_dna, &local_truth, true));
    }
    Some(found)
}

/// One side's local sequence and, per symbol, whether the side wrote it.
#[derive(Debug, Default)]
struct Local {
    text: String,
    written: Vec<bool>,
}

impl Local {
    fn push(&mut self, symbols: &str, written: bool) {
        self.text.push_str(symbols);
        self.written
            .extend(std::iter::repeat_n(written, symbols.len()));
    }
}

impl Block {
    fn of(start: i64, end: i64, dna: &Local, truth: &Local, open: bool) -> Self {
        Self {
            start,
            end,
            dna: dna.text.clone(),
            truth: truth.text.clone(),
            truth_written: truth.written.clone(),
            open,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `ACCCCTGG` at 1..=8: a C run at 2..=5.
    const REFERENCE: &str = "ACCCCTGG";

    /// Edits from tokens: `7A` (substitution), `5-` (deletion), `2.1C`
    /// (insertion), `2CC-` (a deletion of two bases), `*` before a token
    /// for a special edit.
    fn edits(tokens: &[&'static str]) -> Vec<Edit<'static>> {
        tokens
            .iter()
            .map(|token| {
                let (special, token) = token
                    .strip_prefix('*')
                    .map_or((false, *token), |token| (true, token));
                let (position, symbol) = token.split_at(token.len() - 1);
                let digits = position
                    .find(|symbol: char| symbol.is_ascii_alphabetic())
                    .unwrap_or(position.len());
                let (position, span) = position.split_at(digits);
                let (base, ordinal) = position
                    .split_once('.')
                    .map_or((position, None), |(base, ordinal)| {
                        (base, ordinal.parse().ok())
                    });
                Edit {
                    base: base.parse().unwrap_or(0),
                    ordinal,
                    span: i64::try_from(span.len().max(1)).unwrap_or(1),
                    symbols: Some(if symbol == "-" { "" } else { symbol }),
                    special,
                    lenient: true,
                }
            })
            .collect()
    }

    /// IUPAC codes for C, T, Y, A, G, R only.
    fn compatible(left: &str, right: &str) -> bool {
        let bases = |code: &str| {
            match code {
                "Y" => "CT",
                "R" => "AG",
                other => other,
            }
            .to_owned()
        };
        let right = bases(right);
        bases(left).chars().any(|base| right.contains(base))
    }

    fn haplotypes(dna: &[&'static str], truth: &[&'static str]) -> Haplotypes {
        Haplotypes::new(REFERENCE, &[(1, 8)], &edits(dna), &edits(truth))
    }

    /// The block texts of the first way, exactly or IUPAC-compatibly.
    fn texts(dna: &[&'static str], truth: &[&'static str], lenient: bool) -> Vec<String> {
        haplotypes(dna, truth)
            .equivalence(lenient.then_some(compatible as fn(&str, &str) -> bool))
            .blocks()
            .map(Block::text)
            .collect()
    }

    #[test]
    fn an_insertion_anywhere_in_a_run_is_the_same_sequence() {
        assert!(texts(&["2.1C"], &["5.1C"], false).is_empty());
        // A substitution that lengthens the run instead does not: the
        // lengths differ from base 2 on and never meet again.
        assert_eq!(texts(&["2.1C"], &["6C"], false), ["2-8+:CCCCCTGG/CCCCCGG"]);
    }

    #[test]
    fn a_truth_code_admits_the_base_of_dna_but_not_the_reverse() {
        assert_eq!(texts(&["7A"], &["7R"], false), ["7-7:A/R"]);
        assert!(texts(&["7A"], &["7R"], true).is_empty());
        // The reference base where DNA writes nothing, against a truth code
        // that admits it, or does not.
        assert!(texts(&[], &["7R"], true).is_empty());
        assert_eq!(texts(&[], &["7Y"], true), ["7-7:G/Y"]);
        // A DNA code against the truth's reference base is a difference: the
        // truth wrote nothing there.
        assert_eq!(texts(&["7R"], &[], true), ["7-7:R/G"]);
        // Both written: the codes share a base.
        assert!(texts(&["7R"], &["7A"], true).is_empty());
    }

    #[test]
    fn a_deletion_and_an_insertion_balance_into_one_block() {
        // Deleting the first C and inserting one after T: CCCTC against
        // CCCCT, a real difference between equal lengths.
        assert_eq!(texts(&["2-", "6.1C"], &[], false), ["2-6:CCCTC/CCCCT"]);
        // Deleting any C of the run is the same.
        assert!(texts(&["2-"], &["5-"], false).is_empty());
        assert_eq!(texts(&["5-"], &[], false), ["5-8+:TGG/CTGG"]);
    }

    #[test]
    fn edits_of_several_bases_cover_their_span() {
        // A two-base deletion is two one-base deletions, not one.
        assert!(texts(&["2CC-"], &["2-", "3-"], false).is_empty());
        assert_eq!(texts(&["2CC-"], &["2-"], false), ["3-8+:CCTGG/CCCTGG"]);
        // Two substituted bases are two substitutions.
        let mnv = Edit {
            base: 7,
            ordinal: None,
            span: 2,
            symbols: Some("AT"),
            special: false,
            lenient: true,
        };
        let equivalence =
            Haplotypes::new(REFERENCE, &[(1, 8)], &[mnv], &edits(&["7A", "8T"])).equivalence(None);
        assert_eq!(equivalence.blocks().count(), 0);
        assert!(equivalence.explains(7, 8));
        // An edit leaving the interval is not applied.
        let past = Edit { base: 8, ..mnv };
        assert_eq!(
            Haplotypes::new(REFERENCE, &[(1, 8)], &[past], &[]),
            Haplotypes::default()
        );
    }

    #[test]
    fn insertions_follow_their_ordinal_order() {
        assert!(texts(&["6.2A", "6.1G"], &["6.1G", "6.2A"], false).is_empty());
        assert_eq!(
            texts(&["6.2A", "6.1G"], &["6.1A", "6.2G"], false),
            ["6-6:TGA/TAG"]
        );
    }

    #[test]
    fn special_edits_cannot_decide_a_difference() {
        // A special length difference in the run 2-5, and the next C run
        // placement written differently: each side with its own special
        // edits differs to the end, with the same ones it does not, so the
        // shared way decides and leaves no block.
        let dna = ["*5.1C", "6.1G"];
        let truth = ["*5.1C", "*5.2C", "7.1G"];
        let equivalence = haplotypes(&dna, &truth).equivalence(None);
        assert_eq!(equivalence.blocks().count(), 0);
        assert!(equivalence.explains(6, 6) && equivalence.explains(7, 7));
        // A special insertion and a deletion that write the other side's
        // substitution (as `455.1T 459DEL` writes `456T`): equal only with
        // each side's own special edits.
        let equivalence = haplotypes(&["*1.1C", "6-"], &["6C"]).equivalence(None);
        assert!(equivalence.explains(1, 1) && equivalence.explains(6, 6));
        assert_eq!(equivalence.blocks().count(), 0);
        // A real difference stays one under every way.
        let equivalence = haplotypes(&["*5.1C", "7A"], &["7T"]).equivalence(None);
        assert!(!equivalence.explains(7, 7));
    }

    #[test]
    fn one_way_decides_a_whole_interval() {
        let block = |start, end| Block {
            start,
            end,
            dna: "A".into(),
            truth: "C".into(),
            truth_written: vec![true],
            open: false,
        };
        let interval = |ways| Haplotypes {
            intervals: vec![Interval {
                start: 1,
                end: 20,
                edits: vec![(3, 3), (15, 15)],
                ways,
            }],
        };
        // The way with DNA's special edits leaves the edit at 3, the way with
        // the truth's the one at 15, and no way both: one way decides, so
        // the two are not explained together.
        let mixed = interval(vec![
            vec![block(1, 20)],
            vec![block(15, 15)],
            vec![block(3, 3)],
        ])
        .equivalence(None);
        assert!(!(mixed.explains(3, 3) && mixed.explains(15, 15)));
        assert_eq!(
            mixed.blocks().map(Block::text).collect::<Vec<_>>(),
            ["1-20:A/C"]
        );
        // A way without a block over a non-special edit decides, even
        // after a way with fewer blocks elsewhere.
        let clear =
            interval(vec![vec![block(1, 20)], vec![block(8, 8), block(10, 10)]]).equivalence(None);
        assert!(clear.explains(3, 3) && clear.explains(15, 15));
        assert_eq!(clear.blocks().count(), 2);
        // Otherwise the first way with the fewest blocks over one.
        let fewest = interval(vec![
            vec![block(3, 3), block(15, 15)],
            vec![block(15, 15)],
            vec![block(3, 3)],
        ])
        .equivalence(None);
        assert!(fewest.explains(3, 3) && !fewest.explains(15, 15));
    }

    #[test]
    fn sides_that_cannot_be_applied_are_not_decided() {
        assert_eq!(haplotypes(&["3A", "3G"], &[]), Haplotypes::default());
        assert_eq!(haplotypes(&[], &["3.1A", "3.1G"]), Haplotypes::default());
        let unusable = Edit {
            base: 3,
            ordinal: None,
            span: 1,
            symbols: None,
            special: false,
            lenient: true,
        };
        assert_eq!(
            Haplotypes::new(REFERENCE, &[(1, 8)], &[unusable], &[]),
            Haplotypes::default()
        );
        // An interval past the reference end.
        assert_eq!(
            Haplotypes::new(REFERENCE, &[(7, 9)], &[], &[]),
            Haplotypes::default()
        );
        assert!(!Equivalence::default().explains(3, 3));
    }

    #[test]
    fn explains_edits_in_built_intervals_outside_blocks() {
        let equivalence = haplotypes(&["2.1C", "7A"], &["5.1C", "7C"]).equivalence(None);
        assert!(equivalence.explains(2, 2) && equivalence.explains(5, 5));
        assert!(!equivalence.explains(7, 7));
        assert!(!equivalence.explains(9, 9));
    }

    #[test]
    fn long_block_sequences_are_shortened() {
        let block = Block {
            start: 1,
            end: 30,
            dna: "A".repeat(30),
            truth: String::new(),
            truth_written: Vec::new(),
            open: true,
        };
        assert_eq!(block.text(), "1-30+:AAAAAAAAAA..AAAAAAAAAA(30)/-");
    }
}
