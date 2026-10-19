//! Post-calling EMPOP export (PROP-0007): an offline query of interpretation
//! ranges and an rCRS-coded mitotype, built from a sample's called records
//! and their rendered calls.
//!
//! EMPOP reads the omission of sequence inside a declared range as the
//! reference, so a range covers only positions that exactly one record calls
//! callably. The export decides nothing about the sequence: the records'
//! callable segments already leave out unresolved positions and the
//! positions of excluded variants (SRS-VAR-019), and the calls are the
//! sample's represented calls as written. A call belongs to its position,
//! an insertion `X.k` to its anchor `X`, as EMPOP numbers it.

use std::ops::Range;

use dna_kernel::plugin::{Contract, PluginDescriptor, PluginFamily};
use dna_kernel::variant::{Variant, VariantKind};

/// One called record of the sample: a consensus segment or a reviewed
/// sequence.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Record {
    /// Reference intervals the record's alignment covers, 0-based half-open.
    pub reference: Vec<Range<usize>>,
    /// Callable reference intervals, 0-based half-open.
    pub callable: Vec<Range<usize>>,
    /// The record's variants excluded from reporting.
    pub excluded: Vec<Variant>,
}

/// One rendered per-base call of the sample.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Call<'a> {
    /// 1-based position of the changed base, or of the base an insertion
    /// follows.
    pub position: usize,
    /// `0` for a substitution or deletion, `1..` for successive inserted bases.
    pub ordinal: usize,
    /// The rendered call, for example `309.1C`.
    pub text: &'a str,
}

/// Why a position inside a record is not in a range.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum CutReason {
    /// No record calls it: an `N`, or a call the record does not trust.
    Unresolved,
    /// It is a position of a variant excluded from reporting.
    ExcludedVariant,
    /// More than one record calls it.
    Overlap,
}

impl CutReason {
    /// Stable machine-readable label recorded in results.
    #[must_use]
    pub fn label(self) -> &'static str {
        match self {
            Self::Unresolved => "unresolved",
            Self::ExcludedVariant => "excluded_variant",
            Self::Overlap => "overlap",
        }
    }
}

/// Positions inside a profile region and a record that no range covers,
/// 1-based inclusive.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Cut {
    /// First position.
    pub start: usize,
    /// Last position.
    pub end: usize,
    /// Why they are not in a range.
    pub reason: CutReason,
}

/// An EMPOP query.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Export {
    /// Interpretation ranges, 1-based inclusive, in the profile's region
    /// order and within a region in position order.
    pub ranges: Vec<[usize; 2]>,
    /// The calls inside the ranges, an insertion by its anchor, distinct, by
    /// position and then ordinal.
    pub mitotype: Vec<String>,
    /// The positions the ranges leave out, and why, in range order.
    pub cut: Vec<Cut>,
}

impl Export {
    /// The ranges as EMPOP takes them: `16024-16365 73-340`.
    #[must_use]
    pub fn ranges_text(&self) -> String {
        self.ranges
            .iter()
            .map(|[start, end]| format!("{start}-{end}"))
            .collect::<Vec<_>>()
            .join(" ")
    }

    /// The mitotype as EMPOP takes it: `73G 263G 309.1C`.
    #[must_use]
    pub fn mitotype_text(&self) -> String {
        self.mitotype.join(" ")
    }
}

/// What the records say about one position.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum State {
    /// No record covers it.
    Outside,
    /// Exactly one record calls it.
    Called,
    /// A record covers it, and none or more than one calls it.
    Cut(CutReason),
}

/// Builds the EMPOP query of a sample over the profile's `regions` (1-based
/// inclusive), from its called `records` and its rendered `calls`.
#[must_use]
pub fn export(regions: &[[usize; 2]], records: &[Record], calls: &[Call<'_>]) -> Export {
    let excluded = records
        .iter()
        .flat_map(|record| record.excluded.iter().flat_map(positions))
        .collect::<Vec<_>>();
    let state = |position: usize| {
        let index = position - 1;
        let covers =
            |segments: &[Range<usize>]| segments.iter().any(|segment| segment.contains(&index));
        match records
            .iter()
            .filter(|record| covers(&record.callable))
            .count()
        {
            1 => State::Called,
            0 if !records.iter().any(|record| covers(&record.reference)) => State::Outside,
            0 if excluded.contains(&index) => State::Cut(CutReason::ExcludedVariant),
            0 => State::Cut(CutReason::Unresolved),
            _ => State::Cut(CutReason::Overlap),
        }
    };
    let mut ranges: Vec<[usize; 2]> = Vec::new();
    let mut cut: Vec<Cut> = Vec::new();
    for &[start, end] in regions {
        let mut previous = State::Outside;
        for position in start.max(1)..=end {
            let current = state(position);
            let continues = current == previous && position > start;
            match current {
                State::Outside => {}
                State::Called if continues => {
                    if let Some(last) = ranges.last_mut() {
                        last[1] = position;
                    }
                }
                State::Called => match ranges.last_mut() {
                    // Ranges that touch across regions merge.
                    Some(last) if last[1] + 1 == position => last[1] = position,
                    _ => ranges.push([position, position]),
                },
                State::Cut(_) if continues => {
                    if let Some(last) = cut.last_mut() {
                        last.end = position;
                    }
                }
                State::Cut(reason) => cut.push(Cut {
                    start: position,
                    end: position,
                    reason,
                }),
            }
            previous = current;
        }
    }
    // EMPOP numbers an insertion `X.k` at `X`: a range holding `X` claims
    // every insertion after it. An excluded insertion cuts `X` itself
    // (SRS-VAR-019), so where `X` is in a range every insertion DNA observed
    // after it is eligible, and written.
    let inside = |position: usize| {
        ranges
            .iter()
            .any(|&[start, end]| (start..=end).contains(&position))
    };
    let mut kept = calls
        .iter()
        .filter(|call| inside(call.position))
        .map(|call| (call.position, call.ordinal, call.text.to_owned()))
        .collect::<Vec<_>>();
    kept.sort();
    kept.dedup();
    Export {
        ranges,
        mitotype: kept.into_iter().map(|(_, _, text)| text).collect(),
        cut,
    }
}

/// The 0-based reference positions of a variant: a substitution's position,
/// a deletion's deleted positions, and both positions beside an insertion's
/// junction (SRS-VAR-019).
fn positions(variant: &Variant) -> Vec<usize> {
    let anchor = variant.position_1based.saturating_sub(1);
    match variant.kind {
        VariantKind::Del => (anchor + 1..anchor + variant.reference.len()).collect(),
        VariantKind::Ins => vec![anchor, anchor + 1],
        _ => vec![anchor],
    }
}

/// The EMPOP query of the `notation` workflow, when the profile declares
/// `[empop]`.
pub const PLUGIN: PluginDescriptor = PluginDescriptor {
    id: "empop",
    family: PluginFamily::PostCalling,
    version: 1,
    provides: &[Contract::EmpopQuery],
    requires: &[Contract::Nomenclature],
    config_sections: &[],
};

#[cfg(test)]
mod tests {
    use super::*;

    /// HVS-I-like and HVS-II-like regions, in the profile's order.
    const REGIONS: [[usize; 2]; 2] = [[101, 120], [11, 30]];

    /// 0-based half-open segments from `(start, end)` pairs.
    fn segments(pairs: &[(usize, usize)]) -> Vec<Range<usize>> {
        pairs.iter().map(|&(start, end)| start..end).collect()
    }

    fn record(reference: (usize, usize), callable: &[(usize, usize)]) -> Record {
        Record {
            reference: segments(&[reference]),
            callable: segments(callable),
            excluded: Vec::new(),
        }
    }

    fn call(position: usize, ordinal: usize, text: &str) -> Call<'_> {
        Call {
            position,
            ordinal,
            text,
        }
    }

    fn reasons(export: &Export) -> Vec<(usize, usize, &'static str)> {
        export
            .cut
            .iter()
            .map(|cut| (cut.start, cut.end, cut.reason.label()))
            .collect()
    }

    // SRS-NOM-018
    #[test]
    fn an_unresolved_position_splits_a_range_and_regions_keep_their_order() {
        // Record one reads 5-40 (0-based 4..40) with 15-16 unresolved;
        // record two reads 101-120 whole.
        let records = [
            record((4, 40), &[(4, 14), (16, 40)]),
            record((100, 120), &[(100, 120)]),
        ];
        let exported = export(&REGIONS, &records, &[]);
        assert_eq!(exported.ranges, [[101, 120], [11, 14], [17, 30]]);
        assert_eq!(exported.ranges_text(), "101-120 11-14 17-30");
        assert_eq!(reasons(&exported), [(15, 16, "unresolved")]);
    }

    // SRS-NOM-018
    #[test]
    fn an_excluded_variant_and_an_excluded_insertion_cut_their_positions() {
        let excluded = |position: usize, reference: &str, alternate: &str, kind| Variant {
            contig: "toy".into(),
            position_1based: position,
            reference: reference.into(),
            alternate: alternate.into(),
            kind,
        };
        // An excluded substitution at 13 and insertion after 20: their
        // positions are already outside the callable segments.
        let records = [Record {
            reference: segments(&[(10, 30)]),
            callable: segments(&[(10, 12), (13, 19), (21, 30)]),
            excluded: vec![
                excluded(13, "A", "G", VariantKind::Snv),
                excluded(20, "A", "AC", VariantKind::Ins),
            ],
        }];
        let exported = export(&[[11, 30]], &records, &[]);
        assert_eq!(exported.ranges, [[11, 12], [14, 19], [22, 30]]);
        assert_eq!(
            reasons(&exported),
            [(13, 13, "excluded_variant"), (20, 21, "excluded_variant")]
        );
    }

    // SRS-NOM-018
    #[test]
    fn positions_two_records_call_are_cut() {
        let records = [record((10, 20), &[(10, 20)]), record((15, 30), &[(15, 30)])];
        let exported = export(&[[11, 30]], &records, &[]);
        assert_eq!(exported.ranges, [[11, 15], [21, 30]]);
        assert_eq!(reasons(&exported), [(16, 20, "overlap")]);
    }

    // SRS-NOM-018
    #[test]
    fn keeps_the_calls_inside_the_ranges_distinct_and_sorted() {
        let records = [record((10, 30), &[(10, 20), (22, 30)])];
        // Ranges 11-20 and 23-30; 21-22 unresolved.
        let calls = [
            call(25, 0, "25G"),
            call(12, 0, "12DEL"),
            // An insertion after a range's last position: EMPOP numbers it
            // at 20, which the range holds.
            call(20, 1, "20.1C"),
            call(21, 0, "21T"),
            // An insertion after a cut position, though 23 is in a range.
            call(22, 1, "22.1A"),
            call(15, 1, "15.1C"),
            call(15, 2, "15.2C"),
            call(12, 0, "12DEL"),
            // An IUPAC code passes through.
            call(28, 0, "28Y"),
        ];
        let exported = export(&[[11, 30]], &records, &calls);
        assert_eq!(exported.mitotype_text(), "12DEL 15.1C 15.2C 20.1C 25G 28Y");
    }

    // SRS-NOM-018
    #[test]
    fn a_region_without_records_exports_nothing() {
        let exported = export(&REGIONS, &[], &[call(15, 0, "15G")]);
        assert_eq!(
            exported,
            Export {
                ranges: Vec::new(),
                mitotype: Vec::new(),
                cut: Vec::new(),
            }
        );
    }
}
