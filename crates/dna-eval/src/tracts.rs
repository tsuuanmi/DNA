//! Every read's tract decodes (PROP-0012 increment 2): `tracts.tsv`, one row
//! per sample, read, and tract, with where the tract's first call lies on
//! the reference, its chosen structure, and its first-run and total length
//! marginals. Evidence only: nothing is compared with the truth here.

use std::fmt::Write as _;

use dna::sample_analysis::{LengthMarginal, Margin, SampleAnalysis};

/// One tract of one read.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Row {
    pub(crate) read: String,
    pub(crate) forward: bool,
    pub(crate) onset: usize,
    pub(crate) onset_reason: Option<&'static str>,
    pub(crate) to_end: bool,
    pub(crate) start: usize,
    pub(crate) reference: Option<usize>,
    pub(crate) run: char,
    pub(crate) mark: Option<char>,
    pub(crate) extent: usize,
    pub(crate) window: usize,
    pub(crate) residual: f64,
    pub(crate) first: Lengths,
    pub(crate) total: Lengths,
    /// The first-run marginal refitted with the read's own signatures and
    /// leakage; `None` when the refit did not fit.
    pub(crate) cleaned: Option<Lengths>,
}

/// What a row shows of one marginal.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Lengths {
    pub(crate) dominant: Option<usize>,
    pub(crate) share: f64,
    pub(crate) next: Option<usize>,
    pub(crate) z_adjusted: Option<f64>,
    pub(crate) margin: Margin,
}

// The decode gives finite shares and residuals, so equality is reflexive.
impl Eq for Row {}

impl Lengths {
    fn of(marginal: &LengthMarginal) -> Self {
        let share = marginal
            .dominant
            .and_then(|length| marginal.lengths.iter().position(|&kept| kept == length))
            .and_then(|index| marginal.shares.get(index))
            .copied()
            .unwrap_or(0.0);
        Self {
            dominant: marginal.dominant,
            share,
            next: marginal.next,
            z_adjusted: marginal.z_adjusted,
            margin: marginal.margin,
        }
    }
}

/// The tract rows of an analysed sample, by read and onset.
pub(crate) fn rows(analysis: &SampleAnalysis) -> Vec<Row> {
    let mut rows = Vec::new();
    for read in analysis.tracts() {
        for placed in read.tracts {
            let tract = &placed.tract;
            rows.push(Row {
                read: read.read.clone(),
                forward: read.forward,
                onset: tract.onset,
                onset_reason: placed.onset_reason,
                to_end: placed.to_end,
                start: tract.call_start,
                reference: placed.reference,
                run: tract.run,
                mark: tract.mark,
                extent: tract.extent,
                window: tract.window,
                residual: tract.residual,
                first: Lengths::of(&tract.first),
                total: Lengths::of(&tract.total),
                cleaned: tract
                    .cleaned
                    .as_ref()
                    .map(|cleaned| Lengths::of(&cleaned.first)),
            });
        }
    }
    rows
}

/// The table's header.
pub(crate) const HEADER: &str = "sample_id\tbatch_date\tread\tstrand\tonset\tonset_reason\tto_end\tstart\treference\trun\tmark\textent\twindow\tresidual\tfirst\tfirst_share\tfirst_next\tfirst_z\tfirst_margin\ttotal\ttotal_share\ttotal_next\ttotal_z\ttotal_margin\tcleaned_first\tcleaned_first_share\n";

/// Appends `sample`'s rows.
pub(crate) fn append(table: &mut String, sample: (&str, &str), rows: &[Row]) {
    let length = |value: Option<usize>| value.map_or_else(|| ".".to_owned(), |v| v.to_string());
    let margin = |value: Option<f64>| value.map_or_else(|| ".".to_owned(), |v| format!("{v:.2}"));
    let state = |margin: Margin| match margin {
        Margin::Clear => "clear",
        Margin::NearEqual => "near_equal",
        Margin::Undefined => "undefined",
    };
    for row in rows {
        let _ = writeln!(
            table,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{}\t{:.3}\t{}\t{}\t{}\t{}\t{:.3}\t{}\t{}\t{}\t{}\t{}",
            sample.0,
            sample.1,
            row.read,
            if row.forward { "forward" } else { "reverse" },
            row.onset,
            row.onset_reason.unwrap_or("."),
            if row.to_end { "yes" } else { "no" },
            row.start,
            length(row.reference),
            row.run,
            row.mark.map_or_else(|| ".".to_owned(), String::from),
            row.extent,
            row.window,
            row.residual,
            length(row.first.dominant),
            row.first.share,
            length(row.first.next),
            margin(row.first.z_adjusted),
            state(row.first.margin),
            length(row.total.dominant),
            row.total.share,
            length(row.total.next),
            margin(row.total.z_adjusted),
            state(row.total.margin),
            length(row.cleaned.and_then(|cleaned| cleaned.dominant)),
            row.cleaned
                .map_or_else(|| ".".to_owned(), |cleaned| format!("{:.3}", cleaned.share)),
        );
    }
}

#[cfg(test)]
mod tests {
    use super::{HEADER, Lengths, Margin, Row, append};

    fn lengths(dominant: usize, z: Option<f64>, margin: Margin) -> Lengths {
        Lengths {
            dominant: Some(dominant),
            share: 0.625,
            next: z.map(|_| dominant + 1),
            z_adjusted: z,
            margin,
        }
    }

    /// A tract with a mark and a refit, and one with neither and no
    /// reference position, under the contract's columns.
    #[test]
    fn writes_one_row_per_tract_under_the_contract() -> crate::contract::Checked {
        let marked = Row {
            read: "r".into(),
            forward: true,
            onset: 220,
            onset_reason: Some("post_homopolymer"),
            to_end: false,
            start: 200,
            reference: Some(303),
            run: 'C',
            mark: Some('T'),
            extent: 13,
            window: 76,
            residual: 0.0123,
            first: lengths(7, Some(51.6), Margin::Clear),
            total: lengths(13, Some(1.5), Margin::NearEqual),
            cleaned: Some(lengths(7, Some(60.3), Margin::Clear)),
        };
        let single = Row {
            forward: false,
            onset_reason: Some("weak_signal"),
            to_end: true,
            reference: None,
            mark: None,
            cleaned: None,
            first: lengths(6, None, Margin::Clear),
            total: Lengths {
                next: Some(7),
                ..lengths(6, None, Margin::Undefined)
            },
            ..marked.clone()
        };
        let mut table = String::from(HEADER);
        append(&mut table, ("S1", "2026-01-12"), &[marked, single]);
        assert_eq!(
            table.lines().skip(1).collect::<Vec<_>>(),
            [
                "S1\t2026-01-12\tr\tforward\t220\tpost_homopolymer\tno\t200\t303\tC\tT\t13\t76\t0.0123\t7\t0.625\t8\t51.60\tclear\t13\t0.625\t14\t1.50\tnear_equal\t7\t0.625",
                "S1\t2026-01-12\tr\treverse\t220\tweak_signal\tyes\t200\t.\tC\t.\t13\t76\t0.0123\t6\t0.625\t.\t.\tclear\t6\t0.625\t7\t.\tundefined\t.\t.",
            ]
        );
        crate::contract::check_table("tracts.tsv", &table)
    }
}
