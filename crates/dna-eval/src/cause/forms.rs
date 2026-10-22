//! How the reviewers write each sequence: per repeat window and window
//! sequence, the written forms and their samples over the whole truth
//! database, and outside the windows, the samples carrying each call.
//!
//! A representation discrepancy is the reviewers' convention or a minority
//! form in the truth depending on which side writes the majority form of
//! the same window sequence. A run writes the index as `forms.tsv`, keyed by
//! the truth database's SHA-256, and a later run on the same database can
//! read it instead of the whole database.

use std::fmt::Write as _;

use std::collections::{BTreeMap, BTreeSet};

use crate::compare::{self, Allele};
use crate::geometry::{Geometry, Window};
use crate::truth::TruthSample;

/// Written forms of the truth.
#[derive(Debug, Default)]
pub(crate) struct Forms {
    /// Per window and window sequence, the samples per written form.
    windows: BTreeMap<(Window, String), BTreeMap<String, usize>>,
    /// The samples per set of calls outside the windows.
    outside: BTreeMap<BTreeSet<String>, usize>,
    /// The geometry's windows with their bounds.
    spans: Vec<(Window, i64, i64)>,
}

/// Which side writes the majority form.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Majority {
    /// DNA's form is the reviewers' majority form.
    Dna,
    /// The truth's form is.
    Truth,
    /// Neither, or a tie.
    Neither,
}

impl Forms {
    /// Indexes every sample of `samples` with IUPAC-free calls in a window
    /// of `geometry`, on `reference`.
    pub(crate) fn new<'a>(
        samples: impl IntoIterator<Item = &'a TruthSample>,
        reference: &str,
        geometry: &Geometry,
    ) -> Self {
        let mut forms = Self {
            spans: geometry.window_bounds().collect(),
            ..Self::default()
        };
        for sample in samples {
            for (window, first, last) in &forms.spans {
                if let Some(sequence) = window_sequence(&sample.alleles, *first, *last, reference) {
                    *forms
                        .windows
                        .entry((window.clone(), sequence))
                        .or_default()
                        .entry(form(&sample.alleles, *first, *last))
                        .or_default() += 1;
                }
            }
            *forms
                .outside
                .entry(
                    sample
                        .alleles
                        .iter()
                        .filter(|allele| allele.reference != allele.alternate)
                        .filter(|allele| forms.outside_windows(allele))
                        .map(Allele::text)
                        .collect(),
                )
                .or_default() += 1;
        }
        forms
    }

    /// Which side writes the majority form in `window`, for DNA's and the
    /// truth's alleles of the same sequence; `None` when the sides agree
    /// there or the sequence holds a call DNA cannot write.
    pub(crate) fn window(
        &self,
        window: &Window,
        dna: &[Allele],
        truth: &[Allele],
        reference: &str,
    ) -> Option<Majority> {
        let &(_, first, last) = self.spans.iter().find(|(span, _, _)| span == window)?;
        let (dna_form, truth_form) = (form(dna, first, last), form(truth, first, last));
        if dna_form == truth_form {
            return None;
        }
        let sequence = window_sequence(truth, first, last, reference)?;
        let counts = self.windows.get(&(window.clone(), sequence));
        let count = |form: &str| {
            counts
                .and_then(|counts| counts.get(form))
                .copied()
                .unwrap_or(0)
        };
        Some(majority(count(&dna_form), count(&truth_form)))
    }

    /// Outside the windows, which side's calls more truth samples carry
    /// whole: a proxy, since no window sequence groups them; `None` when
    /// neither side has a call there.
    pub(crate) fn outside(&self, dna: &[Allele], truth: &[Allele]) -> Option<Majority> {
        let calls = |alleles: &[Allele]| {
            alleles
                .iter()
                .filter(|allele| self.outside_windows(allele))
                .map(Allele::text)
                .collect::<BTreeSet<_>>()
        };
        let (dna, truth) = (calls(dna), calls(truth));
        if dna.is_empty() && truth.is_empty() {
            return None;
        }
        let carrying = |calls: &BTreeSet<String>| {
            self.outside
                .iter()
                .filter(|(set, _)| !calls.is_empty() && calls.is_subset(set))
                .map(|(_, samples)| samples)
                .sum::<usize>()
        };
        Some(majority(carrying(&dna), carrying(&truth)))
    }
}

const HEADER: &str = "kind\twindow\tsequence\tform\tsamples\n";

impl Forms {
    /// `forms.tsv`: the index, keyed by the truth database's SHA-256. An
    /// empty form or call set is written `-`.
    pub(crate) fn table(&self, truth_sha256: &str) -> String {
        let mut text = format!("# truth_database_sha256\t{truth_sha256}\n{HEADER}");
        let shown = |text: &str| {
            if text.is_empty() {
                "-".to_owned()
            } else {
                text.to_owned()
            }
        };
        for ((window, sequence), forms) in &self.windows {
            for (form, samples) in forms {
                let _ = writeln!(
                    text,
                    "window\t{}\t{sequence}\t{}\t{samples}",
                    window.label(),
                    shown(form)
                );
            }
        }
        for (calls, samples) in &self.outside {
            let calls = calls
                .iter()
                .map(String::as_str)
                .collect::<Vec<_>>()
                .join(" ");
            let _ = writeln!(text, "outside\t-\t-\t{}\t{samples}", shown(&calls));
        }
        text
    }

    /// Reads a `forms.tsv` written for the truth database whose SHA-256 is
    /// `truth_sha256`, over the windows of `geometry`.
    ///
    /// # Errors
    ///
    /// Returns a message when it was written for another database, names a
    /// window `geometry` lacks, or a line is malformed.
    pub(crate) fn parse(
        text: &str,
        truth_sha256: &str,
        geometry: &Geometry,
    ) -> Result<Self, String> {
        let mut lines = text.lines();
        let recorded = lines
            .next()
            .and_then(|line| line.strip_prefix("# truth_database_sha256\t"))
            .ok_or("the first line must record the truth database's SHA-256")?;
        if recorded != truth_sha256 {
            return Err(format!(
                "written for truth database {recorded}, not {truth_sha256}"
            ));
        }
        if lines.next() != Some(HEADER.trim_end()) {
            return Err("the second line must be the column header".to_owned());
        }
        let mut forms = Self {
            spans: geometry.window_bounds().collect(),
            ..Self::default()
        };
        for (index, line) in lines.enumerate() {
            let invalid = || format!("line {} is malformed", index + 3);
            let fields = line.split('\t').collect::<Vec<_>>();
            let [kind, window, sequence, form, samples] = fields[..] else {
                return Err(invalid());
            };
            let samples = samples.parse::<usize>().map_err(|_| invalid())?;
            let form = if form == "-" { "" } else { form };
            match kind {
                "window" => {
                    let window = geometry
                        .window_named(window)
                        .filter(|window| *window != Window::Other)
                        .ok_or_else(|| {
                            format!("line {} names an unknown window {window}", index + 3)
                        })?;
                    forms
                        .windows
                        .entry((window, sequence.to_owned()))
                        .or_default()
                        .insert(form.to_owned(), samples);
                }
                "outside" => {
                    forms.outside.insert(
                        form.split_whitespace().map(str::to_owned).collect(),
                        samples,
                    );
                }
                _ => return Err(invalid()),
            }
        }
        Ok(forms)
    }
}

fn majority(dna: usize, truth: usize) -> Majority {
    match dna.cmp(&truth) {
        std::cmp::Ordering::Greater => Majority::Dna,
        std::cmp::Ordering::Less => Majority::Truth,
        std::cmp::Ordering::Equal => Majority::Neither,
    }
}

impl Forms {
    /// Whether an allele's anchor lies outside every window.
    fn outside_windows(&self, allele: &Allele) -> bool {
        compare::anchor(&allele.position).is_some_and(|anchor| {
            !self
                .spans
                .iter()
                .any(|(_, first, last)| (*first..=*last).contains(&anchor))
        })
    }
}

/// The calls anchored in `first..=last`, sorted, as one text.
fn form(alleles: &[Allele], first: i64, last: i64) -> String {
    let mut calls = alleles
        .iter()
        .filter(|allele| allele.reference != allele.alternate)
        .filter(|allele| {
            compare::anchor(&allele.position).is_some_and(|anchor| (first..=last).contains(&anchor))
        })
        .map(Allele::text)
        .collect::<Vec<_>>();
    calls.sort();
    calls.dedup();
    calls.join(" ")
}

/// The sequence over `first..=last` with the calls anchored there applied,
/// insertions after their anchor; `None` when a call there is not a plain
/// A/C/G/T or deletion call.
fn window_sequence(alleles: &[Allele], first: i64, last: i64, reference: &str) -> Option<String> {
    let mut bases = BTreeMap::new();
    for allele in alleles
        .iter()
        .filter(|allele| allele.reference != allele.alternate)
    {
        let anchor = compare::anchor(&allele.position)?;
        if !(first..=last).contains(&anchor) {
            continue;
        }
        let ordinal = allele
            .position
            .split_once('.')
            .map_or(Some(0), |(_, ordinal)| ordinal.parse::<usize>().ok())?;
        let base = match allele.alternate.as_str() {
            "-" => None,
            base if base.len() == 1 && "ACGT".contains(base) => base.chars().next(),
            _ => return None,
        };
        bases.insert((anchor, ordinal), base);
    }
    let mut sequence = String::new();
    for position in first..=last {
        let index = usize::try_from(position - 1).ok()?;
        match bases.get(&(position, 0)) {
            Some(Some(base)) => sequence.push(*base),
            Some(None) => {}
            None => sequence.push(char::from(*reference.as_bytes().get(index)?)),
        }
        sequence.extend(
            bases
                .range((position, 1)..(position + 1, 0))
                .filter_map(|(_, base)| *base),
        );
    }
    Some(sequence)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample(calls: &[(&str, &str, &str)]) -> TruthSample {
        TruthSample {
            batch: "b".into(),
            status: None,
            alleles: calls
                .iter()
                .map(|(position, reference, alternate)| Allele::new(position, reference, alternate))
                .collect(),
            intervals: Vec::new(),
        }
    }

    #[test]
    fn finds_the_side_that_writes_the_majority_form() -> crate::contract::Checked {
        // A reference with an AC repeat in HVS-III's window, 512-524, as in
        // rCRS.
        let geometry = crate::geometry::tests::rcrs();
        let window = geometry.window_at(514);
        let reference = format!("{}GCACACACACACCG{}", "T".repeat(512), "T".repeat(40));
        let majority = sample(&[("523", "A", "-"), ("524", "C", "-")]);
        let minority = sample(&[("514", "C", "-"), ("515", "A", "-")]);
        let forms = Forms::new([&majority, &majority, &minority], &reference, geometry);
        assert_eq!(
            forms.window(&window, &majority.alleles, &minority.alleles, &reference),
            Some(Majority::Dna)
        );
        assert_eq!(
            forms.window(&window, &minority.alleles, &majority.alleles, &reference),
            Some(Majority::Truth)
        );
        assert_eq!(
            forms.window(&window, &majority.alleles, &majority.alleles, &reference),
            None
        );
        let read = Forms::parse(&forms.table("abc"), "abc", geometry);
        assert_eq!(
            read.as_ref().map(|read| read.table("abc")),
            Ok(forms.table("abc"))
        );
        assert!(Forms::parse(&forms.table("abc"), "other", geometry).is_err());
        let other = sample(&[("73", "A", "G")]);
        let forms = Forms::new([&other, &other, &majority], &reference, geometry);
        assert_eq!(
            forms.outside(&other.alleles, &majority.alleles),
            Some(Majority::Dna)
        );
        crate::contract::check_table("forms.tsv", &forms.table(&"e".repeat(64)))
    }
}
