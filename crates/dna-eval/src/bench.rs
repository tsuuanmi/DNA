//! `dna-eval bench`: the artefact regression bench.
//!
//! A case is one adjudicated Sanger artefact or truth error, checked against
//! the raw signal, in one sample of a development panel, with the outcome
//! DNA must keep at the review level: the alleles it calls over a span, the
//! check items over it, and whether it calls the span at all. The bench runs
//! DNA on every case's sample and tests each expectation, so a signal change
//! that breaks known-good behaviour on an adjudicated case shows before a
//! panel run. It reads no truth: the expectation is the adjudicated outcome.
//!
//! The case list names samples, so it stays with the panel data, outside the
//! repository; only its format and this runner are committed.

use std::collections::BTreeSet;
use std::fmt::Write as _;

use crate::compare::{Allele, Failure, anchor};
use crate::runner::Called;

/// The column header of a case list.
pub(crate) const HEADER: [&str; 10] = [
    "case",
    "sample_id",
    "panel",
    "reads",
    "class",
    "first",
    "last",
    "expect",
    "status",
    "source",
];

/// Most positions a case's span holds: a case is one artefact, not a region.
const MAX_SPAN: i64 = 1_000;

/// The header of `bench.tsv`.
pub(crate) const RESULTS_HEADER: &str = "case\tsample_id\tpanel\treads\tclass\tfirst\tlast\texpect\tstatus\tverdict\talleles\tchecks\tcallable\n";

/// One adjudicated case.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct Case {
    /// Unique name of the case.
    pub(crate) name: String,
    pub(crate) sample: String,
    /// The panel the sample is read from, `<panel>-<version>` of a manifest
    /// header, for example `dev-v5`.
    pub(crate) panel: String,
    /// The reads that show the artefact; descriptive only, since DNA runs on
    /// every trace of the sample.
    pub(crate) reads: String,
    /// The artefact or truth-error class.
    pub(crate) class: String,
    /// The span the expectation reads, 1-based and inclusive.
    pub(crate) first: i64,
    pub(crate) last: i64,
    pub(crate) expect: Expect,
    /// A known failure: reported, not gated, until DNA meets it.
    pub(crate) known: bool,
}

/// What DNA must show over a case's span: any one of the alternatives,
/// each a set of terms that must all hold.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct Expect {
    pub(crate) text: String,
    alternatives: Vec<Vec<Term>>,
}

/// One condition over a case's span.
#[derive(Debug, Clone, PartialEq, Eq)]
enum Term {
    /// The alleles anchored in the span are exactly these (none when empty).
    Alleles(BTreeSet<String>),
    /// None of these alleles is called.
    Absent(BTreeSet<String>),
    /// No insertion or deletion is anchored in the span.
    NoIndel,
    /// A check item overlaps the span, carrying the label when one is given.
    Check(Option<String>),
    /// Every position of the span lies in a callable interval.
    Called,
    /// No position of the span lies in a callable interval: the call is
    /// withheld.
    Uncalled,
}

impl Expect {
    /// Parses `term;term|term`: `|` separates alternatives, `;` the terms
    /// that must all hold.
    pub(crate) fn parse(text: &str) -> Result<Self, String> {
        let alternatives = text
            .split('|')
            .map(|alternative| {
                alternative
                    .split(';')
                    .map(|term| Term::parse(term.trim()))
                    .collect::<Result<Vec<_>, _>>()
            })
            .collect::<Result<Vec<_>, _>>()?;
        Ok(Self {
            text: text.to_owned(),
            alternatives,
        })
    }

    /// Whether DNA's outcome meets the expectation over `first..=last`.
    pub(crate) fn holds(&self, observed: &Observed) -> bool {
        self.alternatives
            .iter()
            .any(|terms| terms.iter().all(|term| term.holds(observed)))
    }
}

impl Term {
    fn parse(text: &str) -> Result<Self, String> {
        let list = |value: &str| -> Result<BTreeSet<String>, String> {
            value
                .split(',')
                .map(str::trim)
                .filter(|allele| !allele.is_empty())
                .map(|allele| {
                    if allele
                        .chars()
                        .all(|c| c.is_ascii_alphanumeric() || c == '.')
                    {
                        Ok(allele.to_owned())
                    } else {
                        Err(format!("invalid allele {allele:?}"))
                    }
                })
                .collect()
        };
        match text.split_once('=') {
            Some(("alleles", value)) => Ok(Self::Alleles(list(value)?)),
            Some(("absent", value)) => {
                let alleles = list(value)?;
                if alleles.is_empty() {
                    return Err("absent= names no allele".to_owned());
                }
                Ok(Self::Absent(alleles))
            }
            Some(("check", label)) if !label.trim().is_empty() => {
                Ok(Self::Check(Some(label.trim().to_owned())))
            }
            None if text == "no_indel" => Ok(Self::NoIndel),
            None if text == "check" => Ok(Self::Check(None)),
            None if text == "called" => Ok(Self::Called),
            None if text == "uncalled" => Ok(Self::Uncalled),
            _ => Err(format!("unknown expectation term {text:?}")),
        }
    }

    fn holds(&self, observed: &Observed) -> bool {
        match self {
            Self::Alleles(alleles) => observed.alleles == *alleles,
            Self::Absent(alleles) => observed.all_alleles.is_disjoint(alleles),
            Self::NoIndel => !observed.indel,
            Self::Check(None) => !observed.checks.is_empty(),
            Self::Check(Some(label)) => observed
                .checks
                .iter()
                .any(|(_, _, labels)| labels.split(' ').any(|name| name == label)),
            Self::Called => observed.callable == observed.positions,
            Self::Uncalled => observed.callable == 0,
        }
    }
}

/// What DNA shows over a case's span.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub(crate) struct Observed {
    /// Allele texts anchored in the span (`573.1C`, `249DEL`).
    alleles: BTreeSet<String>,
    /// Every allele text of the sample.
    all_alleles: BTreeSet<String>,
    /// Whether an insertion or deletion is anchored in the span.
    indel: bool,
    /// Check items overlapping the span, with their labels.
    checks: Vec<(i64, i64, String)>,
    /// Span positions inside a callable interval.
    callable: i64,
    positions: i64,
}

impl Observed {
    /// DNA's outcome over `first..=last`.
    pub(crate) fn of(called: &Called, first: i64, last: i64) -> Self {
        let within = |allele: &&Allele| {
            anchor(&allele.position).is_some_and(|anchor| (first..=last).contains(&anchor))
        };
        let in_span = called.alleles.iter().filter(within).collect::<Vec<_>>();
        Self {
            alleles: in_span.iter().map(|allele| allele.text()).collect(),
            all_alleles: called.alleles.iter().map(Allele::text).collect(),
            indel: in_span
                .iter()
                .any(|allele| allele.reference == "-" || allele.alternate == "-"),
            checks: called
                .check_labels
                .iter()
                .filter(|(start, end, _)| *start <= last && *end >= first)
                .cloned()
                .collect(),
            callable: (first..=last)
                .filter(|position| {
                    called
                        .intervals
                        .iter()
                        .any(|&(start, end)| (start..=end).contains(position))
                })
                .count()
                .try_into()
                .unwrap_or(i64::MAX),
            positions: last - first + 1,
        }
    }

    fn alleles_text(&self) -> String {
        if self.alleles.is_empty() {
            "-".to_owned()
        } else {
            self.alleles.iter().cloned().collect::<Vec<_>>().join(",")
        }
    }

    fn checks_text(&self) -> String {
        if self.checks.is_empty() {
            "-".to_owned()
        } else {
            self.checks
                .iter()
                .map(|(first, last, labels)| format!("{first}-{last}:{}", labels.replace(' ', "+")))
                .collect::<Vec<_>>()
                .join(",")
        }
    }
}

/// A case's verdict.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) enum Verdict {
    /// A guarded case fails: a regression.
    Fail,
    /// DNA failed on the sample, so the case cannot be judged.
    Error,
    /// A known failure still fails.
    Known,
    /// A known failure now passes: promote it to a guard.
    Fixed,
    Pass,
}

impl Verdict {
    const fn label(self) -> &'static str {
        match self {
            Self::Fail => "FAIL",
            Self::Error => "ERROR",
            Self::Known => "KNOWN",
            Self::Fixed => "FIXED",
            Self::Pass => "PASS",
        }
    }

    /// Whether the verdict fails the bench.
    pub(crate) const fn gates(self) -> bool {
        matches!(self, Self::Fail | Self::Error)
    }
}

/// Judges a case against DNA's outcome for its sample.
pub(crate) fn judge(case: &Case, result: &Result<Called, Failure>) -> (Verdict, Option<Observed>) {
    let Ok(called) = result else {
        return (Verdict::Error, None);
    };
    let observed = Observed::of(called, case.first, case.last);
    let verdict = match (case.expect.holds(&observed), case.known) {
        (true, false) => Verdict::Pass,
        (false, false) => Verdict::Fail,
        (false, true) => Verdict::Known,
        (true, true) => Verdict::Fixed,
    };
    (verdict, Some(observed))
}

/// Parses a case list: a TSV with [`HEADER`] and no comment lines.
pub(crate) fn parse(text: &str) -> Result<Vec<Case>, String> {
    let mut lines = text.lines().enumerate();
    let header = lines
        .next()
        .map(|(_, line)| line.split('\t').collect::<Vec<_>>())
        .ok_or("the case list is empty")?;
    if header != HEADER {
        return Err(format!(
            "the column header must be {:?}, found {header:?}",
            HEADER.join("\t")
        ));
    }
    let mut cases = Vec::new();
    let mut names = BTreeSet::new();
    for (index, line) in lines {
        let number = index + 1;
        if line.trim().is_empty() {
            continue;
        }
        let fields = line.split('\t').collect::<Vec<_>>();
        if fields.len() != HEADER.len() {
            return Err(format!(
                "line {number} has {} fields, expected {}",
                fields.len(),
                HEADER.len()
            ));
        }
        let position = |value: &str, column: &str| {
            value
                .parse::<i64>()
                .ok()
                .filter(|position| *position >= 1)
                .ok_or_else(|| format!("line {number} has an invalid {column} {value:?}"))
        };
        let (first, last) = (position(fields[5], "first")?, position(fields[6], "last")?);
        if last < first || last - first >= MAX_SPAN {
            return Err(format!(
                "line {number} has span {first}-{last}: last must not precede first, \
                 and a span holds at most {MAX_SPAN} positions"
            ));
        }
        let known = match fields[8] {
            "guard" => false,
            "known" => true,
            other => {
                return Err(format!(
                    "line {number} has status {other:?}, expected guard or known"
                ));
            }
        };
        for (column, value) in [
            ("case", fields[0]),
            ("sample_id", fields[1]),
            ("panel", fields[2]),
            ("class", fields[4]),
        ] {
            if value.is_empty() {
                return Err(format!("line {number} has an empty {column}"));
            }
        }
        if !names.insert(fields[0]) {
            return Err(format!("line {number} repeats case {:?}", fields[0]));
        }
        let expect =
            Expect::parse(fields[7]).map_err(|message| format!("line {number}: {message}"))?;
        cases.push(Case {
            name: fields[0].to_owned(),
            sample: fields[1].to_owned(),
            panel: fields[2].to_owned(),
            reads: fields[3].to_owned(),
            class: fields[4].to_owned(),
            first,
            last,
            expect,
            known,
        });
    }
    if cases.is_empty() {
        return Err("the case list holds no case".to_owned());
    }
    Ok(cases)
}

/// The bench's report and `bench.tsv` from every case's verdict, and
/// whether it fails.
pub(crate) struct Report {
    pub(crate) text: String,
    pub(crate) table: String,
    pub(crate) failed: bool,
}

/// Builds the report: one line per case that is not a plain pass, the
/// counts per class, and the verdict line.
pub(crate) fn report(
    revision: &str,
    cases: &[Case],
    judged: &[(Verdict, Option<Observed>, Option<String>)],
) -> Report {
    let mut text = format!(
        "# Artefact bench\n\nHarness {revision}; {} cases.\n\n",
        cases.len()
    );
    let mut table = String::from(RESULTS_HEADER);
    let mut classes = std::collections::BTreeMap::<&str, [usize; 5]>::new();
    let mut lines = String::new();
    for (case, (verdict, observed, failure)) in cases.iter().zip(judged) {
        let counts = classes.entry(case.class.as_str()).or_default();
        counts[*verdict as usize] += 1;
        let (alleles, checks, callable) = observed.as_ref().map_or_else(
            || ("-".to_owned(), "-".to_owned(), "-".to_owned()),
            |observed| {
                (
                    observed.alleles_text(),
                    observed.checks_text(),
                    format!("{}/{}", observed.callable, observed.positions),
                )
            },
        );
        let _ = writeln!(
            table,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{alleles}\t{checks}\t{callable}",
            case.name,
            case.sample,
            case.panel,
            case.reads,
            case.class,
            case.first,
            case.last,
            case.expect.text,
            if case.known { "known" } else { "guard" },
            verdict.label(),
        );
        if *verdict != Verdict::Pass {
            let seen = failure.clone().unwrap_or_else(|| {
                format!("alleles {alleles}; checks {checks}; callable {callable}")
            });
            let _ = writeln!(
                lines,
                "| {} | {} | {} {} | {} | {}-{} | `{}` | {seen} |",
                verdict.label(),
                case.name,
                case.sample,
                case.panel,
                case.class,
                case.first,
                case.last,
                case.expect.text,
            );
        }
    }
    if !lines.is_empty() {
        text.push_str("| Verdict | Case | Sample | Class | Span | Expected | Observed |\n");
        text.push_str("|---|---|---|---|---|---|---|\n");
        text.push_str(&lines);
        text.push('\n');
    }
    text.push_str(
        "| Class | Pass | Fail | Error | Known | Fixed |\n|---|---:|---:|---:|---:|---:|\n",
    );
    for (class, counts) in &classes {
        let _ = writeln!(
            text,
            "| {class} | {} | {} | {} | {} | {} |",
            counts[Verdict::Pass as usize],
            counts[Verdict::Fail as usize],
            counts[Verdict::Error as usize],
            counts[Verdict::Known as usize],
            counts[Verdict::Fixed as usize],
        );
    }
    let failing = judged
        .iter()
        .filter(|(verdict, ..)| verdict.gates())
        .count();
    let _ = writeln!(
        text,
        "\n{}",
        if failing == 0 {
            "BENCH PASS".to_owned()
        } else {
            format!("BENCH FAIL: {failing} case(s) regressed or could not be judged")
        }
    );
    Report {
        text,
        table,
        failed: failing > 0,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn called(
        alleles: &[(&str, &str, &str)],
        checks: &[(i64, i64, &str)],
        intervals: &[(i64, i64)],
    ) -> Called {
        Called {
            rejected: 0,
            alleles: alleles
                .iter()
                .map(|(position, reference, alternate)| Allele::new(position, reference, alternate))
                .collect(),
            not_recorded: BTreeSet::new(),
            intervals: intervals.to_vec(),
            signals: crate::lod::Signals::default(),
            loci: crate::cause::loci::Loci::default(),
            identity: crate::runner::Identity {
                configuration_sha256: String::new(),
                reference_name: String::new(),
                reference_sha256: String::new(),
                profile_id: String::new(),
                profile_sha256: String::new(),
            },
            empop: None,
            notes: Vec::new(),
            checks: checks
                .iter()
                .map(|&(first, last, _)| (first, last))
                .collect(),
            check_labels: checks
                .iter()
                .map(|&(first, last, labels)| (first, last, labels.to_owned()))
                .collect(),
            components: Vec::new(),
            tracts: Vec::new(),
        }
    }

    fn case(first: i64, last: i64, expect: &str, known: bool) -> Case {
        Case {
            name: "c".into(),
            sample: "S1".into(),
            panel: "dev-v1".into(),
            reads: "HV3R".into(),
            class: "extra_peak".into(),
            first,
            last,
            expect: Expect::parse(expect).unwrap_or_else(|message| panic!("{message}")),
            known,
        }
    }

    fn verdict(case: &Case, called: &Called) -> Verdict {
        judge(case, &Ok(called.clone())).0
    }

    #[test]
    fn alleles_must_match_exactly_over_the_span() {
        let dna = called(
            &[
                ("573.1", "-", "C"),
                ("573.2", "-", "C"),
                ("16519", "T", "C"),
            ],
            &[],
            &[(438, 576)],
        );
        assert_eq!(
            verdict(&case(568, 574, "alleles=573.1C,573.2C", false), &dna),
            Verdict::Pass
        );
        // One C fewer or more is a regression.
        assert_eq!(
            verdict(&case(568, 574, "alleles=573.1C", false), &dna),
            Verdict::Fail
        );
        // Alleles outside the span do not count.
        assert_eq!(
            verdict(&case(100, 200, "alleles=", false), &dna),
            Verdict::Pass
        );
    }

    #[test]
    fn no_indel_and_absent_read_insertions_and_deletions() {
        let dna = called(&[("560.1", "-", "C"), ("249", "A", "-")], &[], &[(73, 576)]);
        assert_eq!(
            verdict(&case(556, 561, "no_indel", false), &dna),
            Verdict::Fail
        );
        assert_eq!(
            verdict(&case(249, 249, "no_indel", false), &dna),
            Verdict::Fail
        );
        assert_eq!(
            verdict(&case(250, 300, "no_indel", false), &dna),
            Verdict::Pass
        );
        assert_eq!(
            verdict(&case(1, 1, "absent=249DEL", false), &dna),
            Verdict::Fail
        );
        assert_eq!(
            verdict(&case(1, 1, "absent=560.2C", false), &dna),
            Verdict::Pass
        );
    }

    #[test]
    fn checks_match_by_overlap_and_label() {
        let dna = called(
            &[("527", "C", "A")],
            &[(527, 527, "one_read reference_second")],
            &[(438, 576)],
        );
        assert_eq!(
            verdict(&case(527, 527, "check=reference_second", false), &dna),
            Verdict::Pass
        );
        assert_eq!(
            verdict(&case(527, 527, "check=reference", false), &dna),
            Verdict::Fail
        );
        assert_eq!(
            verdict(&case(530, 530, "check", false), &dna),
            Verdict::Fail
        );
    }

    #[test]
    fn callability_and_alternatives() {
        let dna = called(
            &[("16166", "A", "-")],
            &[(16170, 16170, "read_edge")],
            &[(86, 185), (187, 340), (16024, 16365)],
        );
        assert_eq!(
            verdict(&case(186, 186, "uncalled", false), &dna),
            Verdict::Pass
        );
        assert_eq!(
            verdict(&case(185, 186, "called", false), &dna),
            Verdict::Fail
        );
        // A deletion is acceptable only under a check.
        let guarded = case(16166, 16179, "no_indel | check", false);
        assert_eq!(verdict(&guarded, &dna), Verdict::Pass);
        let unchecked = called(&[("16166", "A", "-")], &[], &[(16024, 16365)]);
        assert_eq!(verdict(&guarded, &unchecked), Verdict::Fail);
    }

    #[test]
    fn known_failures_report_without_gating() {
        let dna = called(&[("186", "C", "T")], &[], &[(73, 340)]);
        let known = case(186, 186, "uncalled", true);
        assert_eq!(verdict(&known, &dna), Verdict::Known);
        assert!(!Verdict::Known.gates());
        let fixed = called(&[], &[], &[(73, 185), (187, 340)]);
        assert_eq!(verdict(&known, &fixed), Verdict::Fixed);
        assert_eq!(judge(&known, &Err(Failure::Panic)).0, Verdict::Error);
        assert!(Verdict::Error.gates());
    }

    #[test]
    fn parses_a_case_list_and_rejects_malformed_rows() {
        let header = HEADER.join("\t");
        let row = "c1\tS1\tdev-v1\tHV3R\textra_peak\t556\t561\tno_indel\tguard\tnote";
        let cases =
            parse(&format!("{header}\n{row}\n")).unwrap_or_else(|message| panic!("{message}"));
        assert_eq!(cases.len(), 1);
        assert_eq!(
            (cases[0].first, cases[0].last, cases[0].known),
            (556, 561, false)
        );
        for (bad, message) in [
            (
                "c1\tS1\tdev-v1\tHV3R\textra_peak\t561\t556\tno_indel\tguard\t-",
                "span 561-556",
            ),
            (
                "c1\tS1\tdev-v1\tHV3R\textra_peak\t1\t1001\tno_indel\tguard\t-",
                "at most 1000",
            ),
            (
                "c1\tS1\tdev-v1\tHV3R\textra_peak\t556\t561\tindel\tguard\t-",
                "unknown expectation",
            ),
            (
                "c1\tS1\tdev-v1\tHV3R\textra_peak\t556\t561\tno_indel\tmaybe\t-",
                "status",
            ),
            (
                "c1\tS1\tdev-v1\tHV3R\textra_peak\t0\t561\tno_indel\tguard\t-",
                "invalid first",
            ),
            ("c1\tS1\tdev-v1\tHV3R", "fields"),
        ] {
            let error = parse(&format!("{header}\n{bad}\n"))
                .err()
                .unwrap_or_default();
            assert!(error.contains(message), "{bad}: {error}");
        }
        let repeated = parse(&format!("{header}\n{row}\n{row}\n"))
            .err()
            .unwrap_or_default();
        assert!(repeated.contains("repeats case"));
        assert!(parse("case\tsample_id\n").is_err());
        assert!(parse(&format!("{header}\n")).is_err());
    }

    #[test]
    fn the_report_lists_non_passing_cases_and_fails_on_regressions() {
        let cases = [case(1, 1, "no_indel", false), case(2, 2, "uncalled", true)];
        let judged = [
            (Verdict::Fail, Some(Observed::default()), None),
            (Verdict::Known, Some(Observed::default()), None),
        ];
        let report = report("abc", &cases, &judged);
        assert!(report.failed);
        assert!(report.text.contains("BENCH FAIL: 1 case(s)"));
        assert!(report.text.contains("| extra_peak | 0 | 1 | 0 | 1 | 0 |"));
        assert_eq!(report.table.lines().count(), 3);
        let passing = report_of_passes();
        assert!(!passing.failed && passing.text.contains("BENCH PASS"));
    }

    fn report_of_passes() -> Report {
        report(
            "abc",
            &[case(1, 1, "no_indel", false)],
            &[(Verdict::Pass, Some(Observed::default()), None)],
        )
    }
}
