//! The golden output of a regress check, `golden.tsv`: every line of DNA's
//! output on the regression set that a change can move, one
//! `sample_id kind key value` line each, in a canonical order, so that two
//! builds compare line by line and a file's SHA-256 names its content.
//!
//! The truth-free kinds are DNA's output alone; `block` and `status` read the
//! truth too. A sample DNA could not analyse has only its `failure` line.

use std::collections::{BTreeMap, BTreeSet};

use sha2::{Digest, Sha256};

use crate::compare::Measure;
use crate::manifest::Sample;
use crate::report::Report;
use crate::runner::Outcome;
use crate::workload::{self, Class, Window};
use crate::zones::Zones;

/// The header of `golden.tsv`.
pub(crate) const HEADER: &str = "sample_id\tkind\tkey\tvalue\n";

/// What a golden line records, in the order lines of a sample come.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub(crate) enum Kind {
    /// Why DNA could not analyse the sample; its only line.
    Failure,
    /// A call in the reviewer notation, recorded or not (SRS-NOM-017).
    Allele,
    /// A callable interval.
    Callable,
    /// A check item of the review list, with its labels.
    Check,
    /// A noted tract run of the review list.
    Note,
    /// Whether the sample is clear, its checked zones, and its positions to
    /// review.
    Review,
    /// The EMPOP query's ranges and mitotype.
    Empop,
    /// Traces set aside.
    Rejected,
    /// An unreviewed difference block, with its class and zone.
    Block,
    /// A status under each measure, the discrepancy class, and its cause.
    Status,
}

impl Kind {
    /// Every kind, in file order.
    pub(crate) const ALL: [Self; 10] = [
        Self::Failure,
        Self::Allele,
        Self::Callable,
        Self::Check,
        Self::Note,
        Self::Review,
        Self::Empop,
        Self::Rejected,
        Self::Block,
        Self::Status,
    ];

    /// The kind's name in `golden.tsv`.
    pub(crate) const fn label(self) -> &'static str {
        match self {
            Self::Failure => "failure",
            Self::Allele => "allele",
            Self::Callable => "callable",
            Self::Check => "check",
            Self::Note => "note",
            Self::Review => "review",
            Self::Empop => "empop",
            Self::Rejected => "rejected",
            Self::Block => "block",
            Self::Status => "status",
        }
    }

    /// The kind `label` names.
    pub(crate) fn parse(label: &str) -> Option<Self> {
        Self::ALL.into_iter().find(|kind| kind.label() == label)
    }

    /// Whether the kind is DNA's output alone, read without the truth.
    pub(crate) const fn truth_free(self) -> bool {
        !matches!(self, Self::Block | Self::Status)
    }
}

/// One golden line of the sample at `sample`, an index into the samples.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct Line {
    pub(crate) sample: usize,
    pub(crate) kind: Kind,
    pub(crate) key: String,
    pub(crate) value: String,
}

/// The canonical order of a key within its kind: by the number its key
/// starts with (a position), keys without one last, then by text.
pub(crate) fn key_order(key: &str) -> (bool, u64, &str) {
    let digits = key.bytes().take_while(u8::is_ascii_digit).count();
    let number = key[..digits].parse::<u64>().ok();
    (number.is_none(), number.unwrap_or(0), key)
}

impl Line {
    /// The line's place in the canonical order.
    fn order(&self) -> (usize, Kind, (bool, u64, &str), &str) {
        (self.sample, self.kind, key_order(&self.key), &self.value)
    }

    /// The positions the line's key or value names, 1-based and closed: an
    /// allele's anchor, a callable interval, a check item, a block, or a
    /// note's span. Other kinds name none.
    pub(crate) fn span(kind: Kind, key: &str, value: Option<&str>) -> Option<(i64, i64)> {
        let range = |text: &str| {
            let (first, last) = text.split_once('-')?;
            Some((first.parse().ok()?, last.parse().ok()?))
        };
        match kind {
            Kind::Allele => {
                let digits = key.bytes().take_while(u8::is_ascii_digit).count();
                let position = key[..digits].parse().ok()?;
                Some((position, position))
            }
            Kind::Callable | Kind::Check => range(key),
            Kind::Block => range(key.split_once(':')?.0),
            Kind::Note => range(value?.split(' ').next()?),
            _ => None,
        }
    }
}

/// A golden output: its lines over the samples, in canonical order.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct Golden {
    /// The sample identifiers, in manifest order.
    pub(crate) samples: Vec<String>,
    pub(crate) lines: Vec<Line>,
}

/// `text` as a golden cell: tabs and line ends become spaces, trailing
/// spaces go, and an empty cell is `-`.
pub(crate) fn cell(text: &str) -> String {
    let text = text.replace(['\t', '\n', '\r'], " ").trim_end().to_owned();
    if text.is_empty() {
        "-".to_owned()
    } else {
        text
    }
}

impl Golden {
    /// The golden output of `lines` over `samples`, put in canonical order.
    pub(crate) fn new(samples: Vec<String>, mut lines: Vec<Line>) -> Self {
        lines.sort_by(|left, right| left.order().cmp(&right.order()));
        Self { samples, lines }
    }

    /// One line's text, with its newline.
    pub(crate) fn line_text(&self, line: &Line) -> String {
        format!(
            "{}\t{}\t{}\t{}\n",
            self.samples[line.sample],
            line.kind.label(),
            line.key,
            line.value
        )
    }

    /// `golden.tsv`'s text.
    pub(crate) fn text(&self) -> String {
        let mut text = String::from(HEADER);
        for line in &self.lines {
            text.push_str(&self.line_text(line));
        }
        text
    }

    /// Reads `golden.tsv` over `samples`, the manifest's identifiers in
    /// order: every line must name one of them and a known kind, and the
    /// lines must be in canonical order.
    pub(crate) fn parse(text: &str, samples: &[String]) -> Result<Self, String> {
        let body = text
            .strip_prefix(HEADER)
            .ok_or("golden.tsv does not start with its header")?;
        let index = samples
            .iter()
            .enumerate()
            .map(|(index, sample)| (sample.as_str(), index))
            .collect::<BTreeMap<_, _>>();
        let mut lines = Vec::new();
        for (number, line) in body.lines().enumerate() {
            let at = number + 2;
            let cells = line.split('\t').collect::<Vec<_>>();
            let [sample, kind, key, value] = cells[..] else {
                return Err(format!("golden.tsv line {at} has {} cells", cells.len()));
            };
            let sample = *index
                .get(sample)
                .ok_or_else(|| format!("golden.tsv line {at} names no manifest sample"))?;
            let kind = Kind::parse(kind)
                .ok_or_else(|| format!("golden.tsv line {at} has an unknown kind {kind:?}"))?;
            if cell(key) != key || cell(value) != value {
                return Err(format!("golden.tsv line {at} has a malformed cell"));
            }
            lines.push(Line {
                sample,
                kind,
                key: key.to_owned(),
                value: value.to_owned(),
            });
        }
        let golden = Self::new(samples.to_vec(), lines);
        if golden.text() != text {
            return Err("golden.tsv is not in canonical order or form".to_owned());
        }
        Ok(golden)
    }

    /// The SHA-256 of the whole file.
    pub(crate) fn sha256(&self) -> String {
        sha256(&self.text())
    }

    /// Per kind, in kind order: the SHA-256 of its lines in file order and
    /// their count.
    pub(crate) fn kind_digests(&self) -> Vec<(Kind, String, usize)> {
        Kind::ALL
            .into_iter()
            .map(|kind| {
                let mut hasher = Sha256::new();
                let mut count = 0;
                for line in self.lines.iter().filter(|line| line.kind == kind) {
                    hasher.update(self.line_text(line));
                    count += 1;
                }
                (kind, format!("{:x}", hasher.finalize()), count)
            })
            .collect()
    }

    /// Per stratum, in name order, given each sample's stratum: its samples
    /// and the SHA-256 of their lines in file order. A stratum digest names
    /// no sample.
    pub(crate) fn stratum_digests(&self, strata: &[&str]) -> Vec<(String, usize, String)> {
        let mut members = BTreeMap::<&str, (usize, Sha256)>::new();
        for stratum in strata {
            members
                .entry(stratum)
                .or_insert_with(|| (0, Sha256::new()))
                .0 += 1;
        }
        for line in &self.lines {
            if let Some((_, hasher)) = strata
                .get(line.sample)
                .and_then(|stratum| members.get_mut(stratum))
            {
                hasher.update(self.line_text(line));
            }
        }
        members
            .into_iter()
            .map(|(stratum, (samples, hasher))| {
                (
                    stratum.to_owned(),
                    samples,
                    format!("{:x}", hasher.finalize()),
                )
            })
            .collect()
    }

    /// The counts the ledger records over the samples.
    pub(crate) fn counts(&self) -> Counts {
        let mut counts = Counts {
            samples: self.samples.len(),
            ..Counts::default()
        };
        let mut with_check = BTreeSet::new();
        for line in &self.lines {
            match line.kind {
                Kind::Failure => counts.failed += 1,
                Kind::Allele => counts.alleles += 1,
                Kind::Check => {
                    counts.checks += 1;
                    with_check.insert(line.sample);
                }
                Kind::Note => counts.notes += 1,
                Kind::Review if line.key == "clear" && line.value == "Y" => counts.clear += 1,
                Kind::Block => {
                    let mut cells = line.value.split(' ');
                    let class = cells.next().and_then(Class::parse);
                    let checked = cells.nth(1);
                    match class {
                        Some(Class::NonHet) => {
                            counts.blocks_non_het += 1;
                            if checked == Some("no") {
                                counts.escapes += 1;
                            }
                        }
                        Some(Class::Het) => counts.blocks_het += 1,
                        Some(Class::LengthNote) => counts.blocks_length_note += 1,
                        Some(Class::TruthN) => counts.blocks_truth_n += 1,
                        Some(Class::Recheck) => counts.blocks_recheck += 1,
                        // Omitted truth alleles are golden lines but no
                        // ledger count.
                        Some(
                            Class::OmittedOther
                            | Class::OmittedNotAnalysed
                            | Class::OmittedTract
                            | Class::NoteMisread,
                        )
                        | None => {}
                    }
                }
                _ => {}
            }
        }
        counts.with_check = with_check.len();
        counts
    }
}

/// The SHA-256 of `text`, in hex.
pub(crate) fn sha256(text: &str) -> String {
    format!("{:x}", Sha256::digest(text.as_bytes()))
}

/// What a golden output holds, as the ledger counts it: `clear` and
/// `with_check` count samples, the others lines; `escapes` are non-het
/// blocks in a zone no check item touches.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub(crate) struct Counts {
    pub(crate) samples: usize,
    pub(crate) failed: usize,
    pub(crate) clear: usize,
    pub(crate) with_check: usize,
    pub(crate) checks: usize,
    pub(crate) notes: usize,
    pub(crate) alleles: usize,
    pub(crate) blocks_non_het: usize,
    pub(crate) escapes: usize,
    pub(crate) blocks_het: usize,
    pub(crate) blocks_length_note: usize,
    pub(crate) blocks_truth_n: usize,
    pub(crate) blocks_recheck: usize,
}

/// What the golden output of a check is read from: the samples in manifest
/// order, DNA's outcomes in the same order, the run's report, and the
/// review windows and zones the run counts in.
pub(crate) struct Sources<'a> {
    pub(crate) samples: &'a [Sample],
    pub(crate) outcomes: &'a [Outcome],
    pub(crate) report: &'a Report<'a>,
    pub(crate) windows: &'a [Window],
    pub(crate) zones: &'a Zones,
}

/// The golden output of one check.
pub(crate) fn build(sources: &Sources<'_>) -> Golden {
    let rows = sources.report.rows();
    let row_of = rows
        .iter()
        .enumerate()
        .map(|(index, row)| (row.sample.id.as_str(), index))
        .collect::<BTreeMap<_, _>>();
    let mut blocks = BTreeMap::<usize, Vec<(String, String)>>::new();
    for block in sources.report.unreviewed_blocks() {
        blocks.entry(block.row).or_default().push((
            format!("{}-{}:{}", block.first, block.last, block.text),
            format!("{} {} {}", block.class.label(), block.zone, block.checked),
        ));
    }
    let mut lines = Vec::new();
    for (sample, (manifest, outcome)) in sources.samples.iter().zip(sources.outcomes).enumerate() {
        let mut line = |kind: Kind, key: &str, value: &str| {
            lines.push(Line {
                sample,
                kind,
                key: cell(key),
                value: cell(value),
            });
        };
        let called = match &outcome.result {
            Ok(called) => called,
            Err(failure) => {
                line(Kind::Failure, "-", failure.kind());
                continue;
            }
        };
        for allele in &called.alleles {
            let recorded = if called.not_recorded.contains(allele) {
                "not_recorded"
            } else {
                "recorded"
            };
            line(Kind::Allele, &allele.text(), recorded);
        }
        for (first, last) in &called.intervals {
            line(Kind::Callable, &format!("{first}-{last}"), "-");
        }
        for (first, last, labels) in &called.check_labels {
            let mut labels = labels.split(' ').collect::<Vec<_>>();
            labels.sort_unstable();
            line(Kind::Check, &format!("{first}-{last}"), &labels.join(" "));
        }
        let cell_of =
            |value: Option<usize>| value.map_or_else(|| ".".to_owned(), |value| value.to_string());
        for note in &called.notes {
            for run in &note.runs {
                line(
                    Kind::Note,
                    &format!("{}:run{}", note.tract, run.run),
                    &format!(
                        "{}-{} {} {} {} {}",
                        note.first,
                        note.last,
                        run.status,
                        cell_of(run.called),
                        cell_of(run.reading),
                        cell_of(run.convention),
                    ),
                );
            }
        }
        // The review's fields read the check items alone, so they are
        // truth-free: the same as `scores.tsv` writes for a compared sample,
        // and written for every analysed one.
        let review = workload::of(
            &called.checks,
            &workload::Differences::default(),
            sources.windows,
            &[],
            Some(sources.zones),
            (&[], &[]),
        );
        let mut zones = review
            .zones
            .iter()
            .map(|&instance| sources.zones.label(instance))
            .collect::<Vec<_>>();
        zones.sort_unstable();
        line(
            Kind::Review,
            "clear",
            if review.checks == 0 { "Y" } else { "N" },
        );
        line(Kind::Review, "zones", &zones.join(","));
        line(Kind::Review, "positions", &review.positions.to_string());
        if let Some(empop) = &called.empop {
            line(Kind::Empop, "ranges", &empop.ranges);
            line(Kind::Empop, "mitotype", &empop.mitotype);
        }
        line(Kind::Rejected, "-", &called.rejected.to_string());
        let Some(&index) = row_of.get(manifest.id.as_str()) else {
            continue;
        };
        for (key, value) in blocks.remove(&index).unwrap_or_default() {
            line(Kind::Block, &key, &value);
        }
        let row = &rows[index];
        for measure in Measure::ALL {
            line(Kind::Status, measure.label(), row.status(measure).label());
        }
        line(
            Kind::Status,
            "class",
            &format!("{} {}", row.discrepancy().label(), row.class_regions()),
        );
        line(
            Kind::Status,
            "cause",
            &sources.report.cause_cells(index).join(" "),
        );
    }
    Golden::new(
        sources
            .samples
            .iter()
            .map(|sample| sample.id.clone())
            .collect(),
        lines,
    )
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;

    /// A golden line.
    pub(crate) fn line(sample: usize, kind: Kind, key: &str, value: &str) -> Line {
        Line {
            sample,
            kind,
            key: key.to_owned(),
            value: value.to_owned(),
        }
    }

    #[test]
    fn orders_lines_canonically_and_reads_them_back() -> Result<(), String> {
        let samples = vec!["S2".to_owned(), "S1".to_owned()];
        let golden = Golden::new(
            samples.clone(),
            vec![
                line(1, Kind::Status, "primary", "concordant"),
                line(0, Kind::Allele, "16189C", "recorded"),
                line(0, Kind::Allele, "309.1C", "recorded"),
                line(0, Kind::Allele, "73G", "recorded"),
                line(0, Kind::Callable, "16024-16365", "-"),
                line(0, Kind::Callable, "73-340", "-"),
                line(0, Kind::Failure, "-", "dna"),
                line(1, Kind::Allele, "263G", "recorded"),
            ],
        );
        let text = golden.text();
        assert_eq!(
            text,
            "sample_id\tkind\tkey\tvalue\n\
             S2\tfailure\t-\tdna\n\
             S2\tallele\t73G\trecorded\n\
             S2\tallele\t309.1C\trecorded\n\
             S2\tallele\t16189C\trecorded\n\
             S2\tcallable\t73-340\t-\n\
             S2\tcallable\t16024-16365\t-\n\
             S1\tallele\t263G\trecorded\n\
             S1\tstatus\tprimary\tconcordant\n"
        );
        // The order does not depend on the order the lines came in.
        let mut reversed = golden.lines.clone();
        reversed.reverse();
        assert_eq!(Golden::new(samples.clone(), reversed).text(), text);
        assert_eq!(Golden::parse(&text, &samples)?, golden);
        // A file out of order, with another sample, or another kind is refused.
        let swapped = text.replace(
            "S2\tallele\t73G\trecorded\nS2\tallele\t309.1C\trecorded\n",
            "S2\tallele\t309.1C\trecorded\nS2\tallele\t73G\trecorded\n",
        );
        assert!(Golden::parse(&swapped, &samples).is_err());
        assert!(Golden::parse(&text.replace("S1\t", "S3\t"), &samples).is_err());
        assert!(Golden::parse(&text.replace("callable", "intervals"), &samples).is_err());
        assert!(Golden::parse(&text.replace("\trecorded\n", " \n"), &samples).is_err());
        Ok(())
    }

    #[test]
    fn digests_name_kinds_and_strata_without_samples() {
        let samples = vec!["S1".to_owned(), "S2".to_owned(), "S3".to_owned()];
        let golden = Golden::new(
            samples,
            vec![
                line(0, Kind::Allele, "73G", "recorded"),
                line(1, Kind::Allele, "73G", "recorded"),
                line(2, Kind::Block, "16-16:A/G", "non_het interior@HV1 no"),
                line(2, Kind::Review, "clear", "Y"),
            ],
        );
        let kinds = golden.kind_digests();
        assert_eq!(kinds.len(), Kind::ALL.len());
        let allele = &kinds[1];
        assert_eq!(allele.0, Kind::Allele);
        assert_eq!(
            allele.1,
            sha256("S1\tallele\t73G\trecorded\nS2\tallele\t73G\trecorded\n")
        );
        assert_eq!(allele.2, 2);
        assert_eq!(kinds[0].1, sha256(""));
        let strata = golden.stratum_digests(&["a", "b", "a"]);
        assert_eq!(
            strata,
            vec![
                (
                    "a".to_owned(),
                    2,
                    sha256(
                        "S1\tallele\t73G\trecorded\nS3\treview\tclear\tY\n\
                         S3\tblock\t16-16:A/G\tnon_het interior@HV1 no\n"
                    )
                ),
                ("b".to_owned(), 1, sha256("S2\tallele\t73G\trecorded\n")),
            ]
        );
        let counts = golden.counts();
        assert_eq!(
            (
                counts.samples,
                counts.alleles,
                counts.clear,
                counts.blocks_non_het,
                counts.escapes
            ),
            (3, 2, 1, 1, 1)
        );
    }

    #[test]
    fn cells_hold_no_tab_line_end_or_trailing_space() {
        assert_eq!(cell(""), "-");
        assert_eq!(cell("a\tb\nc "), "a b c");
        assert_eq!(cell("concordant "), "concordant");
        assert_eq!(cell("1 2"), "1 2");
    }

    #[test]
    fn lines_name_their_positions() {
        assert_eq!(Line::span(Kind::Allele, "309.1C", None), Some((309, 309)));
        assert_eq!(Line::span(Kind::Allele, "249DEL", None), Some((249, 249)));
        assert_eq!(Line::span(Kind::Callable, "73-340", None), Some((73, 340)));
        assert_eq!(
            Line::span(Kind::Block, "303-315:CCCCCCCTCCCCC/CCCCCCCCTCCCCC", None),
            Some((303, 315))
        );
        assert_eq!(
            Line::span(
                Kind::Note,
                "hvs2_c_tract:run1",
                Some("303-309 unconfirmed . 8 .")
            ),
            Some((303, 309))
        );
        assert_eq!(Line::span(Kind::Review, "clear", Some("Y")), None);
    }

    /// The golden output of a synthetic panel run with `threads` workers.
    fn synthetic(threads: usize) -> Result<String, Box<dyn std::error::Error>> {
        use std::fmt::Write as _;
        use std::path::PathBuf;

        use crate::cause::forms::Forms;
        use crate::compare::{Allele, Failure};
        use crate::runner::{Called, Identity, Threads};
        use crate::truth::TruthCheck;
        use crate::zones::Widths;

        let directory = tempfile::tempdir()?;
        let trace = directory.path().join("trace.ab1");
        std::fs::write(&trace, b"trace")?;
        let lines = (0..6)
            .map(|index| {
                format!(
                    "S{index}\tb1\t{}\t1\t1\t0\t{}",
                    if index % 2 == 0 { "even" } else { "odd" },
                    trace.display()
                )
            })
            .collect::<Vec<_>>();
        let text =
            crate::manifest::tests::manifest(&lines.iter().map(String::as_str).collect::<Vec<_>>());
        let manifest = crate::manifest::parse(text.as_bytes())?;
        let mut sql = String::from(crate::truth::tests::FIXTURE);
        sql.push_str(
            "DELETE FROM samples; DELETE FROM variants; DELETE FROM sample_status; \
             DELETE FROM intervals;",
        );
        for index in 0..6 {
            let _ = write!(
                sql,
                "INSERT INTO samples VALUES ('S{index}', 'b1', 2, 1, 0, 0, 0);
                 INSERT INTO sample_status VALUES ('b1', 'S{index}', 'x', 'OK', 1, 1, '', 1, 1);
                 INSERT INTO intervals VALUES ('S{index}', 'HV1', 16024, 16365),
                    ('S{index}', 'HV2', 73, 340);
                 INSERT INTO variants VALUES ('S{index}', 'snp', '263', 263, 'A', 'G');"
            );
        }
        let path = crate::truth::tests::database(directory.path(), &sql)?;
        let wanted = manifest
            .samples
            .iter()
            .map(|sample| sample.id.as_str())
            .collect();
        let truth = crate::truth::load(&path, &wanted)?;
        let called = |sample: &Sample, _: &[PathBuf]| {
            let index = sample.id[1..]
                .parse::<usize>()
                .map_err(|_| Failure::Panic)?;
            let (alleles, checks) = match index {
                // A DNA failure.
                0 => return Err(Failure::Dna("no read".to_owned())),
                // Concordant and clear.
                1 => (vec![Allele::new("263", "A", "G")], Vec::new()),
                // A DNA-only allele in HV1's interior, under no check: an
                // escape.
                2 => (
                    vec![Allele::new("263", "A", "G"), Allele::new("16100", "C", "T")],
                    Vec::new(),
                ),
                // The same allele under a check item.
                _ => (
                    vec![Allele::new("263", "A", "G"), Allele::new("16100", "C", "T")],
                    vec![(16_098, 16_102, "one_read".to_owned())],
                ),
            };
            Ok(Called {
                components: Vec::new(),
                tracts: Vec::new(),
                length_context: None,
                not_analysed: Vec::new(),
                rejected: index % 2,
                alleles,
                not_recorded: BTreeSet::new(),
                empop: None,
                checks: checks
                    .iter()
                    .map(|(first, last, _)| (*first, *last))
                    .collect(),
                check_labels: checks,
                notes: Vec::new(),
                signals: crate::lod::Signals::default(),
                loci: crate::cause::loci::Loci::default(),
                intervals: vec![(73, 340), (16_024, 16_365)],
                identity: Identity {
                    configuration_sha256: "c".into(),
                    reference_name: "r".into(),
                    reference_sha256: "s".into(),
                    profile_id: "p".into(),
                    profile_sha256: "q".into(),
                },
            })
        };
        let outcomes = crate::runner::run_with(&manifest.samples, Threads::new(threads)?, called)?;
        let geometry = crate::geometry::tests::rcrs();
        let rows = crate::report::rows(
            &manifest.samples,
            &outcomes,
            &truth.samples,
            crate::compare::tests::rcrs(),
            (geometry, &[]),
        );
        let check = TruthCheck::new(&manifest, &truth);
        let widths = Widths {
            after: 20,
            edge: 40,
            overrun: 2,
        };
        let zones = Zones::new(geometry.regions(), geometry.repeats(), widths);
        let report = Report::new(&rows, (geometry, &check), (&Forms::default(), "")).with_review(
            10,
            None,
            Some(zones.clone()),
            None,
            Vec::new(),
        );
        let windows = workload::windows(10, geometry.regions());
        Ok(build(&Sources {
            samples: &manifest.samples,
            outcomes: &outcomes,
            report: &report,
            windows: &windows,
            zones: &zones,
        })
        .text())
    }

    #[test]
    fn the_golden_output_is_byte_stable_and_reads_dna_and_the_review()
    -> Result<(), Box<dyn std::error::Error>> {
        let text = synthetic(1)?;
        assert_eq!(text, synthetic(5)?);
        let lines = |sample: &str| {
            text.lines()
                .filter(|line| line.starts_with(&format!("{sample}\t")))
                .collect::<Vec<_>>()
        };
        // A failed sample has its failure line only.
        assert_eq!(lines("S0"), ["S0\tfailure\t-\tdna"]);
        // A clear concordant sample: its alleles, intervals, review, traces,
        // and statuses, in kind order.
        let clear = lines("S1");
        assert_eq!(
            clear[..7],
            [
                "S1\tallele\t263G\trecorded",
                "S1\tcallable\t73-340\t-",
                "S1\tcallable\t16024-16365\t-",
                "S1\treview\tclear\tY",
                "S1\treview\tpositions\t0",
                "S1\treview\tzones\t-",
                "S1\trejected\t-\t1",
            ]
        );
        assert!(clear.contains(&"S1\tstatus\thaplotype_lod\tconcordant"));
        assert!(!clear.iter().any(|line| line.contains("\tblock\t")));
        // The unchecked DNA-only allele is a non-het block in an unchecked
        // zone; under a check item it is reviewed, no block.
        assert!(lines("S2").iter().any(|line| {
            line.contains("\tblock\t16100-16100:") && line.ends_with("\tnon_het interior@HV1 no")
        }));
        let checked = lines("S3");
        assert!(checked.contains(&"S3\tcheck\t16098-16102\tone_read"));
        assert!(checked.contains(&"S3\treview\tclear\tN"));
        assert!(checked.contains(&"S3\treview\tzones\tinterior@HV1"));
        assert!(!checked.iter().any(|line| line.contains("\tblock\t")));
        let golden = Golden::parse(
            &text,
            &(0..6).map(|index| format!("S{index}")).collect::<Vec<_>>(),
        )?;
        let counts = golden.counts();
        assert_eq!(
            (
                counts.samples,
                counts.failed,
                counts.clear,
                counts.with_check,
                counts.escapes
            ),
            (6, 1, 2, 3, 1)
        );
        Ok(())
    }
}
