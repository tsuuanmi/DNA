//! The changes of a candidate golden output against the ledger head, the
//! intent each one matches, and the verdict: the forbidden rules first,
//! then the declared intents.

use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write as _;

use super::golden::{Golden, Kind, Line, key_order};
use super::intents::{Change, Intent, Intents, Place, glob};
use crate::geometry::Span;
use crate::workload::Class;
use crate::zones::Zones;

/// The header of `changes.tsv`.
pub(crate) const HEADER: &str = "sample_id\tstratum\tkind\tkey\tchange\tbefore\tafter\tintent\n";

/// A rule of the verdict, in the order the rules are evaluated.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) enum Rule {
    /// A sample's unreviewed non-het blocks in a zone no check item touches
    /// rise in number.
    Escape,
    /// A flagged discordant sample becomes clear while still discordant.
    FlagLost,
    /// The set of failed samples differs.
    Failures,
    /// A truth-free line changed while the version key did not.
    NoBump,
    /// The workload counting version moved and no `workload` intent row
    /// declares the new one.
    WorkloadUndeclared,
    /// An intent row explains more samples than it allows.
    IntentExceeded,
    /// A changed line matches no intent row.
    Unexplained,
    /// With `--rebase`, a truth-free line changed.
    RebaseDnaChange,
}

impl Rule {
    /// Every rule, in evaluation order.
    #[cfg(test)]
    pub(crate) const ALL: [Self; 8] = [
        Self::Escape,
        Self::FlagLost,
        Self::Failures,
        Self::NoBump,
        Self::WorkloadUndeclared,
        Self::IntentExceeded,
        Self::Unexplained,
        Self::RebaseDnaChange,
    ];

    pub(crate) const fn label(self) -> &'static str {
        match self {
            Self::Escape => "forbidden:escape",
            Self::FlagLost => "forbidden:flag_lost",
            Self::Failures => "forbidden:failures",
            Self::NoBump => "forbidden:no_bump",
            Self::WorkloadUndeclared => "workload_undeclared",
            Self::IntentExceeded => "intent_exceeded",
            Self::Unexplained => "unexplained",
            Self::RebaseDnaChange => "rebase_dna_change",
        }
    }
}

/// The verdict of a check.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Verdict {
    /// No head: the golden output seeds the ledger.
    Initial,
    /// No line changed.
    Same,
    /// Every changed line is explained and no rule fired.
    Explained,
    /// With `--rebase`, only evaluation lines changed.
    Rebase,
    /// A rule fired.
    Fail,
}

impl Verdict {
    pub(crate) const fn label(self) -> &'static str {
        match self {
            Self::Initial => "REGRESS INITIAL",
            Self::Same => "REGRESS SAME",
            Self::Explained => "REGRESS EXPLAINED",
            Self::Rebase => "REGRESS REBASE",
            Self::Fail => "REGRESS FAIL",
        }
    }

    /// The verdict a label of `run.json` names.
    pub(crate) fn parse(label: &str) -> Option<Self> {
        [
            Self::Initial,
            Self::Same,
            Self::Explained,
            Self::Rebase,
            Self::Fail,
        ]
        .into_iter()
        .find(|verdict| verdict.label() == label)
    }
}

/// What explains a changed line.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Label {
    /// The intent row of that index.
    Intent(usize),
    /// An evaluation line of a rebase.
    Rebase,
    /// A rule forbids it.
    Forbidden(Rule),
    /// Nothing.
    Unexplained,
}

/// One changed golden line.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct Changed {
    pub(crate) sample: usize,
    pub(crate) kind: Kind,
    pub(crate) key: String,
    pub(crate) change: Change,
    pub(crate) before: Option<String>,
    pub(crate) after: Option<String>,
    pub(crate) label: Label,
}

impl Changed {
    /// The positions the line names, from its key or either value.
    fn span(&self) -> Option<(i64, i64)> {
        Line::span(self.kind, &self.key, self.after.as_deref())
            .or_else(|| Line::span(self.kind, &self.key, self.before.as_deref()))
    }

    /// The values the line had or has.
    fn values(&self) -> impl Iterator<Item = &str> {
        self.before.iter().chain(&self.after).map(String::as_str)
    }

    /// Whether the line is explained: by an intent row or as a rebase.
    pub(crate) const fn explained(&self) -> bool {
        matches!(self.label, Label::Intent(_) | Label::Rebase)
    }
}

/// The changed lines of `candidate` against `head`, over the same samples,
/// in golden order. Lines are keyed by sample, kind, and key; a key a sample
/// repeats is compared occurrence by occurrence.
pub(crate) fn changes(head: &Golden, candidate: &Golden) -> Vec<Changed> {
    type Key<'g> = (usize, Kind, (bool, u64, &'g str), usize);
    fn keyed(golden: &Golden) -> BTreeMap<Key<'_>, &str> {
        let mut seen = BTreeMap::<(usize, Kind, &str), usize>::new();
        golden
            .lines
            .iter()
            .map(|line| {
                let occurrence = seen.entry((line.sample, line.kind, &line.key)).or_default();
                *occurrence += 1;
                (
                    (line.sample, line.kind, key_order(&line.key), *occurrence),
                    line.value.as_str(),
                )
            })
            .collect()
    }
    let (before, after) = (keyed(head), keyed(candidate));
    let keys = before.keys().chain(after.keys()).collect::<BTreeSet<_>>();
    keys.into_iter()
        .filter_map(|key| {
            let (old, new) = (before.get(key).copied(), after.get(key).copied());
            let change = match (old, new) {
                (Some(old), Some(new)) if old == new => return None,
                (Some(_), Some(_)) => Change::Changed,
                (None, Some(_)) => Change::Added,
                (Some(_), None) => Change::Removed,
                (None, None) => return None,
            };
            let (sample, kind, (_, _, text), _) = *key;
            Some(Changed {
                sample,
                kind,
                key: text.to_owned(),
                change,
                before: old.map(str::to_owned),
                after: new.map(str::to_owned),
                label: Label::Unexplained,
            })
        })
        .collect()
}

/// What the verdict reads besides the two golden outputs.
pub(crate) struct Context<'a> {
    /// The zones the intents' `where` names.
    pub(crate) zones: &'a Zones,
    /// The profile's regions.
    pub(crate) regions: &'a [Span],
    /// Whether the check is a rebase.
    pub(crate) rebase: bool,
    /// Whether the candidate's version key differs from the head's.
    pub(crate) version_moved: bool,
    /// The workload counting versions of the head and the candidate.
    pub(crate) workload: (u32, u32),
}

/// A check's judgement.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct Judgement {
    pub(crate) changes: Vec<Changed>,
    /// The rules that fired, in rule order.
    pub(crate) rules: Vec<Rule>,
    pub(crate) verdict: Verdict,
    /// Per intent row: the lines and distinct samples it explains.
    pub(crate) intents: Vec<(usize, usize)>,
    /// With `--rebase`: blocks that became non-het escapes (reported; a
    /// rebase changes the truth they are measured against).
    pub(crate) rebase_escapes: usize,
    /// The changed escape lines, by index into `changes`, of each sample and
    /// zone that gains an escape line while its escape count does not rise
    /// (reported: an escape reshaped, never `forbidden:escape`).
    pub(crate) reshaped_escapes: Vec<usize>,
    /// The `workload` intent row that declares a moved counting version.
    pub(crate) workload_intent: Option<usize>,
}

impl Judgement {
    /// Samples with a changed line.
    pub(crate) fn changed_samples(&self) -> BTreeSet<usize> {
        self.changes.iter().map(|change| change.sample).collect()
    }

    /// Changed lines an intent row or the rebase explains.
    pub(crate) fn explained_lines(&self) -> usize {
        self.changes
            .iter()
            .filter(|change| change.explained())
            .count()
    }
}

/// The zone of a block value that is a non-het block in an unchecked zone.
fn escape(value: Option<&str>) -> Option<&str> {
    let cells = value?.split(' ').collect::<Vec<_>>();
    (cells.first().and_then(|class| Class::parse(class)) == Some(Class::NonHet)
        && cells.get(2) == Some(&"no"))
    .then_some(cells.get(1).copied())
    .flatten()
}

/// The change of each sample's escape count per zone, candidate against
/// head; the unchanged lines count alike in both.
fn escape_counts(changes: &[Changed]) -> BTreeMap<(usize, String), i64> {
    let mut counts = BTreeMap::<(usize, String), i64>::new();
    for change in changes.iter().filter(|change| change.kind == Kind::Block) {
        for (value, step) in [(&change.before, -1), (&change.after, 1)] {
            if let Some(zone) = escape(value.as_deref()) {
                *counts.entry((change.sample, zone.to_owned())).or_default() += step;
            }
        }
    }
    counts
}

/// The value of `kind` `key` of each sample.
fn values<'g>(golden: &'g Golden, kind: Kind, key: &str) -> BTreeMap<usize, &'g str> {
    golden
        .lines
        .iter()
        .filter(|line| line.kind == kind && line.key == key)
        .map(|line| (line.sample, line.value.as_str()))
        .collect()
}

/// What a sample's changed blocks and check items touch: their zones and
/// spans, which a line without a position of its own is placed by.
#[derive(Default)]
struct Touched {
    zones: BTreeSet<String>,
    spans: Vec<(i64, i64)>,
}

/// The zones a changed line lies in: a block's zone, or the zones a check
/// item touches; none for other kinds.
fn zones_of(change: &Changed, zones: &Zones) -> BTreeSet<String> {
    match change.kind {
        Kind::Block => change
            .values()
            .filter_map(|value| value.split(' ').nth(1))
            .map(str::to_owned)
            .collect(),
        Kind::Check => change
            .span()
            .map(|(first, last)| {
                zones
                    .touched(first, last)
                    .into_iter()
                    .map(|instance| zones.label(instance))
                    .collect()
            })
            .unwrap_or_default(),
        _ => BTreeSet::new(),
    }
}

/// Whether `intent` explains `change`, placed by `touched` when it names
/// no position.
fn matches(
    intent: &Intent,
    change: &Changed,
    touched: Option<&Touched>,
    context: &Context<'_>,
) -> bool {
    if intent.workload
        || intent.kind.is_some_and(|kind| kind != change.kind)
        || intent.change.is_some_and(|kind| kind != change.change)
    {
        return false;
    }
    let span = change.span();
    let overlaps = |(first, last): (i64, i64), (from, to): (i64, i64)| first <= to && from <= last;
    let region = |name: &str| {
        context
            .regions
            .iter()
            .find(|region| region.name == name)
            .map(|region| (region.first, region.last))
    };
    let placed = match &intent.place {
        Place::Anywhere => true,
        Place::Positions(first, last) => span.is_some_and(|span| overlaps(span, (*first, *last))),
        Place::Zone(zone) => match change.kind {
            Kind::Block | Kind::Check => zones_of(change, context.zones).contains(zone),
            Kind::Allele | Kind::Callable | Kind::Note => false,
            _ => touched.is_some_and(|touched| touched.zones.contains(zone)),
        },
        Place::Region(name) => region(name).is_some_and(|bounds| match span {
            Some(span) => overlaps(span, bounds),
            None => touched
                .is_some_and(|touched| touched.spans.iter().any(|span| overlaps(*span, bounds))),
        }),
    };
    placed
        && (glob(&intent.pattern, &change.key)
            || change.values().any(|value| {
                glob(&intent.pattern, value)
                    || (change.kind == Kind::Check
                        && value.split(' ').any(|label| glob(&intent.pattern, label)))
            }))
}

/// Judges `candidate` against `head` with the declared `intents`.
pub(crate) fn judge(
    head: &Golden,
    candidate: &Golden,
    intents: Option<&Intents>,
    context: &Context<'_>,
) -> Judgement {
    let mut changes = changes(head, candidate);
    let mut rules = BTreeSet::new();
    let mut rebase_escapes = 0;

    // A sample's escapes in a zone that rise in number; an escape whose
    // text or span changes while the count holds is only reported.
    let counts = escape_counts(&changes);
    let rises = |sample: usize, zone: &str| {
        counts
            .get(&(sample, zone.to_owned()))
            .is_some_and(|step| *step > 0)
    };
    let reshaped = counts
        .keys()
        .filter(|(sample, zone)| {
            !context.rebase
                && !rises(*sample, zone)
                && changes.iter().any(|change| {
                    change.sample == *sample
                        && change.kind == Kind::Block
                        && escape(change.after.as_deref()) == Some(zone.as_str())
                })
        })
        .collect::<BTreeSet<_>>();
    let reshaped_escapes = changes
        .iter()
        .enumerate()
        .filter(|(_, change)| {
            change.kind == Kind::Block
                && [&change.before, &change.after].into_iter().any(|value| {
                    escape(value.as_deref())
                        .is_some_and(|zone| reshaped.contains(&(change.sample, zone.to_owned())))
                })
        })
        .map(|(index, _)| index)
        .collect::<Vec<_>>();

    // A flagged discordant sample that becomes clear while still discordant.
    let lod = |golden| values(golden, Kind::Status, "haplotype_lod");
    let clear = |golden| values(golden, Kind::Review, "clear");
    let (lod_before, lod_after) = (lod(head), lod(candidate));
    let (clear_before, clear_after) = (clear(head), clear(candidate));
    let flag_lost = |sample: usize| {
        lod_before.get(&sample) == Some(&"discordant")
            && lod_after.get(&sample) == Some(&"discordant")
            && clear_before.get(&sample) == Some(&"N")
            && clear_after.get(&sample) == Some(&"Y")
    };
    // Samples whose failure line came or went.
    let failed = changes
        .iter()
        .filter(|change| change.kind == Kind::Failure && change.change != Change::Changed)
        .map(|change| change.sample)
        .collect::<BTreeSet<_>>();

    for change in &mut changes {
        let forbidden = if !context.rebase
            && change.kind == Kind::Block
            && escape(change.after.as_deref()).is_some_and(|zone| rises(change.sample, zone))
        {
            Some(Rule::Escape)
        } else if !context.rebase
            && change.kind == Kind::Review
            && change.key == "clear"
            && flag_lost(change.sample)
        {
            Some(Rule::FlagLost)
        } else if failed.contains(&change.sample) {
            Some(Rule::Failures)
        } else if change.kind.truth_free() && !context.version_moved && !context.rebase {
            Some(Rule::NoBump)
        } else if change.kind.truth_free() && context.rebase {
            Some(Rule::RebaseDnaChange)
        } else {
            None
        };
        if context.rebase
            && change.kind == Kind::Block
            && escape(change.after.as_deref()).is_some()
            && escape(change.before.as_deref()).is_none()
        {
            rebase_escapes += 1;
        }
        if let Some(rule) = forbidden {
            change.label = Label::Forbidden(rule);
            rules.insert(rule);
        } else if context.rebase {
            change.label = Label::Rebase;
        }
    }

    // Lines without a position of their own are placed by their sample's
    // changed blocks and check items.
    let mut touched = BTreeMap::<usize, Touched>::new();
    for change in changes
        .iter()
        .filter(|change| matches!(change.kind, Kind::Block | Kind::Check))
    {
        let entry = touched.entry(change.sample).or_default();
        entry.zones.extend(zones_of(change, context.zones));
        entry.spans.extend(change.span());
    }
    let rows = intents.map_or(&[][..], |intents| intents.rows.as_slice());
    let mut explained = vec![(0, BTreeSet::new()); rows.len()];
    for change in &mut changes {
        if change.label != Label::Unexplained {
            continue;
        }
        if let Some(index) = rows
            .iter()
            .position(|intent| matches(intent, change, touched.get(&change.sample), context))
        {
            change.label = Label::Intent(index);
            explained[index].0 += 1;
            explained[index].1.insert(change.sample);
        } else {
            rules.insert(Rule::Unexplained);
        }
    }
    // A moved counting version reclassifies evaluation lines without any
    // version-key change: a `workload` row must declare the new version.
    let (head_counting, counting) = context.workload;
    let workload_intent = (head_counting != counting)
        .then(|| {
            rows.iter()
                .position(|intent| intent.workload && glob(&intent.pattern, &counting.to_string()))
        })
        .flatten();
    if head_counting != counting && workload_intent.is_none() {
        rules.insert(Rule::WorkloadUndeclared);
    }
    for (intent, (_, samples)) in rows.iter().zip(&explained) {
        if intent
            .max_samples
            .is_some_and(|maximum| samples.len() > maximum)
        {
            rules.insert(Rule::IntentExceeded);
        }
    }
    let verdict = if !rules.is_empty() {
        Verdict::Fail
    } else if context.rebase {
        Verdict::Rebase
    } else if changes.is_empty() {
        Verdict::Same
    } else {
        Verdict::Explained
    };
    Judgement {
        changes,
        rules: rules.into_iter().collect(),
        verdict,
        intents: explained
            .into_iter()
            .map(|(lines, samples)| (lines, samples.len()))
            .collect(),
        rebase_escapes,
        reshaped_escapes,
        workload_intent,
    }
}

/// `changes.tsv`: one row per changed line, in golden order, with each
/// sample's stratum and what explains the line.
pub(crate) fn table(
    judgement: &Judgement,
    samples: &[String],
    strata: &[&str],
    intents: Option<&Intents>,
) -> String {
    let mut text = String::from(HEADER);
    for change in &judgement.changes {
        let label = match change.label {
            Label::Intent(index) => intents
                .and_then(|intents| intents.rows.get(index))
                .map_or("unexplained", |intent| intent.id.as_str()),
            Label::Rebase => "rebase",
            Label::Forbidden(rule) => rule.label(),
            Label::Unexplained => "unexplained",
        };
        let _ = writeln!(
            text,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{label}",
            samples[change.sample],
            strata[change.sample],
            change.kind.label(),
            change.key,
            change.change.label(),
            change.before.as_deref().unwrap_or("-"),
            change.after.as_deref().unwrap_or("-"),
        );
    }
    text
}

#[cfg(test)]
mod tests {
    use super::super::golden::tests::line;
    use super::super::intents;
    use super::*;
    use crate::geometry::tests::rcrs;
    use crate::zones::Widths;

    fn golden(lines: Vec<Line>) -> Golden {
        Golden::new(
            vec!["S1".to_owned(), "S2".to_owned(), "S3".to_owned()],
            lines,
        )
    }

    fn zones() -> Zones {
        Zones::new(
            rcrs().regions(),
            rcrs().repeats(),
            Widths {
                after: 20,
                edge: 40,
                overrun: 2,
            },
        )
    }

    fn judged(
        head: &Golden,
        candidate: &Golden,
        intents: Option<&str>,
        flags: (bool, bool),
    ) -> Judgement {
        counted(head, candidate, intents, flags, (2, 2))
    }

    /// Judges with the head's and the candidate's counting versions.
    fn counted(
        head: &Golden,
        candidate: &Golden,
        intents: Option<&str>,
        (rebase, version_moved): (bool, bool),
        workload: (u32, u32),
    ) -> Judgement {
        let zones = zones();
        let labels = zones
            .instances()
            .into_iter()
            .map(|instance| zones.label(instance))
            .collect();
        let regions = rcrs()
            .regions()
            .iter()
            .map(|region| region.name.clone())
            .collect();
        let intents = intents.map(|rows| {
            intents::parse(
                &format!(
                    "# regress-intent series=regress-v1 base=r0001 pr=#1\n{}\n{rows}",
                    intents::COLUMNS
                ),
                &labels,
                &regions,
            )
            .unwrap_or_else(|error| panic!("{error}"))
        });
        judge(
            head,
            candidate,
            intents.as_ref(),
            &Context {
                zones: &zones,
                regions: rcrs().regions(),
                rebase,
                version_moved,
                workload,
            },
        )
    }

    #[test]
    fn a_moved_counting_version_needs_a_workload_row() {
        let block = |value: &str| golden(vec![line(0, Kind::Block, "16-16:A/G", value)]);
        let (head, candidate) = (
            block("het interior@HV1 no"),
            block("recheck interior@HV1 no"),
        );
        let lines = "blocks\tblock\tchanged\t*\t*\t-\t0\t#1 the new classes\n";
        let declared = format!("{lines}counting\tworkload\tany\t*\t2\t-\t0\t#1 version 2\n");
        // Moved and undeclared: the lines are explained, the version is not.
        let undeclared = counted(&head, &candidate, Some(lines), (false, false), (1, 2));
        assert_eq!(undeclared.rules, vec![Rule::WorkloadUndeclared]);
        assert_eq!(undeclared.verdict, Verdict::Fail);
        // A row of another version does not declare it.
        let other = lines.to_owned() + "counting\tworkload\tany\t*\t3\t-\t0\t#1 version 3\n";
        assert_eq!(
            counted(&head, &candidate, Some(&other), (false, false), (1, 2)).rules,
            vec![Rule::WorkloadUndeclared]
        );
        let explained = counted(&head, &candidate, Some(&declared), (false, false), (1, 2));
        assert_eq!(explained.verdict, Verdict::Explained);
        assert_eq!(explained.workload_intent, Some(1));
        // The workload row explains no line, and a moved version without a
        // changed line is SAME once declared.
        assert_eq!(explained.intents[1], (0, 0));
        let same = counted(&head, &head, Some(&declared), (false, false), (1, 2));
        assert_eq!(same.verdict, Verdict::Same);
        // An unmoved version needs no row; a row then declares nothing.
        let unmoved = counted(&head, &candidate, Some(&declared), (false, false), (2, 2));
        assert_eq!(unmoved.verdict, Verdict::Explained);
        assert_eq!(unmoved.workload_intent, None);
    }

    #[test]
    fn detects_added_removed_and_changed_lines() {
        let head = golden(vec![
            line(0, Kind::Allele, "73G", "recorded"),
            line(0, Kind::Allele, "263G", "recorded"),
            line(1, Kind::Check, "303-315", "tract_length"),
            line(1, Kind::Check, "303-315", "one_read"),
        ]);
        let candidate = golden(vec![
            line(0, Kind::Allele, "73G", "recorded"),
            line(0, Kind::Allele, "263G", "not_recorded"),
            line(0, Kind::Allele, "16189C", "recorded"),
            line(1, Kind::Check, "303-315", "one_read"),
        ]);
        let changes = changes(&head, &candidate);
        let summary = changes
            .iter()
            .map(|change| {
                (
                    change.sample,
                    change.key.as_str(),
                    change.change,
                    change.before.as_deref(),
                    change.after.as_deref(),
                )
            })
            .collect::<Vec<_>>();
        assert_eq!(
            summary,
            vec![
                (
                    0,
                    "263G",
                    Change::Changed,
                    Some("recorded"),
                    Some("not_recorded")
                ),
                (0, "16189C", Change::Added, None, Some("recorded")),
                // A repeated key compares occurrence by occurrence, in
                // value order.
                (1, "303-315", Change::Removed, Some("tract_length"), None),
            ]
        );
        assert!(super::changes(&head, &head).is_empty());
    }

    #[test]
    fn no_change_is_same_and_any_change_without_intents_fails() {
        let head = golden(vec![line(0, Kind::Allele, "73G", "recorded")]);
        let same = judged(&head, &head, None, (false, true));
        assert_eq!((same.verdict, same.rules.len()), (Verdict::Same, 0));
        let candidate = golden(vec![line(0, Kind::Allele, "73G", "not_recorded")]);
        let judgement = judged(&head, &candidate, None, (false, true));
        assert_eq!(judgement.verdict, Verdict::Fail);
        assert_eq!(judgement.rules, vec![Rule::Unexplained]);
    }

    #[test]
    fn intents_match_by_every_where_form_glob_and_bound() {
        let head = golden(vec![
            line(0, Kind::Check, "303-315", "one_read tract_length"),
            line(0, Kind::Review, "clear", "N"),
            line(1, Kind::Allele, "16189C", "recorded"),
            line(
                2,
                Kind::Note,
                "hvs2_c_tract:run1",
                "303-309 unconfirmed . 8 .",
            ),
        ]);
        let candidate = golden(vec![
            line(0, Kind::Review, "clear", "Y"),
            line(1, Kind::Allele, "16189C", "not_recorded"),
            line(
                2,
                Kind::Note,
                "hvs2_c_tract:run1",
                "303-309 confirmed 8 8 .",
            ),
        ]);
        // Zone, region, positions, and glob forms each explain their line;
        // the review line follows its sample's removed check item's zone.
        let judgement = judged(
            &head,
            &candidate,
            Some(
                "check\tcheck\tremoved\twindow@HV2\ttract_*\t1\t0\t#1 checks\n\
                 review\treview\tchanged\twindow@HV2\t*\t-\t0\t#1 clear\n\
                 hv1\tallele\tany\tHV1\t*\t-\t0\t#1 HV1\n\
                 note\tnote\tchanged\t300-305\t*confirmed*\t-\t1\t#1 notes\n\
                 unused\tempop\tany\t*\t*\t-\t0\t#1 nothing\n",
            ),
            (false, true),
        );
        assert_eq!(judgement.verdict, Verdict::Explained, "{judgement:?}");
        assert_eq!(
            judgement
                .changes
                .iter()
                .map(|change| change.label)
                .collect::<Vec<_>>(),
            vec![
                Label::Intent(0),
                Label::Intent(1),
                Label::Intent(2),
                Label::Intent(3)
            ]
        );
        assert_eq!(judgement.intents[4], (0, 0));
        // A label of a check matches on its own; another region, zone, or
        // range does not.
        for (rows, explained) in [
            ("a\tcheck\tremoved\t*\tone_read\t-\t0\tr\n", true),
            ("a\tcheck\tremoved\t*\tone\t-\t0\tr\n", false),
            ("a\tcheck\tremoved\tHV1\t*\t-\t0\tr\n", false),
            ("a\tcheck\tremoved\twindow@HV1\t*\t-\t0\tr\n", false),
            ("a\tcheck\tremoved\t316-400\t*\t-\t0\tr\n", false),
            ("a\tcheck\tadded\t*\t*\t-\t0\tr\n", false),
            ("a\t*\tany\t315-315\t*\t-\t0\tr\n", true),
        ] {
            let head = golden(vec![line(
                0,
                Kind::Check,
                "303-315",
                "one_read tract_length",
            )]);
            let judgement = judged(&head, &golden(Vec::new()), Some(rows), (false, true));
            assert_eq!(judgement.verdict == Verdict::Explained, explained, "{rows}");
        }
        // A row bounded to one sample fails on two.
        let head = golden(vec![
            line(0, Kind::Allele, "73G", "recorded"),
            line(1, Kind::Allele, "73G", "recorded"),
        ]);
        let candidate = golden(Vec::new());
        let judgement = judged(
            &head,
            &candidate,
            Some("a\tallele\tremoved\t*\t73G\t1\t0\tr\n"),
            (false, true),
        );
        assert_eq!(judgement.rules, vec![Rule::IntentExceeded]);
        assert_eq!(judgement.intents, vec![(2, 2)]);
    }

    #[test]
    fn forbidden_rules_fail_whatever_the_intents() {
        let any = Some("all\t*\tany\t*\t*\t-\t0\tr\n");
        // A new non-het block in an unchecked zone, or a block becoming one.
        let head = golden(vec![line(
            0,
            Kind::Block,
            "16-16:A/G",
            "het interior@HV1 no",
        )]);
        let candidate = golden(vec![line(
            0,
            Kind::Block,
            "16-16:A/G",
            "non_het interior@HV1 no",
        )]);
        assert_eq!(
            judged(&head, &candidate, any, (false, true)).rules,
            vec![Rule::Escape]
        );
        let added = golden(vec![
            line(0, Kind::Block, "16-16:A/G", "het interior@HV1 no"),
            line(1, Kind::Block, "20-20:A/G", "non_het interior@HV1 no"),
        ]);
        assert_eq!(
            judged(&head, &added, any, (false, true)).rules,
            vec![Rule::Escape]
        );
        // A checked non-het block is no escape.
        let checked = golden(vec![line(
            0,
            Kind::Block,
            "16-16:A/G",
            "non_het interior@HV1 yes",
        )]);
        assert_eq!(
            judged(&head, &checked, any, (false, true)).verdict,
            Verdict::Explained
        );
        // A discordant sample that loses its flag.
        let flagged = |clear: &str| {
            golden(vec![
                line(0, Kind::Review, "clear", clear),
                line(0, Kind::Status, "haplotype_lod", "discordant"),
            ])
        };
        assert_eq!(
            judged(&flagged("N"), &flagged("Y"), any, (false, true)).rules,
            vec![Rule::FlagLost]
        );
        // A sample that fails, and all its lines with it.
        let analysed = golden(vec![
            line(0, Kind::Allele, "73G", "recorded"),
            line(0, Kind::Status, "primary", "concordant"),
        ]);
        let failed = golden(vec![line(0, Kind::Failure, "-", "dna")]);
        let judgement = judged(&analysed, &failed, any, (false, true));
        assert_eq!(judgement.rules, vec![Rule::Failures]);
        assert!(
            judgement
                .changes
                .iter()
                .all(|change| change.label == Label::Forbidden(Rule::Failures))
        );
        // A truth-free change without a version bump; an evaluation change
        // needs none.
        let head = golden(vec![
            line(0, Kind::Allele, "73G", "recorded"),
            line(0, Kind::Status, "primary", "concordant"),
        ]);
        let candidate = golden(vec![
            line(0, Kind::Allele, "73G", "not_recorded"),
            line(0, Kind::Status, "primary", "discordant"),
        ]);
        let judgement = judged(&head, &candidate, any, (false, false));
        assert_eq!(judgement.rules, vec![Rule::NoBump]);
        assert_eq!(judgement.changes[1].label, Label::Intent(0));
        assert_eq!(
            judged(&head, &candidate, any, (false, true)).verdict,
            Verdict::Explained
        );
    }

    #[test]
    fn escapes_fail_only_when_a_sample_count_rises_in_a_zone() {
        let any = Some("all\t*\tany\t*\t*\t-\t0\tr\n");
        let escape = |sample, key| line(sample, Kind::Block, key, "non_het interior@HV1 no");
        let head = golden(vec![escape(0, "16-16:A/G"), escape(1, "30-30:C/T")]);
        // The same escape with another text or span, count unchanged: no
        // forbidden rule, the intents judge the lines, and the changed
        // escape lines are reported.
        for key in ["16-16:A/R", "16-17:AG/GA"] {
            let candidate = golden(vec![escape(0, key), escape(1, "30-30:C/T")]);
            let judgement = judged(&head, &candidate, any, (false, true));
            assert_eq!(judgement.verdict, Verdict::Explained, "{key}");
            assert_eq!(judgement.reshaped_escapes, vec![0, 1], "{key}");
            let judgement = judged(&head, &candidate, None, (false, true));
            assert_eq!(judgement.rules, vec![Rule::Unexplained], "{key}");
        }
        // A count that rises fails, whatever the intents; only the new
        // escape lines of the risen zone carry the rule.
        let risen = golden(vec![
            escape(0, "16-16:A/G"),
            escape(0, "20-20:A/G"),
            escape(1, "30-30:C/T"),
        ]);
        let judgement = judged(&head, &risen, any, (false, true));
        assert_eq!(judgement.rules, vec![Rule::Escape]);
        assert_eq!(judgement.changes[0].label, Label::Forbidden(Rule::Escape));
        assert!(judgement.reshaped_escapes.is_empty());
        // An escape moving into another zone raises that zone's count.
        let moved = golden(vec![
            line(0, Kind::Block, "16-16:A/G", "non_het after@HV1 no"),
            escape(1, "30-30:C/T"),
        ]);
        assert_eq!(
            judged(&head, &moved, any, (false, true)).rules,
            vec![Rule::Escape]
        );
        // One escape removed in one sample and one added in another: the
        // sample whose count rose fails; the other is no reshape.
        let swapped = golden(vec![escape(0, "16-16:A/G"), escape(2, "40-40:G/A")]);
        let judgement = judged(&head, &swapped, any, (false, true));
        assert_eq!(judgement.rules, vec![Rule::Escape]);
        let labels = judgement
            .changes
            .iter()
            .map(|change| (change.sample, change.label))
            .collect::<Vec<_>>();
        assert_eq!(
            labels,
            vec![(1, Label::Intent(0)), (2, Label::Forbidden(Rule::Escape))]
        );
        assert!(judgement.reshaped_escapes.is_empty());
        // A fixed escape is neither forbidden nor reported.
        let fixed = golden(vec![escape(0, "16-16:A/G")]);
        let judgement = judged(&head, &fixed, any, (false, true));
        assert_eq!(judgement.verdict, Verdict::Explained);
        assert!(judgement.reshaped_escapes.is_empty());
    }

    #[test]
    fn a_rebase_accepts_evaluation_changes_only() {
        let head = golden(vec![
            line(0, Kind::Allele, "73G", "recorded"),
            line(0, Kind::Block, "16-16:A/G", "recheck interior@HV1 no"),
            line(0, Kind::Status, "primary", "discordant"),
        ]);
        // The truth moved: evaluation lines change freely, even into an
        // escape, which is reported.
        let candidate = golden(vec![
            line(0, Kind::Allele, "73G", "recorded"),
            line(0, Kind::Block, "16-16:A/G", "non_het interior@HV1 no"),
            line(0, Kind::Status, "primary", "concordant"),
        ]);
        let judgement = judged(&head, &candidate, None, (true, false));
        assert_eq!(judgement.verdict, Verdict::Rebase);
        assert_eq!(judgement.rebase_escapes, 1);
        assert!(
            judgement
                .changes
                .iter()
                .all(|change| change.label == Label::Rebase)
        );
        assert_eq!(
            judged(&head, &head, None, (true, false)).verdict,
            Verdict::Rebase
        );
        // A truth-free line may not change in a rebase.
        let changed = golden(vec![
            line(0, Kind::Allele, "73G", "not_recorded"),
            line(0, Kind::Block, "16-16:A/G", "recheck interior@HV1 no"),
            line(0, Kind::Status, "primary", "discordant"),
        ]);
        let judgement = judged(&head, &changed, None, (true, false));
        assert_eq!(judgement.rules, vec![Rule::RebaseDnaChange]);
    }

    #[test]
    fn writes_the_change_table() {
        let head = golden(vec![line(1, Kind::Allele, "73G", "recorded")]);
        let candidate = golden(vec![line(1, Kind::Allele, "73G", "not_recorded")]);
        let judgement = judged(&head, &candidate, None, (false, true));
        assert_eq!(
            table(&judgement, &head.samples, &["a", "b", "c"], None),
            format!("{HEADER}S2\tb\tallele\t73G\tchanged\trecorded\tnot_recorded\tunexplained\n")
        );
    }
}
