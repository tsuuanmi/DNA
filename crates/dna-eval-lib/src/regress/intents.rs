//! A regress intent file: the output changes a pull request declares, each
//! row naming what lines it explains (`kind`, `change`, `where`, `match`)
//! and how many samples it may touch (`max_samples`), or, with the kind
//! `workload`, a new workload counting version. It names no sample.

use std::collections::BTreeSet;

use super::golden::Kind;

/// The `kind` of a row that declares a new workload counting version.
pub(crate) const WORKLOAD: &str = "workload";

/// The column header of an intent file.
pub(crate) const COLUMNS: &str = "id\tkind\tchange\twhere\tmatch\tmax_samples\tposthoc\treason";

/// How a golden line changed between the head and the candidate.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) enum Change {
    /// The key is only in the candidate.
    Added,
    /// The key is only in the head.
    Removed,
    /// The same key has another value.
    Changed,
}

impl Change {
    pub(crate) const ALL: [Self; 3] = [Self::Added, Self::Removed, Self::Changed];

    pub(crate) const fn label(self) -> &'static str {
        match self {
            Self::Added => "added",
            Self::Removed => "removed",
            Self::Changed => "changed",
        }
    }
}

/// Where an intent row applies.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) enum Place {
    /// `*`: anywhere.
    Anywhere,
    /// `a-b`: lines whose positions overlap `a..=b`.
    Positions(i64, i64),
    /// A zone instance such as `window@HV2`.
    Zone(String),
    /// A region such as `HV1`.
    Region(String),
}

/// One declared change.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct Intent {
    pub(crate) id: String,
    /// The kind, or any kind (`*`); none for a `workload` row.
    pub(crate) kind: Option<Kind>,
    /// Whether the row declares a new workload counting version (kind
    /// `workload`, `match` a glob on the new version) rather than lines.
    pub(crate) workload: bool,
    /// The change, or any change (`any`).
    pub(crate) change: Option<Change>,
    pub(crate) place: Place,
    /// A glob (`*`, `?`) on the key or a value.
    pub(crate) pattern: String,
    /// The most samples whose lines the row may explain, when bounded.
    pub(crate) max_samples: Option<usize>,
    /// Whether the row was added after a check showed the change.
    pub(crate) posthoc: bool,
    pub(crate) reason: String,
}

/// An intent file: its header fields and rows.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct Intents {
    /// The ledger series, for example `regress-v1`.
    pub(crate) series: String,
    /// The head entry the changes are declared against.
    pub(crate) base: String,
    /// The pull request, for example `#290`.
    pub(crate) pr: String,
    pub(crate) rows: Vec<Intent>,
}

impl Intents {
    /// The rows declared after a check showed their change.
    pub(crate) fn posthoc_rows(&self) -> usize {
        self.rows.iter().filter(|row| row.posthoc).count()
    }
}

/// Reads an intent file whose zones must be among `zones` (labels such as
/// `window@HV2`) and regions among `regions`.
pub(crate) fn parse(
    text: &str,
    zones: &BTreeSet<String>,
    regions: &BTreeSet<String>,
) -> Result<Intents, String> {
    let mut lines = text
        .strip_suffix('\n')
        .ok_or("the file does not end in a newline")?
        .split('\n');
    let header = lines
        .next()
        .and_then(|line| line.strip_prefix("# regress-intent "))
        .ok_or("the first line is not '# regress-intent series=… base=… pr=…'")?;
    let mut fields = [("series", None), ("base", None), ("pr", None)];
    for field in header.split(' ') {
        let (name, value) = field
            .split_once('=')
            .ok_or_else(|| format!("header field {field:?} is not name=value"))?;
        let slot = fields
            .iter_mut()
            .find(|(known, _)| *known == name)
            .ok_or_else(|| format!("unknown header field {name:?}"))?;
        if value.is_empty() || slot.1.replace(value.to_owned()).is_some() {
            return Err(format!("header field {name} is empty or repeated"));
        }
    }
    let [series, base, pr] = fields.map(|(name, value)| value.ok_or(name));
    let missing = |name: &str| format!("the header lacks {name}=");
    let (series, base, pr) = (
        series.map_err(missing)?,
        base.map_err(missing)?,
        pr.map_err(missing)?,
    );
    if lines.next() != Some(COLUMNS) {
        return Err(format!(
            "the second line is not the column header {COLUMNS:?}"
        ));
    }
    let mut rows = Vec::<Intent>::new();
    for (number, line) in lines.enumerate() {
        let at = number + 3;
        let row = row(line, zones, regions).map_err(|message| format!("line {at}: {message}"))?;
        if rows.iter().any(|other| other.id == row.id) {
            return Err(format!("line {at}: id {} is repeated", row.id));
        }
        rows.push(row);
    }
    Ok(Intents {
        series,
        base,
        pr,
        rows,
    })
}

fn row(line: &str, zones: &BTreeSet<String>, regions: &BTreeSet<String>) -> Result<Intent, String> {
    let cells = line.split('\t').collect::<Vec<_>>();
    let [
        id,
        kind,
        change,
        place,
        pattern,
        max_samples,
        posthoc,
        reason,
    ] = cells[..]
    else {
        return Err(format!("{} cells, expected 8", cells.len()));
    };
    if id.is_empty()
        || !id
            .bytes()
            .all(|byte| byte.is_ascii_lowercase() || byte.is_ascii_digit() || b"_-".contains(&byte))
    {
        return Err(format!("id {id:?} is not [a-z0-9_-]+"));
    }
    let workload = kind == WORKLOAD;
    let kind = match kind {
        "*" | WORKLOAD => None,
        label => Some(Kind::parse(label).ok_or_else(|| format!("unknown kind {label:?}"))?),
    };
    let change = match change {
        "any" => None,
        label => Some(
            Change::ALL
                .into_iter()
                .find(|change| change.label() == label)
                .ok_or_else(|| format!("unknown change {label:?}"))?,
        ),
    };
    let place = if place == "*" {
        Place::Anywhere
    } else if let Some((first, last)) = place.split_once('-')
        && let (Ok(first), Ok(last)) = (first.parse::<i64>(), last.parse::<i64>())
    {
        if first < 1 || last < first {
            return Err(format!("position range {place:?} is empty"));
        }
        Place::Positions(first, last)
    } else if zones.contains(place) {
        Place::Zone(place.to_owned())
    } else if regions.contains(place) {
        Place::Region(place.to_owned())
    } else {
        return Err(format!(
            "where {place:?} is no position range, zone, or region of the profile"
        ));
    };
    if pattern.is_empty() {
        return Err("match is empty".to_owned());
    }
    let max_samples = match max_samples {
        "-" => None,
        text => Some(
            text.parse::<usize>()
                .ok()
                .filter(|count| *count > 0)
                .ok_or_else(|| format!("max_samples {text:?} is not a positive integer"))?,
        ),
    };
    let posthoc = match posthoc {
        "0" => false,
        "1" => true,
        other => return Err(format!("posthoc {other:?} is not 0 or 1")),
    };
    if reason.trim().is_empty() {
        return Err("reason is empty".to_owned());
    }
    if workload
        && (matches!(change, Some(Change::Added | Change::Removed))
            || place != Place::Anywhere
            || max_samples.is_some())
    {
        return Err(
            "a workload row declares a version: change changed or any, where *, max_samples -"
                .to_owned(),
        );
    }
    Ok(Intent {
        id: id.to_owned(),
        kind,
        workload,
        change,
        place,
        pattern: pattern.to_owned(),
        max_samples,
        posthoc,
        reason: reason.to_owned(),
    })
}

/// Whether `text` matches the glob `pattern`: `*` matches any run of
/// characters, `?` one character, anything else itself.
pub(crate) fn glob(pattern: &str, text: &str) -> bool {
    let (pattern, text) = (
        pattern.chars().collect::<Vec<_>>(),
        text.chars().collect::<Vec<_>>(),
    );
    let (mut p, mut t) = (0, 0);
    let mut star: Option<(usize, usize)> = None;
    while t < text.len() {
        match pattern.get(p) {
            Some('*') => {
                star = Some((p, t));
                p += 1;
            }
            Some(&symbol) if symbol == '?' || symbol == text[t] => {
                p += 1;
                t += 1;
            }
            _ => match star {
                Some((at, matched)) => {
                    p = at + 1;
                    t = matched + 1;
                    star = Some((at, matched + 1));
                }
                None => return false,
            },
        }
    }
    pattern[p..].iter().all(|symbol| *symbol == '*')
}

#[cfg(test)]
mod tests {
    use super::*;

    fn names(items: &[&str]) -> BTreeSet<String> {
        items.iter().map(|item| (*item).to_owned()).collect()
    }

    #[test]
    fn globs_match_stars_and_single_characters() {
        for (pattern, text, matched) in [
            ("*", "", true),
            ("*", "anything", true),
            ("tract_length", "tract_length", true),
            ("tract_*", "tract_length", true),
            ("*length", "tract_length", true),
            ("309.?C", "309.1C", true),
            ("309.?C", "309.12C", false),
            ("a*b*c", "axxbyyc", true),
            ("a*b*c", "axxbyy", false),
            ("one", "one_read", false),
        ] {
            assert_eq!(glob(pattern, text), matched, "{pattern} {text}");
        }
    }

    #[test]
    fn reads_every_where_form_and_refuses_malformed_rows() -> Result<(), String> {
        let zones = names(&["window@HV2", "start@HV1"]);
        let regions = names(&["HV1", "HV2", "HV3"]);
        let head = format!("# regress-intent series=regress-v1 base=r0001 pr=#290\n{COLUMNS}\n");
        let text = format!(
            "{head}note-hv2\tcheck\tremoved\twindow@HV2\ttract_length\t-\t0\t#222 length notes\n\
             note-hv2b\tnote\tadded\t303-316\t*\t12\t1\t#222 the same samples\n\
             any\t*\tany\tHV1\t*\t-\t0\tADR-0090 reads\n"
        );
        let intents = parse(&text, &zones, &regions)?;
        assert_eq!(
            (
                intents.series.as_str(),
                intents.base.as_str(),
                intents.pr.as_str()
            ),
            ("regress-v1", "r0001", "#290")
        );
        assert_eq!(intents.rows.len(), 3);
        assert_eq!(intents.rows[0].place, Place::Zone("window@HV2".to_owned()));
        assert_eq!(intents.rows[1].place, Place::Positions(303, 316));
        assert_eq!(intents.rows[1].max_samples, Some(12));
        assert_eq!(intents.rows[2].place, Place::Region("HV1".to_owned()));
        assert_eq!((intents.rows[2].kind, intents.rows[2].change), (None, None));
        assert_eq!(intents.posthoc_rows(), 1);
        assert!(!intents.rows[0].workload);
        let counting = parse(
            &format!("{head}counting\tworkload\tchanged\t*\t2\t-\t0\t#282 version 2\n"),
            &zones,
            &regions,
        )?;
        assert!(counting.rows[0].workload && counting.rows[0].kind.is_none());
        for bad in [
            "x\tworkload\tadded\t*\t2\t-\t0\treason\n",
            "x\tworkload\tchanged\tHV1\t2\t-\t0\treason\n",
            "x\tworkload\tany\t*\t2\t3\t0\treason\n",
            "x\tcheck\tremoved\t*\t*\t-\t0\n",
            "X\tcheck\tremoved\t*\t*\t-\t0\treason\n",
            "x\tchecks\tremoved\t*\t*\t-\t0\treason\n",
            "x\tcheck\tgone\t*\t*\t-\t0\treason\n",
            "x\tcheck\tremoved\twindow@HV9\t*\t-\t0\treason\n",
            "x\tcheck\tremoved\t20-10\t*\t-\t0\treason\n",
            "x\tcheck\tremoved\t*\t*\t0\t0\treason\n",
            "x\tcheck\tremoved\t*\t*\t-\t2\treason\n",
            "x\tcheck\tremoved\t*\t*\t-\t0\t \n",
            "x\tcheck\tremoved\t*\t\t-\t0\treason\n",
            "x\tcheck\tremoved\t*\t*\t-\t0\treason\nx\tnote\tadded\t*\t*\t-\t0\treason\n",
        ] {
            assert!(
                parse(&format!("{head}{bad}"), &zones, &regions).is_err(),
                "{bad}"
            );
        }
        for header in [
            "# regress-intent series=regress-v1 base=r0001\n",
            "# regress-intent series=regress-v1 base=r0001 pr=#1 pr=#2\n",
            "# intent series=regress-v1 base=r0001 pr=#1\n",
        ] {
            assert!(parse(&format!("{header}{COLUMNS}\n"), &zones, &regions).is_err());
        }
        Ok(())
    }
}
