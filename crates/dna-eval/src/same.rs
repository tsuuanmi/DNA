//! `dna-eval same`: whether two runs' outputs are byte-identical, and which
//! files differ.
//!
//! Every file of a run but `run.json` is compared by its SHA-256. A run
//! records the digests of its files in `run.json` (`outputs`); they are
//! checked against the files, so a file changed after its run is refused
//! rather than compared. A run of an earlier revision records none, and its
//! files are hashed alone. `run.json` itself records paths, timings, and the
//! harness revision, so it is never compared.

use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::fs::File;
use std::io::{BufRead as _, BufReader};
use std::path::Path;

use crate::error::DiffError;
use crate::rundir::{self, RUN_JSON};

/// The comparison of two runs' outputs.
pub(crate) struct Sameness {
    /// The report, for people.
    pub(crate) report: String,
    /// Whether every output is in both runs and byte-identical.
    pub(crate) same: bool,
}

/// How one file compares.
#[derive(Debug, PartialEq, Eq)]
enum Comparison {
    Identical,
    /// Both runs hold the file with other bytes; the first line that differs,
    /// when the files can be read by line.
    Differs(Option<usize>),
    OnlyIn(char),
}

/// Compares the outputs of the runs in `first` and `second`.
///
/// # Errors
///
/// Returns [`DiffError`] when a run cannot be read, or a file differs from
/// the digest its `run.json` records.
pub(crate) fn compare(first: &Path, second: &Path) -> Result<Sameness, DiffError> {
    let (a, a_recorded) = outputs(first)?;
    let (b, b_recorded) = outputs(second)?;
    let mut names = a.keys().chain(b.keys()).collect::<Vec<_>>();
    names.sort_unstable();
    names.dedup();
    let comparisons = names
        .into_iter()
        .map(|name| {
            let comparison = match (a.get(name), b.get(name)) {
                (Some(left), Some(right)) if left == right => Comparison::Identical,
                (Some(_), Some(_)) => {
                    Comparison::Differs(first_difference(&first.join(name), &second.join(name)))
                }
                (Some(_), None) => Comparison::OnlyIn('a'),
                (None, _) => Comparison::OnlyIn('b'),
            };
            (name.as_str(), comparison)
        })
        .collect::<Vec<_>>();
    let same = comparisons
        .iter()
        .all(|(_, comparison)| *comparison == Comparison::Identical);
    Ok(Sameness {
        report: report(
            (first, a_recorded),
            (second, b_recorded),
            &comparisons,
            same,
        ),
        same,
    })
}

/// The digest of every output of the run in `directory`, and whether its
/// `run.json` recorded them.
fn outputs(directory: &Path) -> Result<(BTreeMap<String, String>, bool), DiffError> {
    let path = directory.join(RUN_JSON);
    let failure = |message: String| DiffError::Run {
        path: path.clone(),
        message,
    };
    let text = std::fs::read(&path).map_err(|error| failure(error.to_string()))?;
    let header = serde_json::from_slice::<serde_json::Value>(&text)
        .map_err(|error| failure(error.to_string()))?;
    let recorded = match header.get("outputs") {
        None => None,
        Some(value) => Some(
            serde_json::from_value::<BTreeMap<String, String>>(value.clone())
                .map_err(|error| failure(format!("outputs: {error}")))?,
        ),
    };
    let found = rundir::digests(directory)?;
    if let Some(recorded) = &recorded
        && *recorded != found
    {
        let changed = recorded
            .keys()
            .chain(found.keys())
            .filter(|name| recorded.get(*name) != found.get(*name))
            .map(String::as_str)
            .collect::<std::collections::BTreeSet<_>>();
        return Err(failure(format!(
            "the files differ from the digests the run recorded: {}",
            changed.into_iter().collect::<Vec<_>>().join(", ")
        )));
    }
    Ok((found, recorded.is_some()))
}

/// The first line, counted from 1, at which two text files differ, or
/// `None` when either cannot be read by line.
fn first_difference(first: &Path, second: &Path) -> Option<usize> {
    let open = |path: &Path| File::open(path).ok().map(BufReader::new);
    let (mut a, mut b) = (open(first)?, open(second)?);
    let (mut left, mut right) = (Vec::new(), Vec::new());
    for number in 1.. {
        left.clear();
        right.clear();
        let read_left = a.read_until(b'\n', &mut left).ok()?;
        let read_right = b.read_until(b'\n', &mut right).ok()?;
        if left != right {
            return Some(number);
        }
        if read_left == 0 && read_right == 0 {
            return None;
        }
    }
    None
}

fn report(
    first: (&Path, bool),
    second: (&Path, bool),
    comparisons: &[(&str, Comparison)],
    same: bool,
) -> String {
    let mut text = String::new();
    for (label, (path, recorded)) in [("a", first), ("b", second)] {
        let _ = writeln!(
            text,
            "{label}: {}{}",
            path.display(),
            if recorded {
                ""
            } else {
                " (run.json records no output digests; the files were hashed)"
            }
        );
    }
    let identical = comparisons
        .iter()
        .filter(|(_, comparison)| *comparison == Comparison::Identical)
        .count();
    let _ = writeln!(text, "identical: {identical} file(s)");
    for (name, comparison) in comparisons {
        match comparison {
            Comparison::Identical => {}
            Comparison::Differs(Some(line)) => {
                let _ = writeln!(text, "differs: {name} (first at line {line})");
            }
            Comparison::Differs(None) => {
                let _ = writeln!(text, "differs: {name}");
            }
            Comparison::OnlyIn(side) => {
                let _ = writeln!(text, "only in {side}: {name}");
            }
        }
    }
    text.push_str(if same {
        "SAME: every output is byte-identical\n"
    } else {
        "DIFFERENT\n"
    });
    text
}

#[cfg(test)]
mod tests {
    use std::fs;

    use super::*;
    use crate::rundir::{CHECKS, SCORES, SIGNALS, write};

    /// A run directory holding `files`, with `run.json` recording their
    /// digests when `record` is set.
    fn run(
        root: &Path,
        name: &str,
        files: &[(&str, &str)],
        record: bool,
    ) -> Result<std::path::PathBuf, Box<dyn std::error::Error>> {
        let directory = root.join(name);
        fs::create_dir(&directory)?;
        for (file, contents) in files {
            write(&directory, file, contents)?;
        }
        let header = if record {
            serde_json::json!({ "outputs": rundir::digests(&directory)? })
        } else {
            serde_json::json!({})
        };
        write(&directory, RUN_JSON, &header.to_string())?;
        Ok(directory)
    }

    #[test]
    fn reports_identical_runs_the_same_whatever_their_headers()
    -> Result<(), Box<dyn std::error::Error>> {
        let root = tempfile::tempdir()?;
        let files = [(SCORES, "a\n1\n"), (CHECKS, "b\n")];
        let a = run(root.path(), "a", &files, true)?;
        let b = run(root.path(), "b", &files, false)?;
        let sameness = compare(&a, &b)?;
        assert!(sameness.same, "{}", sameness.report);
        assert!(sameness.report.contains("identical: 2 file(s)"));
        assert!(sameness.report.contains("records no output digests"));
        Ok(())
    }

    #[test]
    fn names_the_files_that_differ_and_the_first_line() -> Result<(), Box<dyn std::error::Error>> {
        let root = tempfile::tempdir()?;
        let a = run(
            root.path(),
            "a",
            &[(SCORES, "h\n1\n2\n"), (CHECKS, "b\n"), (SIGNALS, "s\n")],
            true,
        )?;
        let b = run(
            root.path(),
            "b",
            &[(SCORES, "h\n1\n3\n"), (CHECKS, "b\n")],
            true,
        )?;
        let sameness = compare(&a, &b)?;
        assert!(!sameness.same);
        assert!(
            sameness
                .report
                .contains("differs: scores.tsv (first at line 3)"),
            "{}",
            sameness.report
        );
        assert!(sameness.report.contains("only in a: signals.tsv"));
        assert!(sameness.report.contains("identical: 1 file(s)"));
        assert!(sameness.report.ends_with("DIFFERENT\n"));
        // A file one line longer differs at the line the other lacks.
        write(&a, CHECKS, "b\nc\n")?;
        assert_eq!(first_difference(&a.join(CHECKS), &b.join(CHECKS)), Some(2));
        assert_eq!(first_difference(&b.join(CHECKS), &b.join(CHECKS)), None);
        Ok(())
    }

    #[test]
    fn refuses_a_file_changed_after_its_run() -> Result<(), Box<dyn std::error::Error>> {
        let root = tempfile::tempdir()?;
        let a = run(root.path(), "a", &[(SCORES, "x\n")], true)?;
        write(&a, SCORES, "y\n")?;
        assert!(matches!(
            compare(&a, &a),
            Err(DiffError::Run { message, .. })
                if message == "the files differ from the digests the run recorded: scores.tsv"
        ));
        Ok(())
    }
}
