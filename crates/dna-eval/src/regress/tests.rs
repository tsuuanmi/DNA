//! Tests of the regress command's comparability, `record`, and the baseline
//! rows, on synthetic run headers and ledgers.

use std::cell::Cell;
use std::collections::BTreeMap;

use super::*;

#[test]
fn the_files_follow_the_contract() -> std::result::Result<(), Box<dyn std::error::Error>> {
    let schema: serde_json::Value = serde_json::from_str(include_str!(
        "../../../../docs/reference/schemas/eval-regress-v1.schema.json"
    ))?;
    let definitions = &schema["$defs"];
    let strings = |value: &serde_json::Value| {
        value
            .as_array()
            .map(|items| {
                items
                    .iter()
                    .filter_map(serde_json::Value::as_str)
                    .map(str::to_owned)
                    .collect::<Vec<_>>()
            })
            .unwrap_or_default()
    };
    let header_of = |name: &str| {
        format!(
            "{}\n",
            strings(&definitions[name]["properties"]["columns"]["const"]).join("\t")
        )
    };
    assert_eq!(header_of("golden"), golden::HEADER);
    assert_eq!(header_of("changes"), judge::HEADER);
    assert_eq!(
        strings(
            &definitions["golden"]["properties"]["rows"]["items"]["properties"]["kind"]["enum"]
        ),
        Kind::ALL.map(|kind| kind.label().to_owned())
    );
    let names = |value: &serde_json::Value| {
        value
            .as_object()
            .map(|object| object.keys().cloned().collect::<BTreeSet<_>>())
            .unwrap_or_default()
    };
    let written = serde_json::to_value(header())?;
    let required = strings(&definitions["run"]["required"])
        .into_iter()
        .collect::<BTreeSet<_>>();
    assert_eq!(names(&written), required);
    assert_eq!(
        names(&written["counts"]),
        names(&definitions["run"]["properties"]["counts"]["properties"])
    );
    Ok(())
}

/// A clean head-like run header.
fn header() -> Header {
    Header {
        harness: "dna-eval".into(),
        harness_version: "0.1.0".into(),
        source: Source {
            revision: Some("a".repeat(40)),
            dirty: Some(false),
        },
        threads: 8,
        series: "regress-v1".into(),
        manifest: ManifestHeader {
            path: "regress-v1.tsv".into(),
            fields: BTreeMap::from([("panel".to_owned(), "regress".to_owned())]),
            body_sha256: "b".repeat(64),
            rows: 2,
        },
        truth: TruthIdentity {
            path: "statistics.sqlite".into(),
            sha256: "c".repeat(64),
            digest: "d".repeat(64),
            recorded_digest: Some("d".repeat(64)),
            state: "unchanged".into(),
            changed_samples: 0,
        },
        recheck: Some(RecheckHeader {
            path: "recheck.tsv".into(),
            sha256: "e".repeat(64),
            rows: 3,
        }),
        configuration_path: "dna.toml".into(),
        reference_path: "rCRS.fasta".into(),
        profile_sha256: "f".repeat(64),
        reference_sha256: "1".repeat(64),
        review_window: 10,
        zones: ZoneHeader {
            after: 20,
            edge: 40,
            overrun: 2,
        },
        run_context: None,
        version_key: VersionKey {
            dna_version: "0.1.0".into(),
            plugins: "consensus=21,sanger=22".into(),
            contracts: "dna.analysis/v10,dna.variants/v1".into(),
            configuration_sha256: "2".repeat(64),
        },
        evidence_cache: None,
        placement_cache: None,
        head: None,
        intent: None,
        rebase: false,
        verdict: Verdict::Initial.label().into(),
        rules: Vec::new(),
        warnings: Vec::new(),
        golden_sha256: golden::sha256(golden::HEADER),
        kinds: vec![KindDigest {
            kind: "failure".into(),
            lines: 0,
            sha256: golden::sha256(""),
        }],
        strata: vec![
            StratumDigest {
                stratum: "a".into(),
                samples: 1,
                golden_sha256: "3".repeat(64),
                changed_samples: 0,
            },
            StratumDigest {
                stratum: "b".into(),
                samples: 1,
                golden_sha256: "4".repeat(64),
                changed_samples: 0,
            },
        ],
        counts: Counts {
            samples: 2,
            ..Counts::default()
        },
        changes: ChangeCounts::default(),
        seconds: Seconds {
            total: 1.0,
            truth: 0.5,
            dna: 0.25,
        },
    }
}

#[test]
fn the_run_header_reads_back() -> Result<()> {
    let header = header();
    let text = header.text()?;
    let read: Header = serde_json::from_str(&text).map_err(|error| refused(error.to_string()))?;
    assert_eq!(read, header);
    Ok(())
}

#[test]
fn refuses_checks_that_measure_something_else() {
    let head = header();
    assert!(comparable(("r0001", &head), &head, false).is_ok());
    let changed = |change: fn(&mut Header)| {
        let mut candidate = head.clone();
        change(&mut candidate);
        candidate
    };
    // A new configuration is a DNA change: accepted, and reported.
    let configured = changed(|header| header.version_key.configuration_sha256 = "9".repeat(64));
    assert!(comparable(("r0001", &head), &configured, false).is_ok());
    for (name, candidate) in [
        (
            "manifest",
            changed(|header| header.manifest.body_sha256 = "9".repeat(64)),
        ),
        (
            "profile",
            changed(|header| header.profile_sha256 = "9".repeat(64)),
        ),
        (
            "reference",
            changed(|header| header.reference_sha256 = "9".repeat(64)),
        ),
        ("window", changed(|header| header.review_window = 20)),
        ("zones", changed(|header| header.zones.edge = 30)),
        (
            "contexts",
            changed(|header| {
                header.run_context = Some(RunContextHeader {
                    directory: "contexts".into(),
                    batches: 1,
                    sha256: "9".repeat(64),
                });
            }),
        ),
        (
            "truth",
            changed(|header| header.truth.digest = "9".repeat(64)),
        ),
        ("re-check", changed(|header| header.recheck = None)),
    ] {
        assert!(
            matches!(
                comparable(("r0001", &head), &candidate, false),
                Err(Error::Regress(_))
            ),
            "{name}"
        );
    }
    // A rebase accepts another truth or re-check list, nothing else, and
    // needs one of them to differ.
    let truth = changed(|header| header.truth.digest = "9".repeat(64));
    assert!(comparable(("r0001", &head), &truth, true).is_ok());
    let recheck = changed(|header| header.recheck = None);
    assert!(comparable(("r0001", &head), &recheck, true).is_ok());
    assert!(comparable(("r0001", &head), &head, true).is_err());
    let both = changed(|header| {
        header.truth.digest = "9".repeat(64);
        header.version_key.configuration_sha256 = "9".repeat(64);
    });
    assert!(comparable(("r0001", &head), &both, true).is_err());
    let profile = changed(|header| {
        header.truth.digest = "9".repeat(64);
        header.profile_sha256 = "9".repeat(64);
    });
    assert!(comparable(("r0001", &head), &profile, true).is_err());
}

#[test]
fn only_regression_manifests_are_checked() -> std::result::Result<(), crate::error::ManifestError> {
    let manifest = |panel: &str| -> std::result::Result<Manifest, crate::error::ManifestError> {
        crate::manifest::parse(
            crate::manifest::tests::manifest(&[])
                .replacen("panel=dev", panel, 1)
                .as_bytes(),
        )
    };
    assert!(regress_panel(&manifest("panel=regress")?).is_ok());
    for panel in [
        "panel=dev",
        "panel=readiness",
        "panel=validation",
        "name=regress",
    ] {
        assert!(
            matches!(regress_panel(&manifest(panel)?), Err(Error::Regress(_))),
            "{panel}"
        );
    }
    Ok(())
}

#[test]
fn an_intent_file_names_no_sample() {
    let samples = BTreeSet::from(["243256", "S1"]);
    assert!(names_no_sample("id\treason\nx\t#222 length notes 303-316\n", &samples).is_ok());
    assert!(names_no_sample("x\tsee 243256\n", &samples).is_err());
    assert!(names_no_sample("x\t*\tS1\n", &samples).is_err());
    assert!(names_no_sample("x\tS12 and 2432567\n", &samples).is_ok());
}

/// A git that answers from fixed values.
struct Fake {
    on_main: bool,
    blob: String,
    asked: Cell<usize>,
}

impl Git for Fake {
    fn blob(&self, _: &Path) -> std::result::Result<String, String> {
        Ok(self.blob.clone())
    }

    fn committed(&self, _: &Path, _: &str) -> std::result::Result<String, String> {
        Ok(self.blob.clone())
    }

    fn on_main(&self, _: &Path, _: &str) -> std::result::Result<bool, String> {
        self.asked.set(self.asked.get() + 1);
        Ok(self.on_main)
    }
}

/// A check output directory of `header`, a ledger directory, and the two
/// baseline files, in a temporary directory.
struct Fixture {
    _directory: tempfile::TempDir,
    args: RecordArgs,
}

fn fixture(header: &Header) -> std::result::Result<Fixture, Box<dyn std::error::Error>> {
    let directory = tempfile::tempdir()?;
    let run = directory.path().join("check");
    fs::create_dir(&run)?;
    fs::write(run.join(GOLDEN), golden::HEADER)?;
    fs::write(run.join(CHANGES), judge::HEADER)?;
    fs::write(run.join(REPORT), "# Regress check\n")?;
    fs::write(run.join(RUN_JSON), header.text()?)?;
    let baseline = directory.path().join("regress-v1-ledger.tsv");
    fs::write(
        &baseline,
        format!(
            "{}\n{}\n",
            ledger::baseline_comment("regress-ledger", "regress-v1"),
            ledger::LEDGER_COLUMNS
        ),
    )?;
    let strata = directory.path().join("regress-v1-strata.tsv");
    fs::write(
        &strata,
        format!(
            "{}\n{}\n",
            ledger::baseline_comment("regress-strata", "regress-v1"),
            ledger::STRATA_COLUMNS
        ),
    )?;
    let args = RecordArgs {
        run,
        ledger: directory.path().join("ledger"),
        entry: "r0001".into(),
        baseline,
        strata,
    };
    Ok(Fixture {
        _directory: directory,
        args,
    })
}

fn git(on_main: bool) -> Fake {
    Fake {
        on_main,
        blob: "5".repeat(40),
        asked: Cell::new(0),
    }
}

#[test]
fn records_a_clean_checked_run_and_appends_the_baselines()
-> std::result::Result<(), Box<dyn std::error::Error>> {
    let header = header();
    let fixture = fixture(&header)?;
    let git = git(true);
    let message = record(&fixture.args, &git, "2026-10-22")?;
    assert!(message.contains("recorded r0001 (initial"), "{message}");
    assert_eq!(git.asked.get(), 1);
    let ledger = &fixture.args.ledger;
    for name in [GOLDEN, CHANGES, REPORT, RUN_JSON] {
        assert!(ledger.join("r0001").join(name).is_file(), "{name}");
    }
    let index = fs::read_to_string(ledger.join(ledger::INDEX))?;
    assert_eq!(
        index,
        format!(
            "{}r0001\tinitial\t-\t{}\t{}\tREGRESS INITIAL\t2026-10-22\n",
            ledger::INDEX_HEADER,
            "a".repeat(40),
            header.golden_sha256
        )
    );
    // The baseline row has one cell per column, in column order.
    let baseline = fs::read_to_string(&fixture.args.baseline)?;
    let row = baseline.lines().last().unwrap_or_default();
    let columns = ledger::LEDGER_COLUMNS.split('\t').collect::<Vec<_>>();
    let cells = row.split('\t').collect::<Vec<_>>();
    assert_eq!(cells.len(), columns.len());
    let cell = |name: &str| {
        columns
            .iter()
            .position(|column| *column == name)
            .map(|index| cells[index])
    };
    assert_eq!(cell("entry"), Some("r0001"));
    assert_eq!(cell("date"), Some("2026-10-22"));
    assert_eq!(cell("kind"), Some("initial"));
    assert_eq!(cell("base"), Some("-"));
    assert_eq!(cell("plugins"), Some("consensus=21,sanger=22"));
    assert_eq!(cell("recheck_sha256"), Some("e".repeat(64).as_str()));
    let kinds = format!("failure={}", golden::sha256(""));
    assert_eq!(cell("kinds_sha256"), Some(kinds.as_str()));
    assert_eq!(cell("samples"), Some("2"));
    assert_eq!(cell("intent_blob"), Some("-"));
    assert_eq!(cell("verdict"), Some("REGRESS INITIAL"));
    let strata = fs::read_to_string(&fixture.args.strata)?;
    assert!(strata.ends_with(&format!(
        "r0001\ta\t1\t{}\t0\nr0001\tb\t1\t{}\t0\n",
        "3".repeat(64),
        "4".repeat(64)
    )));
    // The same entry again is refused.
    assert!(record(&fixture.args, &git, "2026-10-22").is_err());
    Ok(())
}

#[test]
fn record_refuses_dirty_failed_moved_existing_and_unmerged_runs()
-> std::result::Result<(), Box<dyn std::error::Error>> {
    let refuses = |header: Header,
                   git: &Fake,
                   entry: &str|
     -> std::result::Result<bool, Box<dyn std::error::Error>> {
        let mut fixture = fixture(&header)?;
        fixture.args.entry = entry.to_owned();
        Ok(matches!(
            record(&fixture.args, git, "2026-10-22"),
            Err(Error::Regress(_))
        ))
    };
    let main = git(true);
    let mut dirty = header();
    dirty.source.dirty = Some(true);
    assert!(refuses(dirty, &main, "r0001")?);
    let mut unknown = header();
    unknown.source.revision = None;
    assert!(refuses(unknown, &main, "r0001")?);
    let mut failed = header();
    failed.verdict = Verdict::Fail.label().into();
    assert!(refuses(failed, &main, "r0001")?);
    // The run compared with a head the ledger does not have.
    let mut moved = header();
    moved.head = Some("r0001".into());
    assert!(refuses(moved, &main, "r0002")?);
    assert!(refuses(header(), &main, "1")?);
    assert!(refuses(header(), &git(false), "r0001")?);
    // An entry the ledger already holds.
    let existing = fixture(&header())?;
    fs::create_dir_all(existing.args.ledger.join("r0001"))?;
    assert!(matches!(
        record(&existing.args, &main, "2026-10-22"),
        Err(Error::Regress(_))
    ));
    // A golden output that is not the one the header records.
    let other = fixture(&header())?;
    fs::write(
        other.args.run.join(GOLDEN),
        format!("{}S1\tfailure\t-\tdna\n", golden::HEADER),
    )?;
    assert!(matches!(
        record(&other.args, &main, "2026-10-22"),
        Err(Error::Regress(_))
    ));
    Ok(())
}

#[test]
fn a_change_entry_follows_its_head() -> std::result::Result<(), Box<dyn std::error::Error>> {
    let fixture = fixture(&header())?;
    record(&fixture.args, &git(true), "2026-10-22")?;
    let mut next = header();
    next.head = Some("r0001".into());
    next.verdict = Verdict::Explained.label().into();
    next.changes = ChangeCounts {
        changed_samples: 1,
        changed_lines: 2,
        explained_lines: 2,
        posthoc_rows: 1,
    };
    next.intent = Some(IntentHeader {
        path: fixture.args.run.join(REPORT),
        blob: "5".repeat(40),
        pr: "#1".into(),
        rows: 2,
        posthoc_rows: 1,
    });
    fs::write(fixture.args.run.join(RUN_JSON), next.text()?)?;
    let args = RecordArgs {
        entry: "r0002".into(),
        ..fixture.args
    };
    record(&args, &git(true), "2026-10-23")?;
    let baseline = fs::read_to_string(&args.baseline)?;
    let row = baseline.lines().last().unwrap_or_default();
    assert!(
        row.starts_with("r0002\t2026-10-23\tchange\tr0001\t"),
        "{row}"
    );
    assert!(row.ends_with(&format!(
        "\t1\t2\t2\t1\t{}\tREGRESS EXPLAINED",
        "5".repeat(40)
    )));
    assert!(args.ledger.join("r0002").join(INTENT).is_file());
    // An intent file changed since the check is refused.
    let mut later = next.clone();
    later.head = Some("r0002".into());
    fs::write(args.run.join(RUN_JSON), later.text()?)?;
    let changed = Fake {
        blob: "6".repeat(40),
        ..git(true)
    };
    let args = RecordArgs {
        entry: "r0003".into(),
        ..args
    };
    assert!(record(&args, &changed, "2026-10-24").is_err());
    Ok(())
}
