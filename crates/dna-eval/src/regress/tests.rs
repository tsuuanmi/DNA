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
        manifest::parse(
            manifest::tests::manifest(&[])
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
    let samples = BTreeSet::from(["900001", "S1"]);
    assert!(names_no_sample("id\treason\nx\t#222 length notes 303-316\n", &samples).is_ok());
    assert!(names_no_sample("x\tsee 900001\n", &samples).is_err());
    assert!(names_no_sample("x\t*\tS1\n", &samples).is_err());
    assert!(names_no_sample("x\tS12 and 9000012\n", &samples).is_ok());
}

/// A git that answers from fixed values.
struct Fake {
    main: String,
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

    fn main_revision(&self, _: &Path) -> std::result::Result<String, String> {
        self.asked.set(self.asked.get() + 1);
        Ok(self.main.clone())
    }
}

/// A clean build's revision in the synthetic headers.
fn clean() -> String {
    "a".repeat(40)
}

/// A git whose main is at `main`.
fn git(main: &str) -> Fake {
    Fake {
        main: main.to_owned(),
        blob: "5".repeat(40),
        asked: Cell::new(0),
    }
}

/// The two committed baseline files, empty, under `directory`.
fn baselines(directory: &Path) -> std::io::Result<(PathBuf, PathBuf)> {
    let baseline = directory.join("regress-v1-ledger.tsv");
    fs::write(
        &baseline,
        format!(
            "{}\n{}\n",
            ledger::baseline_comment("regress-ledger", "regress-v1"),
            ledger::LEDGER_COLUMNS
        ),
    )?;
    let strata = directory.join("regress-v1-strata.tsv");
    fs::write(
        &strata,
        format!(
            "{}\n{}\n",
            ledger::baseline_comment("regress-strata", "regress-v1"),
            ledger::STRATA_COLUMNS
        ),
    )?;
    Ok((baseline, strata))
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
    let (baseline, strata) = baselines(directory.path())?;
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
    let main = git(&clean());
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
    // A build of an ancestor of main, or of a branch, is not main's commit.
    let ahead = git(&"9".repeat(40));
    assert!(refuses(header(), &ahead, "r0001")?);
    assert_eq!(ahead.asked.get(), 1);
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
    // A strata baseline that does not end at the head.
    let behind = fixture(&header())?;
    fs::write(
        &behind.args.strata,
        format!(
            "{}\n{}\nr0001\ta\t1\t{}\t0\n",
            ledger::baseline_comment("regress-strata", "regress-v1"),
            ledger::STRATA_COLUMNS,
            "3".repeat(64)
        ),
    )?;
    let mut later = behind.args;
    later.entry = "r0002".into();
    assert!(matches!(
        record(&later, &main, "2026-10-22"),
        Err(Error::Regress(_))
    ));
    Ok(())
}

/// A synthetic regression set on disk: a truth database, a `panel=regress`
/// manifest derived from it, traces DNA cannot read, a re-check list, the
/// repository's configuration and reference, a ledger, and the two
/// baselines.
struct Set {
    directory: tempfile::TempDir,
    manifest: PathBuf,
    truth: PathBuf,
    recheck: PathBuf,
    ledger: PathBuf,
    baseline: PathBuf,
    strata: PathBuf,
}

impl Set {
    fn new() -> std::result::Result<Self, Box<dyn std::error::Error>> {
        let directory = tempfile::tempdir()?;
        let truth = crate::truth::tests::database(directory.path(), crate::truth::tests::FIXTURE)?;
        let database = crate::truth::load(&truth, &BTreeSet::from(["S1", "S2", "S3"]))?;
        let unreadable = directory.path().join("unreadable.ab1");
        fs::write(&unreadable, b"not a trace")?;
        let missing = directory.path().join("missing.ab1");
        let rows = [
            format!("S1\tb1\tbench\t1\t1\t0\t{}", unreadable.display()),
            format!("S2\tb2\trecent\t1\t1\t0\t{}", missing.display()),
            format!(
                "S3\tb1\trecent\t1\t1\t0\t{};{}",
                unreadable.display(),
                missing.display()
            ),
        ];
        let seed = manifest::tests::manifest(&rows.each_ref().map(String::as_str)).replacen(
            "panel=dev",
            "panel=regress",
            1,
        );
        let parent = manifest::parse(seed.as_bytes())?;
        let source = crate::panel::Source {
            path: &truth,
            sha256: "abc",
            database: &database,
        };
        let manifest = directory.path().join("regress-v1.tsv");
        fs::write(
            &manifest,
            crate::panel::derive(&parent, &source, "v1", "2026-10-22").text,
        )?;
        let recheck = directory.path().join("recheck.tsv");
        fs::write(&recheck, format!("{}\n", crate::recheck::HEADER))?;
        let (baseline, strata) = baselines(directory.path())?;
        Ok(Self {
            manifest,
            truth,
            recheck,
            ledger: directory.path().join("ledger"),
            baseline,
            strata,
            directory,
        })
    }

    /// `regress check` arguments writing to `output`.
    fn check(&self, output: &str) -> CheckArgs {
        let repository = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
        CheckArgs {
            manifest: self.manifest.clone(),
            truth: self.truth.clone(),
            recheck: self.recheck.clone(),
            dna: DnaOptions {
                config: repository.join("config/dna.toml"),
                reference: repository.join("references/rCRS.fasta"),
                threads: 2,
                evidence_cache: None,
                run_context: None,
            },
            review: ReviewOptions {
                review_window: crate::conventions::REVIEW_WINDOW,
                zone_after: crate::conventions::ZONE_AFTER,
                zone_edge: crate::conventions::ZONE_EDGE,
                zone_overrun: crate::conventions::ZONE_OVERRUN,
            },
            forms: None,
            ledger: self.ledger.clone(),
            intent: None,
            rebase: false,
            output: self.directory.path().join(output),
        }
    }

    /// `regress record` arguments for the check output `run` as `entry`.
    fn record(&self, run: &str, entry: &str) -> RecordArgs {
        RecordArgs {
            run: self.directory.path().join(run),
            ledger: self.ledger.clone(),
            entry: entry.to_owned(),
            baseline: self.baseline.clone(),
            strata: self.strata.clone(),
        }
    }

    /// Marks the check output `run` as a clean build of [`clean`], as a
    /// build of main's commit records it.
    fn as_clean(&self, run: &str) -> std::result::Result<(), Box<dyn std::error::Error>> {
        let path = self.directory.path().join(run).join(RUN_JSON);
        let mut header: Header = serde_json::from_str(&fs::read_to_string(&path)?)?;
        header.source = Source {
            revision: Some(clean()),
            dirty: Some(false),
        };
        fs::write(&path, header.text()?)?;
        Ok(())
    }
}

#[test]
fn checks_and_records_a_synthetic_set_end_to_end()
-> std::result::Result<(), Box<dyn std::error::Error>> {
    let set = Set::new()?;
    let main = git(&clean());

    // The first check seeds the series.
    let (report, verdict) = check(&set.check("first"), &main)?;
    assert_eq!(verdict, Verdict::Initial);
    assert!(report.contains("**Verdict: REGRESS INITIAL**"), "{report}");
    let first = set.directory.path().join("first");
    let golden_text = fs::read_to_string(first.join(GOLDEN))?;
    // Every sample fails: S1's trace is unreadable, S2's is missing.
    assert!(golden_text.starts_with(golden::HEADER));
    assert_eq!(
        golden_text
            .lines()
            .skip(1)
            .map(|line| line.split('\t').take(2).collect::<Vec<_>>().join(" "))
            .collect::<Vec<_>>(),
        ["S1 failure", "S2 failure", "S3 failure"]
    );
    assert_eq!(fs::read_to_string(first.join(CHANGES))?, judge::HEADER);
    // The same check again writes the same golden output.
    check(&set.check("again"), &main)?;
    assert_eq!(
        fs::read_to_string(set.directory.path().join("again").join(GOLDEN))?,
        golden_text
    );

    // A run whose run.json claims another verdict is not recorded.
    set.as_clean("first")?;
    let path = first.join(RUN_JSON);
    let honest = fs::read_to_string(&path)?;
    let mut claimed: Header = serde_json::from_str(&honest)?;
    claimed.counts.failed = 0;
    fs::write(&path, claimed.text()?)?;
    let refused = record(&set.record("first", "r0001"), &main, "2026-10-22");
    assert!(
        matches!(&refused, Err(Error::Regress(message)) if message.contains("re-judgement")),
        "{refused:?}"
    );
    assert!(!set.ledger.join("r0001").exists());
    fs::write(&path, honest)?;

    // The first entry.
    let message = record(&set.record("first", "r0001"), &main, "2026-10-22")?;
    assert!(message.contains("recorded r0001 (initial"), "{message}");
    for name in [GOLDEN, CHANGES, REPORT, RUN_JSON] {
        assert!(set.ledger.join("r0001").join(name).is_file(), "{name}");
    }
    let index = fs::read_to_string(set.ledger.join(ledger::INDEX))?;
    assert_eq!(
        index,
        format!(
            "{}r0001\tinitial\t-\t{}\t{}\tREGRESS INITIAL\t2026-10-22\n",
            ledger::INDEX_HEADER,
            clean(),
            golden::sha256(&golden_text)
        )
    );
    // The baseline row has one cell per column, in column order.
    let baseline = fs::read_to_string(&set.baseline)?;
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
    assert_eq!(cell("kind"), Some("initial"));
    assert_eq!(cell("base"), Some("-"));
    assert_eq!(cell("samples"), Some("3"));
    assert_eq!(cell("failed"), Some("3"));
    assert_eq!(cell("intent_blob"), Some("-"));
    assert_eq!(cell("verdict"), Some("REGRESS INITIAL"));
    let strata = fs::read_to_string(&set.strata)?;
    assert_eq!(
        strata
            .lines()
            .skip(2)
            .map(|line| line.split('\t').take(3).collect::<Vec<_>>().join(" "))
            .collect::<Vec<_>>(),
        ["r0001 bench 1", "r0001 recent 2"]
    );
    // Nothing is left staged.
    assert_eq!(
        fs::read_dir(&set.ledger)?.count(),
        2,
        "the ledger holds index.tsv and r0001 only"
    );
    // The same entry again is refused.
    assert!(record(&set.record("first", "r0001"), &main, "2026-10-22").is_err());

    // Against the head: the same output is SAME, and records as a change.
    let mut same = set.check("same");
    let intent = set.directory.path().join("intent.tsv");
    fs::write(
        &intent,
        format!(
            "# regress-intent series=regress-v1 base=r0001 pr=#1\n{}\n\
             hv1\tallele\tadded\tHV1\t*\t-\t0\t#1 a synthetic intent\n",
            intents::COLUMNS
        ),
    )?;
    // The intent file is checked against the harness's revision, which a
    // build without git does not record.
    let with_intent = crate::source().revision.is_some();
    if with_intent {
        same.intent = Some(intent.clone());
    }
    let (report, verdict) = check(&same, &main)?;
    assert_eq!(verdict, Verdict::Same);
    if with_intent {
        assert!(report.contains("unused intent hv1"), "{report}");
    }
    set.as_clean("same")?;
    record(&set.record("same", "r0002"), &main, "2026-10-23")?;
    let baseline = fs::read_to_string(&set.baseline)?;
    let row = baseline.lines().last().unwrap_or_default();
    assert!(
        row.starts_with("r0002\t2026-10-23\tchange\tr0001\t"),
        "{row}"
    );
    assert!(row.ends_with("\tREGRESS SAME"), "{row}");
    assert_eq!(set.ledger.join("r0002").join(INTENT).is_file(), with_intent);

    // A truth the manifest does not pin is refused: against the head
    // without --rebase, and for a series' first check.
    let refreshed = set.directory.path().join("refreshed");
    fs::create_dir(&refreshed)?;
    let sql = format!(
        "{} INSERT INTO sample_status VALUES ('b1', 'S3', 'x', 'HOLD', 1, 1, '', 1, 1);",
        crate::truth::tests::FIXTURE
    );
    let mut moved = set.check("moved");
    moved.truth = crate::truth::tests::database(&refreshed, &sql)?;
    assert!(matches!(check(&moved, &main), Err(Error::Regress(_))));
    let mut first_moved = set.check("first-moved");
    first_moved.truth = moved.truth.clone();
    first_moved.ledger = set.directory.path().join("empty-ledger");
    let refused = check(&first_moved, &main);
    assert!(
        matches!(&refused, Err(Error::Regress(message)) if message.contains("pins")),
        "{refused:?}"
    );

    // Another series on this ledger names the series.
    let mut v2 = set.check("v2");
    let text = fs::read_to_string(&set.manifest)?.replacen("version=v1", "version=v2", 1);
    v2.manifest = set.directory.path().join("regress-v2.tsv");
    fs::write(&v2.manifest, text)?;
    let refused = check(&v2, &main);
    assert!(
        matches!(&refused, Err(Error::Regress(message)) if message.contains("another series")),
        "{refused:?}"
    );
    Ok(())
}
