//! The two-strand run-length rule (SRS-CONS-013) against the research
//! reference's vectors: twelve HVS-II runs from development panel v5,
//! with every view's top, next, `p`, and standard error, the strands'
//! spread, and the outcome.

use std::path::Path;

use dna_consensus::length::{Outcome, RuleParams, Run, View, decide, derived, direct, in_phase};
use dna_kernel::read_evidence::{LengthMarginal, Margin};
use serde_json::Value;

fn number(value: &Value) -> Option<f64> {
    value.as_f64()
}

fn marginal(value: &Value) -> LengthMarginal {
    let lengths = value["lengths"]
        .as_array()
        .into_iter()
        .flatten()
        .filter_map(Value::as_u64)
        .map(|length| usize::try_from(length).unwrap_or_default())
        .collect::<Vec<_>>();
    let numbers = |key: &str| {
        value[key]
            .as_array()
            .into_iter()
            .flatten()
            .filter_map(Value::as_f64)
            .collect::<Vec<_>>()
    };
    LengthMarginal {
        dominant: lengths.first().copied(),
        next: None,
        lengths,
        shares: numbers("shares"),
        errors: numbers("errors"),
        z: number(&value["z"]),
        z_adjusted: number(&value["z_adjusted"]),
        margin: Margin::Clear,
    }
}

fn close(actual: f64, expected: f64, what: &str, case: &str) {
    assert!(
        (actual - expected).abs() <= 1e-9 * expected.abs().max(1.0),
        "{case}: {what} {actual} != {expected}"
    );
}

// SRS-CONS-013
#[test]
fn decides_the_reference_runs() -> Result<(), Box<dyn std::error::Error>> {
    let directory = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/vectors/two_strand");
    let mut files = std::fs::read_dir(&directory)?
        .filter_map(|entry| entry.ok().map(|entry| entry.path()))
        .collect::<Vec<_>>();
    files.sort();
    assert_eq!(files.len(), 12);
    for file in files {
        let vector: Value = serde_json::from_str(&std::fs::read_to_string(&file)?)?;
        let case = vector["case"].as_str().unwrap_or_default().to_owned();
        let run = if vector["run"] == "run1" {
            Run::Before
        } else {
            Run::After
        };
        let tau = &vector["parameters"]["tau"];
        let params = RuleParams {
            confidence: number(&vector["parameters"]["k"]).ok_or("k")?,
            forward: [
                number(&tau["F/run1"]).ok_or("tau")?,
                number(&tau["F/run2"]).ok_or("tau")?,
            ],
            reverse: [
                number(&tau["R/run1"]).ok_or("tau")?,
                number(&tau["R/run2"]).ok_or("tau")?,
            ],
        };
        // Each read's view of the run, F reads before R reads, as the
        // reference lists them.
        let mut views: Vec<View> = Vec::new();
        for read in vector["reads"].as_array().into_iter().flatten() {
            let forward = read["strand"] == "F";
            let key = if run == Run::Before { "run1" } else { "run2" };
            if read["source"] == "in_phase" {
                let length = read["lengths"][key].as_u64().ok_or("in-phase length")?;
                views.push(in_phase(forward, usize::try_from(length)?));
                continue;
            }
            let expected_mark = if forward { "T" } else { "A" };
            let mark = read["mark"].as_str();
            if mark.is_some_and(|mark| mark != expected_mark) {
                continue;
            }
            let entered = if forward || read["opens_at"] == "run1" {
                Run::Before
            } else {
                Run::After
            };
            let first = marginal(&read["first"]);
            let view = if run == entered {
                direct(forward, &first)
            } else if mark == Some(expected_mark)
                && (forward && entered == Run::Before || !forward && entered == Run::After)
                && direct(forward, &first).is_some()
            {
                derived(forward, &first, &marginal(&read["total"]))
            } else {
                None
            };
            views.extend(view);
        }
        views.sort_by_key(|view| !view.forward);
        let expected_views = vector["intermediate"]["views"].as_array();
        if let Some(expected) = expected_views {
            assert_eq!(views.len(), expected.len(), "{case}: views");
            for (view, wanted) in views.iter().zip(expected) {
                assert_eq!(
                    view.top as u64,
                    wanted["top"].as_u64().unwrap_or_default(),
                    "{case}: top"
                );
                assert_eq!(
                    view.next.map(|next| next as u64),
                    wanted["next"].as_u64(),
                    "{case}: next"
                );
                close(view.p, number(&wanted["p"]).ok_or("p")?, "p", &case);
                match (view.se, number(&wanted["se"])) {
                    (Some(se), Some(wanted)) => close(se, wanted, "se", &case),
                    (None, None) => {}
                    other => panic!("{case}: se {other:?}"),
                }
            }
        }
        let (outcome, spread) = decide(&views, run, &params);
        if let Some(d) = number(&vector["intermediate"]["d"]) {
            close(spread.ok_or("d")?, d, "d", &case);
        }
        let wanted = &vector["outcome"];
        let expected = match wanted["kind"].as_str() {
            Some("decided") => {
                Outcome::Decided(usize::try_from(wanted["length"].as_u64().ok_or("length")?)?)
            }
            Some("near_equal") => {
                Outcome::NearEqual(usize::try_from(wanted["pair"][0].as_u64().ok_or("pair")?)?)
            }
            Some("strands_differ") => Outcome::StrandsDiffer,
            _ => Outcome::OneStrand,
        };
        assert_eq!(outcome, expected, "{case}: outcome");
    }
    Ok(())
}
