//! `dna-eval approve`: the single pre-registered approval run of triage pass
//! rules on the hold-out (PROP-0005 approval protocol).
//!
//! Everything that could change the verdict is frozen in a committed
//! pre-registration ([`prereg`]) and checked before DNA runs: the file is
//! tracked and unmodified, the harness was built from the frozen revision
//! with a clean tree, the manifest, configuration, and truth copy have their
//! frozen checksums, and the panel's truth digest is unchanged. One DNA pass
//! then feeds every hypothesis of the fixed sequence, which are decided and
//! scored in order: a hypothesis is tested only when its condition on the
//! previous one holds.
//!
//! Per tested hypothesis, each gate bounds the leak rate among passed
//! samples three ways ([`stats`]): the count-based Clopper–Pearson bound, the
//! Korn–Graubard bound of the corpus-weighted rate, and the 95th percentile
//! of a batch-cluster bootstrap. The gate passes only when all three are
//! within the bound: at zero leaks the bootstrap is degenerate and bounds
//! nothing by itself. A marginal endpoint bounds the leaks among the samples
//! a hypothesis passes and another does not; a utility floor applies to the
//! corpus-weighted pass rate. Samples DNA fails never pass; retired samples
//! are left out. Every leak is listed for adjudication, which never changes
//! the verdict.

pub(crate) mod power;
pub(crate) mod prereg;
pub(crate) mod stats;

use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write as _;
use std::path::{Component, Path};
use std::process::Command;

use serde::Serialize;

use crate::compare::Measure;
use crate::error::{Error, Result};
use crate::estimate::{self, Observation};
use crate::report::{Row, Status};
use crate::runner::Threads;
use crate::triage::detect::{Footprints, Outcome, decide};
use crate::triage::score::{ALPHA, upper_bound};
use prereg::{Gate, Hypothesis, Prereg, Requires, UtilityForm};
use stats::{Draw, bootstrap, korn_graubard};

/// One-sided normal quantile for the utility floor's lower bound.
const Z_ONE_SIDED_95: f64 = 1.645;

/// The bounded result of one gate.
#[derive(Debug, Clone, Serialize)]
struct GateResult {
    measure: &'static str,
    bound: f64,
    passed_samples: usize,
    leaks: usize,
    count_upper: f64,
    weighted_rate: f64,
    korn_graubard_upper: f64,
    bootstrap_upper: f64,
    bootstrap_empty_draws: usize,
    pass: bool,
}

#[derive(Debug, Clone, Serialize)]
struct MarginalResult {
    against: String,
    measure: &'static str,
    bound: f64,
    /// The pre-registered lower-bound projection of `newly_passed`.
    projected_n: usize,
    projected_from: prereg::Projection,
    newly_passed: usize,
    leaks: usize,
    /// Leaks with a discordant position the one-strand exception lets one
    /// strand read, and the rest; reported apart, not gating.
    leaks_inside_exception: usize,
    leaks_outside_exception: usize,
    count_upper: f64,
    pass: bool,
}

#[derive(Debug, Clone, Serialize)]
struct UtilityResult {
    form: &'static str,
    minimum: f64,
    pass_rate: f64,
    lower_bound: f64,
    pass: bool,
}

#[derive(Debug, Clone, Serialize)]
struct HypothesisResult {
    name: String,
    policy: String,
    policy_sha256: String,
    tested: bool,
    verdict: &'static str,
    gates: Vec<GateResult>,
    marginal: Option<MarginalResult>,
    utility: Option<UtilityResult>,
    /// Leaks among passed samples under every measure, unweighted: reported,
    /// not gating.
    secondary: BTreeMap<&'static str, (usize, usize)>,
}

#[derive(Debug, Serialize)]
struct Approval {
    prereg_id: String,
    prereg_sha256: String,
    dna_revision: String,
    manifest_body_sha256: String,
    panel_truth_digest: String,
    truth_sha256: String,
    config_sha256: String,
    alpha: f64,
    bootstrap_draws: u64,
    bootstrap_seed: u64,
    rng: &'static str,
    samples: usize,
    active_samples: usize,
    batches: usize,
    hypotheses: Vec<HypothesisResult>,
}

/// Runs the pre-registered approval and writes its outputs; returns a
/// one-line summary.
pub(crate) fn approve(prereg_path: &Path, output: &Path, threads: usize) -> Result<String> {
    let threads = Threads::new(threads)?;
    let bytes = crate::read(prereg_path, "pre-registration", 1 << 20)?;
    let prereg = prereg::parse(&bytes).map_err(Error::Prereg)?;
    committed(prereg_path)?;
    let source = crate::source();
    if source.revision.as_deref() != Some(prereg.inputs.dna_revision.as_str())
        || source.dirty != Some(false)
    {
        return Err(mismatch(
            "harness build revision",
            &prereg.inputs.dna_revision,
            &format!(
                "{} (dirty: {:?})",
                source.revision.as_deref().unwrap_or("unknown"),
                source.dirty
            ),
        ));
    }
    under_results(output)?;
    if output.exists() {
        return Err(Error::OutputNotEmpty {
            path: output.to_path_buf(),
        });
    }
    let inputs = &prereg.inputs;
    let manifest = crate::manifest::parse(&crate::read(
        &inputs.manifest,
        "manifest",
        crate::MAX_MANIFEST_BYTES,
    )?)?;
    check(
        "manifest body sha256",
        &inputs.manifest_body_sha256,
        &manifest.body_sha256,
    )?;
    check(
        "config sha256",
        &inputs.config_sha256,
        &crate::runner::sha256_file(&inputs.config).map_err(|source| Error::Read {
            kind: "configuration",
            path: inputs.config.clone(),
            source,
        })?,
    )?;
    let wanted = manifest
        .samples
        .iter()
        .map(|sample| sample.id.as_str())
        .collect::<BTreeSet<_>>();
    let truth = crate::load_truth(&inputs.truth, &wanted)?;
    check("truth sha256", &inputs.truth_sha256, &truth.sha256)?;
    let truth_check = crate::truth::TruthCheck::new(&manifest, &truth.database);
    check(
        "panel truth digest",
        &inputs.panel_truth_digest,
        truth_check.recorded.as_deref().unwrap_or("none"),
    )?;
    check(
        "current panel truth digest",
        &inputs.panel_truth_digest,
        &truth_check.current,
    )?;
    let reference = dna_kernel::reference::load(
        &inputs.reference,
        dna_kernel::model::reference::ReferenceTopology::Linear,
    )
    .map_err(|error| Error::Reference {
        path: inputs.reference.clone(),
        message: error.to_string(),
    })?;

    let outcomes = crate::runner::run(
        &manifest.samples,
        (&inputs.reference, &reference.sequence_sha256),
        &inputs.config,
        threads,
        true,
        None,
        &mut Vec::new(),
    )?;
    let rows = crate::report::rows(
        &manifest.samples,
        &outcomes,
        &truth.database.samples,
        &reference.sequence,
    );
    let (approval, leaks) = score(&prereg, &rows, &reference.sequence);

    crate::prepare_output(output)?;
    let mut json = serde_json::to_string_pretty(&approval).map_err(Error::Serialize)?;
    json.push('\n');
    crate::write(output, "approval.json", &json)?;
    crate::write(output, "approval.md", &markdown(&approval))?;
    crate::write(output, "leaks.tsv", &leaks)?;
    crate::write(
        output,
        "prereg.toml",
        std::str::from_utf8(&bytes).unwrap_or_default(),
    )?;
    let verdicts = approval
        .hypotheses
        .iter()
        .map(|hypothesis| format!("{} {}", hypothesis.name, hypothesis.verdict))
        .collect::<Vec<_>>()
        .join(", ");
    Ok(format!(
        "dna-eval approve: {} ({}); {verdicts}; outputs in {}",
        approval.prereg_id,
        approval.prereg_sha256,
        output.display()
    ))
}

/// Decides and scores every hypothesis of the fixed sequence; returns the
/// approval and the leak list for adjudication.
fn score(prereg: &Prereg, rows: &[Row<'_>], reference: &str) -> (Approval, String) {
    let active = rows
        .iter()
        .filter(|row| !row.sample.retired)
        .collect::<Vec<_>>();
    let footprints = Footprints::learn(active.iter().filter_map(|row| row.outcome.triage.as_ref()))
        .with_reference(reference);
    let batches = active
        .iter()
        .map(|row| row.sample.batch.as_str())
        .collect::<BTreeSet<_>>()
        .into_iter()
        .enumerate()
        .map(|(index, batch)| (batch, index))
        .collect::<BTreeMap<_, _>>();
    // Per hypothesis, whether each active sample passes.
    let mut pass_sets = BTreeMap::<&str, Vec<bool>>::new();
    let mut results = Vec::new();
    let mut leaks = String::from("sample_id\thypothesis\tmeasure\tstratum\tbatch\tverdict\n");
    let mut previous: Option<&HypothesisResult> = None;
    for hypothesis in &prereg.hypotheses {
        let passed = active
            .iter()
            .map(|row| {
                row.outcome.result.is_ok()
                    && row.outcome.triage.as_ref().is_some_and(|evidence| {
                        decide(evidence, &hypothesis.policy, &footprints).outcome == Outcome::Pass
                    })
            })
            .collect::<Vec<_>>();
        let tested = tested(
            hypothesis.requires,
            previous.map(|previous| previous.verdict),
        );
        let result = if tested {
            let result = test(
                prereg,
                hypothesis,
                (&active, &passed),
                &pass_sets,
                &batches,
                &footprints,
            );
            for gate in &prereg.gates {
                for (row, _) in active
                    .iter()
                    .zip(&passed)
                    .filter(|(row, passed)| **passed && leaked(row, gate.measure))
                {
                    let _ = writeln!(
                        leaks,
                        "{}\t{}\t{}\t{}\t{}\t",
                        row.sample.id,
                        hypothesis.name,
                        gate.measure.label(),
                        row.sample.stratum,
                        row.sample.batch
                    );
                }
            }
            if let Some(marginal) = &hypothesis.marginal
                && let Some(against) = pass_sets.get(marginal.against.as_str())
            {
                for ((row, passed), before) in active.iter().zip(&passed).zip(against) {
                    if *passed && !*before && leaked(row, marginal.measure) {
                        let _ = writeln!(
                            leaks,
                            "{}\t{}\tmarginal_{}\t{}\t{}\t",
                            row.sample.id,
                            hypothesis.name,
                            marginal.measure.label(),
                            row.sample.stratum,
                            row.sample.batch
                        );
                    }
                }
            }
            result
        } else {
            HypothesisResult {
                name: hypothesis.name.clone(),
                policy: hypothesis.policy.parameters(),
                policy_sha256: hypothesis.policy_sha256.clone(),
                tested: false,
                verdict: "not_tested",
                gates: Vec::new(),
                marginal: None,
                utility: None,
                secondary: BTreeMap::new(),
            }
        };
        pass_sets.insert(hypothesis.name.as_str(), passed);
        results.push(result);
        previous = results.last();
    }
    let approval = Approval {
        prereg_id: prereg.id.clone(),
        prereg_sha256: prereg.sha256.clone(),
        dna_revision: prereg.inputs.dna_revision.clone(),
        manifest_body_sha256: prereg.inputs.manifest_body_sha256.clone(),
        panel_truth_digest: prereg.inputs.panel_truth_digest.clone(),
        truth_sha256: prereg.inputs.truth_sha256.clone(),
        config_sha256: prereg.inputs.config_sha256.clone(),
        alpha: ALPHA,
        bootstrap_draws: prereg.bootstrap_draws,
        bootstrap_seed: prereg.bootstrap_seed,
        rng: "splitmix64(seed ^ draw * 0x9E3779B97F4A7C15)",
        samples: rows.len(),
        active_samples: active.len(),
        batches: batches.len(),
        hypotheses: results,
    };
    (approval, leaks)
}

/// Whether a hypothesis of the fixed sequence is tested, given the verdict
/// of the previous one (`not_tested` when it was skipped): the sequence
/// stops at the first failure.
fn tested(requires: Requires, previous: Option<&str>) -> bool {
    match requires {
        Requires::First => true,
        Requires::PreviousPass => previous == Some("pass"),
    }
}

/// Whether a leaked sample has a discordant position under `measure` that
/// one strand alone reads and the hypothesis's one-strand exception lets
/// pass.
fn inside_exception(
    row: &Row<'_>,
    measure: Measure,
    hypothesis: &Hypothesis,
    footprints: &Footprints,
) -> bool {
    let (Ok(comparison), Some(evidence)) = (row.comparison(measure), row.outcome.triage.as_ref())
    else {
        return false;
    };
    comparison.positions.iter().any(|position| {
        crate::compare::anchor(position)
            .and_then(|anchor| Some((anchor, evidence.locate(anchor)?)))
            .is_some_and(|(anchor, (region, offset))| {
                let one_strand = evidence.regions[region]
                    .positions
                    .get(offset)
                    .is_some_and(|position| position.strands() == 1);
                let coordinate = usize::try_from(anchor).unwrap_or(0);
                one_strand
                    && !hypothesis
                        .policy
                        .one_strand
                        .applies_at(coordinate, &footprints.repeats)
            })
    })
}

/// Whether a passed sample is a leak under `measure`: anything but
/// concordant.
fn leaked(row: &Row<'_>, measure: Measure) -> bool {
    row.status(measure) != Status::Concordant
}

fn test(
    prereg: &Prereg,
    hypothesis: &Hypothesis,
    (active, passed): (&[&Row<'_>], &[bool]),
    pass_sets: &BTreeMap<&str, Vec<bool>>,
    batches: &BTreeMap<&str, usize>,
    footprints: &Footprints,
) -> HypothesisResult {
    let gates = prereg
        .gates
        .iter()
        .map(|gate| gate_result(*gate, prereg, active, passed, batches))
        .collect::<Vec<_>>();
    let marginal = hypothesis.marginal.as_ref().map(|marginal| {
        let against = pass_sets.get(marginal.against.as_str());
        let newly = active
            .iter()
            .zip(passed)
            .enumerate()
            .filter(|(index, (_, passed))| {
                **passed
                    && !against
                        .and_then(|other| other.get(*index))
                        .copied()
                        .unwrap_or(false)
            })
            .map(|(_, (row, _))| *row)
            .collect::<Vec<_>>();
        let leaking = newly
            .iter()
            .filter(|row| leaked(row, marginal.measure))
            .collect::<Vec<_>>();
        let inside = leaking
            .iter()
            .filter(|row| inside_exception(row, marginal.measure, hypothesis, footprints))
            .count();
        let count_upper = upper_bound(leaking.len(), newly.len());
        MarginalResult {
            against: marginal.against.clone(),
            measure: marginal.measure.label(),
            bound: marginal.upper_bound,
            projected_n: marginal.projected_n,
            projected_from: marginal.projected_from.clone(),
            newly_passed: newly.len(),
            leaks: leaking.len(),
            leaks_inside_exception: inside,
            leaks_outside_exception: leaking.len() - inside,
            count_upper,
            pass: count_upper <= marginal.upper_bound,
        }
    });
    let pass_rate = estimate::estimate(
        &active
            .iter()
            .zip(passed)
            .map(|(row, passed)| Observation {
                stratum: &row.sample.stratum,
                weight: row.sample.weight_corpus.value,
                event: *passed,
            })
            .collect::<Vec<_>>(),
    );
    let utility = pass_rate.map(|estimate| {
        let lower_bound = estimate.rate - Z_ONE_SIDED_95 * estimate.variance.sqrt();
        let (form, value) = match prereg.utility.form {
            UtilityForm::Point => ("point", estimate.rate),
            UtilityForm::LowerBound => ("lower_bound", lower_bound),
        };
        UtilityResult {
            form,
            minimum: prereg.utility.minimum_pass_rate,
            pass_rate: estimate.rate,
            lower_bound,
            pass: value >= prereg.utility.minimum_pass_rate,
        }
    });
    let secondary = Measure::ALL
        .into_iter()
        .map(|measure| {
            let passed_rows = active.iter().zip(passed).filter(|(_, passed)| **passed);
            let (count, leaks) = passed_rows.fold((0, 0), |(count, leaks), (row, _)| {
                (count + 1, leaks + usize::from(leaked(row, measure)))
            });
            (measure.label(), (leaks, count))
        })
        .collect();
    let verdict = if gates.iter().all(|gate| gate.pass)
        && marginal.as_ref().is_none_or(|marginal| marginal.pass)
        && utility.as_ref().is_some_and(|utility| utility.pass)
    {
        "pass"
    } else {
        "fail"
    };
    HypothesisResult {
        name: hypothesis.name.clone(),
        policy: hypothesis.policy.parameters(),
        policy_sha256: hypothesis.policy_sha256.clone(),
        tested: true,
        verdict,
        gates,
        marginal,
        utility,
        secondary,
    }
}

fn gate_result(
    gate: Gate,
    prereg: &Prereg,
    active: &[&Row<'_>],
    passed: &[bool],
    batches: &BTreeMap<&str, usize>,
) -> GateResult {
    let passed_rows = active
        .iter()
        .zip(passed)
        .filter(|(_, passed)| **passed)
        .map(|(row, _)| *row)
        .collect::<Vec<_>>();
    let leaks = passed_rows
        .iter()
        .filter(|row| leaked(row, gate.measure))
        .count();
    let count_upper = upper_bound(leaks, passed_rows.len());
    let weighted = estimate::estimate(
        &passed_rows
            .iter()
            .map(|row| Observation {
                stratum: &row.sample.stratum,
                weight: row.sample.weight_corpus.value,
                event: leaked(row, gate.measure),
            })
            .collect::<Vec<_>>(),
    );
    let (weighted_rate, korn_graubard_upper) = weighted.map_or((0.0, 1.0), |estimate| {
        (
            estimate.rate,
            korn_graubard(estimate.rate, estimate.variance, estimate.samples),
        )
    });
    let draws = active
        .iter()
        .zip(passed)
        .map(|(row, passed)| Draw {
            batch: batches.get(row.sample.batch.as_str()).copied().unwrap_or(0),
            weight: row.sample.weight_corpus.value,
            passed: *passed,
            leaked: leaked(row, gate.measure),
        })
        .collect::<Vec<_>>();
    let cluster = bootstrap(
        &draws,
        batches.len(),
        prereg.bootstrap_draws,
        prereg.bootstrap_seed,
    );
    GateResult {
        measure: gate.measure.label(),
        bound: gate.upper_bound,
        passed_samples: passed_rows.len(),
        leaks,
        count_upper,
        weighted_rate,
        korn_graubard_upper,
        bootstrap_upper: cluster.upper,
        bootstrap_empty_draws: cluster.empty_draws,
        pass: count_upper <= gate.upper_bound
            && korn_graubard_upper <= gate.upper_bound
            && cluster.upper <= gate.upper_bound,
    }
}

/// `approval.md`: the verdicts in aggregate; no sample identifier.
fn markdown(approval: &Approval) -> String {
    let mut text = format!(
        "# Approval run `{}`\n\n\
         Pre-registration sha256 `{}`, DNA revision `{}`, {} samples ({} active) in {} batches.\n\n\
         Each gate passes only when the count-based Clopper–Pearson bound, the Korn–Graubard \
         bound of the corpus-weighted rate, and the batch-bootstrap 95th percentile are all \
         within the bound (one-sided, alpha {}). At zero leaks the bootstrap bound is 0 and is \
         no evidence by itself.\n",
        approval.prereg_id,
        approval.prereg_sha256,
        approval.dna_revision,
        approval.samples,
        approval.active_samples,
        approval.batches,
        approval.alpha
    );
    for hypothesis in &approval.hypotheses {
        let _ = write!(
            text,
            "\n## {}: {}\n\nPolicy `{}` (sha256 `{}`).\n",
            hypothesis.name, hypothesis.verdict, hypothesis.policy, hypothesis.policy_sha256
        );
        if !hypothesis.tested {
            continue;
        }
        text.push_str(
            "\n| Gate | Bound | Passed | Leaks | Count UB | Weighted rate | KG UB | Bootstrap UB | Pass |\n\
             |---|---:|---:|---:|---:|---:|---:|---:|---|\n",
        );
        for gate in &hypothesis.gates {
            let _ = writeln!(
                text,
                "| {} | {:.4} | {} | {} | {:.4} | {:.4} | {:.4} | {:.4} | {} |",
                gate.measure,
                gate.bound,
                gate.passed_samples,
                gate.leaks,
                gate.count_upper,
                gate.weighted_rate,
                gate.korn_graubard_upper,
                gate.bootstrap_upper,
                gate.pass
            );
        }
        if let Some(marginal) = &hypothesis.marginal {
            let _ = writeln!(
                text,
                "\nMarginal against {} ({}): {} newly passed (projected {}), {} leaks ({} at a \
                 position the one-strand exception lets one strand read, {} elsewhere), UB {:.4} \
                 (bound {:.4}): {}.",
                marginal.against,
                marginal.measure,
                marginal.newly_passed,
                marginal.projected_n,
                marginal.leaks,
                marginal.leaks_inside_exception,
                marginal.leaks_outside_exception,
                marginal.count_upper,
                marginal.bound,
                marginal.pass
            );
        }
        if let Some(utility) = &hypothesis.utility {
            let _ = writeln!(
                text,
                "\nUtility ({} ≥ {:.4}): pass rate {:.4}, lower bound {:.4}: {}.",
                utility.form, utility.minimum, utility.pass_rate, utility.lower_bound, utility.pass
            );
        }
        text.push_str("\n| Measure | Leaks / passed |\n|---|---:|\n");
        for (measure, (leaks, passed)) in &hypothesis.secondary {
            let _ = writeln!(text, "| {measure} | {leaks} / {passed} |");
        }
    }
    text
}

fn check(what: &'static str, expected: &str, found: &str) -> Result<()> {
    if expected == found {
        Ok(())
    } else {
        Err(mismatch(what, expected, found))
    }
}

fn mismatch(what: &'static str, expected: &str, found: &str) -> Error {
    Error::Frozen {
        what,
        expected: expected.to_owned(),
        found: found.to_owned(),
    }
}

/// The output directory must lie under a `results` directory, which the
/// repository ignores.
fn under_results(output: &Path) -> Result<()> {
    if output
        .components()
        .any(|component| component == Component::Normal("results".as_ref()))
    {
        Ok(())
    } else {
        Err(Error::OutputOutsideResults {
            path: output.to_path_buf(),
        })
    }
}

/// The pre-registration must be tracked by git and unmodified.
fn committed(path: &Path) -> Result<()> {
    let directory = path
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));
    let name = path.file_name().map_or(path, Path::new);
    let git = |arguments: &[&str]| {
        Command::new("git")
            .arg("-C")
            .arg(directory)
            .args(arguments)
            .arg(name)
            .output()
            .is_ok_and(|output| output.status.success())
    };
    if git(&["ls-files", "--error-unmatch", "--"]) && git(&["diff", "--quiet", "HEAD", "--"]) {
        Ok(())
    } else {
        Err(Error::PreregNotCommitted {
            path: path.to_path_buf(),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::stats::beta_upper;
    use super::*;

    #[test]
    fn outputs_must_lie_under_results() {
        assert!(under_results(Path::new("results/approval/run-1")).is_ok());
        assert!(under_results(Path::new("/tmp/run-1")).is_err());
    }

    #[test]
    fn the_fixed_sequence_tests_in_order() {
        assert!(tested(Requires::First, None));
        assert!(tested(Requires::PreviousPass, Some("pass")));
        assert!(!tested(Requires::PreviousPass, Some("fail")));
        // A skipped hypothesis skips every later one.
        assert!(!tested(Requires::PreviousPass, Some("not_tested")));
    }

    #[test]
    fn frozen_values_must_match() {
        assert!(check("config sha256", "ab", "ab").is_ok());
        assert!(matches!(
            check("config sha256", "ab", "cd"),
            Err(Error::Frozen {
                what: "config sha256",
                ..
            })
        ));
    }

    #[test]
    fn beta_and_count_bounds_agree_at_zero() {
        assert!((beta_upper(0.0, 1600.0) - upper_bound(0, 1600)).abs() < 1e-6);
    }
}
