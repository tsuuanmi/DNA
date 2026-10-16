//! The pre-registration file of an approval run (PROP-0005 approval
//! protocol): what is frozen before the hold-out is read.
//!
//! Every bound, the utility form, and every input identity are written down
//! before the run. A value still written as `"pending"` stops the run, so a
//! file the user has not finished can never be run.

use std::collections::BTreeSet;
use std::path::PathBuf;

use serde::Deserialize;
use sha2::{Digest, Sha256};

use crate::compare::Measure;
use crate::error::PreregError;
use crate::evidence::Cause;
use crate::triage::detect::{Policy, Proxy};
use crate::triage::minor::{Kind, Rule};

/// A number, or the word `"pending"`.
#[derive(Debug, Clone, Deserialize)]
#[serde(untagged)]
enum Pinned<T> {
    Value(T),
    Text(String),
}

impl<T> Pinned<T> {
    fn value(self, field: &str) -> Result<T, PreregError> {
        match self {
            Self::Value(value) => Ok(value),
            Self::Text(text) if text == "pending" => Err(PreregError::Pending(field.to_owned())),
            Self::Text(text) => Err(invalid(field, &format!("{text:?} is not a number"))),
        }
    }
}

fn invalid(field: &str, message: &str) -> PreregError {
    PreregError::Invalid {
        field: field.to_owned(),
        message: message.to_owned(),
    }
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct RawPrereg {
    id: String,
    inputs: RawInputs,
    analysis: RawAnalysis,
    utility: RawUtility,
    gates: Vec<RawGate>,
    hypotheses: Vec<RawHypothesis>,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct RawInputs {
    manifest: PathBuf,
    manifest_body_sha256: String,
    panel_truth_digest: String,
    truth: PathBuf,
    truth_sha256: String,
    config: PathBuf,
    config_sha256: String,
    reference: PathBuf,
    dna_revision: String,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct RawAnalysis {
    bootstrap_draws: u64,
    bootstrap_seed: u64,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct RawUtility {
    form: String,
    minimum_pass_rate: Pinned<f64>,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct RawGate {
    measure: String,
    upper_bound: Pinned<f64>,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct RawHypothesis {
    name: String,
    requires: String,
    policy: RawPolicy,
    #[serde(default)]
    marginal: Option<RawMarginal>,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct RawPolicy {
    short_calls: usize,
    causes: Vec<String>,
    one_strand: bool,
    proxy: String,
    minor: String,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct RawMarginal {
    against: String,
    measure: String,
    upper_bound: Pinned<f64>,
}

/// A pre-registration ready to run.
#[derive(Debug, Clone)]
pub(crate) struct Prereg {
    pub(crate) id: String,
    /// SHA-256 of the file as read.
    pub(crate) sha256: String,
    pub(crate) inputs: Inputs,
    pub(crate) bootstrap_draws: u64,
    pub(crate) bootstrap_seed: u64,
    pub(crate) utility: Utility,
    pub(crate) gates: Vec<Gate>,
    pub(crate) hypotheses: Vec<Hypothesis>,
}

/// The inputs an approval run must find exactly as frozen.
#[derive(Debug, Clone)]
pub(crate) struct Inputs {
    pub(crate) manifest: PathBuf,
    pub(crate) manifest_body_sha256: String,
    pub(crate) panel_truth_digest: String,
    pub(crate) truth: PathBuf,
    pub(crate) truth_sha256: String,
    pub(crate) config: PathBuf,
    pub(crate) config_sha256: String,
    pub(crate) reference: PathBuf,
    pub(crate) dna_revision: String,
}

/// Whether the utility floor applies to the pass-rate point estimate or to
/// its one-sided 95 % lower bound.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum UtilityForm {
    Point,
    LowerBound,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Utility {
    pub(crate) form: UtilityForm,
    pub(crate) minimum_pass_rate: f64,
}

/// A safety endpoint every tested hypothesis must meet: the leak rate among
/// passed samples under `measure`, bounded three ways.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Gate {
    pub(crate) measure: Measure,
    pub(crate) upper_bound: f64,
}

/// When a hypothesis of the fixed sequence is tested. The sequence stops at
/// the first failure: testing on after a failure would spend alpha again
/// and lose the family-wise error control of a fixed-sequence procedure.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Requires {
    /// Always: the first hypothesis.
    First,
    /// Only if the previous hypothesis was tested and passed.
    PreviousPass,
}

/// The risk a hypothesis adds over another: leaks among the samples it
/// passes and `against` does not.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Marginal {
    pub(crate) against: String,
    pub(crate) measure: Measure,
    pub(crate) upper_bound: f64,
}

#[derive(Debug, Clone)]
pub(crate) struct Hypothesis {
    pub(crate) name: String,
    pub(crate) requires: Requires,
    pub(crate) policy: Policy,
    /// SHA-256 of the policy's full definition.
    pub(crate) policy_sha256: String,
    pub(crate) marginal: Option<Marginal>,
}

/// Parses and validates a pre-registration.
pub(crate) fn parse(bytes: &[u8]) -> Result<Prereg, PreregError> {
    let text = std::str::from_utf8(bytes)
        .map_err(|_| PreregError::Format("the file is not UTF-8".to_owned()))?;
    let raw: RawPrereg =
        toml::from_str(text).map_err(|error| PreregError::Format(error.to_string()))?;
    if raw.analysis.bootstrap_draws == 0 {
        return Err(invalid("analysis.bootstrap_draws", "must be at least 1"));
    }
    let form = match raw.utility.form.as_str() {
        "point" => UtilityForm::Point,
        "lower_bound" => UtilityForm::LowerBound,
        "pending" => return Err(PreregError::Pending("utility.form".to_owned())),
        other => {
            return Err(invalid(
                "utility.form",
                &format!("{other:?} is not point or lower_bound"),
            ));
        }
    };
    let utility = Utility {
        form,
        minimum_pass_rate: fraction(raw.utility.minimum_pass_rate, "utility.minimum_pass_rate")?,
    };
    if raw.gates.is_empty() {
        return Err(invalid("gates", "at least one gate is required"));
    }
    let gates = raw
        .gates
        .into_iter()
        .enumerate()
        .map(|(index, gate)| {
            Ok(Gate {
                measure: measure(&gate.measure, &format!("gates[{index}].measure"))?,
                upper_bound: fraction(gate.upper_bound, &format!("gates[{index}].upper_bound"))?,
            })
        })
        .collect::<Result<Vec<_>, PreregError>>()?;
    if raw.hypotheses.is_empty() {
        return Err(invalid("hypotheses", "at least one hypothesis is required"));
    }
    let mut names = BTreeSet::new();
    let mut hypotheses = Vec::new();
    for (index, hypothesis) in raw.hypotheses.into_iter().enumerate() {
        let field = |name: &str| format!("hypotheses[{index}].{name}");
        if !names.insert(hypothesis.name.clone()) {
            return Err(invalid(&field("name"), "names must be unique"));
        }
        let requires = match (index, hypothesis.requires.as_str()) {
            (0, "first") => Requires::First,
            (0, _) => {
                return Err(invalid(
                    &field("requires"),
                    "the first hypothesis requires \"first\"",
                ));
            }
            (_, "previous_pass") => Requires::PreviousPass,
            (_, other) => {
                return Err(invalid(
                    &field("requires"),
                    &format!(
                        "{other:?}: a later hypothesis requires \"previous_pass\"; the fixed \
                         sequence stops at the first failure"
                    ),
                ));
            }
        };
        let policy = policy(&hypothesis.name, &hypothesis.policy, &field("policy"))?;
        let marginal = hypothesis
            .marginal
            .map(|marginal| {
                if !names.contains(&marginal.against) || marginal.against == hypothesis.name {
                    return Err(invalid(
                        &field("marginal.against"),
                        "must name an earlier hypothesis",
                    ));
                }
                Ok(Marginal {
                    against: marginal.against,
                    measure: measure(&marginal.measure, &field("marginal.measure"))?,
                    upper_bound: fraction(marginal.upper_bound, &field("marginal.upper_bound"))?,
                })
            })
            .transpose()?;
        hypotheses.push(Hypothesis {
            name: hypothesis.name,
            requires,
            policy_sha256: hex(&Sha256::digest(definition(&policy).as_bytes())),
            policy,
            marginal,
        });
    }
    let inputs = raw.inputs;
    for (name, value) in [
        ("inputs.manifest_body_sha256", &inputs.manifest_body_sha256),
        ("inputs.panel_truth_digest", &inputs.panel_truth_digest),
        ("inputs.truth_sha256", &inputs.truth_sha256),
        ("inputs.config_sha256", &inputs.config_sha256),
        ("inputs.dna_revision", &inputs.dna_revision),
    ] {
        if value == "pending" {
            return Err(PreregError::Pending(name.to_owned()));
        }
        if value.is_empty() || !value.bytes().all(|byte| byte.is_ascii_hexdigit()) {
            return Err(invalid(name, "must be lowercase hexadecimal"));
        }
    }
    Ok(Prereg {
        id: raw.id,
        sha256: hex(&Sha256::digest(bytes)),
        inputs: Inputs {
            manifest: inputs.manifest,
            manifest_body_sha256: inputs.manifest_body_sha256,
            panel_truth_digest: inputs.panel_truth_digest,
            truth: inputs.truth,
            truth_sha256: inputs.truth_sha256,
            config: inputs.config,
            config_sha256: inputs.config_sha256,
            reference: inputs.reference,
            dna_revision: inputs.dna_revision,
        },
        bootstrap_draws: raw.analysis.bootstrap_draws,
        bootstrap_seed: raw.analysis.bootstrap_seed,
        utility,
        gates,
        hypotheses,
    })
}

/// The canonical text of a policy's full definition.
pub(crate) fn definition(policy: &Policy) -> String {
    format!("{}\n{}", policy.name, policy.parameters())
}

fn fraction(value: Pinned<f64>, field: &str) -> Result<f64, PreregError> {
    let value = value.value(field)?;
    if !(value.is_finite() && (0.0..=1.0).contains(&value)) {
        return Err(invalid(field, "must lie in 0..=1"));
    }
    Ok(value)
}

fn measure(label: &str, field: &str) -> Result<Measure, PreregError> {
    Measure::ALL
        .into_iter()
        .find(|measure| measure.label() == label)
        .ok_or_else(|| invalid(field, &format!("{label:?} is not a measure")))
}

fn policy(name: &str, raw: &RawPolicy, field: &str) -> Result<Policy, PreregError> {
    let causes = raw
        .causes
        .iter()
        .map(|label| {
            Cause::ALL
                .into_iter()
                .find(|cause| cause.label() == label)
                .ok_or_else(|| {
                    invalid(
                        &format!("{field}.causes"),
                        &format!("{label:?} is not a cause"),
                    )
                })
        })
        .collect::<Result<BTreeSet<_>, PreregError>>()?;
    // A passed region must never hold an abstention (PROP-0005 approval
    // protocol): every policy keeps the `n` cause.
    if !causes.contains(&Cause::N) {
        return Err(invalid(&format!("{field}.causes"), "must contain \"n\""));
    }
    let proxy = match raw.proxy.as_str() {
        "primer" => Proxy::Primer,
        "any" => Proxy::Any,
        other => {
            return Err(invalid(
                &format!("{field}.proxy"),
                &format!("{other:?} is not primer or any"),
            ));
        }
    };
    let minor = if raw.minor == "off" {
        None
    } else {
        let (kind, threshold) = raw.minor.split_once('@').ok_or_else(|| {
            invalid(
                &format!("{field}.minor"),
                "must be <kind>@<thousandths> or off",
            )
        })?;
        let kind = Kind::ALL
            .into_iter()
            .find(|candidate| candidate.label() == kind)
            .ok_or_else(|| {
                invalid(
                    &format!("{field}.minor"),
                    &format!("{kind:?} is not a kind"),
                )
            })?;
        let threshold = threshold
            .parse::<u16>()
            .map_err(|_| invalid(&format!("{field}.minor"), "the threshold is not a number"))?;
        Some(Rule { kind, threshold })
    };
    Ok(Policy {
        name: name.to_owned(),
        short_calls: raw.short_calls,
        causes,
        one_strand: raw.one_strand,
        proxy,
        minor,
    })
}

fn hex(bytes: &[u8]) -> String {
    use std::fmt::Write as _;
    bytes.iter().fold(String::new(), |mut text, byte| {
        let _ = write!(text, "{byte:02x}");
        text
    })
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;

    pub(crate) const VALID: &str = r#"
id = "triage-v1"

[inputs]
manifest = "results/test-panel/validation-v3.tsv"
manifest_body_sha256 = "aa"
panel_truth_digest = "bb"
truth = "results/approval/statistics.sqlite"
truth_sha256 = "cc"
config = "config/dna.toml"
config_sha256 = "dd"
reference = "references/rCRS.fasta"
dna_revision = "ee"

[analysis]
bootstrap_draws = 100
bootstrap_seed = 7

[utility]
form = "point"
minimum_pass_rate = 0.10

[[gates]]
measure = "haplotype_iupac_compatible"
upper_bound = 0.01

[[hypotheses]]
name = "H-base"
requires = "first"
policy = { short_calls = 100, causes = ["n", "contested", "masked", "phase_loss", "reference_frame"], one_strand = true, proxy = "primer", minor = "upstream@100" }

[[hypotheses]]
name = "H-os"
requires = "previous_pass"
policy = { short_calls = 100, causes = ["n", "contested"], one_strand = false, proxy = "primer", minor = "upstream@100" }
marginal = { against = "H-base", measure = "primary", upper_bound = 0.02 }
"#;

    #[test]
    fn parses_a_complete_file_and_pins_each_policy() -> Result<(), PreregError> {
        let prereg = parse(VALID.as_bytes())?;
        assert_eq!(prereg.hypotheses.len(), 2);
        assert_eq!(
            prereg.hypotheses[0].policy,
            Policy {
                name: "H-base".to_owned(),
                ..Policy::draft()
            }
        );
        assert_ne!(
            prereg.hypotheses[0].policy_sha256,
            prereg.hypotheses[1].policy_sha256
        );
        assert_eq!(prereg.gates[0].measure, Measure::HaplotypeIupacCompatible);
        Ok(())
    }

    #[test]
    fn refuses_pending_values_and_policies_without_n() {
        for (from, to, error) in [
            (
                "upper_bound = 0.01",
                "upper_bound = \"pending\"",
                "gates[0].upper_bound is still \"pending\"",
            ),
            (
                "form = \"point\"",
                "form = \"pending\"",
                "utility.form is still \"pending\"",
            ),
            (
                "minimum_pass_rate = 0.10",
                "minimum_pass_rate = \"pending\"",
                "utility.minimum_pass_rate is still \"pending\"",
            ),
            (
                "dna_revision = \"ee\"",
                "dna_revision = \"pending\"",
                "inputs.dna_revision is still \"pending\"",
            ),
            (
                "causes = [\"n\", \"contested\"]",
                "causes = [\"contested\"]",
                "hypotheses[1].policy.causes: must contain \"n\"",
            ),
            (
                "against = \"H-base\"",
                "against = \"H-later\"",
                "hypotheses[1].marginal.against: must name an earlier hypothesis",
            ),
            (
                "requires = \"previous_pass\"",
                "requires = \"previous_fail\"",
                "hypotheses[1].requires: \"previous_fail\": a later hypothesis requires \
                 \"previous_pass\"; the fixed sequence stops at the first failure",
            ),
        ] {
            let text = VALID.replacen(from, to, 1);
            assert_eq!(
                parse(text.as_bytes())
                    .map(|_| ())
                    .map_err(|error| error.to_string()),
                Err(error.to_owned()),
                "{to}"
            );
        }
    }
}
