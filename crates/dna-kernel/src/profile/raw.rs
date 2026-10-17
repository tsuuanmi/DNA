//! Strict profile TOML records and their validation.

use std::collections::BTreeSet;

use serde::Deserialize;

use super::window::{Anchor, NomenclatureWindow, WindowRule};
use super::{
    Conformance, ConformanceRule, Convention, ConventionKind, ConventionLabel, IndelPlacement,
    Notation, NotationStyle, Profile, ProfileIdentity,
};
use crate::error::{ProfileError, Result};
use crate::model::nucleotide::{Nucleotide, is_canonical};
use crate::model::reference::ReferenceTopology;
use crate::reference::MAX_REFERENCE_LENGTH;

/// Profile schema version this build accepts.
const SCHEMA_VERSION: u32 = 1;
/// Longest accepted profile identifier.
const MAX_ID_LENGTH: usize = 64;

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub(super) struct RawProfile {
    schema_version: u32,
    id: String,
    reference: RawReference,
    variant_calling: RawVariantCalling,
    normalization: Option<RawNormalization>,
    nomenclature: Option<RawNomenclature>,
    notation: Option<RawNotation>,
    conformance: Option<RawConformance>,
    empop: Option<RawEmpop>,
}

/// `[empop]` takes no keys: declaring it turns the export on.
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct RawEmpop {}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct RawConformance {
    rules: Vec<ConformanceRule>,
    minimum_run_length: usize,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct RawReference {
    sequence_sha256: Option<String>,
    topology: ReferenceTopology,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct RawVariantCalling {
    regions: Vec<[usize; 2]>,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct RawNormalization {
    indel_placement: IndelPlacement,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct RawNotation {
    style: RawNotationStyle,
}

#[derive(Debug, Clone, Copy, Deserialize)]
#[serde(rename_all = "snake_case")]
enum RawNotationStyle {
    PerBaseDecimal,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct RawNomenclature {
    windows: Vec<RawWindow>,
    #[serde(default)]
    conventions: Vec<RawConvention>,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct RawConvention {
    name: String,
    label: ConventionLabel,
    kind: RawConventionKind,
    /// Inclusive 1-based window.
    window: [usize; 2],
    base: char,
    minimum_run: usize,
}

#[derive(Debug, Clone, Copy, Deserialize)]
#[serde(rename_all = "snake_case")]
enum RawConventionKind {
    RunLength,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct RawWindow {
    name: String,
    /// 1-based position of the first window base.
    start: usize,
    sequence: String,
    structure: RawStructure,
    rules: Vec<RawRule>,
    canonical_haplotype: Option<Vec<RawSubstitution>>,
}

#[derive(Debug, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
enum RawStructure {
    AnchoredHomopolymer {
        repeat_base: char,
        /// 1-based position of the anchor base.
        anchor: usize,
        anchor_base: char,
    },
    TandemRepeat {
        motif: String,
    },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Deserialize)]
#[serde(rename_all = "snake_case")]
enum RawRule {
    AnchoredRunLengths,
    AnchoredTract,
    AnchorDuplication,
    AnchorDeletion,
    TerminalRepeatInsertion,
    CanonicalHaplotype,
    MotifShift,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct RawSubstitution {
    /// 1-based substituted position.
    position: usize,
    base: char,
}

/// Parses strict profile TOML without validating its scientific content.
pub(super) fn parse(text: &str) -> std::result::Result<RawProfile, toml::de::Error> {
    toml::from_str(text)
}

impl RawProfile {
    /// Validates every profile rule and records the file identity.
    pub(super) fn validate(self, sha256: String) -> Result<Profile> {
        if self.schema_version != SCHEMA_VERSION {
            return Err(ProfileError::UnsupportedSchemaVersion {
                found: self.schema_version,
                expected: SCHEMA_VERSION,
            }
            .into());
        }
        if !valid_id(&self.id) {
            return Err(ProfileError::Constraint(
                "id must be 1-64 lowercase ASCII letters, digits, '.', '_' or '-', starting with a letter or digit",
            )
            .into());
        }
        if let Some(digest) = &self.reference.sequence_sha256
            && !(digest.len() == 64
                && digest
                    .bytes()
                    .all(|byte| matches!(byte, b'0'..=b'9' | b'a'..=b'f')))
        {
            return Err(ProfileError::Constraint(
                "reference.sequence_sha256 must be 64 lowercase hexadecimal digits",
            )
            .into());
        }
        validate_regions(&self.variant_calling.regions)?;
        let notation = match (self.normalization, self.notation) {
            (Some(normalization), Some(notation)) => Some(Notation {
                indel_placement: normalization.indel_placement,
                style: match notation.style {
                    RawNotationStyle::PerBaseDecimal => NotationStyle::PerBaseDecimal,
                },
            }),
            (None, None) => None,
            _ => {
                return Err(ProfileError::Constraint(
                    "normalization and notation must be declared together",
                )
                .into());
            }
        };
        let conformance = match self.conformance {
            None => None,
            Some(_) if notation.is_none() => {
                return Err(ProfileError::Constraint(
                    "conformance requires the notation it checks",
                )
                .into());
            }
            Some(raw) => {
                let mut sorted = raw.rules.clone();
                sorted.sort_unstable();
                sorted.dedup();
                if raw.rules.is_empty() || sorted.len() != raw.rules.len() {
                    return Err(ProfileError::Constraint(
                        "conformance.rules must name at least one rule, each once",
                    )
                    .into());
                }
                if raw.minimum_run_length < 2 {
                    return Err(ProfileError::Constraint(
                        "conformance.minimum_run_length must be at least 2",
                    )
                    .into());
                }
                Some(Conformance {
                    rules: raw.rules,
                    minimum_run_length: raw.minimum_run_length,
                })
            }
        };
        if self.empop.is_some() && notation.is_none() {
            return Err(ProfileError::Constraint("empop requires the notation it exports").into());
        }
        let (windows, conventions) = match self.nomenclature {
            Some(nomenclature) if nomenclature.windows.is_empty() => {
                return Err(ProfileError::Constraint(
                    "nomenclature.windows must contain at least one window",
                )
                .into());
            }
            Some(nomenclature) => {
                if !nomenclature.conventions.is_empty() && notation.is_none() {
                    return Err(ProfileError::Constraint(
                        "nomenclature.conventions require the notation they label",
                    )
                    .into());
                }
                (
                    validate_windows(nomenclature.windows)?,
                    validate_conventions(nomenclature.conventions)?,
                )
            }
            None => (Vec::new(), Vec::new()),
        };
        Ok(Profile {
            identity: ProfileIdentity {
                id: self.id,
                sha256,
            },
            reference_sha256: self.reference.sequence_sha256,
            topology: self.reference.topology,
            regions: self.variant_calling.regions,
            windows,
            notation,
            conformance,
            conventions,
            empop: self.empop.is_some(),
        })
    }
}

fn valid_id(id: &str) -> bool {
    id.len() <= MAX_ID_LENGTH
        && id
            .bytes()
            .next()
            .is_some_and(|first| first.is_ascii_lowercase() || first.is_ascii_digit())
        && id.bytes().all(|byte| {
            byte.is_ascii_lowercase() || byte.is_ascii_digit() || matches!(byte, b'.' | b'_' | b'-')
        })
}

fn validate_regions(regions: &[[usize; 2]]) -> Result<()> {
    if regions.is_empty() {
        return Err(ProfileError::Constraint(
            "variant_calling.regions must contain at least one inclusive range",
        )
        .into());
    }
    for (index, &[start, end]) in regions.iter().enumerate() {
        if start == 0 || start > end || end > MAX_REFERENCE_LENGTH {
            return Err(ProfileError::RegionOutOfBounds {
                index,
                maximum: MAX_REFERENCE_LENGTH,
            }
            .into());
        }
    }
    Ok(())
}

fn validate_windows(raw: Vec<RawWindow>) -> Result<Vec<NomenclatureWindow>> {
    let mut windows: Vec<NomenclatureWindow> = Vec::with_capacity(raw.len());
    for window in raw {
        let validated = validate_window(window)?;
        if let Some(previous) = windows.last()
            && validated.start < previous.end()
        {
            return Err(invalid(
                &validated.name,
                "windows must be in reference order without overlap",
            ));
        }
        if windows.iter().any(|other| other.name == validated.name) {
            return Err(invalid(&validated.name, "window names must be unique"));
        }
        windows.push(validated);
    }
    Ok(windows)
}

fn validate_conventions(raw: Vec<RawConvention>) -> Result<Vec<Convention>> {
    let mut conventions: Vec<Convention> = Vec::with_capacity(raw.len());
    for convention in raw {
        let RawConvention {
            name,
            label,
            kind,
            window,
            base,
            minimum_run,
        } = convention;
        let invalid = |reason| ProfileError::InvalidConvention {
            convention: name.clone(),
            reason,
        };
        if name.trim().is_empty() {
            return Err(ProfileError::Constraint(
                "nomenclature convention names must not be empty",
            )
            .into());
        }
        let [start, end] = window;
        if start == 0 || start > end || end > MAX_REFERENCE_LENGTH {
            return Err(
                invalid("window must satisfy 1 <= start <= end within the reference").into(),
            );
        }
        if let Some(previous) = conventions.last()
            && start <= previous.window[1]
        {
            return Err(invalid("conventions must be in reference order without overlap").into());
        }
        if conventions.iter().any(|other| other.name == name) {
            return Err(invalid("convention names must be unique").into());
        }
        let kind = match kind {
            RawConventionKind::RunLength => {
                let base = base_byte(base).ok_or_else(|| invalid("base must be A/C/G/T"))?;
                if minimum_run < 2 {
                    return Err(invalid("minimum_run must be at least 2").into());
                }
                ConventionKind::RunLength { base, minimum_run }
            }
        };
        conventions.push(Convention {
            name,
            label,
            kind,
            window,
        });
    }
    Ok(conventions)
}

fn validate_window(window: RawWindow) -> Result<NomenclatureWindow> {
    let RawWindow {
        name,
        start,
        sequence,
        structure,
        rules,
        canonical_haplotype,
    } = window;
    if name.trim().is_empty() {
        return Err(ProfileError::Constraint("nomenclature window names must not be empty").into());
    }
    if start == 0 {
        return Err(invalid(&name, "start must be a 1-based position"));
    }
    if sequence.is_empty() || !sequence.chars().all(is_canonical) {
        return Err(invalid(
            &name,
            "sequence must be non-empty uppercase A/C/G/T",
        ));
    }
    let offset = start - 1;
    if offset + sequence.len() > MAX_REFERENCE_LENGTH {
        return Err(invalid(
            &name,
            "window ends beyond the supported reference length",
        ));
    }
    if rules.is_empty() {
        return Err(invalid(&name, "rules must name at least one rule"));
    }
    if rules.iter().collect::<BTreeSet<_>>().len() != rules.len() {
        return Err(invalid(&name, "rules must not repeat"));
    }
    if canonical_haplotype.is_some() != rules.contains(&RawRule::CanonicalHaplotype) {
        return Err(invalid(
            &name,
            "canonical_haplotype must be given exactly when the canonical_haplotype rule is used",
        ));
    }

    let bytes = sequence.as_bytes();
    let (anchor, motif) = match structure {
        RawStructure::AnchoredHomopolymer {
            repeat_base,
            anchor,
            anchor_base,
        } => {
            let (Some(repeat), Some(base)) = (base_byte(repeat_base), base_byte(anchor_base))
            else {
                return Err(invalid(
                    &name,
                    "repeat_base and anchor_base must be A/C/G/T",
                ));
            };
            let index = anchor
                .checked_sub(start)
                .filter(|&index| index < bytes.len())
                .ok_or_else(|| invalid(&name, "anchor must lie inside the window"))?;
            if repeat == base || bytes[index] != base {
                return Err(invalid(
                    &name,
                    "anchor_base must differ from repeat_base and match the window sequence",
                ));
            }
            (
                Some(Anchor {
                    index,
                    base,
                    repeat,
                }),
                None,
            )
        }
        RawStructure::TandemRepeat { motif } => {
            if motif.is_empty()
                || !motif.chars().all(is_canonical)
                || motif.len() >= sequence.len()
                || !sequence.ends_with(motif.as_str())
            {
                return Err(invalid(
                    &name,
                    "motif must be A/C/G/T, shorter than the window, and end it",
                ));
            }
            (None, Some(motif))
        }
    };

    let mut validated = Vec::with_capacity(rules.len());
    for rule in rules {
        validated.push(match rule {
            RawRule::AnchoredRunLengths => {
                let anchor = require_anchor(&name, anchor)?;
                if sequence.matches(char::from(anchor.base)).count() != 1
                    || bytes
                        .iter()
                        .any(|&byte| byte != anchor.base && byte != anchor.repeat)
                {
                    return Err(invalid(
                        &name,
                        "anchored_run_lengths needs a window of repeat bases around one anchor",
                    ));
                }
                WindowRule::AnchoredRunLengths(anchor)
            }
            RawRule::AnchoredTract => {
                let anchor = require_anchor(&name, anchor)?;
                let lead = bytes
                    .iter()
                    .take_while(|&&byte| byte == bytes[0])
                    .count();
                let runs = &bytes[lead..];
                if bytes[0] == anchor.repeat
                    || bytes[0] == anchor.base
                    || anchor.index <= lead
                    || !repeat_runs_around(runs, anchor)
                {
                    return Err(invalid(
                        &name,
                        "anchored_tract needs a lead run, then repeat bases around one anchor with a run on each side",
                    ));
                }
                WindowRule::AnchoredTract { anchor, lead }
            }
            RawRule::AnchorDuplication => {
                WindowRule::AnchorDuplication(require_run_after_anchor(&name, anchor, bytes)?)
            }
            RawRule::AnchorDeletion => {
                WindowRule::AnchorDeletion(require_run_after_anchor(&name, anchor, bytes)?)
            }
            RawRule::TerminalRepeatInsertion => {
                WindowRule::TerminalRepeatInsertion(require_anchor(&name, anchor)?)
            }
            RawRule::CanonicalHaplotype => WindowRule::CanonicalHaplotype(canonical_substitutions(
                &name,
                start,
                bytes,
                canonical_haplotype.as_deref(),
            )?),
            RawRule::MotifShift => WindowRule::MotifShift {
                motif: motif
                    .clone()
                    .ok_or_else(|| invalid(&name, "motif_shift needs a tandem_repeat structure"))?,
            },
        });
    }

    Ok(NomenclatureWindow {
        name,
        start: offset,
        sequence,
        rules: validated,
    })
}

/// Whether `bytes` are repeat bases around exactly one anchor base, with at
/// least one repeat base on each side.
fn repeat_runs_around(bytes: &[u8], anchor: Anchor) -> bool {
    let first = bytes.iter().position(|&byte| byte == anchor.base);
    first.is_some()
        && first == bytes.iter().rposition(|&byte| byte == anchor.base)
        && bytes
            .iter()
            .all(|&byte| byte == anchor.base || byte == anchor.repeat)
        && bytes.first() == Some(&anchor.repeat)
        && bytes.last() == Some(&anchor.repeat)
}

fn require_anchor(name: &str, anchor: Option<Anchor>) -> Result<Anchor> {
    anchor.ok_or_else(|| invalid(name, "anchor rules need an anchored_homopolymer structure"))
}

fn require_run_after_anchor(name: &str, anchor: Option<Anchor>, bytes: &[u8]) -> Result<Anchor> {
    let anchor = require_anchor(name, anchor)?;
    if bytes.get(anchor.index + 1) != Some(&anchor.repeat) {
        return Err(invalid(
            name,
            "anchor_duplication and anchor_deletion need a repeat run after the anchor",
        ));
    }
    Ok(anchor)
}

fn canonical_substitutions(
    name: &str,
    start: usize,
    bytes: &[u8],
    raw: Option<&[RawSubstitution]>,
) -> Result<Vec<(usize, u8)>> {
    let raw = raw.filter(|raw| !raw.is_empty()).ok_or_else(|| {
        invalid(
            name,
            "canonical_haplotype must list at least one substitution",
        )
    })?;
    let mut substitutions = Vec::with_capacity(raw.len());
    for substitution in raw {
        let index = substitution
            .position
            .checked_sub(start)
            .filter(|&index| index < bytes.len())
            .ok_or_else(|| {
                invalid(
                    name,
                    "canonical_haplotype positions must lie inside the window",
                )
            })?;
        let base = base_byte(substitution.base)
            .filter(|&base| base != bytes[index])
            .ok_or_else(|| {
                invalid(
                    name,
                    "canonical_haplotype bases must be A/C/G/T and differ from the reference",
                )
            })?;
        if substitutions.iter().any(|&(other, _)| other == index) {
            return Err(invalid(
                name,
                "canonical_haplotype positions must be unique",
            ));
        }
        substitutions.push((index, base));
    }
    substitutions.sort_unstable();
    Ok(substitutions)
}

fn base_byte(base: char) -> Option<u8> {
    Nucleotide::from_char(base).map(Nucleotide::as_byte)
}

fn invalid(window: &str, reason: &'static str) -> crate::error::Error {
    ProfileError::InvalidWindow {
        window: window.to_owned(),
        reason,
    }
    .into()
}

#[cfg(test)]
mod tests {
    use crate::error::Error;

    use super::*;

    const VALID: &str = r#"
schema_version = 1
id = "toy"
[reference]
topology = "linear"
[variant_calling]
regions = [[1, 40]]
[normalization]
indel_placement = "right"
[notation]
style = "per_base_decimal"
[[nomenclature.windows]]
name = "poly-G"
start = 11
sequence = "GGGAGG"
structure = { kind = "anchored_homopolymer", repeat_base = "G", anchor = 14, anchor_base = "A" }
rules = ["anchored_run_lengths", "anchor_deletion", "canonical_haplotype"]
canonical_haplotype = [{ position = 12, base = "T" }]
[[nomenclature.windows]]
name = "TG-repeat"
start = 21
sequence = "ATGTGTG"
structure = { kind = "tandem_repeat", motif = "TG" }
rules = ["motif_shift"]
"#;

    fn validated(text: &str) -> std::result::Result<Result<Profile>, toml::de::Error> {
        Ok(parse(text)?.validate("0".repeat(64)))
    }

    fn window_reason(text: &str) -> Option<&'static str> {
        match validated(text) {
            Ok(Err(Error::Profile(ProfileError::InvalidWindow { reason, .. }))) => Some(reason),
            _ => None,
        }
    }

    #[test]
    fn validates_a_target_without_any_mtdna_knowledge()
    -> std::result::Result<(), Box<dyn std::error::Error>> {
        let profile = validated(VALID)??;
        assert_eq!(profile.identity.id, "toy");
        assert_eq!(profile.reference_sha256, None);
        let anchor = Anchor {
            index: 3,
            base: b'A',
            repeat: b'G',
        };
        assert_eq!(
            profile.windows,
            [
                NomenclatureWindow {
                    name: "poly-G".into(),
                    start: 10,
                    sequence: "GGGAGG".into(),
                    rules: vec![
                        WindowRule::AnchoredRunLengths(anchor),
                        WindowRule::AnchorDeletion(anchor),
                        WindowRule::CanonicalHaplotype(vec![(1, b'T')]),
                    ],
                },
                NomenclatureWindow {
                    name: "TG-repeat".into(),
                    start: 20,
                    sequence: "ATGTGTG".into(),
                    rules: vec![WindowRule::MotifShift { motif: "TG".into() }],
                },
            ]
        );
        Ok(())
    }

    #[test]
    fn representation_sections_are_optional() -> std::result::Result<(), Box<dyn std::error::Error>>
    {
        let minimal = "schema_version = 1\nid = \"plain\"\n[reference]\ntopology = \"circular\"\n[variant_calling]\nregions = [[1, 10]]\n";
        let profile = validated(minimal)??;
        assert_eq!(profile.notation, None);
        assert!(profile.windows.is_empty());
        Ok(())
    }

    #[test]
    fn validates_report_conventions() -> std::result::Result<(), Box<dyn std::error::Error>> {
        let convention = |name: &str, window: &str, base: &str, minimum: &str| {
            format!(
                "[[nomenclature.conventions]]\nname = \"{name}\"\nlabel = \"not_recorded\"\n\
                 kind = \"run_length\"\nwindow = {window}\nbase = \"{base}\"\nminimum_run = {minimum}\n"
            )
        };
        let declared = format!(
            "{VALID}{}{}",
            convention("g-length", "[11, 16]", "G", "5"),
            convention("t-length", "[21, 27]", "T", "3")
        );
        let conventions = validated(&declared)??.conventions;
        assert_eq!(
            conventions
                .iter()
                .map(|convention| (convention.name.as_str(), convention.window, convention.kind))
                .collect::<Vec<_>>(),
            [
                (
                    "g-length",
                    [11, 16],
                    ConventionKind::RunLength {
                        base: b'G',
                        minimum_run: 5
                    }
                ),
                (
                    "t-length",
                    [21, 27],
                    ConventionKind::RunLength {
                        base: b'T',
                        minimum_run: 3
                    }
                ),
            ]
        );
        assert_eq!(conventions[0].label, ConventionLabel::NotRecorded);
        assert!(validated(VALID)??.conventions.is_empty());
        let reason = |text: &str| match validated(text) {
            Ok(Err(Error::Profile(ProfileError::InvalidConvention { reason, .. }))) => Some(reason),
            _ => None,
        };
        let overlapping = format!(
            "{VALID}{}{}",
            convention("g-length", "[11, 16]", "G", "5"),
            convention("t-length", "[16, 27]", "T", "3")
        );
        assert_eq!(
            reason(&overlapping),
            Some("conventions must be in reference order without overlap")
        );
        assert_eq!(
            reason(&format!("{VALID}{}", convention("n", "[16, 11]", "G", "5"))),
            Some("window must satisfy 1 <= start <= end within the reference")
        );
        assert_eq!(
            reason(&format!("{VALID}{}", convention("n", "[11, 16]", "N", "5"))),
            Some("base must be A/C/G/T")
        );
        assert_eq!(
            reason(&format!("{VALID}{}", convention("n", "[11, 16]", "G", "1"))),
            Some("minimum_run must be at least 2")
        );
        let without_notation = declared
            .replace("[normalization]\nindel_placement = \"right\"\n", "")
            .replace("[notation]\nstyle = \"per_base_decimal\"\n", "");
        assert!(matches!(
            validated(&without_notation)?,
            Err(Error::Profile(ProfileError::Constraint(_)))
        ));
        Ok(())
    }

    #[test]
    fn declares_the_empop_export_with_notation()
    -> std::result::Result<(), Box<dyn std::error::Error>> {
        assert!(!validated(VALID)??.empop);
        assert!(validated(&format!("{VALID}[empop]\n"))??.empop);
        assert!(parse(&format!("{VALID}[empop]\nranges = 1\n")).is_err());
        let without_notation = format!("{VALID}[empop]\n")
            .replace("[normalization]\nindel_placement = \"right\"\n", "")
            .replace("[notation]\nstyle = \"per_base_decimal\"\n", "");
        assert!(matches!(
            validated(&without_notation)?,
            Err(Error::Profile(ProfileError::Constraint(_)))
        ));
        Ok(())
    }

    #[test]
    fn validates_conformance_rules() -> std::result::Result<(), Box<dyn std::error::Error>> {
        let declared = format!(
            "{VALID}[conformance]\nrules = [\"insertion_matches_run\", \"insertion_at_run_end\"]\nminimum_run_length = 4\n"
        );
        assert_eq!(
            validated(&declared)??.conformance,
            Some(Conformance {
                rules: vec![
                    ConformanceRule::InsertionMatchesRun,
                    ConformanceRule::InsertionAtRunEnd,
                ],
                minimum_run_length: 4,
            })
        );
        assert_eq!(validated(VALID)??.conformance, None);
        for invalid in [
            declared.replace("\"insertion_matches_run\", \"insertion_at_run_end\"", ""),
            declared.replace("\"insertion_at_run_end\"", "\"insertion_matches_run\""),
            declared.replace("minimum_run_length = 4", "minimum_run_length = 1"),
            declared
                .replace("[normalization]\nindel_placement = \"right\"\n", "")
                .replace("[notation]\nstyle = \"per_base_decimal\"\n", ""),
        ] {
            assert!(
                matches!(
                    validated(&invalid)?,
                    Err(Error::Profile(ProfileError::Constraint(_)))
                ),
                "unexpectedly accepted {invalid}"
            );
        }
        assert!(parse(&declared.replace("insertion_at_run_end", "anything")).is_err());
        Ok(())
    }

    #[test]
    fn rejects_unknown_keys_and_old_schemas() -> std::result::Result<(), toml::de::Error> {
        assert!(parse(&VALID.replace("id = \"toy\"", "id = \"toy\"\nextra = 1")).is_err());
        assert!(parse(&VALID.replace("motif = \"TG\" }", "motif = \"TG\", anchor = 1 }")).is_err());
        assert!(matches!(
            validated(&VALID.replace("schema_version = 1", "schema_version = 2"))?,
            Err(Error::Profile(ProfileError::UnsupportedSchemaVersion {
                found: 2,
                expected: 1
            }))
        ));
        Ok(())
    }

    #[test]
    fn rejects_invalid_identity_regions_and_dependencies()
    -> std::result::Result<(), toml::de::Error> {
        for invalid in [
            VALID.replace("id = \"toy\"", "id = \"\""),
            VALID.replace("id = \"toy\"", "id = \"Toy\""),
            VALID.replace("id = \"toy\"", "id = \"-toy\""),
            VALID.replace(
                "topology = \"linear\"",
                "topology = \"linear\"\nsequence_sha256 = \"ABC\"",
            ),
            VALID.replace("regions = [[1, 40]]", "regions = []"),
            VALID.replace("[normalization]\nindel_placement = \"right\"\n", ""),
            VALID.replace("[notation]\nstyle = \"per_base_decimal\"\n", ""),
            VALID.replace("name = \"poly-G\"", "name = \" \""),
            format!(
                "{}[nomenclature]\nwindows = []\n",
                &VALID[..VALID
                    .find("[[nomenclature.windows]]")
                    .unwrap_or(VALID.len())]
            ),
        ] {
            assert!(
                matches!(
                    validated(&invalid)?,
                    Err(Error::Profile(ProfileError::Constraint(_)))
                ),
                "unexpectedly accepted {invalid}"
            );
        }
        for region in ["[[0, 4]]", "[[5, 4]]", "[[1, 50001]]"] {
            assert!(matches!(
                validated(&VALID.replace("[[1, 40]]", region))?,
                Err(Error::Profile(ProfileError::RegionOutOfBounds {
                    index: 0,
                    ..
                }))
            ));
        }
        Ok(())
    }

    #[test]
    fn rejects_inconsistent_windows() {
        let cases = [
            (
                "sequence = \"GGGAGG\"",
                "sequence = \"GGNAGG\"",
                "sequence must be non-empty uppercase A/C/G/T",
            ),
            (
                "anchor = 14",
                "anchor = 20",
                "anchor must lie inside the window",
            ),
            (
                "anchor = 14",
                "anchor = 13",
                "anchor_base must differ from repeat_base and match the window sequence",
            ),
            (
                "motif = \"TG\"",
                "motif = \"CA\"",
                "motif must be A/C/G/T, shorter than the window, and end it",
            ),
            (
                "rules = [\"motif_shift\"]",
                "rules = [\"anchor_deletion\"]",
                "anchor rules need an anchored_homopolymer structure",
            ),
            (
                "rules = [\"anchored_run_lengths\", \"anchor_deletion\", \"canonical_haplotype\"]",
                "rules = [\"motif_shift\", \"canonical_haplotype\"]",
                "motif_shift needs a tandem_repeat structure",
            ),
            (
                "rules = [\"motif_shift\"]",
                "rules = []",
                "rules must name at least one rule",
            ),
            (
                "rules = [\"motif_shift\"]",
                "rules = [\"motif_shift\", \"motif_shift\"]",
                "rules must not repeat",
            ),
            (
                "canonical_haplotype = [{ position = 12, base = \"T\" }]\n",
                "",
                "canonical_haplotype must be given exactly when the canonical_haplotype rule is used",
            ),
            (
                "position = 12, base = \"T\"",
                "position = 12, base = \"G\"",
                "canonical_haplotype bases must be A/C/G/T and differ from the reference",
            ),
            (
                "position = 12, base = \"T\"",
                "position = 30, base = \"T\"",
                "canonical_haplotype positions must lie inside the window",
            ),
            (
                "start = 21",
                "start = 15",
                "windows must be in reference order without overlap",
            ),
            (
                "name = \"TG-repeat\"",
                "name = \"poly-G\"",
                "window names must be unique",
            ),
            (
                "sequence = \"GGGAGG\"",
                "sequence = \"GGGAGT\"",
                "anchored_run_lengths needs a window of repeat bases around one anchor",
            ),
            (
                "start = 11",
                "start = 0",
                "start must be a 1-based position",
            ),
            (
                "start = 21",
                "start = 49999",
                "window ends beyond the supported reference length",
            ),
            (
                "canonical_haplotype = [{ position = 12, base = \"T\" }]",
                "canonical_haplotype = []",
                "canonical_haplotype must list at least one substitution",
            ),
            (
                "[{ position = 12, base = \"T\" }]",
                "[{ position = 12, base = \"T\" }, { position = 12, base = \"C\" }]",
                "canonical_haplotype positions must be unique",
            ),
            (
                "repeat_base = \"G\"",
                "repeat_base = \"N\"",
                "repeat_base and anchor_base must be A/C/G/T",
            ),
            (
                "repeat_base = \"G\"",
                "repeat_base = \"A\"",
                "anchor_base must differ from repeat_base and match the window sequence",
            ),
        ];
        for (from, to, reason) in cases {
            assert_eq!(
                window_reason(&VALID.replace(from, to)),
                Some(reason),
                "{from} -> {to}"
            );
        }
        let no_run = VALID
            .replace("sequence = \"GGGAGG\"", "sequence = \"GGGGGA\"")
            .replace("anchor = 14", "anchor = 16")
            .replace("\"anchored_run_lengths\", ", "");
        assert_eq!(
            window_reason(&no_run),
            Some("anchor_duplication and anchor_deletion need a repeat run after the anchor")
        );
    }
}
