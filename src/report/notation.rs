//! The `per_base_decimal` notation style for represented variants.
//!
//! Rendering is mechanical serialization (ADR-0060 §8): each changed base gets
//! one call. A substitution is `<position><base>` (`73G`), each deleted base is
//! `<position>DEL` (`249DEL`), and each inserted base is
//! `<anchor>.<ordinal><base>` after the preceding reference base (`309.1C`,
//! `309.2C`); an insertion before the first base uses anchor `0`.

use std::collections::BTreeMap;

use crate::model::notation_result::{
    ConformanceResult, FindingResult, NotationProvenanceResult, NotationResult, SourceResult,
};
use crate::model::result::ReferenceResult;
use crate::model::sample_result::{SampleNotationCallResult, SampleNotationResult};
use crate::report::json::{project_plugins, project_profile};
use dna_kernel::error::{ReportError, RepresentationError};
use dna_kernel::model::reference::Reference;
use dna_kernel::plugin::PluginDescriptor;
use dna_kernel::profile::{ConformanceRule, NotationStyle, ProfileIdentity};
use dna_kernel::variant::Variant;
use dna_post::conformance::Finding;
use dna_post::variant_representation::{sort_edits, variants_to_edits};

/// Every read's represented calls and the profile style to render them in.
pub(crate) struct SampleNotation {
    pub(crate) style: NotationStyle,
    pub(crate) reads: Vec<ReadRepresentation>,
}

/// One read's eligible calls in the profile representation.
pub(crate) struct ReadRepresentation {
    pub(crate) input_sha256: String,
    pub(crate) variants: Vec<Variant>,
}

/// One read of a variants document: its content identity and reviewer-facing
/// name, in document order.
pub(crate) struct NamedRead {
    pub(crate) sha256: String,
    pub(crate) name: String,
}

/// Inputs consumed to build one immutable notation document.
pub(crate) struct CompletedNotation {
    pub(crate) sample_id: String,
    pub(crate) reference: Reference,
    pub(crate) profile: ProfileIdentity,
    pub(crate) configuration_sha256: String,
    /// SHA-256 of the source variants document.
    pub(crate) source_sha256: String,
    pub(crate) reads: Vec<NamedRead>,
    pub(crate) notation: SampleNotation,
    /// The profile's rules and each read's findings, in document read order,
    /// when the profile declares conformance.
    pub(crate) conformance: Option<(Vec<ConformanceRule>, Vec<Vec<Finding>>)>,
    /// Plugins of the workflow, in execution order.
    pub(crate) plugins: &'static [&'static PluginDescriptor],
}

/// Builds `dna.notation/v1` without filesystem side effects.
pub(crate) fn build(completed: CompletedNotation) -> dna_kernel::error::Result<NotationResult> {
    let CompletedNotation {
        sample_id,
        reference,
        profile,
        configuration_sha256,
        source_sha256,
        reads,
        notation,
        conformance,
        plugins,
    } = completed;
    let identities = reads
        .iter()
        .map(|read| read.sha256.as_str())
        .collect::<Vec<_>>();
    let names = reads
        .iter()
        .map(|read| read.name.clone())
        .collect::<Vec<_>>();
    let conformance = conformance
        .map(
            |(rules, findings)| -> dna_kernel::error::Result<ConformanceResult> {
                if findings.len() != reads.len() {
                    return Err(ReportError::Inconsistent(
                        "conformance findings do not match the document reads",
                    )
                    .into());
                }
                let findings = findings
                    .into_iter()
                    .zip(&names)
                    .flat_map(|(findings, name)| {
                        findings.into_iter().map(move |finding| (finding, name))
                    })
                    .map(|(finding, name)| {
                        Ok(FindingResult {
                            rule: finding.rule.label(),
                            read: name.clone(),
                            calls: render(&reference.name, &reference.sequence, &[finding.variant])
                                .map_err(ReportError::Representation)?
                                .into_iter()
                                .map(|call| call.text)
                                .collect(),
                        })
                    })
                    .collect::<dna_kernel::error::Result<_>>()?;
                Ok(ConformanceResult {
                    rules: rules.into_iter().map(ConformanceRule::label).collect(),
                    findings,
                })
            },
        )
        .transpose()?;
    let notation = project(&reference, &identities, &names, notation)?;
    Ok(NotationResult {
        schema_version: "dna.notation/v1",
        sample_id,
        provenance: NotationProvenanceResult {
            reference: ReferenceResult {
                name: reference.name,
                topology: reference.topology,
                sha256: reference.sequence_sha256,
            },
            configuration_sha256,
            profile: project_profile(profile),
            plugins: project_plugins(plugins),
            source: SourceResult {
                schema_version: "dna.variants/v1",
                sha256: source_sha256,
            },
        },
        notation,
        conformance,
    })
}

/// Renders each read's represented calls and lists, per distinct call, the reads
/// containing it in the document's read order. This is not a consensus: reads
/// that disagree contribute different calls.
///
/// `identities` and `names` give the content identity and reviewer-facing name
/// of every read, in the document's read order.
pub(super) fn project(
    reference: &Reference,
    identities: &[&str],
    names: &[String],
    notation: SampleNotation,
) -> dna_kernel::error::Result<SampleNotationResult> {
    let SampleNotation { style, reads } = notation;
    let render = match style {
        NotationStyle::PerBaseDecimal => render,
    };
    let mut calls: BTreeMap<NotationCall, Vec<usize>> = BTreeMap::new();
    for representation in reads {
        let read_index = identities
            .iter()
            .position(|identity| *identity == representation.input_sha256)
            .ok_or(ReportError::Inconsistent(
                "notation references a read missing from the document",
            ))?;
        for call in render(
            &reference.name,
            &reference.sequence,
            &representation.variants,
        )
        .map_err(ReportError::Representation)?
        {
            calls.entry(call).or_default().push(read_index);
        }
    }
    let calls = calls
        .into_iter()
        .map(|(call, mut reads)| {
            reads.sort_unstable();
            reads.dedup();
            Ok(SampleNotationCallResult {
                call: call.text,
                reads: reads
                    .into_iter()
                    .map(|index| {
                        names
                            .get(index)
                            .cloned()
                            .ok_or(ReportError::MissingRead { index }.into())
                    })
                    .collect::<dna_kernel::error::Result<_>>()?,
            })
        })
        .collect::<dna_kernel::error::Result<_>>()?;
    Ok(SampleNotationResult {
        style: style.label(),
        calls,
    })
}

/// One rendered per-base call with its reference-order sort key.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub(super) struct NotationCall {
    /// 1-based position of the changed base, or of the base an insertion follows.
    pub(super) position: usize,
    /// `0` for substitutions and deletions, `1..` for successive inserted bases.
    pub(super) ordinal: usize,
    /// The rendered call, for example `309.1C`.
    pub(super) text: String,
}

/// Renders represented variants as per-base calls in reference order.
pub(super) fn render(
    contig: &str,
    reference: &str,
    variants: &[Variant],
) -> Result<Vec<NotationCall>, RepresentationError> {
    let mut edits = variants_to_edits(contig, reference, variants)?;
    sort_edits(&mut edits);
    let mut calls = Vec::new();
    for edit in edits {
        let deleted = edit.end - edit.start;
        match (deleted, edit.alternate.len()) {
            (0, _) => {
                calls.extend(edit.alternate.chars().enumerate().map(|(index, base)| {
                    NotationCall {
                        position: edit.start,
                        ordinal: index + 1,
                        text: format!("{}.{}{base}", edit.start, index + 1),
                    }
                }));
            }
            (_, 0) => calls.extend((edit.start..edit.end).map(|index| NotationCall {
                position: index + 1,
                ordinal: 0,
                text: format!("{}DEL", index + 1),
            })),
            (1, 1) => calls.push(NotationCall {
                position: edit.start + 1,
                ordinal: 0,
                text: format!("{}{}", edit.start + 1, edit.alternate),
            }),
            _ => return Err(RepresentationError::UnsupportedReplacement),
        }
    }
    calls.sort();
    Ok(calls)
}

#[cfg(test)]
mod tests {
    use dna_kernel::variant::{Variant, VariantKind};

    use super::*;

    const REFERENCE: &str = "GATTACA";

    fn variant(position: usize, reference: &str, alternate: &str, kind: VariantKind) -> Variant {
        Variant {
            contig: "synthetic".into(),
            position_1based: position,
            reference: reference.into(),
            alternate: alternate.into(),
            kind,
        }
    }

    fn texts(variants: &[Variant]) -> Result<Vec<String>, RepresentationError> {
        Ok(render("synthetic", REFERENCE, variants)?
            .into_iter()
            .map(|call| call.text)
            .collect())
    }

    #[test]
    fn renders_one_call_per_changed_base_in_position_order() -> Result<(), RepresentationError> {
        let variants = [
            variant(4, "T", "TCC", VariantKind::Ins),
            variant(2, "A", "C", VariantKind::Snv),
            variant(5, "ACA", "A", VariantKind::Del),
        ];

        assert_eq!(texts(&variants)?, ["2C", "4.1C", "4.2C", "6DEL", "7DEL"]);
        Ok(())
    }

    #[test]
    fn orders_an_insertion_after_a_change_at_its_anchor() -> Result<(), RepresentationError> {
        let variants = [
            variant(3, "T", "TG", VariantKind::Ins),
            variant(3, "T", "A", VariantKind::Snv),
        ];

        assert_eq!(texts(&variants)?, ["3A", "3.1G"]);
        Ok(())
    }

    #[test]
    fn renders_a_leading_insertion_before_the_first_base() -> Result<(), RepresentationError> {
        assert_eq!(texts(&[variant(1, "G", "CG", VariantKind::Ins)])?, ["0.1C"]);
        Ok(())
    }

    #[test]
    fn rejects_variants_that_disagree_with_the_reference() {
        assert_eq!(
            texts(&[variant(2, "G", "C", VariantKind::Snv)]),
            Err(RepresentationError::ReferenceAlleleMismatch)
        );
    }
}
