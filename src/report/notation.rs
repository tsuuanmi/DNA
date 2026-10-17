//! The `per_base_decimal` notation style for represented variants.
//!
//! Rendering is mechanical serialization (ADR-0060 §8): each changed base gets
//! one call. A substitution is `<position><base>` (`73G`), each deleted base is
//! `<position>DEL` (`249DEL`), and each inserted base is
//! `<anchor>.<ordinal><base>` after the preceding reference base (`309.1C`,
//! `309.2C`); an insertion before the first base uses anchor `0`.

use std::collections::BTreeMap;

use crate::model::notation_result::{
    ConformanceResult, CutResult, EmpopResult, FindingResult, NotationProvenanceResult,
    NotationResult, SourceResult,
};
use crate::model::result::ReferenceResult;
use crate::model::sample_result::{
    ConventionResult, SampleNotationCallResult, SampleNotationResult,
};
use crate::report::json::{project_plugins, project_profile};
use dna_empop::conformance::Finding;
use dna_empop::export::Export;
use dna_kernel::error::{ReportError, RepresentationError};
use dna_kernel::model::reference::Reference;
use dna_kernel::plugin::PluginDescriptor;
use dna_kernel::profile::{ConformanceRule, Convention, NotationStyle, ProfileIdentity};
use dna_kernel::variant::Variant;
use dna_normalize::variant_nomenclature::conventions::{Call, Labels};
use dna_normalize::variant_representation::{sort_edits, variants_to_edits};

/// Every read's represented calls, the profile style to render them in, and
/// the profile's report conventions that label them.
pub(crate) struct SampleNotation {
    pub(crate) style: NotationStyle,
    pub(crate) conventions: Vec<Convention>,
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
    /// The EMPOP query, when the profile declares `[empop]`.
    pub(crate) empop: Option<Export>,
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
        empop,
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
    let notation = project(&reference, &identities, &names, &notation)?;
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
        empop: empop.map(project_empop),
    })
}

/// Projects an EMPOP query as `notation/v1` and `review/v1` publish it.
pub(super) fn project_empop(export: Export) -> EmpopResult {
    EmpopResult {
        ranges_text: export.ranges_text(),
        mitotype: export.mitotype_text(),
        ranges: export.ranges,
        cut: export
            .cut
            .into_iter()
            .map(|cut| CutResult {
                start: cut.start,
                end: cut.end,
                reason: cut.reason.label(),
            })
            .collect(),
    }
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
    notation: &SampleNotation,
) -> dna_kernel::error::Result<SampleNotationResult> {
    let index = |identity: &str| {
        identities
            .iter()
            .position(|candidate| *candidate == identity)
            .ok_or(ReportError::Inconsistent(
                "notation references a read missing from the document",
            ))
    };
    for representation in &notation.reads {
        index(&representation.input_sha256)?;
    }
    let calls = calls(reference, notation)?
        .into_iter()
        .map(|(call, CallReads { reads, convention })| {
            let mut reads = reads
                .into_iter()
                .map(index)
                .collect::<Result<Vec<_>, _>>()?;
            reads.sort_unstable();
            reads.dedup();
            Ok(SampleNotationCallResult {
                call: call.text,
                convention: convention.map(|convention| ConventionResult {
                    name: convention.name.clone(),
                    label: convention.label.label(),
                }),
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
        style: notation.style.label(),
        calls,
    })
}

/// One rendered per-base call with its reference-order sort key.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) struct NotationCall {
    /// 1-based position of the changed base, or of the base an insertion follows.
    pub(crate) position: usize,
    /// `0` for substitutions and deletions, `1..` for successive inserted bases.
    pub(crate) ordinal: usize,
    /// The rendered call, for example `309.1C`.
    pub(crate) text: String,
    /// The reference base; `None` for an inserted base.
    pub(crate) reference: Option<char>,
    /// The called base; `None` for a deleted base.
    pub(crate) alternate: Option<char>,
}

/// The reads containing one distinct call, and the report convention that
/// labels it.
pub(crate) struct CallReads<'a> {
    /// Content identities of the reads, in the notation's read order.
    pub(crate) reads: Vec<&'a str>,
    /// The convention every read containing the call labels it with
    /// (SRS-NOM-017); `None` when one read labels it otherwise or not at all.
    pub(crate) convention: Option<&'a Convention>,
}

/// The distinct calls of every read's represented variants, in reference
/// order, each with the reads containing it and its convention: the one
/// rendering `project` and the in-process consensus call build on.
pub(crate) fn calls<'a>(
    reference: &Reference,
    notation: &'a SampleNotation,
) -> dna_kernel::error::Result<BTreeMap<NotationCall, CallReads<'a>>> {
    let render = match notation.style {
        NotationStyle::PerBaseDecimal => render,
    };
    let mut calls = BTreeMap::<NotationCall, CallReads<'a>>::new();
    for representation in &notation.reads {
        let labels = Labels::new(&notation.conventions, reference, &representation.variants);
        for call in render(
            &reference.name,
            &reference.sequence,
            &representation.variants,
        )
        .map_err(ReportError::Representation)?
        {
            let convention = labels
                .of(Call {
                    position: call.position,
                    ordinal: call.ordinal,
                    reference: call.reference,
                    alternate: call.alternate,
                })
                .map_err(ReportError::Representation)?;
            let read = representation.input_sha256.as_str();
            match calls.entry(call) {
                std::collections::btree_map::Entry::Vacant(entry) => {
                    entry.insert(CallReads {
                        reads: vec![read],
                        convention,
                    });
                }
                std::collections::btree_map::Entry::Occupied(mut entry) => {
                    let entry = entry.get_mut();
                    entry.reads.push(read);
                    if entry.convention != convention {
                        entry.convention = None;
                    }
                }
            }
        }
    }
    Ok(calls)
}

/// Renders represented variants as per-base calls in reference order.
pub(super) fn render(
    contig: &str,
    reference: &str,
    variants: &[Variant],
) -> Result<Vec<NotationCall>, RepresentationError> {
    let mut edits = variants_to_edits(contig, reference, variants)?;
    sort_edits(&mut edits);
    let base = |index: usize| reference.as_bytes().get(index).copied().map(char::from);
    let mut calls = Vec::new();
    for edit in edits {
        let deleted = edit.end - edit.start;
        match (deleted, edit.alternate.len()) {
            (0, _) => {
                calls.extend(edit.alternate.chars().enumerate().map(|(index, inserted)| {
                    NotationCall {
                        position: edit.start,
                        ordinal: index + 1,
                        text: format!("{}.{}{inserted}", edit.start, index + 1),
                        reference: None,
                        alternate: Some(inserted),
                    }
                }));
            }
            (_, 0) => calls.extend((edit.start..edit.end).map(|index| NotationCall {
                position: index + 1,
                ordinal: 0,
                text: format!("{}DEL", index + 1),
                reference: base(index),
                alternate: None,
            })),
            (1, 1) => calls.push(NotationCall {
                position: edit.start + 1,
                ordinal: 0,
                text: format!("{}{}", edit.start + 1, edit.alternate),
                reference: base(edit.start),
                alternate: edit.alternate.chars().next(),
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
    fn renders_the_alleles_of_each_call() -> Result<(), RepresentationError> {
        let variants = [
            variant(4, "T", "TC", VariantKind::Ins),
            variant(2, "A", "C", VariantKind::Snv),
            variant(5, "AC", "A", VariantKind::Del),
        ];
        let alleles = render("synthetic", REFERENCE, &variants)?
            .into_iter()
            .map(|call| (call.text, call.reference, call.alternate))
            .collect::<Vec<_>>();
        assert_eq!(
            alleles,
            [
                ("2C".to_owned(), Some('A'), Some('C')),
                ("4.1C".to_owned(), None, Some('C')),
                ("6DEL".to_owned(), Some('C'), None),
            ]
        );
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
    fn labels_a_call_only_when_every_read_containing_it_does() -> dna_kernel::error::Result<()> {
        use dna_kernel::checksum::hex_sha256;
        use dna_kernel::model::reference::ReferenceTopology;
        use dna_kernel::profile::{ConventionKind, ConventionLabel};

        // `CCCCTCCC` at 1-8, then `A`: with `5C` the run is C8.
        let sequence = "CCCCTCCCAGG".to_owned();
        let reference = Reference {
            name: "synthetic".into(),
            sequence_sha256: hex_sha256(sequence.as_bytes()),
            sequence,
            topology: ReferenceTopology::Linear,
        };
        let joined = variant(5, "T", "C", VariantKind::Snv);
        let longer = variant(8, "C", "CC", VariantKind::Ins);
        let read = |name: &str, variants: Vec<Variant>| ReadRepresentation {
            input_sha256: name.into(),
            variants,
        };
        let notation = SampleNotation {
            style: NotationStyle::PerBaseDecimal,
            conventions: vec![Convention {
                name: "c_length".into(),
                label: ConventionLabel::NotRecorded,
                kind: ConventionKind::RunLength {
                    base: b'C',
                    minimum_run: 8,
                },
                window: [1, 8],
            }],
            // Read c has no `5C`: its `8.1C` lengthens a C3, not the C8.
            reads: vec![
                read("a", vec![joined.clone(), longer.clone()]),
                read("b", vec![joined, longer.clone()]),
                read("c", vec![longer]),
            ],
        };
        let labelled =
            |notation: &SampleNotation| -> dna_kernel::error::Result<Vec<(String, bool)>> {
                Ok(calls(&reference, notation)?
                    .into_iter()
                    .map(|(call, reads)| (call.text, reads.convention.is_some()))
                    .collect())
            };
        assert_eq!(
            labelled(&notation)?,
            [("5C".to_owned(), false), ("8.1C".to_owned(), false)]
        );
        let agreeing = SampleNotation {
            reads: notation.reads.into_iter().take(2).collect(),
            ..notation
        };
        assert_eq!(
            labelled(&agreeing)?,
            [("5C".to_owned(), false), ("8.1C".to_owned(), true)]
        );
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
