//! The `per_base_decimal` notation style for represented variants.
//!
//! Rendering is mechanical serialization (ADR-0060 §8): each changed base gets
//! one call. A substitution is `<position><base>` (`73G`), each deleted base is
//! `<position>DEL` (`249DEL`), and each inserted base is
//! `<anchor>.<ordinal><base>` after the preceding reference base (`309.1C`,
//! `309.2C`); an insertion before the first base uses anchor `0`.

use crate::error::RepresentationError;
use crate::variant::Variant;
use crate::variant_representation::{sort_edits, variants_to_edits};

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
    use crate::variant::{Variant, VariantKind};

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
