//! Report conventions (SRS-NOM-017): which rendered calls of a read a
//! profile convention covers. A covered call is labelled, never dropped or
//! rewritten, so the haplotype stays complete.
//!
//! A rule matches only on the read's represented calls and the reference,
//! never on read evidence.

use std::cell::OnceCell;

use crate::variant_representation::{sort_edits, variants_to_edits};
use dna_kernel::error::RepresentationError;
use dna_kernel::model::reference::Reference;
use dna_kernel::profile::{Convention, ConventionKind};
use dna_kernel::variant::Variant;

/// One rendered per-base call of a read.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Call {
    /// 1-based position of the changed base, or of the base an insertion
    /// follows.
    pub position: usize,
    /// `0` for a substitution or deletion, `1..` for successive inserted bases.
    pub ordinal: usize,
    /// The reference base; `None` for an inserted base.
    pub reference: Option<char>,
    /// The called base; `None` for a deleted base.
    pub alternate: Option<char>,
}

/// The conventions of a profile over one read's represented calls.
#[derive(Debug)]
pub struct Labels<'a, 'r> {
    conventions: &'a [Convention],
    reference: &'r Reference,
    variants: &'r [Variant],
    /// The read's represented haplotype, built on first need.
    haplotype: OnceCell<Result<Haplotype, RepresentationError>>,
}

impl<'a, 'r> Labels<'a, 'r> {
    /// The conventions over the represented `variants` of one read.
    #[must_use]
    pub fn new(
        conventions: &'a [Convention],
        reference: &'r Reference,
        variants: &'r [Variant],
    ) -> Self {
        Self {
            conventions,
            reference,
            variants,
            haplotype: OnceCell::new(),
        }
    }

    /// The convention covering `call`, if any.
    ///
    /// # Errors
    ///
    /// Returns [`RepresentationError`] when the read's represented variants
    /// do not apply to the reference.
    pub fn of(&self, call: Call) -> Result<Option<&'a Convention>, RepresentationError> {
        for convention in self.conventions {
            let [start, end] = convention.window;
            if !(start..=end).contains(&call.position) {
                continue;
            }
            let ConventionKind::RunLength { base, minimum_run } = convention.kind;
            let base = char::from(base);
            let changed = match (call.reference, call.alternate) {
                (None, Some(inserted)) => inserted == base,
                (Some(deleted), None) => deleted == base,
                _ => false,
            };
            if !changed {
                continue;
            }
            let haplotype = self
                .haplotype
                .get_or_init(|| Haplotype::of(self.reference, self.variants))
                .as_ref()
                .map_err(Clone::clone)?;
            if haplotype.run_at(call, base as u8) >= minimum_run {
                return Ok(Some(convention));
            }
        }
        Ok(None)
    }
}

/// A represented haplotype, each base keyed by its place in reference order:
/// reference position `p` (1-based) is `(p, 0)`, and the `k`-th base
/// inserted after it is `(p, k)`.
#[derive(Debug)]
struct Haplotype {
    bases: Vec<(u8, (usize, usize))>,
}

impl Haplotype {
    fn of(reference: &Reference, variants: &[Variant]) -> Result<Self, RepresentationError> {
        let sequence = reference.sequence.as_bytes();
        let mut edits = variants_to_edits(&reference.name, &reference.sequence, variants)?;
        sort_edits(&mut edits);
        let mut bases = Vec::with_capacity(sequence.len());
        let mut cursor = 0;
        let keep = |bases: &mut Vec<_>, range: std::ops::Range<usize>| {
            bases.extend(range.map(|index| (sequence[index], (index + 1, 0))));
        };
        for edit in edits {
            if edit.start < cursor || edit.end > sequence.len() {
                return Err(RepresentationError::OverlappingEdits);
            }
            keep(&mut bases, cursor..edit.start);
            if edit.start == edit.end {
                bases.extend(
                    edit.alternate
                        .bytes()
                        .enumerate()
                        .map(|(index, inserted)| (inserted, (edit.start, index + 1))),
                );
            } else {
                // A substitution keeps its position; a deletion drops it.
                bases.extend(
                    edit.alternate
                        .bytes()
                        .zip(edit.start..edit.end)
                        .map(|(called, index)| (called, (index + 1, 0))),
                );
            }
            cursor = edit.end;
        }
        keep(&mut bases, cursor..sequence.len());
        Ok(Self { bases })
    }

    /// The length of the uninterrupted run of `base` that `call` changes,
    /// counting the call's own base: the run holding an inserted base, or
    /// the run across the junction a deleted base leaves, plus that base.
    fn run_at(&self, call: Call, base: u8) -> usize {
        let key = (call.position, call.ordinal);
        let (left, right, deleted) = if call.ordinal == 0 {
            let right = self.bases.partition_point(|(_, at)| *at <= key);
            (right.checked_sub(1), right, 1)
        } else {
            match self.bases.binary_search_by_key(&key, |(_, at)| *at) {
                Ok(index) => (Some(index), index + 1, 0),
                Err(_) => return 0,
            }
        };
        let before = left.map_or(0, |left| {
            self.bases[..=left]
                .iter()
                .rev()
                .take_while(|(called, _)| *called == base)
                .count()
        });
        let after = self.bases[right..]
            .iter()
            .take_while(|(called, _)| *called == base)
            .count();
        before + after + deleted
    }
}

#[cfg(test)]
mod tests {
    use dna_kernel::checksum::hex_sha256;
    use dna_kernel::model::reference::ReferenceTopology;
    use dna_kernel::profile::ConventionLabel;
    use dna_kernel::variant::VariantKind;

    use super::*;

    /// HVS-I 16180-16194 shape at 1-based 11-25: `AAAACCCCCTCCCCA`.
    fn reference() -> Reference {
        let sequence = format!("{}AAAACCCCCTCCCCA{}", "G".repeat(10), "G".repeat(10));
        Reference {
            name: "toy".into(),
            sequence_sha256: hex_sha256(sequence.as_bytes()),
            sequence,
            topology: ReferenceTopology::Linear,
        }
    }

    fn conventions() -> Vec<Convention> {
        vec![Convention {
            name: "c_length".into(),
            label: ConventionLabel::NotRecorded,
            kind: ConventionKind::RunLength {
                base: b'C',
                minimum_run: 10,
            },
            window: [11, 24],
        }]
    }

    fn variant(position: usize, reference: &str, alternate: &str, kind: VariantKind) -> Variant {
        Variant {
            contig: "toy".into(),
            position_1based: position,
            reference: reference.into(),
            alternate: alternate.into(),
            kind,
        }
    }

    /// `20C`: the anchor T becomes C and joins both runs into C10 at 15-24.
    fn joined() -> Variant {
        variant(20, "T", "C", VariantKind::Snv)
    }

    fn inserted(position: usize, ordinal: usize) -> Call {
        Call {
            position,
            ordinal,
            reference: None,
            alternate: Some('C'),
        }
    }

    fn deleted(position: usize) -> Call {
        Call {
            position,
            ordinal: 0,
            reference: Some('C'),
            alternate: None,
        }
    }

    fn names(variants: &[Variant], calls: &[Call]) -> Result<Vec<bool>, RepresentationError> {
        let reference = reference();
        let conventions = conventions();
        let labels = Labels::new(&conventions, &reference, variants);
        calls
            .iter()
            .map(|call| Ok(labels.of(*call)?.is_some()))
            .collect()
    }

    #[test]
    fn labels_a_length_call_of_a_long_enough_run() -> Result<(), RepresentationError> {
        // C10 plus one inserted C at 24 is C11; the substitution is content.
        let variants = [joined(), variant(24, "C", "CC", VariantKind::Ins)];
        let substitution = Call {
            position: 20,
            ordinal: 0,
            reference: Some('T'),
            alternate: Some('C'),
        };
        assert_eq!(
            names(&variants, &[substitution, inserted(24, 1)])?,
            [false, true]
        );
        Ok(())
    }

    #[test]
    fn counts_a_deleted_base_in_the_run_it_shortens() -> Result<(), RepresentationError> {
        // C10 less its last C leaves C9: the call changes a C10.
        let variants = [joined(), variant(23, "CC", "C", VariantKind::Del)];
        assert_eq!(names(&variants, &[deleted(24)])?, [true]);
        Ok(())
    }

    #[test]
    fn leaves_a_length_call_of_a_short_run_unlabelled() -> Result<(), RepresentationError> {
        // With the anchor T, the runs stay C5 and C4: lengths reads resolve.
        let insertion = [variant(19, "C", "CC", VariantKind::Ins)];
        assert_eq!(names(&insertion, &[inserted(19, 1)])?, [false]);
        let deletion = [variant(23, "CC", "C", VariantKind::Del)];
        assert_eq!(names(&deletion, &[deleted(24)])?, [false]);
        Ok(())
    }

    #[test]
    fn covers_only_calls_inside_the_window() -> Result<(), RepresentationError> {
        // An insertion anchored after the window is not covered.
        let variants = [joined(), variant(25, "A", "AC", VariantKind::Ins)];
        assert_eq!(names(&variants, &[inserted(25, 1)])?, [false]);
        Ok(())
    }
}
