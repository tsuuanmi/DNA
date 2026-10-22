//! Run lengths a profile convention leaves unrecorded (SRS-NOM-017) are not
//! compared: the reviewers do not write them, so a length DNA writes there and
//! the reference length the truth keeps are not a difference, and the rare
//! length a truth does write is not one either.
//!
//! The same rule applies to both sides, from each side's own alleles: an
//! insertion of the convention's base anchored in its window, or a deletion
//! of that base in it, is left out when the side's haplotype with its
//! substitutions only has an uninterrupted run of that base of at least
//! `minimum_run` bases at the call. Leaving the side's insertions and
//! deletions out of the run keeps an over-estimated length from carrying a
//! run across the minimum.

use std::collections::BTreeMap;

use dna_kernel::profile::{Convention, ConventionKind, ConventionLabel};

use crate::compare::{Allele, anchor};

/// `alleles` without those a convention of `conventions` leaves unrecorded,
/// on the 1-based `reference`.
pub(crate) fn recorded(
    alleles: &[Allele],
    reference: &str,
    conventions: &[Convention],
) -> Vec<Allele> {
    let reference = reference.as_bytes();
    // The side's substitutions, by 1-based position.
    let substituted = alleles
        .iter()
        .filter(|allele| {
            !allele.position.contains('.')
                && allele.reference.len() == 1
                && allele.alternate.len() == 1
                && allele.reference != "-"
                && allele.alternate != "-"
        })
        .filter_map(|allele| Some((anchor(&allele.position)?, allele.alternate.as_bytes()[0])))
        .collect::<BTreeMap<_, _>>();
    let base_at = |position: i64| -> Option<u8> {
        let index = usize::try_from(position.checked_sub(1)?).ok()?;
        substituted
            .get(&position)
            .copied()
            .or_else(|| reference.get(index).copied())
    };
    // The length of the run of `base` holding `position`.
    let run = |position: i64, base: u8| -> usize {
        if base_at(position) != Some(base) {
            return 0;
        }
        let (mut first, mut last) = (position, position);
        while base_at(first - 1) == Some(base) {
            first -= 1;
        }
        while base_at(last + 1) == Some(base) {
            last += 1;
        }
        usize::try_from(last - first + 1).unwrap_or(0)
    };
    let unrecorded = |allele: &Allele| {
        conventions.iter().any(|convention| {
            let ConventionKind::RunLength { base, minimum_run } = convention.kind;
            if convention.label != ConventionLabel::NotRecorded {
                return false;
            }
            let Some(position) = anchor(&allele.position) else {
                return false;
            };
            let window = i64::try_from(convention.window[0]).unwrap_or(i64::MAX)
                ..=i64::try_from(convention.window[1]).unwrap_or(i64::MIN);
            let of_base = |text: &str| !text.is_empty() && text.bytes().all(|byte| byte == base);
            let inserted = allele.reference == "-" && of_base(&allele.alternate);
            let deleted = allele.alternate == "-" && of_base(&allele.reference);
            // An insertion lies in the run its anchor or the next base starts.
            let length = if inserted {
                run(position, base).max(run(position + 1, base))
            } else if deleted {
                run(position, base)
            } else {
                0
            };
            (inserted || deleted) && window.contains(&position) && length >= minimum_run
        })
    };
    alleles
        .iter()
        .filter(|allele| !unrecorded(allele))
        .cloned()
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `AAAACCCCCTCCCCA` from 1, a C run of nine broken by a T at 10.
    const REFERENCE: &str = "AAAACCCCCTCCCCAG";

    fn convention() -> Convention {
        Convention {
            name: "c_length".into(),
            label: ConventionLabel::NotRecorded,
            kind: ConventionKind::RunLength {
                base: b'C',
                minimum_run: 10,
            },
            window: [1, 14],
        }
    }

    fn positions(alleles: &[Allele]) -> Vec<&str> {
        alleles
            .iter()
            .map(|allele| allele.position.as_str())
            .collect()
    }

    #[test]
    fn leaves_out_run_lengths_only_where_the_substituted_run_reaches_the_minimum() {
        let inserted = Allele::new("14.1", "-", "C");
        let deleted = Allele::new("14", "C", "-");
        let closing = Allele::new("10", "T", "C");
        // With T10C the run is ten C: its length is not compared, its
        // substitution is.
        let alleles = [closing.clone(), inserted.clone()];
        assert_eq!(
            positions(&recorded(&alleles, REFERENCE, &[convention()])),
            ["10"]
        );
        let alleles = [closing, deleted.clone()];
        assert_eq!(
            positions(&recorded(&alleles, REFERENCE, &[convention()])),
            ["10"]
        );
        // Without it the runs are five and four C: lengths are compared, even
        // when an insertion would carry the run past ten.
        let longer = [
            Allele::new("14.1", "-", "C"),
            Allele::new("14.2", "-", "C"),
            Allele::new("14.3", "-", "C"),
            Allele::new("14.4", "-", "C"),
            Allele::new("14.5", "-", "C"),
            Allele::new("14.6", "-", "C"),
        ];
        assert_eq!(
            recorded(&longer, REFERENCE, &[convention()]).len(),
            longer.len()
        );
        assert_eq!(
            positions(&recorded(&[deleted], REFERENCE, &[convention()])),
            ["14"]
        );
        // No convention, nothing left out.
        assert_eq!(positions(&recorded(&[inserted], REFERENCE, &[])), ["14.1"]);
    }
}
