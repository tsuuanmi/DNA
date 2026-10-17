//! 3'/right-most sequence-equivalent indel placement.

use crate::variant_representation::{RepresentationResult, SequenceEdit, apply_edits, sort_edits};

pub(super) fn right_align(
    reference: &str,
    alternate_sequence: &str,
    source_edits: &[SequenceEdit],
) -> RepresentationResult<Vec<SequenceEdit>> {
    let mut normalized = source_edits.to_vec();

    loop {
        let mut order = (0..normalized.len()).collect::<Vec<_>>();
        order.sort_by_key(|index| std::cmp::Reverse(normalized[*index].start));

        let mut moved = false;
        for index in order {
            let Some(shifted) = shift_right_once(reference, &normalized, index) else {
                continue;
            };
            let mut candidate = normalized.clone();
            candidate[index] = shifted;
            if apply_edits(reference, &candidate)? == alternate_sequence {
                normalized = candidate;
                moved = true;
                break;
            }
        }
        if !moved {
            break;
        }
    }

    sort_edits(&mut normalized);
    Ok(normalized)
}

/// The edit at `index` moved one base to the right, when the base it crosses
/// allows it. A deletion crosses a reference base equal to its first deleted
/// base. An insertion crosses a reference base equal to its first inserted
/// base; while it does not yet follow a reference base equal to its first
/// inserted base, it also crosses a base that a substitution makes equal to
/// it, so it reaches the 3' end of a reference run of its base (`16172C` takes
/// an inserted C from 16171 to 16174), yet one at the end of a run stays there
/// (`16262.1C 16263C`).
fn shift_right_once(reference: &str, edits: &[SequenceEdit], index: usize) -> Option<SequenceEdit> {
    let reference_bytes = reference.as_bytes();
    let edit = &edits[index];

    if edit.start == edit.end {
        let first_inserted = *edit.alternate.as_bytes().first()?;
        let attached = edit
            .start
            .checked_sub(1)
            .and_then(|before| reference_bytes.get(before))
            == Some(&first_inserted);
        let substituted = edits
            .iter()
            .find(|other| {
                other.start == edit.start
                    && other.end == edit.start + 1
                    && other.alternate.len() == 1
            })
            .and_then(|substitution| substitution.alternate.as_bytes().first().copied());
        let crossed = match substituted {
            Some(base) if !attached => base,
            _ => *reference_bytes.get(edit.start)?,
        };
        if first_inserted != crossed {
            return None;
        }
        let mut alternate = edit.alternate.as_bytes()[1..].to_vec();
        alternate.push(crossed);
        return Some(SequenceEdit {
            start: edit.start + 1,
            end: edit.end + 1,
            alternate: String::from_utf8(alternate).ok()?,
        });
    }

    if edit.alternate.is_empty() {
        let first_deleted = *reference_bytes.get(edit.start)?;
        let following = *reference_bytes.get(edit.end)?;
        if first_deleted != following {
            return None;
        }
        return Some(SequenceEdit {
            start: edit.start + 1,
            end: edit.end + 1,
            alternate: String::new(),
        });
    }

    None
}
