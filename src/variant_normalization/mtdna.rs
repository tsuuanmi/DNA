//! Human-mtDNA 3'/right-most sequence-equivalent indel placement.

use crate::error::Result;

use super::edit::{SequenceEdit, apply_edits, sort_edits};

pub(super) fn right_align(
    reference: &str,
    alternate_sequence: &str,
    source_edits: &[SequenceEdit],
) -> Result<Vec<SequenceEdit>> {
    let mut normalized = source_edits.to_vec();

    loop {
        let mut order = (0..normalized.len()).collect::<Vec<_>>();
        order.sort_by_key(|index| std::cmp::Reverse(normalized[*index].start));

        let mut moved = false;
        for index in order {
            let Some(shifted) = shift_right_once(reference, &normalized[index]) else {
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

fn shift_right_once(reference: &str, edit: &SequenceEdit) -> Option<SequenceEdit> {
    let reference_bytes = reference.as_bytes();

    if edit.start == edit.end {
        let first_inserted = *edit.alternate.as_bytes().first()?;
        let crossed = *reference_bytes.get(edit.start)?;
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
