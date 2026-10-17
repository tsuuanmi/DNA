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
            let candidate = match shift_right_once(reference, &normalized, index) {
                Some(shifted) => {
                    let mut candidate = normalized.clone();
                    candidate[index] = shifted;
                    candidate
                }
                None => match trade_deletion(reference, &normalized, index)? {
                    Some(candidate) => candidate,
                    None => continue,
                },
            };
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

/// The farthest a deletion is moved past substitutions in one step.
const MAX_TRADE: usize = 20;

/// The deletion at `index` placed as far 3' as the same haplotype allows
/// without more edits: deleting reference bases further right and writing
/// the bases in between as substitutions, which may replace substitutions
/// the deletion passes (`16293DEL 16298C` becomes `16293C 16298DEL`, as the
/// lab's reviewers write it). Only substitutions may lie in the stretch it
/// crosses; `None` when no farther placement keeps the edit count.
fn trade_deletion(
    reference: &str,
    edits: &[SequenceEdit],
    index: usize,
) -> RepresentationResult<Option<Vec<SequenceEdit>>> {
    let deletion = &edits[index];
    if !deletion.alternate.is_empty() || deletion.end <= deletion.start {
        return Ok(None);
    }
    let (start, length) = (deletion.start, deletion.end - deletion.start);
    let mut best = None;
    for target in start + 1..=(start + MAX_TRADE).min(reference.len().saturating_sub(length)) {
        let end = target + length;
        let inside = |edit: &SequenceEdit| {
            edit.start < end && (edit.end > start || (edit.start == edit.end && edit.start > start))
        };
        let crossed = edits
            .iter()
            .enumerate()
            .filter(|(other, edit)| *other != index && inside(edit))
            .map(|(_, edit)| edit)
            .collect::<Vec<_>>();
        if crossed
            .iter()
            .any(|edit| edit.end != edit.start + 1 || edit.alternate.len() != 1)
        {
            break;
        }
        let local = std::iter::once(deletion)
            .chain(crossed.iter().copied())
            .map(|edit| SequenceEdit {
                start: edit.start - start,
                end: edit.end - start,
                alternate: edit.alternate.clone(),
            })
            .collect::<Vec<_>>();
        let shown = apply_edits(&reference[start..end], &local)?;
        let substitutions = shown
            .bytes()
            .zip(reference[start..target].bytes())
            .enumerate()
            .filter(|(_, (observed, expected))| observed != expected)
            .map(|(offset, (observed, _))| SequenceEdit {
                start: start + offset,
                end: start + offset + 1,
                alternate: char::from(observed).to_string(),
            })
            .collect::<Vec<_>>();
        let mut candidate = edits
            .iter()
            .enumerate()
            .filter(|(other, edit)| *other != index && !inside(edit))
            .map(|(_, edit)| edit.clone())
            .collect::<Vec<_>>();
        candidate.extend(substitutions);
        candidate.push(SequenceEdit {
            start: target,
            end,
            alternate: String::new(),
        });
        if candidate.len() <= edits.len() {
            best = Some(candidate);
        }
    }
    Ok(best)
}
