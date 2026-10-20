//! 3'/right-most sequence-equivalent indel placement.

use super::TradeBounds;
use crate::variant_representation::{RepresentationResult, SequenceEdit, apply_edits, sort_edits};
use dna_kernel::model::placement;

pub(super) fn right_align(
    reference: &str,
    alternate_sequence: &str,
    source_edits: &[SequenceEdit],
    bounds: TradeBounds<'_>,
    max_trade: usize,
) -> RepresentationResult<Vec<SequenceEdit>> {
    let mut normalized = source_edits.to_vec();

    loop {
        let mut order = (0..normalized.len()).collect::<Vec<_>>();
        order.sort_by_key(|index| std::cmp::Reverse(normalized[*index].start));

        let mut moved = false;
        for index in order {
            // A shift whose edits overlap describes no haplotype; the deletion
            // may still trade past the substitution it runs into.
            let describes = |candidate: &[SequenceEdit]| {
                apply_edits(reference, candidate)
                    .is_ok_and(|sequence| sequence == alternate_sequence)
            };
            let shifted = shift_right_once(reference, &normalized, index)
                .map(|shifted| {
                    let mut candidate = normalized.clone();
                    candidate[index] = shifted;
                    candidate
                })
                .filter(|candidate| describes(candidate));
            let candidate = match shifted {
                Some(candidate) => Some(candidate),
                None => trade_deletion(reference, &normalized, index, bounds, max_trade)?
                    .filter(|candidate| describes(candidate)),
            };
            if let Some(candidate) = candidate {
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
        if !placement::crosses(first_inserted, crossed) {
            return None;
        }
        return Some(SequenceEdit {
            start: edit.start + 1,
            end: edit.end + 1,
            alternate: String::from_utf8(placement::rotated(edit.alternate.as_bytes(), crossed))
                .ok()?,
        });
    }

    if edit.alternate.is_empty() {
        let first_deleted = *reference_bytes.get(edit.start)?;
        let following = *reference_bytes.get(edit.end)?;
        if !placement::crosses(first_deleted, following) {
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

/// The deletion at `index` placed as far 3' as the same haplotype allows,
/// within `max_trade` bases,
/// without more edits: deleting reference bases further right and writing
/// the bases in between as substitutions, which may replace substitutions
/// the deletion passes (`16293DEL 16298C` becomes `16293C 16298DEL`, as the
/// lab's reviewers write it). Only substitutions may lie in the stretch it
/// crosses, every base of the stretch must be read, and the stretch may not
/// cross a boundary of a bounding window; `None` when no farther placement
/// keeps the edit count.
fn trade_deletion(
    reference: &str,
    edits: &[SequenceEdit],
    index: usize,
    bounds: TradeBounds<'_>,
    max_trade: usize,
) -> RepresentationResult<Option<Vec<SequenceEdit>>> {
    let deletion = &edits[index];
    if !deletion.alternate.is_empty() || deletion.end <= deletion.start {
        return Ok(None);
    }
    let (start, length) = (deletion.start, deletion.end - deletion.start);
    let mut best = None;
    for target in start + 1..=(start + max_trade).min(reference.len().saturating_sub(length)) {
        let end = target + length;
        // A longer stretch keeps any boundary or unread base of a shorter one.
        let split = bounds.windows.iter().any(|window| {
            window.start < end
                && start < window.end
                && !(window.start <= start && end <= window.end)
        });
        let unread = bounds.readable.is_some_and(|readable| {
            !(start..end).all(|base| readable.iter().any(|range| range.contains(&base)))
        });
        if split || unread {
            break;
        }
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

#[cfg(test)]
mod tests {
    use super::*;

    fn edit(start: usize, end: usize, alternate: &str) -> SequenceEdit {
        SequenceEdit {
            start,
            end,
            alternate: alternate.to_owned(),
        }
    }

    /// As 16293DEL 16297C 16298C over CCCACCCTT: the deletion of A at 3
    /// trades to 7, and 3 becomes C.
    const REFERENCE: &str = "CCCACCCTTAA";

    fn trade(bounds: TradeBounds<'_>) -> RepresentationResult<Vec<SequenceEdit>> {
        let source = [edit(3, 4, ""), edit(7, 8, "C"), edit(8, 9, "C")];
        let alternate = apply_edits(REFERENCE, &source)?;
        right_align(REFERENCE, &alternate, &source, bounds, 20)
    }

    // SRS-VN-008
    #[test]
    fn a_deletion_that_runs_into_a_substitution_trades_past_it() -> RepresentationResult<()> {
        // `16297DEL 16298C` over `TT`: the shift would delete the substituted
        // base, so the deletion trades instead and gives `16297C 16298DEL`.
        let reference = "CCCTTA";
        let source = [edit(3, 4, ""), edit(4, 5, "C")];
        let alternate = apply_edits(reference, &source)?;
        assert_eq!(
            right_align(reference, &alternate, &source, TradeBounds::default(), 20)?,
            [edit(3, 4, "C"), edit(4, 5, "")]
        );
        Ok(())
    }

    // SRS-VN-008
    #[test]
    fn a_deletion_trades_within_its_bounds() -> RepresentationResult<()> {
        let traded = vec![edit(3, 4, "C"), edit(7, 8, "C"), edit(8, 9, "")];
        assert_eq!(trade(TradeBounds::default())?, traded);
        assert_eq!(
            trade(TradeBounds {
                windows: &[0..9, 9..11],
                readable: Some(&[0..5, 5..11]),
            })?,
            traded
        );
        Ok(())
    }

    // SRS-VN-008
    #[test]
    fn a_deletion_trades_neither_across_a_window_boundary_nor_an_unread_base()
    -> RepresentationResult<()> {
        let kept = vec![edit(3, 4, ""), edit(7, 8, "C"), edit(8, 9, "C")];
        assert_eq!(
            trade(TradeBounds {
                windows: &[0..4, 4..11],
                readable: None,
            })?,
            kept
        );
        // Windows meeting at 8 keep the deletion at the first one's last base.
        assert_eq!(
            trade(TradeBounds {
                windows: &[0..8, 8..11],
                readable: None,
            })?,
            vec![edit(3, 4, "C"), edit(7, 8, ""), edit(8, 9, "C")]
        );
        assert_eq!(
            trade(TradeBounds {
                windows: &[],
                readable: Some(&[0..5, 6..11]),
            })?,
            kept
        );
        Ok(())
    }

    // SRS-VN-004
    #[test]
    fn every_accepted_movement_reconstructs_the_source_haplotype() -> RepresentationResult<()> {
        let source = [edit(3, 4, ""), edit(7, 8, "C"), edit(8, 9, "C")];
        let alternate = apply_edits(REFERENCE, &source)?;
        for bounds in [
            TradeBounds::default(),
            TradeBounds {
                windows: &[0..4, 4..11],
                readable: None,
            },
        ] {
            let normalized = trade(bounds)?;
            assert_eq!(apply_edits(REFERENCE, &normalized)?, alternate);
        }
        // A deletion that would run onto a substituted base is rejected as a
        // plain shift: its result would delete the substituted base instead.
        let reference = "CCCTTA";
        let source = [edit(3, 4, ""), edit(4, 5, "C")];
        let alternate = apply_edits(reference, &source)?;
        let plain = [edit(4, 5, ""), edit(4, 5, "C")];
        assert!(apply_edits(reference, &plain).map_or(true, |sequence| sequence != alternate));
        let normalized = right_align(reference, &alternate, &source, TradeBounds::default(), 20)?;
        assert_eq!(apply_edits(reference, &normalized)?, alternate);
        Ok(())
    }
}
