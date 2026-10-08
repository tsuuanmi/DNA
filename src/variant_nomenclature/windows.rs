//! Profile-driven window representation engine.
//!
//! A profile declares reference windows, each with an ordered list of
//! representation rules; the first rule whose candidate reconstructs exactly
//! the same window haplotype wins, and a window no rule matches keeps its
//! normalized representation.

use crate::error::{NomenclatureError, Result};
use crate::model::nucleotide::is_canonical;
use crate::profile::{Anchor, NomenclatureWindow, WindowRule};
use crate::variant_representation::{SequenceEdit, apply_edits, sort_edits};

/// Represents every window; edits outside all windows pass through. The
/// reference must carry every window sequence (checked again here so the
/// engine never rewrites edits against a different sequence).
pub(super) fn represent_windows(
    windows: &[NomenclatureWindow],
    reference: &str,
    normalized: &[SequenceEdit],
) -> Result<Vec<SequenceEdit>> {
    let mut remaining = normalized.to_vec();
    let mut output = Vec::with_capacity(normalized.len());
    for window in windows {
        let end = window.end();
        if reference.get(window.start..end) != Some(window.sequence.as_str()) {
            return Err(NomenclatureError::ReferenceMotifMismatch {
                window: window.name.clone(),
            }
            .into());
        }
        let (inside, outside): (Vec<_>, Vec<_>) = remaining
            .into_iter()
            .partition(|edit| intersects(edit, window.start, end));
        remaining = outside;
        if inside.is_empty() {
            continue;
        }
        if !inside.iter().all(|edit| contained(edit, window.start, end)) {
            return Err(NomenclatureError::WindowCrossing {
                window: window.name.clone(),
            }
            .into());
        }
        output.extend(
            represent(window, &shift_left(&inside, window.start))?
                .into_iter()
                .map(|edit| SequenceEdit {
                    start: edit.start + window.start,
                    end: edit.end + window.start,
                    alternate: edit.alternate,
                }),
        );
    }
    output.extend(remaining);
    sort_edits(&mut output);
    Ok(output)
}

/// Applies the first rule that reproduces the local haplotype, or keeps
/// `source` when none does.
fn represent(window: &NomenclatureWindow, source: &[SequenceEdit]) -> Result<Vec<SequenceEdit>> {
    let reference = window.sequence.as_str();
    let alternate = apply_edits(reference, source).map_err(NomenclatureError::from)?;
    for rule in &window.rules {
        let Some(mut candidate) = candidate(reference, rule, source, &alternate) else {
            continue;
        };
        sort_edits(&mut candidate);
        return if apply_edits(reference, &candidate).ok().as_deref() == Some(alternate.as_str()) {
            Ok(candidate)
        } else {
            Err(NomenclatureError::LocalHaplotypeChanged {
                window: window.name.clone(),
            }
            .into())
        };
    }
    Ok(source.to_vec())
}

fn candidate(
    reference: &str,
    rule: &WindowRule,
    source: &[SequenceEdit],
    alternate: &str,
) -> Option<Vec<SequenceEdit>> {
    let candidate = match rule {
        WindowRule::AnchoredRunLengths(anchor) => {
            anchored_run_lengths(reference, *anchor, alternate)
        }
        WindowRule::AnchorDuplication(anchor) => anchor_duplication(reference, *anchor, alternate),
        WindowRule::AnchorDeletion(anchor) => anchor_deletion(reference, *anchor, source),
        WindowRule::TerminalRepeatInsertion(anchor) => {
            terminal_repeat_insertion(reference, anchor.repeat, alternate, source.len())
        }
        WindowRule::CanonicalHaplotype(substitutions) => Some(
            substitutions
                .iter()
                .map(|&(index, base)| substitution(index, base))
                .collect(),
        ),
        WindowRule::MotifShift { motif } => motif_shift(reference, motif, alternate),
    }?;
    // A rule proposes; only a candidate that reproduces the haplotype counts.
    (apply_edits(reference, &candidate).ok().as_deref() == Some(alternate)).then_some(candidate)
}

fn anchored_run_lengths(
    reference: &str,
    anchor: Anchor,
    alternate: &str,
) -> Option<Vec<SequenceEdit>> {
    let bytes = alternate.as_bytes();
    if bytes
        .iter()
        .any(|&base| base != anchor.repeat && base != anchor.base)
        || alternate.matches(char::from(anchor.base)).count() != 1
    {
        return None;
    }
    let alternate_anchor = bytes.iter().position(|&base| base == anchor.base)?;
    let right_run = reference.len() - anchor.index - 1;
    let alternate_right_run = alternate.len() - alternate_anchor - 1;
    let mut edits = Vec::new();
    if alternate_anchor < anchor.index {
        edits.push(edit(alternate_anchor, anchor.index, ""));
    } else if alternate_anchor > anchor.index {
        edits.push(insertion(
            anchor.index,
            anchor.repeat,
            alternate_anchor - anchor.index,
        ));
    }
    if alternate_right_run < right_run {
        edits.push(edit(
            reference.len() - (right_run - alternate_right_run),
            reference.len(),
            "",
        ));
    } else if alternate_right_run > right_run {
        edits.push(insertion(
            reference.len(),
            anchor.repeat,
            alternate_right_run - right_run,
        ));
    }
    (!edits.is_empty()).then_some(edits)
}

fn anchor_duplication(
    reference: &str,
    anchor: Anchor,
    alternate: &str,
) -> Option<Vec<SequenceEdit>> {
    let (head, tail) = reference.split_at(anchor.index + 1);
    let expected = format!("{head}{}{tail}", char::from(anchor.base));
    (alternate == expected && tail.as_bytes().first() == Some(&anchor.repeat)).then(|| {
        vec![
            substitution(anchor.index + 1, anchor.base),
            insertion(reference.len(), anchor.repeat, 1),
        ]
    })
}

fn anchor_deletion(
    reference: &str,
    anchor: Anchor,
    source: &[SequenceEdit],
) -> Option<Vec<SequenceEdit>> {
    let deletion = edit(anchor.index, anchor.index + 1, "");
    if !source.contains(&deletion) {
        return None;
    }
    let run_end = anchor.index
        + 1
        + reference.as_bytes()[anchor.index + 1..]
            .iter()
            .take_while(|&&base| base == anchor.repeat)
            .count();
    let terminal = run_end
        .checked_sub(1)
        .filter(|&terminal| terminal > anchor.index)?;
    let mut candidate: Vec<_> = source
        .iter()
        .filter(|&item| *item != deletion)
        .cloned()
        .collect();
    candidate.push(substitution(anchor.index, anchor.repeat));
    candidate.push(edit(terminal, terminal + 1, ""));
    Some(candidate)
}

fn terminal_repeat_insertion(
    reference: &str,
    repeat: u8,
    alternate: &str,
    source_edits: usize,
) -> Option<Vec<SequenceEdit>> {
    let remainder = alternate.strip_suffix(char::from(repeat))?;
    if remainder.len() != reference.len() || !alternate.chars().all(is_canonical) {
        return None;
    }
    let mut candidate = substitutions(reference, remainder);
    if candidate.len() > source_edits {
        return None;
    }
    candidate.push(insertion(reference.len(), repeat, 1));
    Some(candidate)
}

fn motif_shift(reference: &str, motif: &str, alternate: &str) -> Option<Vec<SequenceEdit>> {
    let terminal = reference.len().checked_sub(motif.len())?;
    if !reference.ends_with(motif) || alternate.len() != terminal {
        return None;
    }
    let mut candidate = substitutions(&reference[..terminal], alternate);
    candidate.push(edit(terminal, reference.len(), ""));
    Some(candidate)
}

/// Position-wise substitutions between two equal-length sequences.
fn substitutions(reference: &str, alternate: &str) -> Vec<SequenceEdit> {
    reference
        .char_indices()
        .zip(alternate.chars())
        .filter(|((_, expected), observed)| expected != observed)
        .map(|((index, _), observed)| SequenceEdit {
            start: index,
            end: index + 1,
            alternate: observed.to_string(),
        })
        .collect()
}

fn edit(start: usize, end: usize, alternate: &str) -> SequenceEdit {
    SequenceEdit {
        start,
        end,
        alternate: alternate.to_owned(),
    }
}

fn substitution(index: usize, base: u8) -> SequenceEdit {
    edit(index, index + 1, &char::from(base).to_string())
}

fn insertion(at: usize, base: u8, repeats: usize) -> SequenceEdit {
    edit(at, at, &char::from(base).to_string().repeat(repeats))
}

fn shift_left(edits: &[SequenceEdit], offset: usize) -> Vec<SequenceEdit> {
    edits
        .iter()
        .map(|item| SequenceEdit {
            start: item.start - offset,
            end: item.end - offset,
            alternate: item.alternate.clone(),
        })
        .collect()
}

/// An insertion belongs to a window when it lies after the window's first base
/// and no later than its end; other edits when they overlap it.
fn intersects(item: &SequenceEdit, start: usize, end: usize) -> bool {
    if item.start == item.end {
        item.start > start && item.start <= end
    } else {
        item.start < end && item.end > start
    }
}

fn contained(item: &SequenceEdit, start: usize, end: usize) -> bool {
    item.start == item.end || (item.start >= start && item.end <= end)
}

#[cfg(test)]
mod tests {
    use std::path::Path;

    use crate::error::{Error, NomenclatureError};
    use crate::model::reference::ReferenceTopology;
    use crate::profile::tests::human_mtdna;
    use crate::reference;

    use super::*;

    type TestResult = std::result::Result<(), Box<dyn std::error::Error>>;

    fn rcrs() -> Result<String> {
        Ok(reference::load(
            &Path::new(env!("CARGO_MANIFEST_DIR")).join("references/rCRS.fasta"),
            ReferenceTopology::Circular,
        )?
        .sequence)
    }

    fn edit(start: usize, end: usize, alternate: &str) -> SequenceEdit {
        SequenceEdit {
            start,
            end,
            alternate: alternate.into(),
        }
    }

    /// Represents `normalized` under the shipped human-mtDNA windows and checks
    /// the haplotype is unchanged.
    fn represented(reference: &str, normalized: &[SequenceEdit]) -> Result<Vec<SequenceEdit>> {
        represented_by(&human_mtdna()?.windows, reference, normalized)
    }

    fn represented_by(
        windows: &[NomenclatureWindow],
        reference: &str,
        normalized: &[SequenceEdit],
    ) -> Result<Vec<SequenceEdit>> {
        let output = represent_windows(windows, reference, normalized)?;
        assert_eq!(
            apply_edits(reference, &output).map_err(NomenclatureError::from)?,
            apply_edits(reference, normalized).map_err(NomenclatureError::from)?,
            "nomenclature must preserve the haplotype"
        );
        Ok(output)
    }

    #[test]
    fn hv2_duplicated_anchor_is_311t_315_1c() -> TestResult {
        let rcrs = rcrs()?;
        assert_eq!(
            represented(&rcrs, &[edit(310, 310, "T")])?,
            [edit(310, 311, "T"), edit(315, 315, "C")]
        );
        Ok(())
    }

    #[test]
    fn hv2_anchor_deletion_is_310c_315del() -> TestResult {
        let rcrs = rcrs()?;
        assert_eq!(
            represented(&rcrs, &[edit(309, 310, "")])?,
            [edit(309, 310, "C"), edit(314, 315, "")]
        );
        Ok(())
    }

    #[test]
    fn hv2_terminal_gain_keeps_empop_315_1c() -> TestResult {
        let rcrs = rcrs()?;
        for inserted in ["G", "T"] {
            assert_eq!(
                represented(&rcrs, &[edit(313, 313, inserted)])?,
                [edit(313, 314, inserted), edit(315, 315, "C")],
                "313.1{inserted} is 314{inserted} 315.1C"
            );
        }
        Ok(())
    }

    #[test]
    fn hv2_substitution_without_length_gain_does_not_add_315_1c() -> TestResult {
        let rcrs = rcrs()?;
        assert_eq!(
            represented(&rcrs, &[edit(313, 314, "G")])?,
            [edit(313, 314, "G")]
        );
        Ok(())
    }

    #[test]
    fn hv3_flank_loss_is_513a_523del_524del() -> TestResult {
        let rcrs = rcrs()?;
        assert_eq!(
            represented(&rcrs, &[edit(512, 514, "")])?,
            [edit(512, 513, "A"), edit(522, 524, "")]
        );
        Ok(())
    }

    #[test]
    fn hv3_motif_loss_keeps_an_independent_substitution() -> TestResult {
        let rcrs = rcrs()?;
        assert_eq!(
            represented(&rcrs, &[edit(520, 522, ""), edit(522, 523, "G")])?,
            [edit(520, 521, "G"), edit(522, 524, "")]
        );
        Ok(())
    }

    #[test]
    fn hv1_jumping_alignment_converges_on_16183c_16184a_16189c() -> TestResult {
        let rcrs = rcrs()?;
        assert_eq!(
            represented(
                &rcrs,
                &[
                    edit(16182, 16182, "C"),
                    edit(16188, 16189, "C"),
                    edit(16192, 16193, "")
                ],
            )?,
            [
                edit(16182, 16183, "C"),
                edit(16183, 16184, "A"),
                edit(16188, 16189, "C")
            ]
        );
        Ok(())
    }

    #[test]
    fn hv1_anchor_deletion_is_16189c_16193del() -> TestResult {
        let rcrs = rcrs()?;
        assert_eq!(
            represented(&rcrs, &[edit(16188, 16189, "")])?,
            [edit(16188, 16189, "C"), edit(16192, 16193, "")]
        );
        Ok(())
    }

    #[test]
    fn edits_outside_every_window_are_unchanged() -> TestResult {
        let rcrs = rcrs()?;
        let outside = [edit(72, 73, "G"), edit(16303, 16304, "C")];
        assert_eq!(represented(&rcrs, &outside)?, outside);
        Ok(())
    }

    #[test]
    fn rejects_a_reference_without_the_validated_window_motif() -> TestResult {
        let mut reference = rcrs()?;
        reference.replace_range(512..524, "GTGTGTGTGTGT");
        assert!(matches!(
            represent_windows(&human_mtdna()?.windows, &reference, &[edit(72, 73, "G")]),
            Err(Error::VariantNomenclature(
                NomenclatureError::ReferenceMotifMismatch { window }
            )) if window == "HVS-III"
        ));
        Ok(())
    }

    #[test]
    fn applies_rules_from_any_target_profile() -> TestResult {
        // A poly-G window around one A anchor, unrelated to mtDNA.
        let anchor = Anchor {
            index: 3,
            base: b'A',
            repeat: b'G',
        };
        let windows = [NomenclatureWindow {
            name: "poly-G".into(),
            start: 2,
            sequence: "GGGAGG".into(),
            rules: vec![
                WindowRule::AnchoredRunLengths(anchor),
                WindowRule::AnchorDeletion(anchor),
            ],
        }];
        let reference = "TTGGGAGGTT";
        // The anchor moved one base left: one G shorter left run, one G longer right run.
        assert_eq!(
            represented_by(&windows, reference, &[edit(4, 5, "A"), edit(5, 6, "G")])?,
            [edit(4, 5, ""), edit(8, 8, "G")]
        );
        // A deleted anchor becomes a repeat base with the last run base deleted.
        assert_eq!(
            represented_by(&windows, reference, &[edit(5, 6, "")])?,
            [edit(5, 6, "G"), edit(7, 8, "")]
        );
        Ok(())
    }

    #[test]
    fn rejects_an_edit_crossing_a_window_boundary() -> TestResult {
        let rcrs = rcrs()?;
        assert!(matches!(
            represent_windows(&human_mtdna()?.windows, &rcrs, &[edit(300, 304, "")]),
            Err(Error::VariantNomenclature(NomenclatureError::WindowCrossing { window }))
                if window == "HVS-II"
        ));
        Ok(())
    }
}
