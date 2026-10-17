//! Profile-driven window representation engine.
//!
//! A profile declares reference windows, each with an ordered list of
//! representation rules; the first rule whose candidate reconstructs exactly
//! the same window haplotype wins, and a window no rule matches keeps its
//! normalized representation.

use std::ops::Range;

use crate::variant_representation::{SequenceEdit, apply_edits, sort_edits};
use dna_kernel::error::{NomenclatureError, Result};
use dna_kernel::model::nucleotide::is_canonical;
use dna_kernel::profile::{Anchor, NomenclatureWindow, WindowRule};

/// Represents every window; edits outside all windows pass through, and so
/// do the edits of a window that an edit crosses: that window keeps its
/// normalized form while the others are still represented. The reference
/// must carry every window sequence (checked again here so the engine never
/// rewrites edits against a different sequence).
pub(super) fn represent_windows(
    windows: &[NomenclatureWindow],
    reference: &str,
    normalized: &[SequenceEdit],
    readable: Option<&[Range<usize>]>,
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
            output.extend(inside);
            continue;
        }
        let local = shift_left(&inside, window.start);
        let (represented, rule) = represent(window, &local)?;
        let represented = represented
            .into_iter()
            .map(|edit| shifted(edit, window.start))
            .collect::<Vec<_>>();
        // A rule may rely only on bases the read reads: the bases its choice
        // of form rests on, and every base it calls. Otherwise the window
        // keeps its normalized form.
        let shown = readable.is_none_or(|readable| {
            let base = |index: usize| {
                readable
                    .iter()
                    .any(|range| range.contains(&(window.start + index)))
            };
            rule.is_none_or(|rule| {
                decisive(window, rule, &local)
                    .iter()
                    .all(|&index| base(index))
            }) && represented.iter().all(|edit| read_at(edit, readable))
        });
        output.extend(if shown { represented } else { inside });
    }
    output.extend(remaining);
    sort_edits(&mut output);
    Ok(output)
}

/// Whether the read reads every reference base `edit` calls: its substituted
/// or deleted bases, or an insertion's anchor base (none before position 1).
fn read_at(edit: &SequenceEdit, readable: &[Range<usize>]) -> bool {
    let (first, last) = if edit.start == edit.end {
        let Some(anchor) = edit.start.checked_sub(1) else {
            return true;
        };
        (anchor, anchor + 1)
    } else {
        (edit.start, edit.end)
    };
    (first..last).all(|base| readable.iter().any(|range| range.contains(&base)))
}

/// The window-local bases a rule's choice of form rests on, beyond the bases
/// it calls: the anchor and the run after it, plus the base before the anchor
/// when the left run's length changes (`anchored_run_lengths`); the end of the
/// lead run, the start of the left run, the anchor, and the run after it
/// (`anchored_tract`); the two terminal motif copies when the length
/// changes (`motif_shift`); and the whole window for `canonical_haplotype`,
/// whose choice compares the whole window with its declared haplotype.
fn decisive(window: &NomenclatureWindow, rule: &WindowRule, source: &[SequenceEdit]) -> Vec<usize> {
    let length = window.sequence.len();
    let alternate = apply_edits(&window.sequence, source).unwrap_or_default();
    match rule {
        WindowRule::AnchoredRunLengths(anchor) => {
            let mut bases = (anchor.index..length).collect::<Vec<_>>();
            if alternate
                .as_bytes()
                .iter()
                .position(|&base| base == anchor.base)
                != Some(anchor.index)
            {
                bases.push(anchor.index - 1);
            }
            bases
        }
        WindowRule::AnchoredTract { anchor, lead } => {
            let mut bases = vec![lead - 1, *lead];
            bases.extend(anchor.index..length);
            bases
        }
        WindowRule::AnchorDuplication(anchor)
        | WindowRule::AnchorDeletion(anchor)
        | WindowRule::TerminalRepeatInsertion(anchor) => (anchor.index..length).collect(),
        WindowRule::MotifShift { motif } => {
            if alternate.len() == length {
                Vec::new()
            } else {
                (length.saturating_sub(2 * motif.len())..length).collect()
            }
        }
        WindowRule::CanonicalHaplotype(_) => (0..length).collect(),
    }
}

/// Applies the first rule that reproduces the local haplotype, with that
/// rule, or keeps `source` when none does.
fn represent<'a>(
    window: &'a NomenclatureWindow,
    source: &[SequenceEdit],
) -> Result<(Vec<SequenceEdit>, Option<&'a WindowRule>)> {
    let reference = window.sequence.as_str();
    let alternate = apply_edits(reference, source).map_err(NomenclatureError::from)?;
    for rule in &window.rules {
        let Some(mut candidate) = candidate(reference, rule, source, &alternate) else {
            continue;
        };
        sort_edits(&mut candidate);
        return if apply_edits(reference, &candidate).ok().as_deref() == Some(alternate.as_str()) {
            Ok((candidate, Some(rule)))
        } else {
            Err(NomenclatureError::LocalHaplotypeChanged {
                window: window.name.clone(),
            }
            .into())
        };
    }
    Ok((source.to_vec(), None))
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
        WindowRule::AnchoredTract { anchor, lead } => {
            anchored_tract(reference, *anchor, *lead, alternate)
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
    {
        return None;
    }
    // The first anchor base splits the window; the left run is repeat bases.
    let alternate_anchor = bytes.iter().position(|&base| base == anchor.base)?;
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
    let right = anchor.index + 1;
    edits.extend(
        left_anchored(
            &reference[right..],
            &alternate[alternate_anchor + 1..],
            anchor.repeat,
        )?
        .into_iter()
        .map(|item| shifted(item, right)),
    );
    (!edits.is_empty()).then_some(edits)
}

fn anchored_tract(
    reference: &str,
    anchor: Anchor,
    lead: usize,
    alternate: &str,
) -> Option<Vec<SequenceEdit>> {
    let lead_base = reference.as_bytes()[0];
    let bytes = alternate.as_bytes();
    let alternate_lead = bytes.iter().take_while(|&&base| base == lead_base).count();
    let rest = &bytes[alternate_lead..];
    let right_run = reference.len() - anchor.index - 1;
    let left_run = anchor.index - lead;
    // The anchored form: the anchor sits after a left run and before a right
    // run of its reference length, every other base is a repeat or anchor
    // base, and the lead run did not grow.
    let anchored = rest.len() > right_run + 1
        && rest[rest.len() - right_run - 1] == anchor.base
        && rest
            .iter()
            .all(|&base| base == anchor.repeat || base == anchor.base)
        && alternate_lead <= lead;
    let whole = || left_anchored(reference, alternate, anchor.repeat);
    if !anchored {
        return whole();
    }
    let alternate_left = rest.len() - right_run - 1;
    let left = &rest[..alternate_left];
    if alternate_left > left_run && left[left_run..].iter().any(|&base| base != anchor.repeat) {
        return whole();
    }
    let mut edits = Vec::new();
    if alternate_lead < lead {
        edits.push(edit(alternate_lead, lead, ""));
    }
    for (offset, &base) in left.iter().take(left_run).enumerate() {
        if base != anchor.repeat {
            edits.push(substitution(lead + offset, base));
        }
    }
    if alternate_left > left_run {
        edits.push(insertion(
            anchor.index,
            anchor.repeat,
            alternate_left - left_run,
        ));
    } else if alternate_left < left_run {
        edits.push(edit(lead + alternate_left, anchor.index, ""));
    }
    let right = anchor.index + 1;
    for (offset, &base) in rest[alternate_left + 1..].iter().enumerate() {
        if base != reference.as_bytes()[right + offset] {
            edits.push(substitution(right + offset, base));
        }
    }
    (!edits.is_empty()).then_some(edits)
}

/// `alternate` written left-anchored against `reference`: position-wise
/// substitutions, then the extra bases as insertions after the last reference
/// base when they are all `repeat`, or the missing bases deleted from the 3'
/// end; `None` when the extra bases are not all `repeat`.
fn left_anchored(reference: &str, alternate: &str, repeat: u8) -> Option<Vec<SequenceEdit>> {
    let shared = reference.len().min(alternate.len());
    let mut edits = substitutions(&reference[..shared], &alternate[..shared]);
    if alternate.len() > reference.len() {
        let extra = &alternate.as_bytes()[reference.len()..];
        if extra.iter().any(|&base| base != repeat) {
            return None;
        }
        edits.push(insertion(reference.len(), repeat, extra.len()));
    } else if alternate.len() < reference.len() {
        edits.push(edit(alternate.len(), reference.len(), ""));
    }
    Some(edits)
}

fn shifted(item: SequenceEdit, offset: usize) -> SequenceEdit {
    SequenceEdit {
        start: item.start + offset,
        end: item.end + offset,
        alternate: item.alternate,
    }
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
    if !reference.ends_with(motif) || motif.is_empty() {
        return None;
    }
    if alternate.len() <= reference.len() {
        let lost = reference.len() - alternate.len();
        if !lost.is_multiple_of(motif.len()) {
            return None;
        }
        let kept = alternate.len();
        let mut candidate = substitutions(&reference[..kept], alternate);
        if lost > 0 {
            candidate.push(edit(kept, reference.len(), ""));
        }
        return Some(candidate);
    }
    let (head, extra) = alternate.split_at(reference.len());
    if !extra.len().is_multiple_of(motif.len()) || extra != motif.repeat(extra.len() / motif.len())
    {
        return None;
    }
    let mut candidate = substitutions(reference, head);
    candidate.push(edit(reference.len(), reference.len(), extra));
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

    use dna_kernel::error::{Error, NomenclatureError};
    use dna_kernel::model::reference::ReferenceTopology;
    use dna_kernel::profile::fixtures::human_mtdna;
    use dna_kernel::reference;

    use super::*;

    type TestResult = std::result::Result<(), Box<dyn std::error::Error>>;

    fn rcrs() -> Result<String> {
        Ok(reference::load(
            &Path::new(env!("CARGO_MANIFEST_DIR")).join("../../references/rCRS.fasta"),
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
        let output = represent_windows(windows, reference, normalized, None)?;
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
    fn the_phylogenetic_form_needs_the_whole_window_read() -> TestResult {
        let rcrs = rcrs()?;
        let source = [
            edit(16182, 16182, "C"),
            edit(16188, 16189, "C"),
            edit(16192, 16193, ""),
        ];
        // 16186 (0-based 16185) is unread: the window keeps its source edits.
        let gapped = [16_000..16_185, 16_186..16_300];
        assert_eq!(
            represent_windows(&human_mtdna()?.windows, &rcrs, &source, Some(&gapped))?,
            source
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
            represent_windows(&human_mtdna()?.windows, &reference, &[edit(72, 73, "G")], None),
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
    fn a_window_an_edit_crosses_keeps_its_form_and_the_others_are_represented() -> TestResult {
        let rcrs = rcrs()?;
        // A deletion of 511-513 crosses HVS-III; a C inserted at the start of
        // the left HVS-II run still converges on 309.1C.
        let crossing = edit(510, 513, "");
        let represented = represent_windows(
            &human_mtdna()?.windows,
            &rcrs,
            &[edit(303, 303, "C"), crossing.clone()],
            None,
        )?;
        assert_eq!(represented, [edit(309, 309, "C"), crossing]);
        Ok(())
    }

    /// The rendered `per_base_decimal` calls of window-local `edits`.
    fn calls(window: &NomenclatureWindow, edits: &[SequenceEdit]) -> Vec<String> {
        let mut calls = Vec::new();
        for item in edits {
            let first = window.start + item.start;
            if item.start == item.end {
                for (ordinal, base) in item.alternate.chars().enumerate() {
                    calls.push((
                        (first, ordinal + 1),
                        format!("{first}.{}{base}", ordinal + 1),
                    ));
                }
            } else if item.alternate.is_empty() {
                for position in first + 1..=window.start + item.end {
                    calls.push(((position, 0), format!("{position}DEL")));
                }
            } else {
                for (offset, base) in item.alternate.chars().enumerate() {
                    let position = first + offset + 1;
                    calls.push(((position, 0), format!("{position}{base}")));
                }
            }
        }
        calls.sort();
        calls.into_iter().map(|(_, call)| call).collect()
    }

    #[test]
    fn writes_the_corpus_sequences_in_the_lab_convention() -> TestResult {
        let profile = human_mtdna()?;
        let table = include_str!("../../tests/fixtures/lab-convention.tsv");
        let mut checked = 0;
        for line in table.lines().filter(|line| !line.starts_with('#')).skip(1) {
            let fields = line.split('\t').collect::<Vec<_>>();
            let [name, sequence, form, ..] = fields[..] else {
                return Err(format!("malformed row {line:?}").into());
            };
            let window = profile
                .windows
                .iter()
                .find(|window| window.name == name)
                .ok_or("unknown window")?;
            if sequence == window.sequence {
                assert_eq!(form, "reference");
                checked += 1;
                continue;
            }
            let source = [edit(0, window.sequence.len(), sequence)];
            let (represented, _) = represent(window, &source)?;
            assert_eq!(
                calls(window, &represented).join(" "),
                form,
                "{name} {sequence}"
            );
            checked += 1;
        }
        assert_eq!(checked, 90);
        Ok(())
    }
}
