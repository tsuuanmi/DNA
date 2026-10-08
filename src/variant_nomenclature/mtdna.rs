//! Human-mtDNA control-region representation policy.
//!
//! The policy is a list of declarative rCRS windows. Each window names its
//! repeat structure and an ordered list of representation rules; the first rule
//! whose candidate reconstructs exactly the same window haplotype wins, and a
//! window no rule matches keeps its normalized representation. Rules follow the
//! forensic conventions also implemented by the legacy `mtdna_raw` pipeline.

use std::path::Path;

use crate::error::{NomenclatureError, Result};
use crate::model::reference::{Reference, ReferenceTopology};
use crate::reference;
use crate::variant_representation::{
    SequenceEdit, apply_edits, render_edits, sort_edits, variants_to_edits,
};

use super::{NomenclatureInput, VariantNomenclatureResult};

/// SHA-256 of the normalized rCRS (`NC_012920.1`) sequence, the only reference
/// the human-mtDNA policies are validated against.
const RCRS_SEQUENCE_SHA256: &str =
    "f156ff3f65bbcc80c7ebb9936dceb96b1477b4f8f535c4e1dbe7baea225cbc66";

/// Repeated base of the anchored poly-C windows.
const REPEAT: u8 = b'C';
/// Anchor base splitting an anchored poly-C window into two runs.
const ANCHOR: u8 = b'T';

/// Repeat structure of one nomenclature window.
#[derive(Clone, Copy)]
enum Structure {
    /// Two `C` runs around one `T` anchor at a window-local index.
    AnchoredPolyC { anchor: usize },
    /// A tandem repeat whose window ends with one copy of `motif`.
    TandemRepeat { motif: &'static str },
}

/// Sequence-preserving representation rules, tried in window order.
#[derive(Clone, Copy)]
enum Rule {
    /// Express anchor movement as `C`-run length change at the run ends
    /// (`309.1C`, `315.1C`, `309DEL`).
    AnchoredRunLengths,
    /// A duplicated anchor is the base after it plus a terminal `C`
    /// (`311T 315.1C`).
    AnchorDuplication,
    /// A deleted anchor becomes a `C` with the last run base deleted
    /// (`310C 315DEL`, `16189C 16193DEL`).
    AnchorDeletion,
    /// A one-base gain ending in `C` keeps the phylogenetic terminal insertion
    /// (EMPOP `315.1C`) and describes the rest as substitutions.
    TerminalRepeatInsertion,
    /// A validated phylogenetic representation, as window-local edits.
    CanonicalHaplotype(&'static [(usize, usize, &'static str)]),
    /// Loss of one motif copy is the terminal copy's deletion plus
    /// position-wise substitutions (`513A 523DEL 524DEL`).
    MotifShift,
}

/// One validated rCRS window of the control-region policy.
struct Window {
    name: &'static str,
    start: usize,
    reference: &'static str,
    structure: Structure,
    rules: &'static [Rule],
}

/// Windows in reference order.
const WINDOWS: [Window; 3] = [
    Window {
        name: "HVS-II",
        start: 302,
        reference: "CCCCCCCTCCCCC",
        structure: Structure::AnchoredPolyC { anchor: 7 },
        rules: &[
            Rule::AnchoredRunLengths,
            Rule::AnchorDuplication,
            Rule::AnchorDeletion,
            Rule::TerminalRepeatInsertion,
        ],
    },
    Window {
        name: "HVS-III",
        start: 512,
        reference: "GCACACACACAC",
        structure: Structure::TandemRepeat { motif: "AC" },
        rules: &[Rule::MotifShift],
    },
    Window {
        name: "HVS-I",
        start: 16180,
        reference: "AAACCCCCTCCCC",
        structure: Structure::AnchoredPolyC { anchor: 8 },
        rules: &[
            Rule::CanonicalHaplotype(&[(2, 3, "C"), (3, 4, "A"), (8, 9, "C")]),
            Rule::AnchorDeletion,
        ],
    },
];

/// Whether `reference` is the rCRS sequence the human-mtDNA policies target.
pub(crate) fn is_rcrs(reference: &Reference) -> bool {
    reference.sequence_sha256 == RCRS_SEQUENCE_SHA256
}

/// Applies the validated human-mtDNA control-region representation policy.
///
/// Inside the HVS-II 303-315 and HVS-I 16181-16193 poly-C windows and the
/// HVS-III 513-524 AC repeat, sequence-equivalent haplotypes are expressed in
/// forensic notation (`309.1C 315.1C`, `513A 523DEL 524DEL`,
/// `16183C 16184A 16189C`). Normalized variants outside the windows, and window
/// haplotypes no rule represents, are retained unchanged.
///
/// # Errors
///
/// Returns [`Error`](crate::error::Error) when the reference cannot be loaded,
/// does not match the input reference identity or a validated window motif, an
/// edit crosses a window boundary, or the input variants do not reconstruct the
/// supplied haplotype.
pub fn apply_control_region(
    reference_path: &Path,
    input: NomenclatureInput<'_>,
) -> Result<VariantNomenclatureResult> {
    let reference = reference::load(reference_path, ReferenceTopology::Circular)?;
    control_region_with(&reference, input)
}

/// Applies the control-region policy against an already loaded reference; see
/// [`apply_control_region`].
pub(crate) fn control_region_with(
    reference: &Reference,
    input: NomenclatureInput<'_>,
) -> Result<VariantNomenclatureResult> {
    if input.reference.name != reference.name || input.reference.sha256 != reference.sequence_sha256
    {
        return Err(NomenclatureError::ReferenceIdentityMismatch.into());
    }

    let normalized_edits = variants_to_edits(
        &reference.name,
        &reference.sequence,
        input.normalized_variants,
    )
    .map_err(NomenclatureError::from)?;
    if apply_edits(&reference.sequence, &normalized_edits).map_err(NomenclatureError::from)?
        != input.alternate_sequence
    {
        return Err(NomenclatureError::InconsistentInput.into());
    }

    let represented_edits = represent_control_region(&reference.sequence, &normalized_edits)?;
    if apply_edits(&reference.sequence, &represented_edits).map_err(NomenclatureError::from)?
        != input.alternate_sequence
    {
        return Err(NomenclatureError::HaplotypeChanged.into());
    }

    let represented_variants =
        render_edits(&reference.name, &reference.sequence, &represented_edits)
            .map_err(NomenclatureError::from)?;

    Ok(VariantNomenclatureResult {
        reference: input.reference.clone(),
        source_variants: input.source_variants.to_vec(),
        alternate_sequence: input.alternate_sequence.to_owned(),
        normalized_variants: input.normalized_variants.to_vec(),
        represented_variants,
    })
}

/// Represents every window that lies inside `reference`; edits outside all
/// windows pass through. A window beyond the reference end cannot hold edits.
fn represent_control_region(
    reference: &str,
    normalized: &[SequenceEdit],
) -> Result<Vec<SequenceEdit>> {
    let mut remaining = normalized.to_vec();
    let mut output = Vec::with_capacity(normalized.len());
    for window in &WINDOWS {
        let end = window.start + window.reference.len();
        if end > reference.len() {
            continue;
        }
        if reference.get(window.start..end) != Some(window.reference) {
            return Err(NomenclatureError::ReferenceMotifMismatch {
                window: window.name,
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
                window: window.name,
            }
            .into());
        }
        output.extend(
            window
                .represent(&shift_left(&inside, window.start))?
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

impl Window {
    /// Applies the first rule that reproduces the local haplotype, or keeps
    /// `source` when none does.
    fn represent(&self, source: &[SequenceEdit]) -> Result<Vec<SequenceEdit>> {
        let alternate = apply_edits(self.reference, source).map_err(NomenclatureError::from)?;
        for &rule in self.rules {
            let Some(mut candidate) = self.candidate(rule, source, &alternate) else {
                continue;
            };
            sort_edits(&mut candidate);
            return if apply_edits(self.reference, &candidate).ok().as_deref()
                == Some(alternate.as_str())
            {
                Ok(candidate)
            } else {
                Err(NomenclatureError::LocalHaplotypeChanged { window: self.name }.into())
            };
        }
        Ok(source.to_vec())
    }

    fn candidate(
        &self,
        rule: Rule,
        source: &[SequenceEdit],
        alternate: &str,
    ) -> Option<Vec<SequenceEdit>> {
        let reference = self.reference;
        let candidate = match (rule, self.structure) {
            (Rule::AnchoredRunLengths, Structure::AnchoredPolyC { anchor }) => {
                anchored_run_lengths(reference, anchor, alternate)
            }
            (Rule::AnchorDuplication, Structure::AnchoredPolyC { anchor }) => {
                anchor_duplication(reference, anchor, alternate)
            }
            (Rule::AnchorDeletion, Structure::AnchoredPolyC { anchor }) => {
                anchor_deletion(reference, anchor, source)
            }
            (Rule::TerminalRepeatInsertion, Structure::AnchoredPolyC { .. }) => {
                terminal_repeat_insertion(reference, alternate, source.len())
            }
            (Rule::CanonicalHaplotype(edits), _) => Some(
                edits
                    .iter()
                    .map(|&(start, end, base)| edit(start, end, base))
                    .collect(),
            ),
            (Rule::MotifShift, Structure::TandemRepeat { motif }) => {
                motif_shift(reference, motif, alternate)
            }
            _ => None,
        }?;
        // A rule proposes; only a candidate that reproduces the haplotype counts.
        (apply_edits(reference, &candidate).ok().as_deref() == Some(alternate)).then_some(candidate)
    }
}

fn anchored_run_lengths(
    reference: &str,
    anchor: usize,
    alternate: &str,
) -> Option<Vec<SequenceEdit>> {
    if alternate
        .bytes()
        .any(|base| base != REPEAT && base != ANCHOR)
        || alternate.matches(char::from(ANCHOR)).count() != 1
    {
        return None;
    }
    let alternate_anchor = alternate.find(char::from(ANCHOR))?;
    let right_run = reference.len() - anchor - 1;
    let alternate_right_run = alternate.len() - alternate_anchor - 1;
    let mut edits = Vec::new();
    if alternate_anchor < anchor {
        edits.push(edit(alternate_anchor, anchor, ""));
    } else if alternate_anchor > anchor {
        edits.push(insertion(anchor, alternate_anchor - anchor));
    }
    if alternate_right_run < right_run {
        edits.push(edit(
            reference.len() - (right_run - alternate_right_run),
            reference.len(),
            "",
        ));
    } else if alternate_right_run > right_run {
        edits.push(insertion(reference.len(), alternate_right_run - right_run));
    }
    (!edits.is_empty()).then_some(edits)
}

fn anchor_duplication(
    reference: &str,
    anchor: usize,
    alternate: &str,
) -> Option<Vec<SequenceEdit>> {
    let (head, tail) = reference.split_at(anchor + 1);
    let expected = format!("{head}{}{tail}", char::from(ANCHOR));
    (alternate == expected && tail.as_bytes().first() == Some(&REPEAT)).then(|| {
        vec![
            edit(anchor + 1, anchor + 2, "T"),
            insertion(reference.len(), 1),
        ]
    })
}

fn anchor_deletion(
    reference: &str,
    anchor: usize,
    source: &[SequenceEdit],
) -> Option<Vec<SequenceEdit>> {
    let deletion = edit(anchor, anchor + 1, "");
    if !source.contains(&deletion) {
        return None;
    }
    let run_end = anchor
        + 1
        + reference.as_bytes()[anchor + 1..]
            .iter()
            .take_while(|&&base| base == REPEAT)
            .count();
    let terminal = run_end
        .checked_sub(1)
        .filter(|&terminal| terminal > anchor)?;
    let mut candidate: Vec<_> = source
        .iter()
        .filter(|&item| *item != deletion)
        .cloned()
        .collect();
    candidate.push(edit(anchor, anchor + 1, "C"));
    candidate.push(edit(terminal, terminal + 1, ""));
    Some(candidate)
}

fn terminal_repeat_insertion(
    reference: &str,
    alternate: &str,
    source_edits: usize,
) -> Option<Vec<SequenceEdit>> {
    let remainder = alternate.strip_suffix(char::from(REPEAT))?;
    if remainder.len() != reference.len()
        || !alternate
            .bytes()
            .all(|base| matches!(base, b'A' | b'C' | b'G' | b'T'))
    {
        return None;
    }
    let mut candidate = substitutions(reference, remainder);
    if candidate.len() > source_edits {
        return None;
    }
    candidate.push(insertion(reference.len(), 1));
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

fn insertion(at: usize, repeats: usize) -> SequenceEdit {
    edit(at, at, &char::from(REPEAT).to_string().repeat(repeats))
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

    /// Represents `normalized` and checks the haplotype is unchanged.
    fn represented(reference: &str, normalized: &[SequenceEdit]) -> Result<Vec<SequenceEdit>> {
        let output = represent_control_region(reference, normalized)?;
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
            represent_control_region(&reference, &[edit(72, 73, "G")]),
            Err(Error::VariantNomenclature(
                NomenclatureError::ReferenceMotifMismatch { window: "HVS-III" }
            ))
        ));
        Ok(())
    }
}
