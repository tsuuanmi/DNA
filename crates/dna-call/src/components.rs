//! The sequences of a read's phase mixtures, each compared with the reference
//! on its own (ADR-0082, SRS-VAR-022).
//!
//! A phase mixture asserts that the read's signal after a run is its own
//! dominant sequence superposed with copies offset by whole calls. A
//! component offset by `d` is the molecule whose run is `d` calls longer
//! (shorter when negative): the read's calls with the run lengthened by its
//! own last `d` calls, or shortened by them. Each component's sequence, the
//! read's calls over its informative interval with that run length (the same
//! sequence on either side of the run), is placed and called as a read of its own, by the same aligner and
//! record rules, on the strand the read itself was placed on, and its observed
//! edits within the mixture's reference span are its own. The span covers the
//! tract and the stretch: the opening run, the runs of its base one other
//! call apart before it in read order (the HVS-II C-tract read from its 3'
//! end opens with `C9` after `C6 T`), and the calls after it. A component that
//! cannot be placed is left out and logged; the mixture and the sample stay.
//! A mixture its components do not explain (`PhaseMixture::explained`) is
//! not decomposed: its stretch is the fit's best guess, not sequences the
//! read holds, so it keeps its offsets and fractions and gets no calls.
//! Nothing here decides a call: the read's own placement and variants are
//! unchanged.

use std::ops::Range;

use dna_align::AlignmentConfig;
use dna_align::align_best;
use dna_align::model::alignment::AlignmentColumn;
use dna_align::model::alignment::Orientation;
use dna_kernel::error::Result;
use dna_kernel::model::nucleotide::Nucleotide;
use dna_kernel::model::nucleotide::reverse_complement;
use dna_kernel::model::reference::Reference;
use dna_kernel::read_evidence::{CallEvidence, EvidenceProfile, PhaseMixture, ReadEvidence};

use crate::config::VariantCallingConfig;
use crate::model::called_read::CalledRead;
use crate::model::variant::ObservedVariant;

/// One sequence of a phase mixture, compared with the reference.
#[derive(Debug, Clone)]
pub struct ComponentCall {
    /// Calls by which this sequence lags the read's dominant one.
    pub offset: i32,
    /// Its share of the fitted signal.
    pub fraction: f64,
    /// The reference interval its mixture's tract and stretch cover: its first and last
    /// positions in read order on the reference strand (1-based, inclusive);
    /// the first exceeds the last when the stretch crosses the origin.
    pub span: [usize; 2],
    /// Its bases over that interval, on the reference strand.
    pub sequence: String,
    /// Its edits within the interval and on the bases either side of it,
    /// which bound its runs, with the record rules' eligibility.
    pub variants: Vec<ObservedVariant>,
}

/// One phase mixture of a read with its components compared.
#[derive(Debug, Clone)]
pub struct MixtureCalls {
    /// The modality's assertion.
    pub mixture: PhaseMixture,
    /// Its components, in the assertion's order; none when the mixture is
    /// not explained.
    pub components: Vec<ComponentCall>,
}

/// The phase mixtures of `read`, each component of an explained one placed
/// on `reference` and called on its own; `regions` are the reportable
/// regions. A component that cannot be placed or called is left out of its
/// mixture and logged.
#[must_use]
pub fn component_calls(
    read: &CalledRead,
    reference: &Reference,
    alignment: &AlignmentConfig,
    calling: &VariantCallingConfig,
    regions: &[[usize; 2]],
) -> Vec<MixtureCalls> {
    let reverse = read.alignment.orientation() == Orientation::Reverse;
    read.evidence
        .mixtures()
        .iter()
        .map(|mixture| {
            if !mixture.explained {
                return MixtureCalls {
                    mixture: mixture.clone(),
                    components: Vec::new(),
                };
            }
            let components = mixture
                .components
                .iter()
                .filter_map(|component| {
                    let (sequence, stretch) =
                        component_sequence(&read.evidence, mixture, component.offset)?;
                    match called(
                        &sequence, stretch, reverse, reference, alignment, calling, regions,
                    ) {
                        Ok(Some((span, bases, variants))) => Some(ComponentCall {
                            offset: component.offset,
                            fraction: component.fraction,
                            span,
                            sequence: bases,
                            variants,
                        }),
                        Ok(None) => {
                            tracing::warn!(
                                event = "mixture_component_unplaced",
                                read = %read.input_name,
                                offset = component.offset,
                                reason = "no_reference_position",
                            );
                            None
                        }
                        Err(error) => {
                            tracing::warn!(
                                event = "mixture_component_unplaced",
                                read = %read.input_name,
                                offset = component.offset,
                                reason = "placement_failed",
                                error = %error,
                            );
                            None
                        }
                    }
                })
                .collect();
            MixtureCalls {
                mixture: mixture.clone(),
                components,
            }
        })
        .collect()
}

/// The read's calls over its informative interval, widened to cover
/// `mixture`, with the run that opens the mixture lengthened or shortened by
/// `offset` of its own last calls: the molecule whose run differs, with the
/// same sequence on either side. Also the calls of the mixture's tract and
/// stretch in it. `None` when the offset would remove more than the run.
fn component_sequence(
    evidence: &ReadEvidence,
    mixture: &PhaseMixture,
    offset: i32,
) -> Option<(String, Range<usize>)> {
    let calls = evidence.calls();
    let start = evidence.informative().start.min(mixture.call_start);
    let end = evidence.informative().end.max(mixture.call_end);
    let tract = tract_start(evidence, mixture.call_start, start);
    let run_end = mixture.call_start + mixture.run_calls;
    let base =
        |range: Range<usize>| -> String { calls[range].iter().map(|call| call.base).collect() };
    let shift = usize::try_from(offset.unsigned_abs()).ok()?;
    let run = base(mixture.call_start..run_end);
    let run = if offset >= 0 {
        let repeated = run.get(run.len().checked_sub(shift)?..)?.to_owned();
        format!("{run}{repeated}")
    } else {
        run.get(..run.len().checked_sub(shift)?)?.to_owned()
    };
    let stretch = tract - start..mixture.call_end - start + run.len() - mixture.run_calls;
    let sequence = format!(
        "{}{run}{}",
        base(start..mixture.call_start),
        base(run_end..end)
    );
    Some((sequence, stretch))
}

/// The calls of `mixture`'s tract and stretch in `read` (ADR-0082): from the
/// first call of the tract its opening run closes (see the span of
/// [`ComponentCall`]) to the stretch's end. The components show these calls'
/// sequence.
#[must_use]
pub fn mixture_calls(read: &ReadEvidence, mixture: &PhaseMixture) -> Range<usize> {
    let start = read.informative().start.min(mixture.call_start);
    tract_start(read, mixture.call_start, start)..mixture.call_end
}

/// The first call of the tract the run at `run_start` closes: that run and
/// the runs of its base, two calls or more, one other call apart before it
/// in read order (`C6 T` before `C9` read from the 3' end), no earlier than
/// `first`.
fn tract_start(evidence: &ReadEvidence, run_start: usize, first: usize) -> usize {
    let base = |index: usize| evidence.calls()[index].base;
    let run = base(run_start);
    let mut tract = run_start;
    while tract >= first + 3
        && base(tract - 1) != run
        && base(tract - 2) == run
        && base(tract - 3) == run
    {
        tract -= 2;
        while tract > first && base(tract - 1) == run {
            tract -= 1;
        }
    }
    tract
}

/// Places `sequence`, read-strand, and calls it; the span, bases on the
/// reference strand, and edits of its calls `stretch`. A `reverse` read's
/// component is placed reverse-complemented, on the strand the read was
/// placed on. `None` when the stretch places on no reference position.
#[allow(clippy::type_complexity)]
fn called(
    sequence: &str,
    stretch: Range<usize>,
    reverse: bool,
    reference: &Reference,
    alignment: &AlignmentConfig,
    calling: &VariantCallingConfig,
    regions: &[[usize; 2]],
) -> Result<Option<([usize; 2], String, Vec<ObservedVariant>)>> {
    let (stranded, stretch) = if reverse {
        let length = sequence.len();
        (
            reverse_complement(sequence),
            length - stretch.end..length - stretch.start,
        )
    } else {
        (sequence.to_owned(), stretch)
    };
    let calls = stranded
        .chars()
        .map(|base| {
            let profile = Nucleotide::from_char(base).map(|base| {
                let mut weights = [0.0; 4];
                weights[base.channel_index()] = 1.0;
                EvidenceProfile { weights }
            });
            CallEvidence::scored(base, profile)
        })
        .collect::<Vec<_>>();
    let length = calls.len();
    let evidence = ReadEvidence::new(calls, 0..length, Vec::new())?;
    let placed = align_best(&evidence, reference, alignment)?;
    if placed.orientation() != Orientation::Forward {
        return Ok(None);
    }
    let positions = placed
        .columns()
        .iter()
        .filter(|column| {
            column
                .original_call_index_0based()
                .is_some_and(|call| stretch.contains(&call))
        })
        .filter_map(AlignmentColumn::reference_index_0based)
        .collect::<Vec<_>>();
    let (Some(&first), Some(&last)) = (positions.first(), positions.last()) else {
        return Ok(None);
    };
    let span = [first + 1, last + 1];
    let bases = stranded[stretch].to_owned();
    let result = crate::call(&placed, reference, &evidence, calling, regions)?;
    // The bases either side of the interval bound its runs, and an insertion
    // or deletion at its first base is anchored on the base before it, so
    // their edits belong to the component too.
    let variants = result
        .observed
        .into_iter()
        .filter(|observed| {
            let position = observed.variant.position_1based;
            within(span, position) || position + 1 == span[0] || position == span[1] + 1
        })
        .collect();
    Ok(Some((span, bases, variants)))
}

/// Whether 1-based `position` lies in `span`, which crosses the origin when
/// its first position exceeds its last.
#[must_use]
pub fn within(span: [usize; 2], position: usize) -> bool {
    if span[0] <= span[1] {
        (span[0]..=span[1]).contains(&position)
    } else {
        position >= span[0] || position <= span[1]
    }
}

// SRS-VAR-022
#[cfg(test)]
mod tests {
    use dna_align::RawAlignmentConfig;
    use dna_kernel::model::reference::ReferenceTopology;
    use dna_kernel::read_evidence::MixtureComponent;

    use super::*;
    use crate::model::variant::VariantCallingResult;

    /// A reference with a run of seven C at 14-20 (1-based).
    const REFERENCE: &str = "TTGACGTCAGTACCCCCCCTGAGTACGAGGACTAGCATGCAAGTCGATGCA";
    /// The run's first call (0-based) and length.
    const RUN: (usize, usize) = (12, 7);

    fn reference() -> Reference {
        Reference {
            name: "ref".into(),
            sequence: REFERENCE.into(),
            topology: ReferenceTopology::Linear,
            sequence_sha256: String::new(),
        }
    }

    fn mixture(offsets: &[i32], explained: bool) -> PhaseMixture {
        PhaseMixture {
            call_start: RUN.0,
            call_end: REFERENCE.len(),
            run_calls: RUN.1,
            components: offsets
                .iter()
                .map(|&offset| MixtureComponent {
                    offset,
                    fraction: 1.0 / offsets.len() as f64,
                })
                .collect(),
            residual: 0.0,
            unexplained: Vec::new(),
            explained,
        }
    }

    fn read(mixture: PhaseMixture) -> Result<CalledRead> {
        let alignment = RawAlignmentConfig {
            match_score: 3,
            mismatch_score: -5,
            ambiguous_score: 0,
            gap_open_score: -10,
            gap_extension_score: -4,
            minimum_callable_bases: 1,
            minimum_identity: 0.8,
            max_stretch: 64,
        }
        .validate()?;
        let evidence = ReadEvidence::clean(REFERENCE).with_mixtures(vec![mixture])?;
        let placed = align_best(&evidence, &reference(), &alignment)?;
        Ok(CalledRead {
            input_name: "r".into(),
            input_sha256: "r".into(),
            reference_sha256: String::new(),
            configuration_sha256: String::new(),
            evidence,
            alignment: placed,
            variants: VariantCallingResult {
                reported: Vec::new(),
                observed: Vec::new(),
                excluded: Vec::new(),
            },
        })
    }

    fn calls(read: &CalledRead) -> Result<Vec<MixtureCalls>> {
        let alignment = RawAlignmentConfig {
            match_score: 3,
            mismatch_score: -5,
            ambiguous_score: 0,
            gap_open_score: -10,
            gap_extension_score: -4,
            minimum_callable_bases: 1,
            minimum_identity: 0.8,
            max_stretch: 64,
        }
        .validate()?;
        let calling = VariantCallingConfig {
            max_indel_length: 50,
            read_end_margin: 0,
        };
        Ok(component_calls(
            read,
            &reference(),
            &alignment,
            &calling,
            &[[1, REFERENCE.len()]],
        ))
    }

    #[test]
    fn a_tract_reaches_back_over_runs_of_its_base_one_call_apart() {
        // C6 T C9 read from the 3' end: the run of nine opens the mixture.
        let evidence = ReadEvidence::clean("GACCCCCCTCCCCCCCCCAG");
        assert_eq!(tract_start(&evidence, 9, 0), 2);
        // A single call of the base before the other call is no run.
        let single = ReadEvidence::clean("GACTCCCCCCCCCAG");
        assert_eq!(tract_start(&single, 4, 0), 4);
        // Never before the first call given.
        assert_eq!(tract_start(&evidence, 9, 5), 5);
    }

    #[test]
    fn a_span_may_cross_the_origin() {
        assert!(within([10, 20], 15));
        assert!(!within([10, 20], 21));
        assert!(within([16500, 40], 16569));
        assert!(within([16500, 40], 1));
        assert!(!within([16500, 40], 41));
    }

    #[test]
    fn a_component_lengthens_or_shortens_the_opening_run() {
        let evidence = ReadEvidence::clean(REFERENCE);
        let mixture = mixture(&[0], true);
        let (longer, stretch) = component_sequence(&evidence, &mixture, 2).unwrap_or_default();
        assert_eq!(stretch, RUN.0..REFERENCE.len() + 2);
        assert_eq!(&longer[RUN.0..RUN.0 + 10], "CCCCCCCCCT");
        assert_eq!(longer.len(), REFERENCE.len() + 2);
        let (shorter, _) = component_sequence(&evidence, &mixture, -3).unwrap_or_default();
        assert_eq!(&shorter[RUN.0..RUN.0 + 5], "CCCCT");
        assert_eq!(component_sequence(&evidence, &mixture, -8), None);
    }

    #[test]
    fn each_component_is_compared_with_the_reference_on_its_own() -> Result<()> {
        let found = calls(&read(mixture(&[0, 1], true))?)?;
        let components = &found[0].components;
        assert_eq!(components.len(), 2);
        // The dominant component is the reference; the lagging one has one C
        // more in the run.
        assert!(components[0].variants.is_empty());
        let inserted = components[1]
            .variants
            .iter()
            .map(|observed| {
                (
                    observed.variant.reference.clone(),
                    observed.variant.alternate.clone(),
                )
            })
            .collect::<Vec<_>>();
        assert_eq!(inserted.len(), 1);
        assert_eq!(inserted[0].1.len(), inserted[0].0.len() + 1);
        assert!(inserted[0].1.ends_with('C'));
        assert_eq!(components[0].span, [RUN.0 + 1, REFERENCE.len()]);
        Ok(())
    }

    #[test]
    fn an_unexplained_mixture_is_not_decomposed() -> Result<()> {
        let found = calls(&read(mixture(&[0, 1], false))?)?;
        assert_eq!(found[0].mixture.components.len(), 2);
        assert!(found[0].components.is_empty());
        Ok(())
    }
}
