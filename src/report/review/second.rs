//! Second peaks, the Sanger part of the review: a call's second peak judged
//! against its run's signal (the sample's own calls taken out) and the
//! configured floor, its cause, a difference whose call shows the reference
//! base as a second peak (`reference_second`), and a possible heteroplasmy
//! that both strands confirm, `explained` when the modality explains every
//! qualifying read's second peak.

use super::Sample;
use super::{Read, round};
use crate::error::{ReportError, Result, RunSignalError};
use crate::model::review_result::{LabelResult, ShareResult};
use crate::sample_analysis::{CallTrust, SampleAnalysis};
use dna_align::model::alignment::Orientation;
use dna_kernel::model::nucleotide::Nucleotide;
use dna_kernel::read_evidence::Placement;
use dna_sample::model::sample_evidence::LocusState;
use dna_sanger::run_signal::{self, Key, RunSignal, Threshold};

/// The run signal without the sample's own admitted reads, so the sample is
/// never judged against itself. A run signal that holds none of the sample's
/// traces is another run's.
pub(super) fn own(analysis: &SampleAnalysis, signal: &RunSignal) -> Result<RunSignal> {
    if !analysis
        .called
        .iter()
        .any(|read| signal.traces.contains(&read.input_sha256))
    {
        return Err(RunSignalError::OtherRun {
            run_id: signal.run_id.clone(),
        }
        .into());
    }
    let parameters = &analysis.loaded.config.sanger.signal_processing.run_signal;
    signal.check(parameters)?;
    let mut own = signal.clone();
    for read in &analysis.called {
        if !signal.traces.contains(&read.input_sha256) {
            continue;
        }
        let sanger = analysis
            .sanger
            .get(&read.input_sha256)
            .ok_or(ReportError::Inconsistent("review read lacks its evidence"))?;
        own = own.without(&run_signal::contributions(
            &sanger.heights(),
            &read.evidence,
            parameters,
        ));
    }
    Ok(own)
}

/// One read that reads a position cleanly: trusted or exempt, unmasked, with
/// a second base.
#[derive(Debug, Clone, Copy, PartialEq)]
struct Clean {
    /// Read index.
    read: usize,
    strand: Orientation,
    /// The call's second base on the reference strand, and its peak share.
    second: char,
    share: f64,
    /// The share the second base must exceed.
    threshold: Threshold,
    /// The bases the same read shows at the previous and the next reference
    /// position, when it shows one.
    neighbours: [Option<char>; 2],
    /// Whether the modality explains the call's second peak by the signal's
    /// own shape (SRS-CALL-022: not `second_unexplained`).
    explained: bool,
}

/// Whether a call reading another base than `reference` shows the
/// reference base as a second peak above the call's threshold
/// (SRS-CALL-020): `reference_second` when the difference rests on that
/// read's strand.
fn reference_second(call: &Clean, reference: char) -> bool {
    call.second == reference && call.share > call.threshold.share
}

/// Whether a `possible_heteroplasmy`'s qualifying reads leave one strand
/// without a second peak the modality leaves unexplained: every qualifying
/// read of that strand shows a neighbour's tail or a dye blob, not a possible
/// second sequence. The label is then context, not a check (on the
/// development panel, none of 23 such flags was a heteroplasmy).
fn explained(qualifying: &[Clean]) -> bool {
    [Orientation::Forward, Orientation::Reverse]
        .into_iter()
        .any(|strand| {
            qualifying
                .iter()
                .filter(|read| read.strand == strand)
                .all(|read| read.explained)
        })
}

/// `possible_heteroplasmy` over the clean reads of one position: a read
/// qualifies when its second peak share exceeds its threshold (SRS-CALL-020)
/// and its second base is not a base it shows next to the position (a
/// neighbour shadow). When both strands
/// read the position cleanly, a qualifying read on each strand must show the
/// same second base; when one strand does, one qualifying read suffices. The
/// base most qualifying reads show wins, then the first in A/C/G/T order. It
/// returns that base and its qualifying reads, in input order.
fn heteroplasmy(clean: &[Clean]) -> Option<(char, Vec<Clean>)> {
    let qualifying = clean
        .iter()
        .filter(|read| {
            read.share > read.threshold.share && !read.neighbours.contains(&Some(read.second))
        })
        .collect::<Vec<_>>();
    let strands = clean.iter().map(|read| read.strand).collect::<Vec<_>>();
    let both = strands.contains(&Orientation::Forward) && strands.contains(&Orientation::Reverse);
    let shows = |base: char, strand: Orientation| {
        qualifying
            .iter()
            .any(|read| read.second == base && read.strand == strand)
    };
    let mut best: Option<(usize, char)> = None;
    for base in Nucleotide::ALL.map(Nucleotide::as_char) {
        let count = qualifying.iter().filter(|read| read.second == base).count();
        // Both strands must confirm the base: a site only one strand reads
        // cleanly gives no label (the one-strand flags were 3 % precise).
        let fires = both && shows(base, Orientation::Forward) && shows(base, Orientation::Reverse);
        if fires && best.is_none_or(|(most, _)| count > most) {
            best = Some((count, base));
        }
    }
    let (_, base) = best?;
    Some((
        base,
        qualifying
            .into_iter()
            .filter(|read| read.second == base)
            .copied()
            .collect(),
    ))
}

impl Sample<'_> {
    /// Why call `call` of `read` shows a second peak: the modality's cause
    /// (`SangerCallSignal::second`), refined with the review's run signal,
    /// the sample's own reads taken out, so that a share its key's run level
    /// holds reads `context`. `None` when the call's second share does not
    /// reach `signal_processing.second_floor` and it is not masked as mixed.
    pub(super) fn second_cause(&self, read: &Read<'_>, call: usize) -> Option<&'static str> {
        let peak = read.sanger.calls.get(call)?.second.as_ref()?;
        let key = Key::of(
            read.evidence,
            call,
            peak.second.base,
            self.second.run_signal.flank,
        );
        Some(run_signal::refine(peak, key.as_ref(), self.signal.as_ref(), self.second).label())
    }

    /// The second base of a call at `position`, on the reference strand, its
    /// peak share (SRS-CALL-015), its threshold (SRS-CALL-020), and the bases
    /// the same read shows at the neighbouring positions.
    fn judged(
        &self,
        position: usize,
        read: usize,
        call: usize,
        trust: &CallTrust,
    ) -> Option<Clean> {
        let read_index = read;
        let read = self.reads.get(read)?;
        let second = trust.second?;
        let key = Key::of(
            read.evidence,
            call,
            second.base,
            self.second.run_signal.flank,
        );
        let neighbour = |other: Option<usize>| {
            other
                .and_then(|other| self.observation(read_index, other))
                .and_then(|other| other.base)
        };
        Some(Clean {
            read: read_index,
            strand: read.strand,
            second: read.strand.reference_base(second.base),
            share: second.share,
            threshold: run_signal::threshold(key.as_ref(), self.signal.as_ref(), self.second),
            neighbours: [
                neighbour(position.checked_sub(1)),
                neighbour(position.checked_add(1)),
            ],
            explained: read
                .evidence
                .calls()
                .get(call)
                .is_some_and(|evidence| !evidence.second_unexplained),
        })
    }

    /// `reference_second` at each position of `first..=last` where read
    /// `index`, one of the reads of the one strand a difference rests on,
    /// reads another base than the reference while its call shows the
    /// reference base as a second peak above the call's threshold
    /// (SRS-CALL-020). One strand cannot tell that base from a broad peak of
    /// another dye over the reference's, nor from a mixture, and its reruns
    /// repeat a dye blob at its primer's start; on development panel v5 most
    /// such bases at reverse read starts were dye blobs over the reference
    /// (SRS-CALL-021).
    pub(super) fn reference_second(
        &self,
        index: usize,
        first: usize,
        last: usize,
    ) -> Result<Vec<LabelResult>> {
        let read = self.read(index)?;
        let mut labels = Vec::new();
        for position in first..=last {
            let Some(observation) = self.observation(index, position) else {
                continue;
            };
            let reference = self.base(position);
            if observation.state != LocusState::Alternate
                || observation.base.is_none_or(|base| base == reference)
            {
                continue;
            }
            let Some(call) = observation.call_index_0based else {
                continue;
            };
            let trust = read.trust.get(call).ok_or(ReportError::Inconsistent(
                "review observation lacks the trust of its call",
            ))?;
            let Some(judged) = self.judged(position, index, call, trust) else {
                continue;
            };
            if reference_second(&judged, reference) {
                labels.push(LabelResult::ReferenceSecond {
                    position,
                    read: read.name.to_owned(),
                    strand: read.strand,
                    share: round(judged.share),
                    threshold: round(judged.threshold.share),
                    calls: judged.threshold.calls,
                });
            }
        }
        Ok(labels)
    }

    /// `possible_heteroplasmy` at one position, when it fires.
    pub(super) fn heteroplasmy(&self, position: usize) -> Option<LabelResult> {
        let locus = self.loci.get(&position)?;
        let clean = locus
            .observations
            .iter()
            .filter_map(|observation| {
                if !matches!(
                    observation.state,
                    LocusState::Reference | LocusState::Alternate
                ) {
                    return None;
                }
                let read = self.reads.get(observation.read_index)?;
                let trust = read.trust.get(observation.call_index_0based?)?;
                if trust.placement != Placement::Scored || !(trust.trusted || trust.exempt) {
                    return None;
                }
                self.judged(
                    position,
                    observation.read_index,
                    observation.call_index_0based?,
                    trust,
                )
            })
            .collect::<Vec<_>>();
        let (base, reads) = heteroplasmy(&clean)?;
        Some(LabelResult::PossibleHeteroplasmy {
            position,
            base,
            explained: explained(&reads),
            reads: reads
                .into_iter()
                .map(|clean| ShareResult {
                    read: self.reads[clean.read].name.to_owned(),
                    strand: clean.strand,
                    share: round(clean.share),
                    threshold: round(clean.threshold.share),
                    calls: clean.threshold.calls,
                })
                .collect(),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::report::review::test_support::{F, R};

    fn clean(read: usize, strand: Orientation, second: char, share: f64) -> Clean {
        Clean {
            read,
            strand,
            second,
            share,
            threshold: Threshold {
                share: 0.12,
                calls: None,
            },
            neighbours: [None, None],
            explained: false,
        }
    }

    #[test]
    fn heteroplasmy_on_both_strands_needs_the_same_base_on_each() {
        let reads = [clean(0, F, 'T', 0.2), clean(1, R, 'T', 0.15)];
        let (base, qualifying) = heteroplasmy(&reads).unwrap_or(('-', Vec::new()));
        assert_eq!(base, 'T');
        assert_eq!(qualifying, reads);

        let other = [clean(0, F, 'T', 0.2), clean(1, R, 'G', 0.3)];
        assert_eq!(heteroplasmy(&other), None);
        let one_side = [clean(0, F, 'T', 0.2), clean(1, R, 'T', 0.05)];
        assert_eq!(heteroplasmy(&one_side), None);
    }

    #[test]
    fn heteroplasmy_needs_both_strands_to_confirm() {
        // One strand alone never fires, however strong its share.
        assert_eq!(heteroplasmy(&[clean(0, R, 'A', 0.4)]), None);
        let one_strand = [clean(0, F, 'T', 0.3), clean(1, F, 'T', 0.3)];
        assert_eq!(heteroplasmy(&one_strand), None);
        // Both strands with the same second base fire, and the reads that
        // show it are reported.
        let both = [
            clean(0, F, 'T', 0.3),
            clean(1, R, 'T', 0.25),
            clean(2, R, 'C', 0.01),
        ];
        let (base, qualifying) = heteroplasmy(&both).unwrap_or(('-', Vec::new()));
        assert_eq!(base, 'T');
        assert_eq!(qualifying, [both[0], both[1]]);
        // The base most reads show wins, then A/C/G/T order.
        let tied = [
            clean(0, F, 'T', 0.2),
            clean(1, R, 'T', 0.2),
            clean(2, F, 'C', 0.2),
            clean(3, R, 'C', 0.2),
            clean(4, R, 'T', 0.2),
        ];
        assert_eq!(heteroplasmy(&tied).map(|(base, _)| base), Some('T'));
        let even = [
            clean(0, F, 'T', 0.2),
            clean(1, R, 'T', 0.2),
            clean(2, F, 'C', 0.2),
            clean(3, R, 'C', 0.2),
        ];
        assert_eq!(heteroplasmy(&even).map(|(base, _)| base), Some('C'));
    }

    #[test]
    fn a_reference_second_needs_the_reference_base_above_the_threshold() {
        // One reverse read reads A where the reference is C, a C under it.
        assert!(reference_second(&clean(0, R, 'C', 0.39), 'C'));
        // Another second base, or one at the threshold, is no such label;
        // neighbours do not matter: a blob covers them too.
        assert!(!reference_second(&clean(0, R, 'G', 0.39), 'C'));
        assert!(!reference_second(&clean(0, R, 'C', 0.12), 'C'));
        let mut beside = clean(0, F, 'C', 0.4);
        beside.neighbours = [Some('C'), None];
        assert!(reference_second(&beside, 'C'));
        // A run level above the share holds it back.
        let mut raised = clean(0, R, 'C', 0.3);
        raised.threshold = Threshold {
            share: 0.35,
            calls: Some(40),
        };
        assert!(!reference_second(&raised, 'C'));
    }

    #[test]
    fn heteroplasmy_needs_a_share_above_the_threshold() {
        let pair = |share| [clean(0, F, 'G', share), clean(1, R, 'G', share)];
        assert_eq!(heteroplasmy(&pair(0.12)), None);
        assert!(heteroplasmy(&pair(0.121)).is_some());
        // A run level above the share on one strand holds the site back.
        let mut raised = pair(0.3);
        raised[0].threshold = Threshold {
            share: 0.35,
            calls: Some(40),
        };
        assert_eq!(heteroplasmy(&raised), None);
        assert_eq!(heteroplasmy(&[]), None);
    }

    #[test]
    fn heteroplasmy_excludes_a_neighbour_shadow() {
        let mut pair = [clean(0, F, 'G', 0.4), clean(1, R, 'G', 0.4)];
        pair[0].neighbours = [Some('A'), Some('G')];
        assert_eq!(heteroplasmy(&pair), None);
        pair[0].neighbours = [Some('G'), None];
        assert_eq!(heteroplasmy(&pair), None);
        pair[0].neighbours = [Some('A'), Some('C')];
        assert!(heteroplasmy(&pair).is_some());
    }

    #[test]
    fn a_strand_whose_second_peaks_the_modality_explains_makes_heteroplasmy_context() {
        let shaped = |mut read: Clean| {
            read.explained = true;
            read
        };
        // Centred, unexplained second peaks on both strands: a check.
        assert!(!explained(&[clean(0, F, 'G', 0.3), clean(1, R, 'G', 0.3)]));
        // One strand's only qualifying read is a neighbour's tail or a blob.
        assert!(explained(&[
            clean(0, F, 'G', 0.3),
            shaped(clean(1, R, 'G', 0.3))
        ]));
        // Another read on that strand still shows an unexplained one.
        assert!(!explained(&[
            clean(0, F, 'G', 0.3),
            shaped(clean(1, R, 'G', 0.3)),
            clean(2, R, 'G', 0.25),
        ]));
        let label = |explained| LabelResult::PossibleHeteroplasmy {
            position: 1,
            base: 'G',
            explained,
            reads: Vec::new(),
        };
        assert!(label(false).unassured() && !label(true).unassured());
        assert_eq!((label(false).rank(), label(true).rank()), (1, 4));
    }
}
