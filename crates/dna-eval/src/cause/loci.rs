//! What the admitted reads show at each position of the profile's regions,
//! as compact flags, for the cause of a sample's discrepancy.

use std::collections::BTreeSet;

use dna::sample_analysis::SampleAnalysis;
use dna::sample_analysis::SangerAttachment;
use dna_core::model::alignment::Orientation;
use dna_core::model::consensus::SiteState;
use dna_core::model::sample_evidence::LocusState;
use dna_kernel::model::nucleotide::Nucleotide;
use dna_sanger::model::callability::{PhaseSegment, PhaseState};

/// A second channel at or above this share of a read's profile, in
/// thousandths, is a faint second base; below it, none (3 %).
pub(crate) const FAINT: u16 = 30;

/// Some admitted read covers the position.
pub(crate) const COVERED: u16 = 1;
/// A forward read observes it cleanly: a base or a gap, unmasked.
pub(crate) const FORWARD: u16 = 2;
/// A reverse read observes it cleanly.
pub(crate) const REVERSE: u16 = 4;
/// A clean read observes something other than the reference there.
pub(crate) const OTHER: u16 = 8;
/// A read covers it with a call masked for phase loss: dephased or mixed
/// signal, or a segment starting after a repeat run.
pub(crate) const PHASE: u16 = 16;
/// A read covers it with a call masked for weak signal or irregular spacing.
pub(crate) const SIGNAL: u16 = 32;
/// A read covers it with an unresolved call (`N`).
pub(crate) const UNRESOLVED: u16 = 64;
/// The consensus decided the position by a tie of clean reads and wrote `N`.
pub(crate) const CONTESTED: u16 = 128;
/// The position lies in a stretch of whole runs the consensus left undecided:
/// no read shows both ends of the run.
pub(crate) const RUN_UNREAD: u16 = 256;
/// A clean read observes it within the read-end margin of its own first or
/// last callable call.
pub(crate) const EDGE: u16 = 512;
/// A clean read observes it away from its own ends.
pub(crate) const INSIDE: u16 = 1024;
/// A clean read observes it within the read-end margin of a call masked for
/// phase loss, or of its own end where it stops at a repeat run.
pub(crate) const BESIDE_PHASE_LOSS: u16 = 2048;
/// A read's call there carries a weak-signal support veto: its highest peak
/// below the minimum, or its relative quality not above the threshold.
pub(crate) const WEAK_CALL: u16 = 4096;
/// A read's call there has more than one qualifying channel.
pub(crate) const MIXED_CALL: u16 = 8192;

/// What the harness reads from DNA's configuration to place a read's trust.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct Margins {
    /// `variant_calling.read_end_margin`: calls within it of an
    /// uninformative call or a read end do not support a variant.
    pub(crate) read_end: usize,
    /// `callability.repeat_min_length`: the shortest repeat run.
    pub(crate) repeat: usize,
    /// `sanger_evidence.minimum_peak_height`.
    pub(crate) minimum_peak: usize,
    /// `sanger_evidence.relative_quality_threshold`.
    pub(crate) relative_quality: usize,
}

/// The flags of every region position, and the faint second bases.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub(crate) struct Loci {
    /// The regions, 1-based inclusive, with the flags of each position.
    regions: Vec<([usize; 2], Vec<u16>)>,
    /// Positions and bases at which some clean read has a channel, other
    /// than its own call, at [`FAINT`] or more.
    pub(crate) faint: BTreeSet<(i64, char)>,
}

impl Loci {
    /// The flags at a 1-based position; 0 outside the regions.
    pub(crate) fn at(&self, position: i64) -> u16 {
        let Ok(position) = usize::try_from(position) else {
            return 0;
        };
        self.regions
            .iter()
            .find(|([start, end], _)| (*start..=*end).contains(&position))
            .and_then(|([start, _], flags)| flags.get(position - start).copied())
            .unwrap_or(0)
    }

    /// Loci over one region from `start` with the given flags, for tests.
    #[cfg(test)]
    pub(crate) fn of(start: usize, flags: Vec<u16>, faint: &[(i64, char)]) -> Self {
        let end = start + flags.len() - 1;
        Self {
            regions: vec![([start, end], flags)],
            faint: faint.iter().copied().collect(),
        }
    }

    fn flags(&mut self, position: usize) -> Option<&mut u16> {
        self.regions.iter_mut().find_map(|([start, end], flags)| {
            (*start..=*end)
                .contains(&position)
                .then(|| flags.get_mut(position - *start))
                .flatten()
        })
    }

    /// Whether clean reads on both strands observe the position and every
    /// clean read shows the reference.
    pub(crate) fn reference_on_both_strands(&self, position: i64) -> bool {
        let flags = self.at(position);
        flags & FORWARD != 0 && flags & REVERSE != 0 && flags & OTHER == 0
    }
}

/// Takes the flags of every region of `analysis`, placing a clean call's
/// trust with `margins`.
///
/// # Errors
///
/// Returns DNA's error when the observations cannot be read.
pub(crate) fn extract(analysis: &SampleAnalysis, margins: Margins) -> dna::error::Result<Loci> {
    let regions = analysis.regions();
    let mut loci = Loci {
        regions: regions
            .iter()
            .map(|&[start, end]| ([start, end], vec![0; end.saturating_sub(start) + 1]))
            .collect(),
        faint: BTreeSet::new(),
    };
    let (Some(first), Some(last)) = (
        regions.iter().map(|[start, _]| *start).min(),
        regions.iter().map(|[_, end]| *end).max(),
    ) else {
        return Ok(loci);
    };
    for site in &analysis.consensus.sites {
        let flag = match site.state {
            SiteState::Contested if site.call.is_none() => CONTESTED,
            SiteState::Unresolved if site.last_0based > site.first_0based => RUN_UNREAD,
            SiteState::Called
                if site.last_0based > site.first_0based
                    && site.call.as_deref().is_some_and(|call| call.contains('N')) =>
            {
                RUN_UNREAD
            }
            _ => continue,
        };
        for position in site.first_0based + 1..=site.last_0based + 1 {
            if let Some(flags) = loci.flags(position) {
                *flags |= flag;
            }
        }
    }
    let reads = analysis
        .evidence
        .reads
        .iter()
        .map(|read| {
            (
                read.alignment.orientation,
                analysis.reads.get(&read.input_sha256),
            )
        })
        .collect::<Vec<_>>();
    let mut faints = BTreeSet::new();
    for locus in analysis.observations(first, last)? {
        let position = locus.position_1based;
        let Some(flags) = loci.flags(position) else {
            continue;
        };
        for observation in &locus.observations {
            let Some(&(orientation, sanger)) = reads.get(observation.read_index) else {
                continue;
            };
            *flags |= COVERED;
            match observation.state {
                LocusState::Reference | LocusState::Alternate | LocusState::Deletion => {
                    *flags |= match orientation {
                        Orientation::Forward => FORWARD,
                        Orientation::Reverse => REVERSE,
                    };
                    if observation.state != LocusState::Reference {
                        *flags |= OTHER;
                    }
                    if let (Some(sanger), Some(call)) = (sanger, observation.call_index_0based) {
                        *flags |= trust(sanger, call, margins);
                    }
                }
                LocusState::Unresolved => *flags |= UNRESOLVED,
                LocusState::Masked => {
                    let segment =
                        sanger
                            .zip(observation.call_index_0based)
                            .and_then(|(sanger, call)| {
                                sanger.callability.segments.iter().find(|segment| {
                                    (segment.call_start_0based..segment.call_end_0based_exclusive)
                                        .contains(&call)
                                })
                            });
                    *flags |= match segment {
                        Some(segment)
                            if segment.after_repeat
                                || matches!(
                                    segment.state,
                                    PhaseState::Dephased | PhaseState::Mixed
                                ) =>
                        {
                            PHASE
                        }
                        _ => SIGNAL,
                    };
                }
            }
            if observation.state == LocusState::Masked {
                continue;
            }
            let (Some(call), Some(profile)) = (observation.base, observation.profile) else {
                continue;
            };
            let Ok(position) = i64::try_from(position) else {
                continue;
            };
            for base in Nucleotide::ALL {
                let weight = profile.weights[base.channel_index()];
                if base.as_char() != call && faint(weight) {
                    faints.insert((position, base.as_char()));
                }
            }
        }
    }
    loci.faint = faints;
    Ok(loci)
}

/// Where a read's clean `call` stands against its uninformative calls: beside
/// a call masked for phase loss, or its own end where the read stops at a
/// repeat run ([`BESIDE_PHASE_LOSS`]); a call with more than one qualifying
/// channel ([`MIXED_CALL`]) or a weak-signal veto ([`WEAK_CALL`]); within the
/// margin of its first or last
/// callable call inside its trim interval ([`EDGE`]); or neither
/// ([`INSIDE`]).
fn trust(sanger: &SangerAttachment, call: usize, margins: Margins) -> u16 {
    let callability = &sanger.callability;
    // The read's own ends: its callable span within its trim interval.
    let (start, end) = (
        callability
            .callable_start_0based
            .max(sanger.quality.trim_start_0based),
        callability
            .callable_end_0based_exclusive
            .min(sanger.quality.trim_end_0based_exclusive),
    );
    let margin = margins.read_end;
    let vetoed = sanger.calls.calls.get(call).map(|base_call| {
        let highest = base_call
            .peaks
            .iter()
            .map(|peak| usize::try_from(peak.height).unwrap_or(0))
            .max()
            .unwrap_or(0);
        let quality = sanger
            .quality
            .per_call
            .get(call)
            .map_or(0, |quality| usize::from(quality.relative_quality_score));
        (
            highest < margins.minimum_peak || quality <= margins.relative_quality,
            base_call.qualifying_channels.len() > 1,
        )
    });
    let near = |segment: &PhaseSegment| {
        segment.call_start_0based <= call + margin
            && call < segment.call_end_0based_exclusive + margin
    };
    let phase_loss = |segment: &PhaseSegment| {
        segment.after_repeat || matches!(segment.state, PhaseState::Dephased | PhaseState::Mixed)
    };
    let interior = callability.segments.iter().filter(|segment| {
        segment.state != PhaseState::InPhase
            && segment.call_start_0based >= start
            && segment.call_end_0based_exclusive <= end
    });
    if interior
        .clone()
        .any(|segment| phase_loss(segment) && near(segment))
    {
        return BESIDE_PHASE_LOSS;
    }
    match vetoed {
        Some((_, true)) => return MIXED_CALL,
        Some((true, false)) => return WEAK_CALL,
        _ => {}
    }
    if call + margin >= end {
        let tail = callability.segments.iter().find(|segment| {
            (segment.call_start_0based..segment.call_end_0based_exclusive).contains(&end)
        });
        let bases = sanger.calls.primary_sequence.as_bytes();
        let window = bases
            .get(end.saturating_sub(margin)..(end + margins.repeat).min(bases.len()))
            .unwrap_or_default();
        if tail.is_some_and(phase_loss) || longest_run(window) >= margins.repeat {
            return BESIDE_PHASE_LOSS;
        }
        return EDGE;
    }
    if call < start + margin { EDGE } else { INSIDE }
}

/// The longest run of one base in `bases`.
fn longest_run(bases: &[u8]) -> usize {
    bases
        .chunk_by(|left, right| left == right)
        .map(<[u8]>::len)
        .max()
        .unwrap_or(0)
}

/// Whether a profile weight reaches [`FAINT`].
fn faint(weight: f64) -> bool {
    weight * 1000.0 >= f64::from(FAINT)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn finds_the_longest_run_of_one_base() {
        assert_eq!(longest_run(b"ACCCCCCCCCTA"), 9);
        assert_eq!(longest_run(b"ACGT"), 1);
        assert_eq!(longest_run(b""), 0);
    }
}
