//! What the admitted reads show at each position of the profile's regions,
//! as compact flags, for the cause of a sample's discrepancy.

use std::collections::BTreeSet;

use dna::sample_analysis::SampleAnalysis;
use dna_core::model::alignment::Orientation;
use dna_core::model::consensus::SiteState;
use dna_core::model::sample_evidence::LocusState;
use dna_kernel::model::nucleotide::Nucleotide;
use dna_sanger::model::callability::PhaseState;

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

/// Takes the flags of every region of `analysis`.
///
/// # Errors
///
/// Returns DNA's error when the observations cannot be read.
pub(crate) fn extract(analysis: &SampleAnalysis) -> dna::error::Result<Loci> {
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

/// Whether a profile weight reaches [`FAINT`].
fn faint(weight: f64) -> bool {
    weight * 1000.0 >= f64::from(FAINT)
}
