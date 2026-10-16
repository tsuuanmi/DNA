//! Minor-channel evidence for the draft `region.minor_signal` detectors: per
//! region position, each clean read's strongest channel other than its own
//! call, with whether it is a neighbour's shadow and the read's local noise.
//!
//! It reads every covered position of the regions in one pass through
//! `SampleAnalysis::observations`, so positions where every read matches the
//! reference are included, and the read's own calls
//! (`SampleAnalysis::reads`) for its neighbours and noise.

use dna::sample_analysis::{SampleAnalysis, SangerAttachment};
use dna_core::model::alignment::Orientation;
use dna_core::model::sample_evidence::{LocusState, SampleLocusObservation};
use dna_kernel::model::nucleotide::Nucleotide;

use super::Strand;

/// Second-channel weights below this are not kept, in thousandths.
pub(crate) const FLOOR: u16 = 50;

/// Calls on each side of a call whose second channels give its local noise.
pub(crate) const NOISE_CALLS: usize = 10;

/// One clean read's second channel at one position.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) struct Minor {
    pub(crate) strand: Strand,
    /// The read's call and its strongest other channel, reference strand.
    pub(crate) pair: [Nucleotide; 2],
    /// The second channel's weight in the read's profile, in thousandths.
    pub(crate) weight: u16,
    /// The second channel's base is the base the read calls just before
    /// (`before`) or just after (`after`) this call, in sequencing order: a
    /// neighbour's shadow rather than a second allele.
    pub(crate) before: bool,
    pub(crate) after: bool,
    /// The median second-channel weight of the read's calls within
    /// [`NOISE_CALLS`] calls on either side, in thousandths.
    pub(crate) noise: u16,
}

impl Minor {
    /// Whether the second channel is a neighbour's base.
    pub(crate) const fn shadow(self) -> bool {
        self.before || self.after
    }

    /// Whether the second channel is the base of the next reference
    /// position: the call after this one on a forward read, before it on a
    /// reverse read.
    pub(crate) const fn downstream(self) -> bool {
        match self.strand {
            Strand::Forward => self.after,
            Strand::Reverse => self.before,
        }
    }

    /// The pair as an unordered key, so reads calling either base agree.
    pub(crate) fn alleles(self) -> [Nucleotide; 2] {
        let [first, second] = self.pair;
        [first.min(second), first.max(second)]
    }
}

/// Per region, per position from its start, the kept [`Minor`]s in read
/// registry order.
pub(crate) type Signals = Vec<Vec<Vec<Minor>>>;

/// Takes the minor-channel evidence of every region of `analysis`, reading
/// the span from the first region's start to the last region's end once.
///
/// # Errors
///
/// Returns DNA's error when the observations cannot be read.
pub(crate) fn extract(analysis: &SampleAnalysis) -> dna::error::Result<Signals> {
    let regions = analysis.regions();
    let mut signals: Signals = regions
        .iter()
        .map(|[start, end]| vec![Vec::new(); end.saturating_sub(*start) + 1])
        .collect();
    let first = regions.iter().map(|[start, _]| *start).min();
    let last = regions.iter().map(|[_, end]| *end).max();
    let (Some(first), Some(last)) = (first, last) else {
        return Ok(signals);
    };
    let reads = analysis
        .evidence
        .reads
        .iter()
        .map(|read| {
            let orientation = read.alignment.orientation;
            let calls = analysis
                .reads
                .get(&read.input_sha256)
                .map(|sanger| Calls::of(sanger, orientation));
            (orientation, calls)
        })
        .collect::<Vec<_>>();
    for locus in analysis.observations(first, last)? {
        let Some((region, offset)) =
            regions
                .iter()
                .enumerate()
                .find_map(|(index, [start, end])| {
                    (*start..=*end)
                        .contains(&locus.position_1based)
                        .then(|| (index, locus.position_1based - start))
                })
        else {
            continue;
        };
        let Some(kept) = signals
            .get_mut(region)
            .and_then(|positions| positions.get_mut(offset))
        else {
            continue;
        };
        for observation in &locus.observations {
            let Some((orientation, Some(calls))) = reads.get(observation.read_index) else {
                continue;
            };
            if let Some(minor) = minor(observation, *orientation, calls) {
                kept.push(minor);
            }
        }
    }
    Ok(signals)
}

/// A read's calls on the reference strand, and each call's second-channel
/// weight in thousandths.
struct Calls {
    bases: Vec<Option<Nucleotide>>,
    seconds: Vec<Option<u16>>,
}

impl Calls {
    fn of(sanger: &SangerAttachment, orientation: Orientation) -> Self {
        let bases = sanger
            .calls
            .calls
            .iter()
            .map(|call| Nucleotide::from_char(orientation.reference_base(call.primary)))
            .collect();
        // The profile and call share a strand, so the read strand serves.
        let seconds = sanger
            .calls
            .calls
            .iter()
            .enumerate()
            .map(|(index, call)| {
                let profile = sanger.locus(index)?.profile?;
                second(Nucleotide::from_char(call.primary)?, profile.weights)
                    .map(|(_, weight)| weight)
            })
            .collect();
        Self { bases, seconds }
    }
}

/// The median of `series` over the calls within [`NOISE_CALLS`] of
/// `index`, itself left out; the lower middle of an even count, 0 when none.
fn local(series: &[Option<u16>], index: usize) -> u16 {
    let first = index.saturating_sub(NOISE_CALLS);
    let last = (index + NOISE_CALLS).min(series.len().saturating_sub(1));
    let mut window = (first..=last)
        .filter(|other| *other != index)
        .filter_map(|other| series.get(other).copied().flatten())
        .collect::<Vec<_>>();
    window.sort_unstable();
    window
        .get(window.len().saturating_sub(1) / 2)
        .copied()
        .unwrap_or(0)
}

/// The minor channel of one clean, profiled observation, when it reaches
/// [`FLOOR`].
fn minor(
    observation: &SampleLocusObservation,
    orientation: Orientation,
    calls: &Calls,
) -> Option<Minor> {
    if !matches!(
        observation.state,
        LocusState::Reference | LocusState::Alternate
    ) {
        return None;
    }
    let (base, index, profile) = (
        Nucleotide::from_char(observation.base?)?,
        observation.call_index_0based?,
        observation.profile?,
    );
    let (secondary, weight) = second(base, profile.weights)?;
    if weight < FLOOR {
        return None;
    }
    let neighbour = |call: Option<usize>| {
        call.and_then(|call| calls.bases.get(call))
            .is_some_and(|base| *base == Some(secondary))
    };
    Some(Minor {
        strand: Strand::of(orientation),
        pair: [base, secondary],
        weight,
        before: neighbour(index.checked_sub(1)),
        after: neighbour(index.checked_add(1)),
        noise: local(&calls.seconds, index),
    })
}

/// The strongest finite channel other than `base`'s, the lowest on a tie,
/// and its weight in thousandths; `None` when no other channel is finite.
fn second(base: Nucleotide, weights: [f64; 4]) -> Option<(Nucleotide, u16)> {
    let (other, weight) = Nucleotide::ALL
        .into_iter()
        .filter(|other| *other != base)
        .map(|other| (other, weights[other.channel_index()]))
        .filter(|(_, weight)| weight.is_finite())
        .fold(
            None::<(Nucleotide, f64)>,
            |best, (other, weight)| match best {
                Some((_, top)) if top.total_cmp(&weight).is_ge() => best,
                _ => Some((other, weight)),
            },
        )?;
    Some((other, thousandths(weight)))
}

/// A finite weight, clamped to `[0, 1]`, in thousandths, rounded half away
/// from zero.
fn thousandths(weight: f64) -> u16 {
    let scaled = (weight.clamp(0.0, 1.0) * 1000.0).round();
    // The least integer not below `scaled`, which is a whole number.
    let (mut low, mut high) = (0_u16, 1000_u16);
    while low < high {
        let middle = low + (high - low) / 2;
        if f64::from(middle) < scaled {
            low = middle + 1;
        } else {
            high = middle;
        }
    }
    low
}

#[cfg(test)]
mod tests {
    use dna_kernel::read_evidence::EvidenceProfile;

    use super::*;

    const C: Nucleotide = Nucleotide::C;
    const T: Nucleotide = Nucleotide::T;

    /// A read's reference-strand calls and second-channel weights.
    fn calls(bases: &str, seconds: &[Option<u16>]) -> Calls {
        Calls {
            bases: bases.chars().map(Nucleotide::from_char).collect(),
            seconds: seconds.to_vec(),
        }
    }

    fn observation(
        state: LocusState,
        base: char,
        call: usize,
        weights: [f64; 4],
    ) -> SampleLocusObservation {
        SampleLocusObservation {
            read_index: 0,
            state,
            base: Some(base),
            call_index_0based: Some(call),
            profile: Some(EvidenceProfile { weights }),
            nucleotide_contribution:
                dna_core::model::sample_evidence::NucleotideContribution::Eligible,
        }
    }

    #[test]
    fn keeps_the_second_channel_with_its_shadow_and_noise() {
        // Reference strand: G, C (the call), A.
        let read = calls("GCA", &[Some(10), Some(300), Some(30)]);
        let at = |weights| {
            minor(
                &observation(LocusState::Reference, 'C', 1, weights),
                Orientation::Forward,
                &read,
            )
        };
        let Some(t) = at([0.0, 0.7, 0.0, 0.3]) else {
            panic!("a minor channel of 0.3 is kept");
        };
        assert_eq!((t.pair, t.weight, t.shadow()), ([C, T], 300, false));
        // The median of 10 and 30 is the lower middle.
        assert_eq!(t.noise, 10);
        assert_eq!(t.alleles(), [C, T]);
        let Some(g) = at([0.0, 0.7, 0.3, 0.0]) else {
            panic!("kept");
        };
        assert!(g.before && !g.after);
        let Some(a) = at([0.3, 0.7, 0.0, 0.0]) else {
            panic!("kept");
        };
        assert!(!a.before && a.after);
        // The floor is inclusive; below it, and masked calls, are dropped.
        assert!(at([0.0, 0.95, 0.0, 0.05]).is_some());
        assert!(at([0.0, 0.951, 0.0, 0.049]).is_none());
        assert!(
            minor(
                &observation(LocusState::Masked, 'C', 1, [0.0, 0.5, 0.0, 0.5]),
                Orientation::Forward,
                &read,
            )
            .is_none()
        );
    }

    #[test]
    fn reads_a_reverse_read_in_sequencing_order() {
        // Read strand `TGC` is reference `GCA` reversed: the call before the
        // middle one in sequencing order is reference `A`.
        let read = Calls {
            bases: "TGC"
                .chars()
                .map(|base| Nucleotide::from_char(Orientation::Reverse.reference_base(base)))
                .collect(),
            seconds: vec![None; 3],
        };
        assert_eq!(
            read.bases,
            [Some(Nucleotide::A), Some(C), Some(Nucleotide::G)]
        );
        let Some(minor) = minor(
            &observation(LocusState::Alternate, 'C', 1, [0.3, 0.7, 0.0, 0.0]),
            Orientation::Reverse,
            &read,
        ) else {
            panic!("kept");
        };
        assert_eq!(minor.strand, Strand::Reverse);
        assert!(minor.before && !minor.after);
        assert_eq!(minor.noise, 0);
    }

    #[test]
    fn picks_the_lowest_finite_channel_on_a_tie() {
        let a = Nucleotide::A;
        assert_eq!(second(a, [0.5, 0.25, 0.25, 0.0]), Some((C, 250)));
        // A NaN channel neither wins nor hides a real minor channel.
        assert_eq!(
            second(a, [0.5, f64::NAN, 0.2, 0.0]),
            Some((Nucleotide::G, 200))
        );
        assert_eq!(second(a, [0.5, 0.1, f64::NAN, 0.2]), Some((T, 200)));
        assert_eq!(
            second(a, [0.5, f64::INFINITY, 0.2, 0.0]),
            Some((Nucleotide::G, 200))
        );
        assert_eq!(second(a, [1.0, f64::NAN, f64::NAN, f64::NAN]), None);
        assert_eq!(thousandths(0.0995), 100);
        assert_eq!(thousandths(1.5), 1000);
    }
}
