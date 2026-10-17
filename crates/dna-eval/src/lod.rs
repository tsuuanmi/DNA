//! The limit of detection of point heteroplasmy for the `haplotype_lod`
//! measure (PROP-0005 approval protocol).
//!
//! A truth IUPAC code of two bases is a heteroplasmy the traces show when the
//! covering admitted reads observe both of its bases there, each as some
//! read's call or as some read's channel other than its own call at [`LOD`]
//! or more of its normalised profile, or when one of them shows as such a
//! channel. Below that limit the measure compares the position
//! IUPAC-compatibly: DNA's base there is concordant when it is one of the
//! code's bases. At or above it, the code must match exactly, as under
//! `haplotype`.

use std::collections::BTreeSet;

use dna::sample_analysis::SampleAnalysis;
use dna_kernel::model::nucleotide::Nucleotide;

use crate::compare::{self, Allele};

/// The limit of detection, in thousandths of a read's normalised A/C/G/T
/// profile: 10 % (user, 2026-10-16).
pub(crate) const LOD: u16 = 100;

/// Per sample, what the covering admitted reads show at 1-based positions.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub(crate) struct Signals {
    /// The positions and bases some read calls.
    pub(crate) calls: BTreeSet<(i64, char)>,
    /// The positions and bases at which some read has a channel, other than
    /// its own call, at [`LOD`] or more.
    pub(crate) seconds: BTreeSet<(i64, char)>,
}

impl Signals {
    /// Whether the reads show the heteroplasmy of the two-base `code` at
    /// `position`.
    pub(crate) fn shown(&self, position: i64, code: &[char]) -> bool {
        let second = |base: &char| self.seconds.contains(&(position, *base));
        code.iter().any(second)
            || code
                .iter()
                .all(|base| second(base) || self.calls.contains(&(position, *base)))
    }
}

/// The signals of `analysis` over the span of its regions.
///
/// # Errors
///
/// Returns DNA's error when the observations cannot be read.
pub(crate) fn signals(analysis: &SampleAnalysis) -> dna::error::Result<Signals> {
    let regions = analysis.regions();
    let (Some(first), Some(last)) = (
        regions.iter().map(|[start, _]| *start).min(),
        regions.iter().map(|[_, end]| *end).max(),
    ) else {
        return Ok(Signals::default());
    };
    let mut signals = Signals::default();
    for locus in analysis.observations(first, last)? {
        let Ok(position) = i64::try_from(locus.position_1based) else {
            continue;
        };
        for observation in &locus.observations {
            let Some(call) = observation.base else {
                continue;
            };
            signals.calls.insert((position, call));
            let Some(profile) = observation.profile else {
                continue;
            };
            for base in Nucleotide::ALL {
                let weight = profile.weights[base.channel_index()];
                if base.as_char() != call && reaches(weight) {
                    signals.seconds.insert((position, base.as_char()));
                }
            }
        }
    }
    Ok(signals)
}

/// Whether a profile weight, rounded to thousandths, reaches [`LOD`].
fn reaches(weight: f64) -> bool {
    weight.is_finite() && (weight * 1000.0).round() >= f64::from(LOD)
}

/// The positions of the truth's two-base IUPAC substitutions below the
/// limit: the reads do not show their heteroplasmy.
pub(crate) fn below(truth: &[Allele], signals: &Signals) -> BTreeSet<String> {
    truth
        .iter()
        .filter(|allele| allele.reference.len() == 1)
        .filter_map(|allele| {
            let code = compare::bases(&allele.alternate);
            let position = allele.position.parse::<i64>().ok()?;
            (code.len() == 2 && !signals.shown(position, code)).then(|| allele.position.clone())
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::compare::{Concordance, Measure, Side};
    use crate::haplotype::Haplotypes;

    /// `haplotype_lod` of DNA's `dna` against the truth's `truth`, both read
    /// over HVS-I, with `signals`.
    fn lod(dna: &[Allele], truth: &[Allele], signals: &Signals) -> Concordance {
        let intervals = [(16024, 16365)];
        let dna = Side {
            alleles: dna,
            intervals: &intervals,
        };
        let truth_side = Side {
            alleles: truth,
            intervals: &intervals,
        };
        let below = below(truth, signals);
        let haplotypes: Haplotypes =
            compare::haplotypes(dna, truth_side, compare::tests::rcrs(), Some(&below));
        compare::compare(dna, truth_side, Measure::HaplotypeLod, &haplotypes, &below).concordance
    }

    fn allele(position: &str, reference: &str, alternate: &str) -> Allele {
        Allele::new(position, reference, alternate)
    }

    fn seconds(bases: &[(i64, char)]) -> Signals {
        Signals {
            calls: BTreeSet::new(),
            seconds: bases.iter().copied().collect(),
        }
    }

    #[test]
    fn reads_that_call_both_bases_show_the_heteroplasmy() {
        // 16129R (G/A): one read calls G and another A, neither with a second
        // channel at the limit.
        let truth = [allele("16129", "G", "R")];
        let both = Signals {
            calls: BTreeSet::from([(16129, 'G'), (16129, 'A')]),
            seconds: BTreeSet::new(),
        };
        assert_eq!(below(&truth, &both), BTreeSet::new());
        // Every read calls G: the heteroplasmy is not shown.
        let one = Signals {
            calls: BTreeSet::from([(16129, 'G')]),
            seconds: BTreeSet::new(),
        };
        assert_eq!(below(&truth, &one), BTreeSet::from(["16129".to_owned()]));
    }

    #[test]
    fn a_heteroplasmy_the_reads_show_must_be_called() {
        // 16129Y (G/A): a read shows A at the limit, and DNA calls G.
        let truth = [allele("16129", "G", "R"), allele("16223", "C", "T")];
        let dna = [allele("16223", "C", "T")];
        let shown = seconds(&[(16129, 'A')]);
        assert_eq!(below(&truth, &shown), BTreeSet::new());
        assert_eq!(lod(&dna, &truth, &shown), Concordance::Discordant);
    }

    #[test]
    fn below_the_limit_one_of_the_codes_bases_is_concordant() {
        let truth = [allele("16129", "G", "R"), allele("16223", "C", "T")];
        // DNA reads the reference G there.
        let dna = [allele("16223", "C", "T")];
        assert_eq!(
            below(&truth, &Signals::default()),
            BTreeSet::from(["16129".to_owned()])
        );
        assert_eq!(
            lod(&dna, &truth, &Signals::default()),
            Concordance::Concordant
        );
        // DNA calls A, the code's other base.
        let alternate = [allele("16129", "G", "A"), allele("16223", "C", "T")];
        assert_eq!(
            lod(&alternate, &truth, &Signals::default()),
            Concordance::Concordant
        );
        // DNA calls C, which the code does not admit.
        let third = [allele("16129", "G", "C"), allele("16223", "C", "T")];
        assert_eq!(
            lod(&third, &truth, &Signals::default()),
            Concordance::Discordant
        );
    }

    #[test]
    fn a_minor_reference_base_shows_the_heteroplasmy_where_dna_calls_the_alternate() {
        // 16093Y, reference T: DNA's reads call C and show T at the limit.
        let truth = [allele("16093", "T", "Y")];
        let dna = [allele("16093", "T", "C")];
        let shown = seconds(&[(16093, 'T')]);
        assert_eq!(lod(&dna, &truth, &shown), Concordance::Discordant);
        assert_eq!(
            lod(&dna, &truth, &Signals::default()),
            Concordance::Concordant
        );
    }

    #[test]
    fn a_weight_reaches_the_limit_in_rounded_thousandths() {
        assert!(reaches(0.0995) && reaches(0.1));
        assert!(!reaches(0.0994) && !reaches(f64::NAN));
    }

    #[test]
    fn only_two_base_substitution_codes_can_fall_below_the_limit() {
        let truth = [
            allele("16189", "T", "C"),
            allele("16189", "T", "B"),
            allele("309.1", "-", "C"),
        ];
        assert_eq!(below(&truth, &Signals::default()), BTreeSet::new());
    }
}
