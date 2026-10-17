//! Discrepancy classes: which part of DNA, if any, owns a sample's
//! difference from the truth.
//!
//! A sample that failed, or with a side that has no allele, is not compared.
//! Any other takes the first class whose measure is not concordant, from the
//! measures that compare the most content to those that compare the least:
//! a sequence that differs even IUPAC-compatibly, then one that differs only
//! at a truth heteroplasmy, then DNA calls only where the reviewers did not
//! analyse, then the same sequence written differently, then a truth allele
//! where DNA has no callable base. The class follows from the statuses and,
//! for `unanalysed`, from where the primary differences lie; a baseline run's
//! classes come from its scores.

use std::collections::{BTreeMap, BTreeSet};

use crate::baseline::Class;
use crate::compare::{self, Measure, Reason};

/// What a sample's difference from the truth is.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) enum Discrepancy {
    /// Every measure concordant.
    None,
    /// Strict discordant only: a truth allele inside the truth's analysed
    /// intervals where DNA has no callable base.
    Coverage,
    /// Primary discordant, haplotype concordant: the same sequence written
    /// differently.
    Representation,
    /// Primary discordant, haplotype concordant, every primary difference a
    /// DNA call outside the truth's analysed intervals: DNA reports where the
    /// reviewers did not look.
    Unanalysed,
    /// Haplotype discordant, haplotype IUPAC-compatible concordant: the
    /// sequences differ only at a truth heteroplasmy.
    Heteroplasmy,
    /// Haplotype IUPAC-compatible discordant: the sequences differ.
    Sequence,
    /// Not applicable under a measure: DNA or the truth has no allele.
    NotComparable,
    /// DNA failed on the sample, or the truth has no record of it.
    Failed,
}

impl Discrepancy {
    pub(crate) const ALL: [Self; 8] = [
        Self::None,
        Self::Coverage,
        Self::Representation,
        Self::Unanalysed,
        Self::Heteroplasmy,
        Self::Sequence,
        Self::NotComparable,
        Self::Failed,
    ];

    /// The measures a class is decided on, in the order they are read.
    pub(crate) const MEASURES: [Measure; 4] = [
        Measure::HaplotypeIupacCompatible,
        Measure::Haplotype,
        Measure::Primary,
        Measure::Strict,
    ];

    pub(crate) const fn label(self) -> &'static str {
        match self {
            Self::None => "none",
            Self::Coverage => "coverage",
            Self::Representation => "representation",
            Self::Unanalysed => "unanalysed",
            Self::Heteroplasmy => "heteroplasmy",
            Self::Sequence => "sequence",
            Self::NotComparable => "not_comparable",
            Self::Failed => "failed",
        }
    }

    /// The part of DNA (or the lab) a difference of this class falls to.
    pub(crate) const fn owner(self) -> &'static str {
        match self {
            Self::None => "-",
            Self::Coverage => "core callability, or the lab's assay",
            Self::Representation => "dna-normalize, or a truth minority form",
            Self::Unanalysed => "a review where the truth has no analysis",
            Self::Heteroplasmy => "core minor signal (deferred)",
            Self::Sequence => "core, or the truth (adjudication)",
            Self::NotComparable => "the input, or coverage",
            Self::Failed => "core or the input",
        }
    }

    /// The class of a sample whose status class under each measure is
    /// `class`, and whose primary differences are all DNA calls outside the
    /// truth's analysed intervals when `unanalysed`; `None` when a deciding
    /// measure has no status (a baseline that does not report it).
    pub(crate) fn of(class: impl Fn(Measure) -> Option<Class>, unanalysed: bool) -> Option<Self> {
        let classes = Self::MEASURES.map(&class);
        for (excluded, discrepancy) in [
            (Class::Failed, Self::Failed),
            (Class::NotApplicable, Self::NotComparable),
        ] {
            if classes.iter().flatten().any(|class| *class == excluded) {
                return Some(discrepancy);
            }
        }
        let mut decided = Self::MEASURES.into_iter().zip(classes);
        let first = decided.try_fold(None, |found, (measure, class)| {
            let class = class?;
            Some(found.or((class != Class::Good).then_some(measure)))
        })?;
        Some(match first {
            Some(Measure::HaplotypeIupacCompatible) => Self::Sequence,
            Some(Measure::Haplotype) => Self::Heteroplasmy,
            Some(Measure::Primary) if unanalysed => Self::Unanalysed,
            Some(Measure::Primary) => Self::Representation,
            Some(_) => Self::Coverage,
            None => Self::None,
        })
    }

    /// The measure whose discordant positions locate the class.
    pub(crate) const fn measure(self) -> Option<Measure> {
        match self {
            Self::Coverage => Some(Measure::Strict),
            Self::Representation | Self::Unanalysed => Some(Measure::Primary),
            Self::Heteroplasmy => Some(Measure::Haplotype),
            Self::Sequence => Some(Measure::HaplotypeIupacCompatible),
            Self::None | Self::NotComparable | Self::Failed => None,
        }
    }

    /// A class from the statuses of one baseline sample, which cannot tell
    /// `unanalysed` from `representation`.
    pub(crate) fn of_statuses(statuses: &BTreeMap<Measure, String>) -> Option<Self> {
        Self::of(
            |measure| statuses.get(&measure).map(|status| Class::of(status)),
            false,
        )
    }

    /// Whether a primary comparison's differences, its `reasons` at its
    /// `positions`, are all DNA calls outside the truth's `analysed`
    /// intervals: a sample with any difference inside them is not
    /// `unanalysed`, since that part can be judged.
    pub(crate) fn only_unanalysed(
        reasons: &BTreeSet<Reason>,
        positions: &BTreeSet<String>,
        analysed: &[(i64, i64)],
    ) -> bool {
        !analysed.is_empty()
            && !positions.is_empty()
            && reasons.iter().all(|reason| *reason == Reason::DnaOnly)
            && positions
                .iter()
                .all(|position| !compare::inside(position, analysed))
    }

    /// The class a label names.
    pub(crate) fn from_label(label: &str) -> Option<Self> {
        Self::ALL.into_iter().find(|class| class.label() == label)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn classify(statuses: [&str; 4]) -> Option<Discrepancy> {
        Discrepancy::of(
            |measure| {
                Discrepancy::MEASURES
                    .iter()
                    .position(|decided| *decided == measure)
                    .map(|index| Class::of(statuses[index]))
            },
            false,
        )
    }

    #[test]
    fn takes_the_first_measure_that_is_not_concordant() {
        // Haplotype IUPAC-compatible, haplotype, primary, strict.
        let good = "concordant";
        let bad = "discordant";
        assert_eq!(classify([good; 4]), Some(Discrepancy::None));
        assert_eq!(
            classify([good, good, good, bad]),
            Some(Discrepancy::Coverage)
        );
        assert_eq!(
            classify([good, good, bad, bad]),
            Some(Discrepancy::Representation)
        );
        assert_eq!(
            classify([good, bad, bad, good]),
            Some(Discrepancy::Heteroplasmy)
        );
        assert_eq!(classify([bad, bad, bad, bad]), Some(Discrepancy::Sequence));
        assert_eq!(
            classify(["not_applicable", good, good, good]),
            Some(Discrepancy::NotComparable)
        );
        assert_eq!(classify(["failed"; 4]), Some(Discrepancy::Failed));
    }

    #[test]
    fn dna_calls_only_outside_the_truths_intervals_are_unanalysed() {
        let representation = |unanalysed| {
            Discrepancy::of(
                |measure| {
                    Some(Class::of(if measure == Measure::Primary {
                        "discordant"
                    } else {
                        "concordant"
                    }))
                },
                unanalysed,
            )
        };
        assert_eq!(representation(true), Some(Discrepancy::Unanalysed));
        assert_eq!(representation(false), Some(Discrepancy::Representation));
        // The truth analysed HVS-II only to 302.
        let analysed = [(73, 302), (438, 576)];
        let dna_only = BTreeSet::from([Reason::DnaOnly]);
        let outside = BTreeSet::from(["310".to_owned(), "315.1".to_owned()]);
        assert!(Discrepancy::only_unanalysed(&dna_only, &outside, &analysed));
        // One difference inside the intervals keeps the sample judged there.
        let both = BTreeSet::from(["263".to_owned(), "310".to_owned()]);
        assert!(!Discrepancy::only_unanalysed(&dna_only, &both, &analysed));
        let shared = BTreeSet::from([Reason::DnaOnly, Reason::SharedPosition]);
        assert!(!Discrepancy::only_unanalysed(&shared, &outside, &analysed));
        // A truth without intervals analysed everything.
        assert!(!Discrepancy::only_unanalysed(&dna_only, &outside, &[]));
    }

    #[test]
    fn needs_every_deciding_measure() {
        let primary_only = BTreeMap::from([(Measure::Primary, "discordant".to_owned())]);
        assert_eq!(Discrepancy::of_statuses(&primary_only), None);
        let failed = BTreeMap::from([(Measure::Primary, "failed".to_owned())]);
        assert_eq!(Discrepancy::of_statuses(&failed), Some(Discrepancy::Failed));
    }
}
