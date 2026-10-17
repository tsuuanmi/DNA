//! Discrepancy classes: which part of DNA, if any, owns a sample's
//! difference from the truth.
//!
//! A sample takes the first class whose measure is not concordant, from the
//! measures that compare the most content to those that compare the least:
//! a sequence that differs even IUPAC-compatibly, then one that differs only
//! at a truth heteroplasmy, then the same sequence written differently, then
//! a truth allele where DNA has no callable base. The class is a function of
//! the statuses alone, so a baseline run's classes come from its scores.

use std::collections::BTreeMap;

use crate::baseline::Class;
use crate::compare::Measure;

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
    /// Haplotype discordant, haplotype IUPAC-compatible concordant: the
    /// sequences differ only at a truth heteroplasmy.
    Heteroplasmy,
    /// Haplotype IUPAC-compatible not concordant: the sequences differ.
    Sequence,
    /// DNA or the truth gave nothing to compare.
    Failed,
}

impl Discrepancy {
    pub(crate) const ALL: [Self; 6] = [
        Self::None,
        Self::Coverage,
        Self::Representation,
        Self::Heteroplasmy,
        Self::Sequence,
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
            Self::Heteroplasmy => "heteroplasmy",
            Self::Sequence => "sequence",
            Self::Failed => "failed",
        }
    }

    /// The part of DNA (or the lab) a difference of this class falls to.
    pub(crate) const fn owner(self) -> &'static str {
        match self {
            Self::None => "-",
            Self::Coverage => "core callability, or the lab's assay",
            Self::Representation => "dna-normalize, or a truth minority form",
            Self::Heteroplasmy => "core minor signal (deferred)",
            Self::Sequence => "core, or the truth (adjudication)",
            Self::Failed => "core or the input",
        }
    }

    /// The class of a sample whose status class under each measure is
    /// `class`; `None` when a deciding measure has no status (a baseline
    /// that does not report it).
    pub(crate) fn of(class: impl Fn(Measure) -> Option<Class>) -> Option<Self> {
        let classes = Self::MEASURES.map(&class);
        if classes
            .iter()
            .flatten()
            .any(|class| *class == Class::Failed)
        {
            return Some(Self::Failed);
        }
        let mut decided = Self::MEASURES.into_iter().zip(classes);
        let first = decided.try_fold(None, |found, (measure, class)| {
            let class = class?;
            Some(found.or((class != Class::Good).then_some(measure)))
        })?;
        Some(match first {
            Some(Measure::HaplotypeIupacCompatible) => Self::Sequence,
            Some(Measure::Haplotype) => Self::Heteroplasmy,
            Some(Measure::Primary) => Self::Representation,
            Some(_) => Self::Coverage,
            None => Self::None,
        })
    }

    /// The measure whose discordant positions locate the class.
    pub(crate) const fn measure(self) -> Option<Measure> {
        match self {
            Self::Coverage => Some(Measure::Strict),
            Self::Representation => Some(Measure::Primary),
            Self::Heteroplasmy => Some(Measure::Haplotype),
            Self::Sequence => Some(Measure::HaplotypeIupacCompatible),
            Self::None | Self::Failed => None,
        }
    }

    /// A class from the statuses of one baseline sample.
    pub(crate) fn of_statuses(statuses: &BTreeMap<Measure, String>) -> Option<Self> {
        Self::of(|measure| statuses.get(&measure).map(|status| Class::of(status)))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn classify(statuses: [&str; 4]) -> Option<Discrepancy> {
        Discrepancy::of(|measure| {
            Discrepancy::MEASURES
                .iter()
                .position(|decided| *decided == measure)
                .map(|index| Class::of(statuses[index]))
        })
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
            Some(Discrepancy::Sequence)
        );
        assert_eq!(classify(["failed"; 4]), Some(Discrepancy::Failed));
    }

    #[test]
    fn needs_every_deciding_measure() {
        let primary_only = BTreeMap::from([(Measure::Primary, "discordant".to_owned())]);
        assert_eq!(Discrepancy::of_statuses(&primary_only), None);
        let failed = BTreeMap::from([(Measure::Primary, "failed".to_owned())]);
        assert_eq!(Discrepancy::of_statuses(&failed), Some(Discrepancy::Failed));
    }
}
