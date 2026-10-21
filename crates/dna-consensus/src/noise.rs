//! A run's own length noise for the two-strand rule (PROP-0014 increment 1,
//! SRS-CONS-014): the run's samples give cross-strand pairs of the rule's
//! views of an anchored window's runs, and their dispersion, shrunk towards
//! the corpus's by empirical-Bayes pseudo-pairs, rescales the rule's
//! replicate noise τ for that run.
//!
//! The statistics are additive (pairs, the sum of squared normalised
//! differences, and a histogram of the normalised differences), so a run's
//! statistics are the sum of its samples', and a sample takes its own out by
//! subtracting what it gives (leave-one-out). Without statistics, or below
//! the minimum pairs, τ is the configured one.

use crate::length::{RuleParams, Run, Source, View};

/// One sample's strand difference over one run of a window: each strand
/// decodes the run once (direct or derived), with the same top, a next
/// length on at least one strand, and both standard errors defined.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CrossStrand {
    /// The run.
    pub run: Run,
    /// `p_F − p_R`.
    pub difference: f64,
    /// `SE_F² + SE_R²`.
    pub variance: f64,
}

/// The cross-strand pair a run's `views` give, `None` without one.
#[must_use]
pub fn cross_strand(views: &[View], run: Run) -> Option<CrossStrand> {
    let decoded = |forward: bool| {
        let mut strand = views
            .iter()
            .filter(|view| view.forward == forward && view.source != Source::InPhase);
        let view = strand.next()?;
        strand.next().is_none().then_some(view)
    };
    let (forward, reverse) = (decoded(true)?, decoded(false)?);
    if forward.top != reverse.top || (forward.next.is_none() && reverse.next.is_none()) {
        return None;
    }
    let (Some(se_forward), Some(se_reverse)) = (forward.se, reverse.se) else {
        return None;
    };
    Some(CrossStrand {
        run,
        difference: forward.p - reverse.p,
        variance: se_forward.mul_add(se_forward, se_reverse * se_reverse),
    })
}

/// Whether pairs are added to statistics or taken out of them.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Change {
    Add,
    TakeOut,
}

/// One run's additive statistics.
#[derive(Debug, Clone, PartialEq)]
pub struct RunStatistics {
    /// Cross-strand pairs.
    pub pairs: usize,
    /// Sum of their squared normalised differences.
    pub sum_z2: f64,
    /// Sum of their fourth powers: with `sum_z2`, the spread of the squared
    /// differences, which sets how far a batch's dispersion varies by
    /// sampling alone.
    pub sum_z4: f64,
    /// Pairs by normalised difference: below the first edge, between each
    /// two edges, and from the last edge on.
    pub histogram: Vec<usize>,
}

/// A run's length noise statistics for run 1 and run 2, with the τ and the
/// histogram edges the differences were normalised and binned with.
#[derive(Debug, Clone, PartialEq)]
pub struct LengthStatistics {
    /// The forward strand's τ for run 1 and run 2.
    pub forward: [f64; 2],
    /// The reverse strand's.
    pub reverse: [f64; 2],
    /// The normalised differences' bin edges, ascending.
    pub edges: Vec<f64>,
    /// Run 1's and run 2's statistics.
    pub runs: [RunStatistics; 2],
}

/// How a run's statistics rescale τ (`consensus.length_context_*`).
#[derive(Debug, Clone, PartialEq)]
pub struct Calibration {
    /// Fewest pairs of a run for its own estimate to count.
    pub minimum_pairs: usize,
    /// The prior's weight for run 1 and run 2, in pairs.
    pub prior_pairs: [f64; 2],
    /// The corpus's mean squared normalised difference for run 1 and run 2:
    /// the dispersion at which τ is the configured one.
    pub reference: [f64; 2],
    /// The smallest scale τ may take.
    pub minimum_scale: f64,
    /// The normalised differences' bin edges, ascending.
    pub edges: Vec<f64>,
}

impl LengthStatistics {
    /// Empty statistics under `params`' τ and `edges`.
    #[must_use]
    pub fn empty(params: &RuleParams, edges: &[f64]) -> Self {
        let run = RunStatistics {
            pairs: 0,
            sum_z2: 0.0,
            sum_z4: 0.0,
            histogram: vec![0; edges.len() + 1],
        };
        Self {
            forward: params.forward,
            reverse: params.reverse,
            edges: edges.to_vec(),
            runs: [run.clone(), run],
        }
    }

    /// The bin of a normalised difference.
    fn bin(&self, normalised: f64) -> usize {
        self.edges.partition_point(|&edge| edge <= normalised)
    }

    /// Adds `pairs`, or takes them out.
    fn apply(&mut self, pairs: &[CrossStrand], change: Change) {
        for pair in pairs {
            let index = pair.run.index();
            let prior = self.forward[index].mul_add(
                self.forward[index],
                self.reverse[index] * self.reverse[index],
            );
            let scale = (pair.variance + prior).sqrt();
            if scale <= 0.0 {
                continue;
            }
            let normalised = pair.difference / scale;
            let bin = self.bin(normalised);
            let run = &mut self.runs[index];
            if change == Change::Add {
                run.pairs += 1;
                run.sum_z2 += normalised * normalised;
                run.sum_z4 += normalised.powi(4);
                run.histogram[bin] += 1;
            } else {
                run.pairs = run.pairs.saturating_sub(1);
                run.sum_z2 = (run.sum_z2 - normalised * normalised).max(0.0);
                run.sum_z4 = (run.sum_z4 - normalised.powi(4)).max(0.0);
                run.histogram[bin] = run.histogram[bin].saturating_sub(1);
            }
        }
    }

    /// These statistics with `pairs` added.
    #[must_use]
    pub fn with(mut self, pairs: &[CrossStrand]) -> Self {
        self.apply(pairs, Change::Add);
        self
    }

    /// These statistics without `pairs`, a sample's own (leave-one-out).
    #[must_use]
    pub fn without(mut self, pairs: &[CrossStrand]) -> Self {
        self.apply(pairs, Change::TakeOut);
        self
    }

    /// Whether the statistics were made with `params`' τ and `edges`, so
    /// they can be compared and combined with a configuration's.
    #[must_use]
    pub fn made_with(&self, params: &RuleParams, edges: &[f64]) -> bool {
        self.forward == params.forward && self.reverse == params.reverse && self.edges == edges
    }

    /// The scale of τ for run 1 and run 2: `√(φ̃ / φ_ref)`, with `φ̃` the
    /// run's mean squared normalised difference shrunk towards the
    /// reference by the prior's pseudo-pairs, at least the minimum scale;
    /// `None` for a run with fewer pairs than the minimum.
    #[must_use]
    pub fn scale(&self, calibration: &Calibration) -> [Option<f64>; 2] {
        std::array::from_fn(|index| {
            let run = &self.runs[index];
            let reference = calibration.reference[index];
            if run.pairs < calibration.minimum_pairs || reference <= 0.0 {
                return None;
            }
            #[allow(clippy::cast_precision_loss)]
            let pairs = run.pairs as f64;
            let prior = calibration.prior_pairs[index];
            let shrunk = prior.mul_add(reference, run.sum_z2) / (prior + pairs);
            Some((shrunk / reference).sqrt().max(calibration.minimum_scale))
        })
    }
}

impl RuleParams {
    /// The rule's constants for a run whose statistics give `scale`: τ of
    /// each run with a scale multiplied by it, at most one half.
    #[must_use]
    pub fn calibrated(&self, scale: [Option<f64>; 2]) -> Self {
        let mut params = *self;
        for (index, scale) in scale.iter().enumerate() {
            if let Some(scale) = scale {
                params.forward[index] = (self.forward[index] * scale).min(0.5);
                params.reverse[index] = (self.reverse[index] * scale).min(0.5);
            }
        }
        params
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::length::in_phase;

    const PARAMS: RuleParams = RuleParams {
        confidence: 2.326,
        forward: [0.065, 0.111],
        reverse: [0.111, 0.128],
    };

    fn calibration() -> Calibration {
        Calibration {
            minimum_pairs: 2,
            prior_pairs: [2.0, 2.0],
            reference: [0.5, 0.5],
            minimum_scale: 1.0,
            edges: vec![-1.0, 0.0, 1.0],
        }
    }

    fn view(forward: bool, source: Source, p: f64, se: f64) -> View {
        View {
            source,
            forward,
            top: 8,
            next: Some(9),
            p,
            se: Some(se),
        }
    }

    // SRS-CONS-014
    #[test]
    fn pairs_need_one_decoded_view_per_strand_with_the_same_top() {
        let views = [
            view(true, Source::Direct, 0.8, 0.03),
            view(false, Source::Derived, 0.7, 0.04),
        ];
        let pair = cross_strand(&views, Run::Before);
        assert!(
            pair.is_some_and(|pair| (pair.difference - 0.1).abs() < 1e-12
                && (pair.variance - 0.0025).abs() < 1e-12)
        );
        // An in-phase view is no decode.
        assert!(cross_strand(&[views[0].clone(), in_phase(false, 8)], Run::Before).is_none());
        // Two decoded views on one strand, or different tops, give none.
        assert!(
            cross_strand(
                &[views[0].clone(), views[0].clone(), views[1].clone()],
                Run::Before
            )
            .is_none()
        );
        let other = View {
            top: 9,
            ..views[1].clone()
        };
        assert!(cross_strand(&[views[0].clone(), other], Run::Before).is_none());
    }

    fn pair(run: Run, normalised: f64) -> CrossStrand {
        let index = run.index();
        let prior = PARAMS.forward[index].powi(2) + PARAMS.reverse[index].powi(2);
        CrossStrand {
            run,
            difference: normalised * prior.sqrt(),
            variance: 0.0,
        }
    }

    // SRS-CONS-014
    #[test]
    fn statistics_add_and_leave_a_sample_out() {
        let edges = calibration().edges;
        let empty = LengthStatistics::empty(&PARAMS, &edges);
        let first = [pair(Run::Before, 2.0), pair(Run::Before, -0.5)];
        let second = [pair(Run::Before, 0.5), pair(Run::After, 3.0)];
        let both = empty.clone().with(&first).with(&second);
        assert_eq!(both.runs[0].pairs, 3);
        assert!((both.runs[0].sum_z2 - (4.0 + 0.25 + 0.25)).abs() < 1e-9);
        assert!((both.runs[0].sum_z4 - (16.0 + 0.0625 + 0.0625)).abs() < 1e-9);
        assert_eq!(both.runs[0].histogram, [0, 1, 1, 1]);
        assert_eq!(both.runs[1].histogram, [0, 0, 0, 1]);
        // Taking a sample out equals never adding it.
        let out = both.without(&second);
        let alone = empty.with(&first);
        assert_eq!(out.runs[0].pairs, alone.runs[0].pairs);
        assert!((out.runs[0].sum_z2 - alone.runs[0].sum_z2).abs() < 1e-9);
        assert!((out.runs[0].sum_z4 - alone.runs[0].sum_z4).abs() < 1e-9);
        assert_eq!(out.runs[0].histogram, alone.runs[0].histogram);
        assert_eq!(out.runs[1], alone.runs[1]);
        assert!(out.made_with(&PARAMS, &edges));
        assert!(!out.made_with(&PARAMS.calibrated([Some(2.0), None]), &edges));
    }

    // SRS-CONS-014
    #[test]
    fn a_noisy_run_raises_tau_and_a_quiet_one_keeps_it() {
        let calibration = calibration();
        let empty = LengthStatistics::empty(&PARAMS, &calibration.edges);
        // Two pairs at a mean squared normalised difference of 2.0, four
        // times the reference: shrunk by two pseudo-pairs to 1.25, a scale
        // of √2.5.
        let noisy = empty.clone().with(&[pair(Run::Before, 2.0_f64.sqrt()); 2]);
        let scale = noisy.scale(&calibration);
        assert!(scale[0].is_some_and(|scale| (scale - 2.5_f64.sqrt()).abs() < 1e-9));
        assert_eq!(scale[1], None);
        let calibrated = PARAMS.calibrated(scale);
        assert!((calibrated.forward[0] - 0.065 * 2.5_f64.sqrt()).abs() < 1e-9);
        assert_eq!(calibrated.forward[1], PARAMS.forward[1]);
        // A quiet run keeps the configured τ: the scale's floor.
        let quiet = empty.clone().with(&[pair(Run::Before, 0.0); 2]);
        assert_eq!(PARAMS.calibrated(quiet.scale(&calibration)), PARAMS);
        // Too few pairs: nothing changes.
        let few = empty.with(&[pair(Run::Before, 5.0)]);
        assert_eq!(few.scale(&calibration), [None, None]);
    }
}
