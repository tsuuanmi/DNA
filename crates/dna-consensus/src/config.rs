//! The `[consensus]` configuration section, owned by the consensus plugin and
//! read by the decided-edits call of its segments.

use serde::Deserialize;

use dna_kernel::bounds;
use dna_kernel::error::{ConfigError, Result};

use crate::length::RuleParams;
use crate::noise::Calibration;

/// The consensus's own settings (SRS-CFG-009).
#[derive(Debug, Clone)]
pub struct ConsensusConfig {
    /// Shortest reference homopolymer that a disputed, unbounded site grows
    /// into as one stretch (SRS-CONS-004).
    pub(crate) seeded_run: usize,
    /// Most undecided positions between two decided sites of one segment.
    pub(crate) maximum_undecided_gap: usize,
    /// Fewest resolved bases a consensus segment needs.
    pub(crate) minimum_resolved_bases: usize,
    /// Cost of a substitution when a site's call is aligned to its reference
    /// interval.
    pub(crate) site_mismatch_cost: u32,
    /// Cost of opening an insertion or a deletion in that alignment.
    pub(crate) site_gap_open_cost: u32,
    /// Cost of extending one.
    pub(crate) site_gap_extension_cost: u32,
    /// Calls from a read's ends within which its insertion or deletion
    /// against a reference decision is set aside, when a read of the other
    /// strand reads the site cleanly at least this far from its own ends
    /// (SRS-CONS-011).
    pub(crate) dissent_edge: usize,
    /// The two-strand run-length rule's confidence level and replicate
    /// noise (SRS-CONS-013).
    pub(crate) length: RuleParams,
    /// How a run context's length statistics rescale the rule's noise
    /// (PROP-0014, SRS-CONS-014).
    pub(crate) calibration: Calibration,
}

impl ConsensusConfig {
    /// The two-strand run-length rule's constants.
    #[must_use]
    pub const fn length_rule(&self) -> &RuleParams {
        &self.length
    }

    /// How a run context's length statistics rescale the rule's noise.
    #[must_use]
    pub const fn length_calibration(&self) -> &Calibration {
        &self.calibration
    }

    /// This configuration with the two-strand rule's constants replaced:
    /// the rule as a run context calibrates it.
    #[must_use]
    pub fn with_length_rule(&self, length: RuleParams) -> Self {
        Self {
            length,
            ..self.clone()
        }
    }
}

/// The `[consensus]` section as written in the configuration.
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RawConsensusConfig {
    seeded_run: usize,
    maximum_undecided_gap: usize,
    minimum_resolved_bases: usize,
    site_mismatch_cost: u32,
    site_gap_open_cost: u32,
    site_gap_extension_cost: u32,
    dissent_edge: usize,
    length_confidence: f64,
    length_noise_forward: [f64; 2],
    length_noise_reverse: [f64; 2],
    length_context_minimum_pairs: usize,
    length_context_prior_pairs: [f64; 2],
    length_context_reference: [f64; 2],
    length_context_minimum_scale: f64,
    length_context_edges: Vec<f64>,
}

impl RawConsensusConfig {
    /// Validates the section. The site costs must keep a length change in
    /// each of two runs two insertions or deletions rather than one gap and a
    /// substitution (`mismatch + extension > open`), and a substitution
    /// cheaper than an insertion plus a deletion (`mismatch < 2 * open`).
    ///
    /// # Errors
    ///
    /// Returns `ConfigError` naming the first value that violates its rule.
    pub fn validate(self) -> Result<ConsensusConfig> {
        if self.seeded_run < 2 {
            return Err(ConfigError::Constraint("consensus.seeded_run must be at least 2").into());
        }
        if self.dissent_edge == 0 {
            return Err(ConfigError::Constraint("consensus.dissent_edge must be positive").into());
        }
        if self.minimum_resolved_bases == 0 {
            return Err(ConfigError::Constraint(
                "consensus.minimum_resolved_bases must be positive",
            )
            .into());
        }
        if self.site_mismatch_cost == 0
            || self.site_gap_open_cost == 0
            || self.site_gap_extension_cost == 0
        {
            return Err(ConfigError::Constraint("consensus site costs must be positive").into());
        }
        if self.site_mismatch_cost + self.site_gap_extension_cost <= self.site_gap_open_cost
            || self.site_mismatch_cost >= 2 * self.site_gap_open_cost
        {
            return Err(ConfigError::Constraint(
                "consensus site costs must satisfy mismatch + extension > open and mismatch < 2 * open",
            )
            .into());
        }
        bounds::positive_finite("consensus.length_confidence", self.length_confidence)?;
        for noise in self
            .length_noise_forward
            .iter()
            .chain(&self.length_noise_reverse)
        {
            bounds::finite_range("consensus.length_noise", *noise, 0.0, 0.5)?;
        }
        if self.length_context_minimum_pairs == 0 {
            return Err(ConfigError::Constraint(
                "consensus.length_context_minimum_pairs must be positive",
            )
            .into());
        }
        for prior in self.length_context_prior_pairs {
            bounds::positive_finite("consensus.length_context_prior_pairs", prior)?;
        }
        for reference in self.length_context_reference {
            bounds::positive_finite("consensus.length_context_reference", reference)?;
        }
        bounds::positive_finite(
            "consensus.length_context_minimum_scale",
            self.length_context_minimum_scale,
        )?;
        if self.length_context_edges.is_empty()
            || !self
                .length_context_edges
                .iter()
                .all(|edge| edge.is_finite())
            || !self
                .length_context_edges
                .windows(2)
                .all(|pair| pair[0] < pair[1])
        {
            return Err(ConfigError::Constraint(
                "consensus.length_context_edges must be finite and strictly ascending",
            )
            .into());
        }
        Ok(ConsensusConfig {
            calibration: Calibration {
                minimum_pairs: self.length_context_minimum_pairs,
                prior_pairs: self.length_context_prior_pairs,
                reference: self.length_context_reference,
                minimum_scale: self.length_context_minimum_scale,
                edges: self.length_context_edges,
            },
            length: RuleParams {
                confidence: self.length_confidence,
                forward: self.length_noise_forward,
                reverse: self.length_noise_reverse,
            },
            seeded_run: self.seeded_run,
            maximum_undecided_gap: self.maximum_undecided_gap,
            minimum_resolved_bases: self.minimum_resolved_bases,
            site_mismatch_cost: self.site_mismatch_cost,
            site_gap_open_cost: self.site_gap_open_cost,
            site_gap_extension_cost: self.site_gap_extension_cost,
            dissent_edge: self.dissent_edge,
        })
    }
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;

    /// The shipped site costs and seeded run, with segments of any length:
    /// the tests' alignment settings take one callable base.
    pub(crate) fn config() -> ConsensusConfig {
        ConsensusConfig {
            seeded_run: 4,
            maximum_undecided_gap: 1,
            minimum_resolved_bases: 1,
            site_mismatch_cost: 2,
            site_gap_open_cost: 3,
            site_gap_extension_cost: 2,
            dissent_edge: 40,
            length: RuleParams {
                confidence: 2.326,
                forward: [0.065, 0.111],
                reverse: [0.111, 0.128],
            },
            calibration: Calibration {
                minimum_pairs: 10,
                prior_pairs: [6.0, 64.0],
                reference: [0.4643, 0.0563],
                minimum_scale: 1.0,
                edges: vec![-2.0, -1.0, 0.0, 1.0, 2.0],
            },
        }
    }

    fn raw() -> RawConsensusConfig {
        RawConsensusConfig {
            seeded_run: 4,
            maximum_undecided_gap: 20,
            minimum_resolved_bases: 20,
            site_mismatch_cost: 2,
            site_gap_open_cost: 3,
            site_gap_extension_cost: 2,
            dissent_edge: 40,
            length_confidence: 2.326,
            length_noise_forward: [0.065, 0.111],
            length_noise_reverse: [0.111, 0.128],
            length_context_minimum_pairs: 10,
            length_context_prior_pairs: [6.0, 64.0],
            length_context_reference: [0.4643, 0.0563],
            length_context_minimum_scale: 1.0,
            length_context_edges: vec![-2.0, -1.0, 0.0, 1.0, 2.0],
        }
    }

    // SRS-CFG-009
    #[test]
    fn validates_the_site_costs_and_counts() {
        assert!(raw().validate().is_ok());
        for invalid in [
            RawConsensusConfig {
                seeded_run: 1,
                ..raw()
            },
            RawConsensusConfig {
                minimum_resolved_bases: 0,
                ..raw()
            },
            RawConsensusConfig {
                site_gap_extension_cost: 0,
                ..raw()
            },
            RawConsensusConfig {
                dissent_edge: 0,
                ..raw()
            },
            RawConsensusConfig {
                length_context_minimum_pairs: 0,
                ..raw()
            },
            RawConsensusConfig {
                length_context_edges: vec![1.0, 0.0],
                ..raw()
            },
            // A two-base gap and a substitution would cost no more than two
            // one-base gaps.
            RawConsensusConfig {
                site_mismatch_cost: 1,
                site_gap_extension_cost: 1,
                ..raw()
            },
            // A substitution would cost as much as an insertion and a deletion.
            RawConsensusConfig {
                site_mismatch_cost: 6,
                ..raw()
            },
        ] {
            assert!(invalid.validate().is_err());
        }
    }
}
