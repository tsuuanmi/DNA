//! The `[consensus]` configuration section, owned by the consensus plugin and
//! read by the decided-edits call of its segments.

use serde::Deserialize;

use dna_kernel::error::{ConfigError, Result};

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
        Ok(ConsensusConfig {
            seeded_run: self.seeded_run,
            maximum_undecided_gap: self.maximum_undecided_gap,
            minimum_resolved_bases: self.minimum_resolved_bases,
            site_mismatch_cost: self.site_mismatch_cost,
            site_gap_open_cost: self.site_gap_open_cost,
            site_gap_extension_cost: self.site_gap_extension_cost,
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
