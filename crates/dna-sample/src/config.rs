//! The `[sample_reconciliation]` configuration section.

use serde::Deserialize;

use dna_kernel::bounds;
use dna_kernel::error::{ConfigError, Result};

/// Cross-read overlap admission settings used before sample consensus.
#[derive(Debug, Clone)]
pub struct SampleReconciliationConfig {
    /// Fewest comparable canonical bases for an admitted overlap.
    pub minimum_comparable_bases: usize,
    /// Smallest canonical-base agreement for an admitted overlap.
    pub minimum_overlap_agreement: f64,
}

/// The `[sample_reconciliation]` section as written in the configuration.
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RawSampleReconciliationConfig {
    /// Fewest comparable bases an overlap needs.
    pub minimum_comparable_bases: usize,
    /// Lowest agreement an overlap needs.
    pub minimum_overlap_agreement: f64,
}

impl RawSampleReconciliationConfig {
    /// Validates the section.
    ///
    /// # Errors
    ///
    /// Returns `ConfigError` naming the first value that violates its rule.
    pub fn validate(self) -> Result<SampleReconciliationConfig> {
        if self.minimum_comparable_bases == 0 {
            return Err(ConfigError::Constraint(
                "sample_reconciliation.minimum_comparable_bases must be positive",
            )
            .into());
        }
        bounds::fraction(
            "sample_reconciliation.minimum_overlap_agreement",
            self.minimum_overlap_agreement,
        )?;
        Ok(SampleReconciliationConfig {
            minimum_comparable_bases: self.minimum_comparable_bases,
            minimum_overlap_agreement: self.minimum_overlap_agreement,
        })
    }
}
