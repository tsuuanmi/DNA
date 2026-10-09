//! The `[alignment]` configuration section and the alignment size cap.

use serde::Deserialize;

use dna_kernel::bounds;
use dna_kernel::error::{ConfigError, Result};

/// Maximum number of traceback cells allocated by Gotoh.
pub const MAX_ALIGNMENT_CELLS: usize = 100_000_000;

/// Pairwise alignment settings.
#[derive(Debug, Clone)]
pub struct AlignmentConfig {
    pub(crate) match_score: i32,
    pub(crate) mismatch_score: i32,
    pub(crate) ambiguous_score: i32,
    pub(crate) gap_open_score: i32,
    pub(crate) gap_extension_score: i32,
    pub(crate) minimum_callable_bases: usize,
    pub(crate) minimum_identity: f64,
}

/// The `[alignment]` section as written in the configuration.
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RawAlignmentConfig {
    pub(crate) match_score: i32,
    pub(crate) mismatch_score: i32,
    pub(crate) ambiguous_score: i32,
    pub(crate) gap_open_score: i32,
    pub(crate) gap_extension_score: i32,
    pub(crate) minimum_callable_bases: usize,
    pub(crate) minimum_identity: f64,
}

impl RawAlignmentConfig {
    /// Validates the section.
    pub(crate) fn validate(self) -> Result<AlignmentConfig> {
        if self.match_score <= 0 {
            return Err(ConfigError::Constraint("alignment.match_score must be positive").into());
        }
        if self.mismatch_score >= 0 || self.gap_open_score >= 0 || self.gap_extension_score >= 0 {
            return Err(ConfigError::Constraint(
                "alignment mismatch and gap scores must be negative",
            )
            .into());
        }
        if self.minimum_callable_bases == 0 {
            return Err(ConfigError::Constraint(
                "alignment.minimum_callable_bases must be positive",
            )
            .into());
        }
        bounds::fraction("alignment.minimum_identity", self.minimum_identity)?;
        Ok(AlignmentConfig {
            match_score: self.match_score,
            mismatch_score: self.mismatch_score,
            ambiguous_score: self.ambiguous_score,
            gap_open_score: self.gap_open_score,
            gap_extension_score: self.gap_extension_score,
            minimum_callable_bases: self.minimum_callable_bases,
            minimum_identity: self.minimum_identity,
        })
    }
}
