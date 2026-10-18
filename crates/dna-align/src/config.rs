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
    /// Longest reference stretch, in bases, a run-structure re-expression may
    /// span (SRS-ALN-016).
    pub(crate) max_stretch: usize,
}

impl AlignmentConfig {
    /// Fewest callable columns a placement needs.
    #[must_use]
    pub const fn minimum_callable_bases(&self) -> usize {
        self.minimum_callable_bases
    }
}

/// The `[alignment]` section as written in the configuration.
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RawAlignmentConfig {
    /// Score of a matching base.
    pub match_score: i32,
    /// Score of a mismatching base.
    pub mismatch_score: i32,
    /// Score of a base against an ambiguous one.
    pub ambiguous_score: i32,
    /// Score of opening a gap.
    pub gap_open_score: i32,
    /// Score of extending a gap.
    pub gap_extension_score: i32,
    /// Fewest callable columns a placement needs.
    pub minimum_callable_bases: usize,
    /// Lowest callable identity a placement needs.
    pub minimum_identity: f64,
    /// Longest reference stretch a run-structure re-expression may span.
    pub max_stretch: usize,
}

impl RawAlignmentConfig {
    /// Validates the section.
    ///
    /// # Errors
    ///
    /// Returns `ConfigError` naming the first value that violates its rule.
    pub fn validate(self) -> Result<AlignmentConfig> {
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
        if self.max_stretch == 0 {
            return Err(ConfigError::Constraint("alignment.max_stretch must be positive").into());
        }
        Ok(AlignmentConfig {
            match_score: self.match_score,
            mismatch_score: self.mismatch_score,
            ambiguous_score: self.ambiguous_score,
            gap_open_score: self.gap_open_score,
            gap_extension_score: self.gap_extension_score,
            minimum_callable_bases: self.minimum_callable_bases,
            minimum_identity: self.minimum_identity,
            max_stretch: self.max_stretch,
        })
    }
}
