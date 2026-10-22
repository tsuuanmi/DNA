//! The `[variant_calling]` configuration section.

use serde::Deserialize;

use dna_kernel::error::{ConfigError, Result};
use dna_kernel::reference::MAX_REFERENCE_LENGTH;

/// Largest supported primary-sequence indel.
pub(crate) const MAX_INDEL_LENGTH: usize = 50;

/// Modality-neutral primary-difference calling settings.
#[derive(Debug, Clone)]
pub struct VariantCallingConfig {
    pub(crate) max_indel_length: usize,
    pub(crate) read_end_margin: usize,
}

impl VariantCallingConfig {
    /// Longest insertion or deletion called.
    #[must_use]
    pub const fn max_indel_length(&self) -> usize {
        self.max_indel_length
    }

    /// Calls this close to an uninformative call (beyond the trim interval or
    /// masked as unresolved) cannot support a variant.
    #[must_use]
    pub const fn read_end_margin(&self) -> usize {
        self.read_end_margin
    }
}

/// The `[variant_calling]` section as written in the configuration.
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RawVariantCallingConfig {
    /// Longest insertion or deletion called.
    pub max_indel_length: usize,
    /// Calls this close to an uninformative call cannot support a variant.
    pub read_end_margin: usize,
}

impl RawVariantCallingConfig {
    /// Validates the section: the indel cap within the supported one, and a
    /// read-end margin no longer than the longest supported reference, past
    /// which no margin changes what a read may support.
    ///
    /// # Errors
    ///
    /// Returns `ConfigError` naming the first value that violates its rule.
    pub fn validate(self) -> Result<VariantCallingConfig> {
        if self.max_indel_length == 0 || self.max_indel_length > MAX_INDEL_LENGTH {
            return Err(ConfigError::MaxIndelLength {
                maximum: MAX_INDEL_LENGTH,
            }
            .into());
        }
        if self.read_end_margin > MAX_REFERENCE_LENGTH {
            return Err(ConfigError::Constraint(
                "variant_calling.read_end_margin must be at most the longest supported reference",
            )
            .into());
        }
        Ok(VariantCallingConfig {
            max_indel_length: self.max_indel_length,
            read_end_margin: self.read_end_margin,
        })
    }
}
