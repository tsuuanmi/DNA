//! The `[variant_calling]` configuration section.

use serde::Deserialize;

use dna_kernel::error::{ConfigError, Result};

/// Largest supported primary-sequence indel.
pub(crate) const MAX_INDEL_LENGTH: usize = 50;

/// Modality-neutral primary-difference calling settings.
#[derive(Debug, Clone)]
pub struct VariantCallingConfig {
    pub(crate) max_indel_length: usize,
    /// Calls this close to an uninformative call (beyond the trim interval or
    /// masked as unresolved) cannot support a variant.
    pub read_end_margin: usize,
}

impl VariantCallingConfig {
    /// Longest insertion or deletion called.
    #[must_use]
    pub const fn max_indel_length(&self) -> usize {
        self.max_indel_length
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
    /// Validates the section.
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
        Ok(VariantCallingConfig {
            max_indel_length: self.max_indel_length,
            read_end_margin: self.read_end_margin,
        })
    }
}
