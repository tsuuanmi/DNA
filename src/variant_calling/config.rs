//! The `[variant_calling]` configuration section.

use serde::Deserialize;

use crate::error::{ConfigError, Result};

/// Largest supported primary-sequence indel.
pub(crate) const MAX_INDEL_LENGTH: usize = 50;

/// Modality-neutral primary-difference calling settings.
#[derive(Debug, Clone)]
pub(crate) struct VariantCallingConfig {
    pub(crate) max_indel_length: usize,
    /// Calls this close to an uninformative call (beyond the trim interval or
    /// masked as unresolved) cannot support a variant.
    pub(crate) read_end_margin: usize,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct RawVariantCallingConfig {
    pub(crate) max_indel_length: usize,
    pub(crate) read_end_margin: usize,
}

impl RawVariantCallingConfig {
    /// Validates the section.
    pub(crate) fn validate(self) -> Result<VariantCallingConfig> {
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
