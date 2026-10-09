//! Shared range checks for configuration values.

use crate::error::{ConfigError, Result};

/// Requires a finite value in `(0, 1]`.
pub(crate) fn fraction(key: &'static str, value: f64) -> Result<()> {
    finite_range(key, value, f64::MIN_POSITIVE, 1.0)
}

/// Requires a finite positive value.
pub(crate) fn positive_finite(key: &'static str, value: f64) -> Result<()> {
    if !value.is_finite() || value <= 0.0 {
        return Err(ConfigError::NotFinitePositive { key }.into());
    }
    Ok(())
}

/// Requires a finite value in `[minimum, maximum]`.
pub(crate) fn finite_range(
    key: &'static str,
    value: f64,
    minimum: f64,
    maximum: f64,
) -> Result<()> {
    if !value.is_finite() || value < minimum || value > maximum {
        return Err(ConfigError::NotFiniteInRange {
            key,
            minimum,
            maximum,
        }
        .into());
    }
    Ok(())
}
