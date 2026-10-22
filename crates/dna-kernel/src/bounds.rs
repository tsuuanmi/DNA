//! Shared range checks for configuration values.

use std::ops::RangeInclusive;

use crate::error::{ConfigError, Result};

/// Requires a count in `range`.
///
/// # Errors
///
/// Returns `ConfigError::CountOutOfRange` naming `key` and the range when
/// `value` lies outside it.
pub fn count_range(key: &'static str, value: usize, range: RangeInclusive<usize>) -> Result<()> {
    if !range.contains(&value) {
        return Err(ConfigError::CountOutOfRange {
            key,
            minimum: *range.start(),
            maximum: *range.end(),
        }
        .into());
    }
    Ok(())
}

/// Requires a count of at least `minimum`.
///
/// # Errors
///
/// Returns `ConfigError::CountBelow` naming `key` and `minimum` when `value`
/// is below it.
pub fn count_at_least(key: &'static str, value: usize, minimum: usize) -> Result<()> {
    if value < minimum {
        return Err(ConfigError::CountBelow { key, minimum }.into());
    }
    Ok(())
}

/// Requires a positive count.
///
/// # Errors
///
/// Returns `ConfigError::NotPositive` naming `key` when `value` is zero.
pub fn positive_count(key: &'static str, value: usize) -> Result<()> {
    if value == 0 {
        return Err(ConfigError::NotPositive { key }.into());
    }
    Ok(())
}

/// Requires a finite value in `(0, 1]`.
///
/// # Errors
///
/// Returns `ConfigError::NotFiniteInRange` naming `key` when `value` is not
/// finite or outside `(0, 1]`.
pub fn fraction(key: &'static str, value: f64) -> Result<()> {
    finite_range(key, value, f64::MIN_POSITIVE, 1.0)
}

/// Requires a finite positive value.
///
/// # Errors
///
/// Returns `ConfigError::NotFinitePositive` naming `key` when `value` is not
/// finite and positive.
pub fn positive_finite(key: &'static str, value: f64) -> Result<()> {
    if !value.is_finite() || value <= 0.0 {
        return Err(ConfigError::NotFinitePositive { key }.into());
    }
    Ok(())
}

/// Requires a finite value in `[minimum, maximum]`.
///
/// # Errors
///
/// Returns `ConfigError::NotFiniteInRange` naming `key` when `value` is not
/// finite or outside `[minimum, maximum]`.
pub fn finite_range(key: &'static str, value: f64, minimum: f64, maximum: f64) -> Result<()> {
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

#[cfg(test)]
mod tests {
    use super::{count_at_least, count_range, positive_count};

    fn message(result: crate::error::Result<()>) -> String {
        result
            .err()
            .map_or_else(|| "ok".to_owned(), |error| error.to_string())
    }

    #[test]
    fn checks_counts_and_names_the_rule() {
        assert_eq!(message(count_range("a.b", 5, 5..=10)), "ok");
        assert_eq!(message(count_range("a.b", 10, 5..=10)), "ok");
        assert_eq!(
            message(count_range("a.b", 11, 5..=10)),
            "invalid configuration value: a.b must be in 5..=10"
        );
        assert_eq!(message(count_at_least("a.b", 2, 2)), "ok");
        assert_eq!(
            message(count_at_least("a.b", 1, 2)),
            "invalid configuration value: a.b must be at least 2"
        );
        assert_eq!(message(positive_count("a.b", 1)), "ok");
        assert_eq!(
            message(positive_count("a.b", 0)),
            "invalid configuration value: a.b must be positive"
        );
    }
}
