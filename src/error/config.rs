//! Scientific configuration failures.

use std::str::Utf8Error;

/// Why a configuration file was rejected after it was read.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum ConfigError {
    /// The configuration file exceeds the accepted size.
    #[error("configuration file exceeds {maximum} bytes")]
    TooLarge {
        /// Largest accepted file size.
        maximum: usize,
    },
    /// The configuration file is not UTF-8.
    #[error("configuration must be UTF-8: {0}")]
    NotUtf8(Utf8Error),
    /// The declared schema version is not the supported one.
    #[error("unsupported schema_version {found}; expected {expected}")]
    UnsupportedSchemaVersion {
        /// Declared schema version.
        found: u32,
        /// Supported schema version.
        expected: u32,
    },
    /// A fixed scientific constraint is violated; the text names the keys.
    #[error("{0}")]
    Constraint(&'static str),
    /// A floating-point value is not finite and positive.
    #[error("{key} must be finite and positive")]
    NotFinitePositive {
        /// Dotted configuration key.
        key: &'static str,
    },
    /// A floating-point value is not finite or lies outside its closed range.
    #[error("{key} must be finite and in [{minimum}, {maximum}]")]
    NotFiniteInRange {
        /// Dotted configuration key.
        key: &'static str,
        /// Smallest accepted value.
        minimum: f64,
        /// Largest accepted value.
        maximum: f64,
    },
    /// `variant_calling.max_indel_length` is outside its supported range.
    #[error("variant_calling.max_indel_length must be in 1..={maximum}")]
    MaxIndelLength {
        /// Largest supported indel length.
        maximum: usize,
    },
    /// `variant_calling.minimum_peak_height` is outside the ABIF peak range.
    #[error("variant_calling.minimum_peak_height must be in 1..={maximum}")]
    MinimumPeakHeight {
        /// Largest representable peak height.
        maximum: i32,
    },
}
