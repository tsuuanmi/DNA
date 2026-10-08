//! Target-profile failures.

use std::str::Utf8Error;

/// Why a target profile was rejected after it was read, or does not apply.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
#[non_exhaustive]
pub enum ProfileError {
    /// The profile file exceeds the accepted size.
    #[error("profile file exceeds {maximum} bytes")]
    TooLarge {
        /// Largest accepted file size.
        maximum: usize,
    },
    /// The profile file is not UTF-8.
    #[error("profile must be UTF-8: {0}")]
    NotUtf8(Utf8Error),
    /// The declared schema version is not the supported one.
    #[error("unsupported profile schema_version {found}; expected {expected}")]
    UnsupportedSchemaVersion {
        /// Declared schema version.
        found: u32,
        /// Supported schema version.
        expected: u32,
    },
    /// A fixed profile rule is violated; the text names the keys.
    #[error("{0}")]
    Constraint(&'static str),
    /// A `variant_calling.regions` entry is empty, reversed, or out of range.
    #[error("variant_calling.regions[{index}] must satisfy 1 <= start <= end <= {maximum}")]
    RegionOutOfBounds {
        /// Region index in the profile list.
        index: usize,
        /// Largest supported reference coordinate.
        maximum: usize,
    },
    /// A nomenclature window definition is inconsistent.
    #[error("nomenclature window {window}: {reason}")]
    InvalidWindow {
        /// Window name as declared by the profile.
        window: String,
        /// The violated rule.
        reason: &'static str,
    },
    /// The supplied reference does not carry a nomenclature window's sequence
    /// at the window position.
    #[error("reference does not carry nomenclature window {window}")]
    WindowNotInReference {
        /// Window name as declared by the profile.
        window: String,
    },
    /// The supplied reference is not the one the profile is validated against.
    #[error("reference sequence does not match profile {profile}")]
    ReferenceMismatch {
        /// Profile identifier.
        profile: String,
    },
}
