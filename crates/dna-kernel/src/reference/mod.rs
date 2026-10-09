//! Strict loading and identity for one required reference FASTA.

mod fasta;

/// Largest accepted FASTA source before sequence normalization.
pub(crate) const MAX_REFERENCE_BYTES: usize = 4 * 1024 * 1024;
/// Largest accepted direct-alignment reference.
pub const MAX_REFERENCE_LENGTH: usize = 50_000;

pub use fasta::load;
