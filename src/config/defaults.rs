//! Non-scientific default paths and hard resource caps.

/// Default strict TOML path when `DNA_CONFIG` is unset.
pub(super) const DEFAULT_CONFIG_PATH: &str = "config/dna.toml";
/// Largest accepted strict TOML file.
pub(super) const MAX_CONFIG_BYTES: usize = 1024 * 1024;
/// Largest accepted ABIF container for Sanger input.
pub(crate) const MAX_ABIF_BYTES: usize = 64 * 1024 * 1024;
/// Largest accepted FASTA source before sequence normalization.
pub(crate) const MAX_REFERENCE_BYTES: usize = 4 * 1024 * 1024;
/// Largest accepted direct-alignment reference.
pub(crate) const MAX_REFERENCE_LENGTH: usize = 50_000;
/// Largest supported primary-sequence indel.
pub(super) const MAX_INDEL_LENGTH: usize = 50;
/// Largest peak height representable by an ABIF signed short.
pub(super) const MAX_PEAK_HEIGHT: i32 = i16::MAX as i32;
/// Maximum number of traceback cells allocated by Gotoh.
pub(crate) const MAX_ALIGNMENT_CELLS: usize = 100_000_000;
