//! Non-scientific default paths and the configuration size cap.

/// Default strict TOML path when `DNA_CONFIG` is unset.
pub(super) const DEFAULT_CONFIG_PATH: &str = "config/dna.toml";
/// Largest accepted strict TOML file.
pub(super) const MAX_CONFIG_BYTES: usize = 1024 * 1024;
