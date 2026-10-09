//! Explicit coordinate conversion at external reporting boundaries.

use dna_kernel::error::{Result, VariantError};

/// Converts a zero-based reference index to a checked one-based position.
pub(crate) fn reference_one_based(position_0based: usize) -> Result<usize> {
    Ok(position_0based
        .checked_add(1)
        .ok_or(VariantError::Overflow(
            "reference position conversion overflow",
        ))?)
}
