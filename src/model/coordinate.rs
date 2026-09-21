//! Explicit coordinate conversion at external reporting boundaries.

use crate::error::{Error, Result};

/// Converts a zero-based reference index to a checked one-based position.
pub(crate) fn reference_one_based(position_0based: usize) -> Result<usize> {
    position_0based
        .checked_add(1)
        .ok_or_else(|| Error::Variant("reference position conversion overflow".into()))
}
