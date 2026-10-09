//! Validated single-record reference sequence.

use serde::{Deserialize, Serialize};

/// Reference topology used by alignment and normalization.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ReferenceTopology {
    /// Reference has distinct ends.
    Linear,
    /// Reference wraps from its last base to its first base.
    Circular,
}

/// Returns the base byte of a reference sequence at a 0-based index, if any.
#[must_use]
pub fn base_at(sequence: &str, index: usize) -> Option<u8> {
    sequence.as_bytes().get(index).copied()
}

/// One normalized FASTA record and its identities.
#[derive(Debug, Clone)]
pub struct Reference {
    /// FASTA record identifier.
    pub name: String,
    /// Normalized upper-case sequence.
    pub sequence: String,
    /// Linear or circular topology.
    pub topology: ReferenceTopology,
    /// SHA-256 of the normalized sequence.
    pub sequence_sha256: String,
}

impl Reference {
    /// Returns the reference length in bases.
    #[must_use]
    pub fn len(&self) -> usize {
        self.sequence.len()
    }

    /// Whether the reference has no bases; a loaded reference never does.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.sequence.is_empty()
    }

    /// Returns the base byte at a 0-based index, if it lies inside the sequence.
    #[must_use]
    pub fn base(&self, index: usize) -> Option<u8> {
        base_at(&self.sequence, index)
    }
}
