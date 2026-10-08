//! Validated single-record reference sequence.

use serde::{Deserialize, Serialize};

/// Reference topology used by alignment and normalization.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub(crate) enum ReferenceTopology {
    /// Reference has distinct ends.
    Linear,
    /// Reference wraps from its last base to its first base.
    Circular,
}

/// Returns the base byte of a reference sequence at a 0-based index, if any.
pub(crate) fn base_at(sequence: &str, index: usize) -> Option<u8> {
    sequence.as_bytes().get(index).copied()
}

/// One normalized FASTA record and its identities.
#[derive(Debug, Clone)]
pub(crate) struct Reference {
    pub(crate) name: String,
    pub(crate) sequence: String,
    pub(crate) topology: ReferenceTopology,
    pub(crate) sequence_sha256: String,
}

impl Reference {
    /// Returns the reference length in bases.
    pub(crate) fn len(&self) -> usize {
        self.sequence.len()
    }

    /// Returns the base byte at a 0-based index, if it lies inside the sequence.
    pub(crate) fn base(&self, index: usize) -> Option<u8> {
        base_at(&self.sequence, index)
    }
}
