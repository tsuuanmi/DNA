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

/// One normalized FASTA record and its identities.
#[derive(Debug, Clone)]
pub struct Reference {
    pub(crate) name: String,
    pub(crate) sequence: String,
    pub(crate) topology: ReferenceTopology,
    pub(crate) sequence_sha256: String,
}

impl Reference {
    /// Returns the reference length in bases.
    pub fn len(&self) -> usize {
        self.sequence.len()
    }
}
