//! Failures of a run signal document (PROP-0009).

use std::path::PathBuf;

use dna_kernel::error::ForeignError;

/// Why a run signal cannot be used.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum RunSignalError {
    /// The document bytes are not JSON of the run signal's shape.
    #[error("{path}: {source}")]
    Parse {
        /// The document.
        path: PathBuf,
        /// Underlying JSON parse failure.
        #[source]
        source: ForeignError,
    },
    /// The document is not a `dna.signal/v1` run signal.
    #[error("unsupported run signal schema {found:?}")]
    UnsupportedSchema {
        /// The schema the document names.
        found: String,
    },
    /// The document was made with another call selection.
    #[error("run signal selection {found} differs from this build's {expected}")]
    Selection {
        /// The document's selection.
        found: u32,
        /// This build's selection.
        expected: u32,
    },
    /// A key or share is outside the run signal's vocabulary.
    #[error("invalid run signal key {context:?}")]
    Key {
        /// The key's context as written.
        context: String,
    },
    /// The keys' counts disagree with the document's call total.
    #[error("run signal keys hold {total} counts for {calls} calls")]
    Total {
        /// The document's selected calls.
        calls: u64,
        /// The keys' counts.
        total: u64,
    },
    /// The run signal shares no trace with the sample: it is another run's.
    #[error("run signal {run_id:?} holds none of the sample's traces")]
    OtherRun {
        /// The run signal's run.
        run_id: String,
    },
}
