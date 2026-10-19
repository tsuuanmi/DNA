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
    /// The run signal was made with another context flank or share resolution
    /// than the configuration's.
    #[error(
        "run signal flank {flank} and resolution {resolution} differ from the configuration's {expected_flank} and {expected_resolution}"
    )]
    Parameters {
        /// The document's flank.
        flank: usize,
        /// The document's resolution.
        resolution: u16,
        /// The configuration's flank.
        expected_flank: usize,
        /// The configuration's resolution.
        expected_resolution: u16,
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
