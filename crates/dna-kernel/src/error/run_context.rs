//! Run-context failures (PROP-0014).

/// Why a run context (`dna.run_context/v1`) cannot be made or used.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum RunContextError {
    /// The file size is zero or exceeds the accepted maximum.
    #[error("run context file size {bytes} is outside 1..={maximum} bytes")]
    FileSize {
        /// Observed file size.
        bytes: u64,
        /// Largest accepted file size.
        maximum: u64,
    },
    /// The document is not a supported run context.
    #[error("unsupported schema version {found:?}; expected \"dna.run_context/v1\"")]
    UnsupportedSchema {
        /// The declared schema version.
        found: String,
    },
    /// A section was made with another method version.
    #[error("run context {section} method {found} differs from this build's {expected}")]
    Method {
        /// The section.
        section: &'static str,
        /// The document's method.
        found: u32,
        /// This build's method.
        expected: u32,
    },
    /// The length statistics were made with another τ, other bin edges, or
    /// over other windows than the configuration's and the profile's.
    #[error(
        "run context length statistics were made with other noise, bin edges, or windows than this analysis's"
    )]
    LengthParameters,
    /// The run context was made against another reference.
    #[error("run context was made against reference {found}, not {expected}")]
    Reference {
        /// The context's reference SHA-256.
        found: String,
        /// The analysis's reference SHA-256.
        expected: String,
    },
    /// The sample's reads are only partly among the run context's reads, so
    /// its own contribution cannot be taken out.
    #[error("the sample's reads are only partly among the run context's reads")]
    PartialReads,
    /// The run's read list names one sample twice or no read.
    #[error("run reads list {line}: {reason}")]
    ReadList {
        /// The list's line, 1-based.
        line: usize,
        /// What is wrong with it.
        reason: &'static str,
    },
    /// A count or statistic is inconsistent.
    #[error("{0}")]
    Inconsistent(&'static str),
}
