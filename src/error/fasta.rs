//! Reference FASTA failures.

use std::str::Utf8Error;

/// Why a reference FASTA was rejected.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum FastaError {
    /// The file size is zero or exceeds the accepted maximum.
    #[error("reference file size {bytes} is outside 1..={maximum} bytes")]
    FileSize {
        /// Observed file size.
        bytes: u64,
        /// Largest accepted file size.
        maximum: usize,
    },
    /// Part of the FASTA is not UTF-8.
    #[error("{part} must be UTF-8: {error}")]
    NotUtf8 {
        /// The offending part: `"reference"`, `"reference identifier"`, or
        /// `"reference sequence"`.
        part: &'static str,
        /// Decoding failure.
        error: Utf8Error,
    },
    /// The FASTA record syntax is invalid.
    #[error("failed to parse FASTA record: {0}")]
    Parse(std::io::Error),
    /// The file contains no record.
    #[error("reference is empty")]
    Empty,
    /// The file contains more than one record.
    #[error("reference must contain exactly one record")]
    MultipleRecords,
    /// The record has a blank identifier.
    #[error("first line must contain a FASTA identifier")]
    MissingIdentifier,
    /// The sequence contains a symbol other than A, C, G, T, or N.
    #[error("unsupported reference base {base:?}")]
    UnsupportedBase {
        /// The offending symbol as written.
        base: char,
    },
    /// The sequence is empty or longer than the accepted maximum.
    #[error("reference length {length} is outside 1..={maximum}")]
    Length {
        /// Observed sequence length.
        length: usize,
        /// Largest accepted sequence length.
        maximum: usize,
    },
}
