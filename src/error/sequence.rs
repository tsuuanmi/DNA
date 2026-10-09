//! Consensus-sequence FASTA failures.

use std::str::Utf8Error;

/// Why a consensus-sequence FASTA was rejected.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum SequenceError {
    /// The call command received no sequence file.
    #[error("no sequence file was supplied")]
    NoSequences,
    /// The file size is zero or exceeds the accepted maximum.
    #[error("sequence file size {bytes} is outside 1..={maximum} bytes")]
    FileSize {
        /// Observed file size.
        bytes: u64,
        /// Largest accepted file size.
        maximum: usize,
    },
    /// Part of the FASTA is not UTF-8.
    #[error("{part} must be UTF-8: {error}")]
    NotUtf8 {
        /// The offending part: `"sequence file"`, `"sequence identifier"`, or
        /// `"sequence"`.
        part: &'static str,
        /// Decoding failure.
        error: Utf8Error,
    },
    /// The FASTA record syntax is invalid.
    #[error("failed to parse FASTA record: {0}")]
    Parse(std::io::Error),
    /// The file contains no record.
    #[error("sequence file contains no record")]
    Empty,
    /// The sample has more sequence records than accepted.
    #[error("more than {maximum} sequence records")]
    TooManyRecords {
        /// Largest accepted record count per sample.
        maximum: usize,
    },
    /// A record has a blank identifier.
    #[error("every FASTA record must have an identifier")]
    MissingIdentifier,
    /// Two records share one identifier.
    #[error("duplicate sequence identifier {name:?}")]
    DuplicateName {
        /// The repeated identifier.
        name: String,
    },
    /// Two records share one sequence.
    #[error("sequence {name:?} repeats an earlier record")]
    DuplicateSequence {
        /// Identifier of the repeated record.
        name: String,
    },
    /// A sequence contains a symbol that is neither a base nor an IUPAC code.
    #[error("unsupported sequence symbol {symbol:?} in {name:?}")]
    UnsupportedSymbol {
        /// Identifier of the record.
        name: String,
        /// The offending symbol as written.
        symbol: char,
    },
    /// A sequence is empty or longer than the accepted maximum.
    #[error("sequence {name:?} length {length} is outside 1..={maximum}")]
    Length {
        /// Identifier of the record.
        name: String,
        /// Observed sequence length.
        length: usize,
        /// Largest accepted sequence length.
        maximum: usize,
    },
}
