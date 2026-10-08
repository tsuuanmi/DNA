//! ABIF container and Sanger trace record failures.

use std::fmt;
use std::str::Utf8Error;

/// Four-byte ABIF directory tag name, for example `PLOC` or `DATA`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Tag(pub [u8; 4]);

impl fmt::Display for Tag {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(&String::from_utf8_lossy(&self.0))
    }
}

/// Why an ABIF container or one of its required Sanger records was rejected.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum AbifError {
    /// The file size is zero or exceeds the accepted maximum.
    #[error("file size {bytes} is outside 1..={maximum} bytes")]
    FileSize {
        /// Observed file size.
        bytes: u64,
        /// Largest accepted file size.
        maximum: usize,
    },
    /// The content does not start with the `ABIF` magic bytes.
    #[error("missing ABIF signature")]
    MissingSignature,
    /// The root directory entry is not a `tdir` directory.
    #[error("root directory tag is not tdir")]
    RootDirectoryTag,
    /// The root directory declares a non-standard entry size.
    #[error("root directory entry size is {size}; expected {expected}")]
    RootDirectoryEntrySize {
        /// Declared entry size.
        size: usize,
        /// Standard ABIF directory entry size.
        expected: usize,
    },
    /// A byte range lies beyond the end of the file.
    #[error("byte range {offset}..{end} exceeds file length {length}")]
    ByteRange {
        /// First requested byte.
        offset: usize,
        /// End of the requested range (exclusive).
        end: usize,
        /// File length.
        length: usize,
    },
    /// A size or offset computation overflowed the platform's address space.
    #[error("{0}")]
    Overflow(&'static str),
    /// A required directory entry is absent.
    #[error("missing required tag {tag}.{number}")]
    MissingTag {
        /// Entry tag.
        tag: Tag,
        /// Entry number.
        number: u32,
    },
    /// A directory entry that must be unique occurs more than once.
    #[error("duplicate tag {tag}.{number}")]
    DuplicateTag {
        /// Entry tag.
        tag: Tag,
        /// Entry number.
        number: u32,
    },
    /// A directory entry declares no elements or zero-sized elements.
    #[error("tag {tag}.{number} has zero element size or count")]
    EmptyEntry {
        /// Entry tag.
        tag: Tag,
        /// Entry number.
        number: u32,
    },
    /// A directory entry's data is smaller than its declared element payload.
    #[error(
        "tag {tag}.{number} data size {data_size} is smaller than element size product {expected_size}"
    )]
    TruncatedEntry {
        /// Entry tag.
        tag: Tag,
        /// Entry number.
        number: u32,
        /// Declared data size.
        data_size: usize,
        /// Element size multiplied by element count.
        expected_size: usize,
    },
    /// A directory entry has an element type or size the record does not allow.
    #[error("tag {tag}.{number} has unsupported element type/size {element_type}/{element_size}")]
    UnsupportedLayout {
        /// Entry tag.
        tag: Tag,
        /// Entry number.
        number: u32,
        /// Declared ABIF element type code.
        element_type: u16,
        /// Declared element size in bytes.
        element_size: usize,
    },
    /// `FWO_.1` does not define a permutation of the four DNA channels.
    #[error("{0}")]
    ChannelOrder(&'static str),
    /// A text record is not ASCII.
    #[error("{record} is not ASCII: {error}")]
    NonAscii {
        /// Human-readable record name.
        record: &'static str,
        /// Decoding failure.
        error: Utf8Error,
    },
    /// The `DATA.9`-`DATA.12` analyzed channels are empty or unequal in length.
    #[error("DATA.9-12 channels must be non-empty and equally sized")]
    UnequalChannels,
    /// `PLOC.2` peak locations are not a valid increasing in-trace sequence.
    #[error("{0}")]
    PeakLocations(&'static str),
    /// The vendor primary base string contains a non-IUPAC symbol.
    #[error("vendor base string contains a non-IUPAC symbol")]
    NonIupacVendorBase,
    /// The trace file name cannot be represented as UTF-8 provenance.
    #[error("AB1 file name is not valid UTF-8")]
    NonUtf8FileName,
}
