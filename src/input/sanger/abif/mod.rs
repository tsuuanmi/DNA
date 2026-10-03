//! ABIF format parsing and Sanger chromatogram decoding.

mod container;
mod decode;
mod reader;

#[cfg(feature = "fuzzing")]
pub(crate) use container::parse as parse_container;
pub(crate) use decode::load;
