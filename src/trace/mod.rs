//! Bounds-checked decoder for one canonical analyzed ABIF/AB1 file.

mod abif;
mod decode;
mod reader;

pub(crate) use decode::load;
#[cfg(feature = "fuzzing")]
pub(crate) use abif::parse as parse_abif;
