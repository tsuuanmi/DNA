//! ABIF format parsing and Sanger chromatogram decoding.

mod container;
mod decode;
mod reader;

/// Runs the bounds-checked ABIF container parser on untrusted bytes, for the
/// fuzz harness; returns whether the bytes parsed.
#[cfg(feature = "fuzzing")]
#[must_use]
pub fn parse_container(bytes: Vec<u8>) -> bool {
    container::parse(bytes).is_ok()
}
pub use decode::load;
