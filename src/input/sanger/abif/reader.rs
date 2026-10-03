//! Checked big-endian reads over untrusted ABIF bytes.

use crate::error::{Error, Result};

/// Bounds-checked view over a binary input.
pub(crate) struct Reader<'a> {
    bytes: &'a [u8],
}

impl<'a> Reader<'a> {
    /// Wraps immutable input bytes.
    pub(crate) const fn new(bytes: &'a [u8]) -> Self {
        Self { bytes }
    }

    /// Returns a checked byte slice.
    pub(crate) fn slice(&self, offset: usize, length: usize) -> Result<&'a [u8]> {
        let end = offset
            .checked_add(length)
            .ok_or_else(|| Error::Abif("byte range overflow".into()))?;
        self.bytes.get(offset..end).ok_or_else(|| {
            Error::Abif(format!(
                "byte range {offset}..{end} exceeds file length {}",
                self.bytes.len()
            ))
        })
    }

    /// Reads a big-endian unsigned 16-bit integer.
    pub(crate) fn u16(&self, offset: usize) -> Result<u16> {
        let bytes = self.slice(offset, 2)?;
        Ok(u16::from_be_bytes([bytes[0], bytes[1]]))
    }

    /// Reads a big-endian signed 16-bit integer.
    pub(crate) fn i16(&self, offset: usize) -> Result<i16> {
        let bytes = self.slice(offset, 2)?;
        Ok(i16::from_be_bytes([bytes[0], bytes[1]]))
    }

    /// Reads a big-endian unsigned 32-bit integer.
    pub(crate) fn u32(&self, offset: usize) -> Result<u32> {
        let bytes = self.slice(offset, 4)?;
        Ok(u32::from_be_bytes([bytes[0], bytes[1], bytes[2], bytes[3]]))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reads_big_endian_signed_and_unsigned_values() -> Result<()> {
        let reader = Reader::new(&[0xFF, 0xFE, 0x01, 0x02, 0x03, 0x04]);
        assert_eq!(reader.i16(0)?, -2);
        assert_eq!(reader.u16(2)?, 0x0102);
        assert_eq!(reader.u32(2)?, 0x01020304);
        Ok(())
    }

    #[test]
    fn rejects_out_of_bounds_and_overflowing_ranges() {
        let reader = Reader::new(&[0, 1]);
        assert!(reader.slice(1, 2).is_err());
        assert!(reader.slice(usize::MAX, 2).is_err());
    }
}
