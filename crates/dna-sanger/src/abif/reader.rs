//! Checked big-endian reads over untrusted ABIF bytes.

use dna_kernel::error::AbifError;

/// Result of a bounds-checked read.
type ReadResult<T> = Result<T, AbifError>;

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
    pub(crate) fn slice(&self, offset: usize, length: usize) -> ReadResult<&'a [u8]> {
        let end = offset
            .checked_add(length)
            .ok_or(AbifError::Overflow("byte range overflow"))?;
        self.bytes.get(offset..end).ok_or(AbifError::ByteRange {
            offset,
            end,
            length: self.bytes.len(),
        })
    }

    /// Reads a big-endian unsigned 16-bit integer.
    pub(crate) fn u16(&self, offset: usize) -> ReadResult<u16> {
        let bytes = self.slice(offset, 2)?;
        Ok(u16::from_be_bytes([bytes[0], bytes[1]]))
    }

    /// Reads a big-endian signed 16-bit integer.
    pub(crate) fn i16(&self, offset: usize) -> ReadResult<i16> {
        let bytes = self.slice(offset, 2)?;
        Ok(i16::from_be_bytes([bytes[0], bytes[1]]))
    }

    /// Reads a big-endian unsigned 32-bit integer.
    pub(crate) fn u32(&self, offset: usize) -> ReadResult<u32> {
        let bytes = self.slice(offset, 4)?;
        Ok(u32::from_be_bytes([bytes[0], bytes[1], bytes[2], bytes[3]]))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reads_big_endian_signed_and_unsigned_values() -> ReadResult<()> {
        let reader = Reader::new(&[0xFF, 0xFE, 0x01, 0x02, 0x03, 0x04]);
        assert_eq!(reader.i16(0)?, -2);
        assert_eq!(reader.u16(2)?, 0x0102);
        assert_eq!(reader.u32(2)?, 0x0102_0304);
        Ok(())
    }

    #[test]
    fn rejects_out_of_bounds_and_overflowing_ranges() {
        let reader = Reader::new(&[0, 1]);
        assert!(matches!(
            reader.slice(1, 2),
            Err(AbifError::ByteRange {
                offset: 1,
                end: 3,
                length: 2
            })
        ));
        assert!(matches!(
            reader.slice(usize::MAX, 2),
            Err(AbifError::Overflow(_))
        ));
    }
}
