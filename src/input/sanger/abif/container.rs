//! Strict ABIF directory parsing and exact tag lookup.

use crate::error::{AbifError, Result, Tag};
use crate::input::sanger::abif::reader::Reader;

const DIRECTORY_ENTRY_SIZE: usize = 28;
const ROOT_ENTRY_OFFSET: usize = 6;

/// One validated ABIF directory entry.
#[derive(Debug, Clone)]
pub(crate) struct AbifEntry {
    pub(crate) tag: [u8; 4],
    pub(crate) number: u32,
    pub(crate) element_type: u16,
    pub(crate) element_size: usize,
    pub(crate) element_count: usize,
    pub(crate) data_size: usize,
    pub(crate) data_offset: usize,
    pub(crate) entry_offset: usize,
}

/// Parsed ABIF directory over owned input bytes.
#[derive(Debug, Clone)]
pub(crate) struct AbifFile {
    pub(crate) bytes: Vec<u8>,
    pub(crate) entries: Vec<AbifEntry>,
}

impl AbifFile {
    /// Returns the unique requested entry.
    pub(crate) fn required(&self, tag: [u8; 4], number: u32) -> Result<&AbifEntry> {
        Ok(self.optional(tag, number)?.ok_or(AbifError::MissingTag {
            tag: Tag(tag),
            number,
        })?)
    }

    /// Returns an optional unique entry and rejects duplicates.
    pub(crate) fn optional(&self, tag: [u8; 4], number: u32) -> Result<Option<&AbifEntry>> {
        let mut matches = self
            .entries
            .iter()
            .filter(|entry| entry.tag == tag && entry.number == number);
        let first = matches.next();
        if matches.next().is_some() {
            return Err(AbifError::DuplicateTag {
                tag: Tag(tag),
                number,
            }
            .into());
        }
        Ok(first)
    }

    /// Returns the logical element payload, excluding reserved allocation padding.
    pub(crate) fn payload(&self, entry: &AbifEntry) -> Result<&[u8]> {
        let offset = if entry.data_size <= 4 {
            entry.entry_offset + 20
        } else {
            entry.data_offset
        };
        let logical_size = element_payload_size(entry.element_size, entry.element_count)?;
        Ok(Reader::new(&self.bytes).slice(offset, logical_size)?)
    }
}

/// Parses the ABIF header and full root directory.
pub(crate) fn parse(bytes: Vec<u8>) -> Result<AbifFile> {
    let reader = Reader::new(&bytes);
    if reader.slice(0, 4)? != b"ABIF" {
        return Err(AbifError::MissingSignature.into());
    }
    let _version = reader.u16(4)?;
    let root = parse_entry(&reader, ROOT_ENTRY_OFFSET)?;
    if &root.tag != b"tdir" {
        return Err(AbifError::RootDirectoryTag.into());
    }
    if root.element_size != DIRECTORY_ENTRY_SIZE {
        return Err(AbifError::RootDirectoryEntrySize {
            size: root.element_size,
            expected: DIRECTORY_ENTRY_SIZE,
        }
        .into());
    }
    let directory_bytes = root
        .element_count
        .checked_mul(DIRECTORY_ENTRY_SIZE)
        .ok_or(AbifError::Overflow("root directory size overflow"))?;
    reader.slice(root.data_offset, root.data_size)?;
    reader.slice(root.data_offset, directory_bytes)?;

    let mut entries = Vec::with_capacity(root.element_count);
    for index in 0..root.element_count {
        let delta = index
            .checked_mul(DIRECTORY_ENTRY_SIZE)
            .ok_or(AbifError::Overflow("directory offset overflow"))?;
        let offset = root
            .data_offset
            .checked_add(delta)
            .ok_or(AbifError::Overflow("directory offset overflow"))?;
        let entry = parse_entry(&reader, offset)?;
        let payload_offset = if entry.data_size <= 4 {
            entry.entry_offset + 20
        } else {
            entry.data_offset
        };
        reader.slice(payload_offset, entry.data_size)?;
        entries.push(entry);
    }

    Ok(AbifFile { bytes, entries })
}

fn parse_entry(reader: &Reader<'_>, offset: usize) -> Result<AbifEntry> {
    let tag_bytes = reader.slice(offset, 4)?;
    let tag = [tag_bytes[0], tag_bytes[1], tag_bytes[2], tag_bytes[3]];
    let number = reader.u32(offset + 4)?;
    let element_type = reader.u16(offset + 8)?;
    let element_size = usize::from(reader.u16(offset + 10)?);
    let element_count = usize::try_from(reader.u32(offset + 12)?)
        .map_err(|_| AbifError::Overflow("element count does not fit memory size"))?;
    let data_size = usize::try_from(reader.u32(offset + 16)?)
        .map_err(|_| AbifError::Overflow("data size does not fit memory size"))?;
    let data_offset = usize::try_from(reader.u32(offset + 20)?)
        .map_err(|_| AbifError::Overflow("data offset does not fit memory size"))?;
    if element_size == 0 || element_count == 0 {
        return Err(AbifError::EmptyEntry {
            tag: Tag(tag),
            number,
        }
        .into());
    }
    let expected_size = element_payload_size(element_size, element_count)?;
    if data_size < expected_size {
        return Err(AbifError::TruncatedEntry {
            tag: Tag(tag),
            number,
            data_size,
            expected_size,
        }
        .into());
    }
    Ok(AbifEntry {
        tag,
        number,
        element_type,
        element_size,
        element_count,
        data_size,
        data_offset,
        entry_offset: offset,
    })
}

fn element_payload_size(
    element_size: usize,
    element_count: usize,
) -> std::result::Result<usize, AbifError> {
    element_size
        .checked_mul(element_count)
        .ok_or(AbifError::Overflow("element payload size overflow"))
}

#[cfg(test)]
mod tests {
    use crate::error::Error;

    use super::*;

    const DIRECTORY_OFFSET: usize = 64;
    const DIRECTORY_OFFSET_U32: u32 = 64;
    const DIRECTORY_ENTRY_SIZE_U16: u16 = 28;
    const DIRECTORY_ENTRY_SIZE_U32: u32 = 28;
    const PAYLOAD_OFFSET_U32: u32 = 92;

    #[test]
    fn accepts_oversized_root_directory_allocation() -> Result<()> {
        let mut bytes = vec![0_u8; DIRECTORY_OFFSET + 56];
        write_header(&mut bytes, 56);
        write_entry(
            &mut bytes,
            DIRECTORY_OFFSET,
            *b"FWO_",
            1,
            2,
            1,
            4,
            4,
            0,
            b"ACGT",
        );

        let abif = parse(bytes)?;
        let entry = abif.required(*b"FWO_", 1)?;
        assert_eq!(abif.payload(entry)?, b"ACGT");
        Ok(())
    }

    #[test]
    fn returns_logical_payload_from_oversized_external_allocation() -> Result<()> {
        let payload_offset = DIRECTORY_OFFSET + DIRECTORY_ENTRY_SIZE;
        let mut bytes = vec![0_u8; payload_offset + 8];
        write_header(&mut bytes, DIRECTORY_ENTRY_SIZE_U32);
        write_entry(
            &mut bytes,
            DIRECTORY_OFFSET,
            *b"TEST",
            1,
            2,
            1,
            4,
            8,
            PAYLOAD_OFFSET_U32,
            &[],
        );
        bytes[payload_offset..payload_offset + 8].copy_from_slice(b"ACGTpad!");

        let abif = parse(bytes)?;
        let entry = abif.required(*b"TEST", 1)?;
        assert_eq!(abif.payload(entry)?, b"ACGT");
        Ok(())
    }

    #[test]
    fn rejects_allocation_smaller_than_logical_payload() {
        let mut bytes = vec![0_u8; DIRECTORY_OFFSET + DIRECTORY_ENTRY_SIZE];
        write_header(&mut bytes, DIRECTORY_ENTRY_SIZE_U32);
        write_entry(
            &mut bytes,
            DIRECTORY_OFFSET,
            *b"TEST",
            1,
            2,
            1,
            4,
            3,
            0,
            b"ACG",
        );

        let result = parse(bytes);
        assert!(matches!(
            result,
            Err(Error::Abif(AbifError::TruncatedEntry {
                tag: Tag(tag),
                number: 1,
                data_size: 3,
                expected_size: 4,
            })) if tag == *b"TEST"
        ));
    }

    #[test]
    fn rejects_content_without_abif_signature() {
        assert!(matches!(
            parse(b"not an abif container".to_vec()),
            Err(Error::Abif(AbifError::MissingSignature))
        ));
    }

    fn write_header(bytes: &mut [u8], root_data_size: u32) {
        bytes[0..4].copy_from_slice(b"ABIF");
        bytes[4..6].copy_from_slice(&101_u16.to_be_bytes());
        write_entry(
            bytes,
            ROOT_ENTRY_OFFSET,
            *b"tdir",
            1,
            1023,
            DIRECTORY_ENTRY_SIZE_U16,
            1,
            root_data_size,
            DIRECTORY_OFFSET_U32,
            &[],
        );
    }

    #[allow(clippy::too_many_arguments)]
    fn write_entry(
        bytes: &mut [u8],
        offset: usize,
        tag: [u8; 4],
        number: u32,
        element_type: u16,
        element_size: u16,
        element_count: u32,
        data_size: u32,
        data_offset: u32,
        inline: &[u8],
    ) {
        bytes[offset..offset + 4].copy_from_slice(&tag);
        bytes[offset + 4..offset + 8].copy_from_slice(&number.to_be_bytes());
        bytes[offset + 8..offset + 10].copy_from_slice(&element_type.to_be_bytes());
        bytes[offset + 10..offset + 12].copy_from_slice(&element_size.to_be_bytes());
        bytes[offset + 12..offset + 16].copy_from_slice(&element_count.to_be_bytes());
        bytes[offset + 16..offset + 20].copy_from_slice(&data_size.to_be_bytes());
        if data_size <= 4 {
            bytes[offset + 20..offset + 24].fill(0);
            bytes[offset + 20..offset + 20 + inline.len()].copy_from_slice(inline);
        } else {
            bytes[offset + 20..offset + 24].copy_from_slice(&data_offset.to_be_bytes());
        }
    }
}
