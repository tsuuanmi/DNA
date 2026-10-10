//! Canonical ABIF tag decoding into a validated chromatogram.

use std::fs;
use std::path::Path;

use dna_kernel::checksum::hex_sha256;

use crate::abif::container::{AbifEntry, AbifFile, parse};
use crate::abif::reader::Reader;
use crate::model::sanger::{Chromatogram, VendorEvidence};
use dna_kernel::error::{AbifError, Error, Result, Tag};
use dna_kernel::model::nucleotide::Nucleotide;

/// Largest accepted ABIF container for Sanger input.
const MAX_ABIF_BYTES: usize = 64 * 1024 * 1024;

const TYPE_BYTE: u16 = 1;
const TYPE_CHAR: u16 = 2;
const TYPE_SHORT: u16 = 4;

/// Reads and decodes one analyzed Sanger trace stored in an ABIF container.
///
/// # Errors
///
/// Returns `Error::Read` when the file cannot be read and `AbifError` when it
/// is not a valid analyzed ABIF trace.
pub fn load(path: &Path) -> Result<Chromatogram> {
    let metadata = fs::metadata(path).map_err(|source| Error::Read {
        kind: "AB1",
        path: path.to_path_buf(),
        source,
    })?;
    if metadata.len() == 0 || metadata.len() > MAX_ABIF_BYTES as u64 {
        return Err(AbifError::FileSize {
            bytes: metadata.len(),
            maximum: MAX_ABIF_BYTES,
        }
        .into());
    }
    let bytes = fs::read(path).map_err(|source| Error::Read {
        kind: "AB1",
        path: path.to_path_buf(),
        source,
    })?;
    if bytes.is_empty() || bytes.len() > MAX_ABIF_BYTES {
        return Err(AbifError::FileSize {
            bytes: u64::try_from(bytes.len()).unwrap_or(u64::MAX),
            maximum: MAX_ABIF_BYTES,
        }
        .into());
    }
    let source_sha256 = hex_sha256(&bytes);
    let abif = parse(bytes)?;
    decode(path, &abif, source_sha256)
}

fn decode(path: &Path, abif: &AbifFile, source_sha256: String) -> Result<Chromatogram> {
    let order_entry = abif.required(*b"FWO_", 1)?;
    require_layout(order_entry, TYPE_CHAR, 1)?;
    let order_bytes = abif.payload(order_entry)?;
    if order_bytes.len() != 4 {
        return Err(AbifError::ChannelOrder("FWO_.1 must contain exactly four bases").into());
    }
    let channel_order = std::str::from_utf8(order_bytes)
        .map_err(|error| AbifError::NonAscii {
            record: "FWO_.1",
            error,
        })?
        .to_owned();
    let mut seen = [false; 4];
    for base in channel_order.chars() {
        let index = channel_index(base).ok_or(AbifError::ChannelOrder(
            "FWO_.1 is not an A/C/G/T permutation",
        ))?;
        if seen[index] {
            return Err(AbifError::ChannelOrder("FWO_.1 repeats a channel").into());
        }
        seen[index] = true;
    }

    let mut raw_channels: [Vec<i32>; 4] = std::array::from_fn(|_| Vec::new());
    for (index, number) in (9_u32..=12).enumerate() {
        let entry = abif.required(*b"DATA", number)?;
        require_layout(entry, TYPE_SHORT, 2)?;
        raw_channels[index] = decode_i16(abif, entry)?;
    }
    let sample_count = raw_channels[0].len();
    if sample_count == 0
        || raw_channels
            .iter()
            .any(|channel| channel.len() != sample_count)
    {
        return Err(AbifError::UnequalChannels.into());
    }
    let mut channels: [Vec<i32>; 4] = std::array::from_fn(|_| Vec::new());
    for (source_index, base) in channel_order.chars().enumerate() {
        let target_index =
            channel_index(base).ok_or(AbifError::ChannelOrder("invalid channel order"))?;
        channels[target_index] = std::mem::take(&mut raw_channels[source_index]);
    }

    let ploc_entry = abif.required(*b"PLOC", 2)?;
    require_layout(ploc_entry, TYPE_SHORT, 2)?;
    let peak_locations: Vec<usize> = decode_i16(abif, ploc_entry)?
        .into_iter()
        .map(|value| {
            usize::try_from(value)
                .map_err(|_| AbifError::PeakLocations("PLOC.2 contains a negative position"))
        })
        .collect::<std::result::Result<_, AbifError>>()?;
    if peak_locations.is_empty() {
        return Err(AbifError::PeakLocations("PLOC.2 is empty").into());
    }
    // A position equal to its predecessor names the same signal event twice:
    // it is merged into one locus, and the vendor entries at the same index
    // are dropped so the vendor series stay aligned with the loci.
    let mut repeated = vec![false; peak_locations.len()];
    for (index, pair) in peak_locations.windows(2).enumerate() {
        if pair[0] > pair[1] {
            return Err(AbifError::PeakLocations("PLOC.2 positions must not decrease").into());
        }
        repeated[index + 1] = pair[0] == pair[1];
    }
    let duplicate_loci = repeated.iter().filter(|&&dropped| dropped).count();
    let locus_positions = kept(peak_locations, &repeated);
    if locus_positions
        .iter()
        .any(|position| *position >= sample_count)
    {
        return Err(
            AbifError::PeakLocations("PLOC.2 position lies outside channel samples").into(),
        );
    }

    let primary = decode_optional_string(abif, *b"PBAS", 2)?.map(|bases| {
        kept(bases.chars().collect(), &repeated)
            .into_iter()
            .collect()
    });
    let qualities = decode_optional_bytes(abif, *b"PCON", 2)?.map(|values| kept(values, &repeated));
    let source_name = path
        .file_name()
        .and_then(|name| name.to_str())
        .ok_or(AbifError::NonUtf8FileName)?
        .to_owned();

    Ok(Chromatogram {
        source_name,
        source_sha256,
        channels,
        locus_positions,
        duplicate_loci,
        vendor: VendorEvidence { primary, qualities },
    })
}

/// The elements whose index is not marked as a repeated peak location; a
/// series longer than the marks keeps its unmarked tail.
fn kept<T>(values: Vec<T>, repeated: &[bool]) -> Vec<T> {
    values
        .into_iter()
        .enumerate()
        .filter(|(index, _)| !repeated.get(*index).copied().unwrap_or(false))
        .map(|(_, value)| value)
        .collect()
}

fn require_layout(entry: &AbifEntry, element_type: u16, element_size: usize) -> Result<()> {
    if entry.element_type != element_type || entry.element_size != element_size {
        return Err(unsupported_layout(entry).into());
    }
    Ok(())
}

fn decode_i16(abif: &AbifFile, entry: &AbifEntry) -> Result<Vec<i32>> {
    let payload = abif.payload(entry)?;
    let reader = Reader::new(payload);
    Ok((0..entry.element_count)
        .map(|index| reader.i16(index * 2).map(i32::from))
        .collect::<std::result::Result<_, AbifError>>()?)
}

fn decode_optional_string(abif: &AbifFile, tag: [u8; 4], number: u32) -> Result<Option<String>> {
    let Some(entry) = abif.optional(tag, number)? else {
        return Ok(None);
    };
    require_layout(entry, TYPE_CHAR, 1)?;
    let payload = abif.payload(entry)?;
    let end = payload
        .iter()
        .position(|byte| *byte == 0)
        .unwrap_or(payload.len());
    let text = std::str::from_utf8(&payload[..end]).map_err(|error| AbifError::NonAscii {
        record: "vendor base string",
        error,
    })?;
    if !text.chars().all(|base| {
        matches!(
            base,
            'A' | 'C'
                | 'G'
                | 'T'
                | 'U'
                | 'R'
                | 'Y'
                | 'S'
                | 'W'
                | 'K'
                | 'M'
                | 'B'
                | 'D'
                | 'H'
                | 'V'
                | 'N'
        )
    }) {
        return Err(AbifError::NonIupacVendorBase.into());
    }
    Ok(Some(text.to_owned()))
}

fn decode_optional_bytes(abif: &AbifFile, tag: [u8; 4], number: u32) -> Result<Option<Vec<u8>>> {
    let Some(entry) = abif.optional(tag, number)? else {
        return Ok(None);
    };
    if entry.element_size != 1 || !matches!(entry.element_type, TYPE_BYTE | TYPE_CHAR) {
        return Err(unsupported_layout(entry).into());
    }
    Ok(Some(abif.payload(entry)?.to_vec()))
}

fn unsupported_layout(entry: &AbifEntry) -> AbifError {
    AbifError::UnsupportedLayout {
        tag: Tag(entry.tag),
        number: entry.number,
        element_type: entry.element_type,
        element_size: entry.element_size,
    }
}

fn channel_index(base: char) -> Option<usize> {
    Nucleotide::from_char(base).map(Nucleotide::channel_index)
}
