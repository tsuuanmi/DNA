//! Canonical ABIF tag decoding into a validated chromatogram.

use std::fs;
use std::path::Path;

use crate::checksum::hex_sha256;
use crate::config::MAX_AB1_BYTES;
use crate::error::{Error, Result};
use crate::model::trace::{Chromatogram, VendorEvidence};
use crate::trace::abif::{AbifEntry, AbifFile, parse};
use crate::trace::reader::Reader;

const TYPE_BYTE: u16 = 1;
const TYPE_CHAR: u16 = 2;
const TYPE_SHORT: u16 = 4;

/// Reads and decodes one canonical analyzed ABIF/AB1 file.
pub(crate) fn load(path: &Path) -> Result<Chromatogram> {
    let metadata = fs::metadata(path).map_err(|source| Error::Read {
        kind: "AB1",
        path: path.to_path_buf(),
        source,
    })?;
    if metadata.len() == 0 || metadata.len() > MAX_AB1_BYTES as u64 {
        return Err(Error::Abif(format!(
            "file size {} is outside 1..={MAX_AB1_BYTES} bytes",
            metadata.len()
        )));
    }
    let bytes = fs::read(path).map_err(|source| Error::Read {
        kind: "AB1",
        path: path.to_path_buf(),
        source,
    })?;
    if bytes.is_empty() || bytes.len() > MAX_AB1_BYTES {
        return Err(Error::Abif(format!(
            "file size {} is outside 1..={MAX_AB1_BYTES} bytes",
            bytes.len()
        )));
    }
    let source_sha256 = hex_sha256(&bytes);
    let abif = parse(bytes)?;
    decode(path, abif, source_sha256)
}

fn decode(path: &Path, abif: AbifFile, source_sha256: String) -> Result<Chromatogram> {
    let order_entry = abif.required(b"FWO_", 1)?;
    require_layout(order_entry, TYPE_CHAR, 1)?;
    let order_bytes = abif.payload(order_entry)?;
    if order_bytes.len() != 4 {
        return Err(Error::Abif("FWO_.1 must contain exactly four bases".into()));
    }
    let channel_order = std::str::from_utf8(order_bytes)
        .map_err(|error| Error::Abif(format!("FWO_.1 is not ASCII: {error}")))?
        .to_owned();
    let mut seen = [false; 4];
    for base in channel_order.chars() {
        let index = channel_index(base)
            .ok_or_else(|| Error::Abif("FWO_.1 is not an A/C/G/T permutation".into()))?;
        if seen[index] {
            return Err(Error::Abif("FWO_.1 repeats a channel".into()));
        }
        seen[index] = true;
    }

    let mut raw_channels: [Vec<i32>; 4] = std::array::from_fn(|_| Vec::new());
    for (index, number) in (9_u32..=12).enumerate() {
        let entry = abif.required(b"DATA", number)?;
        require_layout(entry, TYPE_SHORT, 2)?;
        raw_channels[index] = decode_i16(&abif, entry)?;
    }
    let sample_count = raw_channels[0].len();
    if sample_count == 0
        || raw_channels
            .iter()
            .any(|channel| channel.len() != sample_count)
    {
        return Err(Error::Abif(
            "DATA.9-12 channels must be non-empty and equally sized".into(),
        ));
    }
    let mut channels: [Vec<i32>; 4] = std::array::from_fn(|_| Vec::new());
    for (source_index, base) in channel_order.chars().enumerate() {
        let target_index =
            channel_index(base).ok_or_else(|| Error::Abif("invalid channel order".into()))?;
        channels[target_index] = std::mem::take(&mut raw_channels[source_index]);
    }

    let ploc_entry = abif.required(b"PLOC", 2)?;
    require_layout(ploc_entry, TYPE_SHORT, 2)?;
    let base_locations: Vec<usize> = decode_i16(&abif, ploc_entry)?
        .into_iter()
        .map(|value| {
            usize::try_from(value)
                .map_err(|_| Error::Abif("PLOC.2 contains a negative position".into()))
        })
        .collect::<Result<_>>()?;
    if base_locations.is_empty() {
        return Err(Error::Abif("PLOC.2 is empty".into()));
    }
    for pair in base_locations.windows(2) {
        if pair[0] >= pair[1] {
            return Err(Error::Abif(
                "PLOC.2 positions must be strictly increasing".into(),
            ));
        }
    }
    if base_locations
        .iter()
        .any(|position| *position >= sample_count)
    {
        return Err(Error::Abif(
            "PLOC.2 position lies outside channel samples".into(),
        ));
    }

    let primary = decode_optional_string(&abif, b"PBAS", 2)?;
    let qualities = decode_optional_bytes(&abif, b"PCON", 2)?;
    let source_name = path
        .file_name()
        .and_then(|name| name.to_str())
        .ok_or_else(|| Error::Abif("AB1 file name is not valid UTF-8".into()))?
        .to_owned();

    Ok(Chromatogram {
        source_name,
        source_sha256,
        channels,
        base_locations,
        vendor: VendorEvidence { primary, qualities },
    })
}

fn require_layout(entry: &AbifEntry, element_type: u16, element_size: usize) -> Result<()> {
    if entry.element_type != element_type || entry.element_size != element_size {
        return Err(Error::Abif(format!(
            "tag {}.{} has unsupported element type/size {}/{}",
            String::from_utf8_lossy(&entry.tag),
            entry.number,
            entry.element_type,
            entry.element_size
        )));
    }
    Ok(())
}

fn decode_i16(abif: &AbifFile, entry: &AbifEntry) -> Result<Vec<i32>> {
    let payload = abif.payload(entry)?;
    let reader = Reader::new(payload);
    (0..entry.element_count)
        .map(|index| reader.i16(index * 2).map(i32::from))
        .collect()
}

fn decode_optional_string(abif: &AbifFile, tag: &[u8; 4], number: u32) -> Result<Option<String>> {
    let Some(entry) = abif.optional(tag, number)? else {
        return Ok(None);
    };
    require_layout(entry, TYPE_CHAR, 1)?;
    let payload = abif.payload(entry)?;
    let end = payload
        .iter()
        .position(|byte| *byte == 0)
        .unwrap_or(payload.len());
    let text = std::str::from_utf8(&payload[..end])
        .map_err(|error| Error::Abif(format!("vendor base string is not ASCII: {error}")))?;
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
        return Err(Error::Abif(
            "vendor base string contains a non-IUPAC symbol".into(),
        ));
    }
    Ok(Some(text.to_owned()))
}

fn decode_optional_bytes(abif: &AbifFile, tag: &[u8; 4], number: u32) -> Result<Option<Vec<u8>>> {
    let Some(entry) = abif.optional(tag, number)? else {
        return Ok(None);
    };
    if entry.element_size != 1 || !matches!(entry.element_type, TYPE_BYTE | TYPE_CHAR) {
        return Err(Error::Abif(format!(
            "tag {}.{} has unsupported element type/size {}/{}",
            String::from_utf8_lossy(&entry.tag),
            entry.number,
            entry.element_type,
            entry.element_size
        )));
    }
    Ok(Some(abif.payload(entry)?.to_vec()))
}

const fn channel_index(base: char) -> Option<usize> {
    match base {
        'A' => Some(0),
        'C' => Some(1),
        'G' => Some(2),
        'T' => Some(3),
        _ => None,
    }
}
