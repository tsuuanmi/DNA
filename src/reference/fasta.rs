//! Strict one-record plain FASTA loading.

use std::fs;
use std::path::Path;

use noodles_fasta as fasta;

use super::{MAX_REFERENCE_BYTES, MAX_REFERENCE_LENGTH};
use crate::checksum::hex_sha256;
use crate::error::{Error, FastaError, Result};
use crate::model::reference::{Reference, ReferenceTopology};

/// Loads one normalized reference record.
pub(crate) fn load(path: &Path, topology: ReferenceTopology) -> Result<Reference> {
    let metadata = fs::metadata(path).map_err(|source| Error::Read {
        kind: "reference",
        path: path.to_path_buf(),
        source,
    })?;
    if metadata.len() == 0 || metadata.len() > MAX_REFERENCE_BYTES as u64 {
        return Err(FastaError::FileSize {
            bytes: metadata.len(),
            maximum: MAX_REFERENCE_BYTES,
        }
        .into());
    }

    let bytes = fs::read(path).map_err(|source| Error::Read {
        kind: "reference",
        path: path.to_path_buf(),
        source,
    })?;
    if bytes.is_empty() || bytes.len() > MAX_REFERENCE_BYTES {
        return Err(FastaError::FileSize {
            bytes: u64::try_from(bytes.len()).unwrap_or(u64::MAX),
            maximum: MAX_REFERENCE_BYTES,
        }
        .into());
    }

    std::str::from_utf8(&bytes).map_err(|error| FastaError::NotUtf8 {
        part: "reference",
        error,
    })?;

    let mut reader = fasta::io::Reader::new(bytes.as_slice());
    let mut records = reader.records();
    let record = records
        .next()
        .transpose()
        .map_err(FastaError::Parse)?
        .ok_or(FastaError::Empty)?;

    if records.next().is_some() {
        return Err(FastaError::MultipleRecords.into());
    }

    let name = std::str::from_utf8(record.name())
        .map_err(|error| FastaError::NotUtf8 {
            part: "reference identifier",
            error,
        })?
        .trim();
    if name.is_empty() {
        return Err(FastaError::MissingIdentifier.into());
    }

    let raw_sequence =
        std::str::from_utf8(record.sequence().as_ref()).map_err(|error| FastaError::NotUtf8 {
            part: "reference sequence",
            error,
        })?;
    let mut sequence = String::new();
    for character in raw_sequence
        .chars()
        .filter(|character| !character.is_whitespace())
    {
        let base = character.to_ascii_uppercase();
        if !matches!(base, 'A' | 'C' | 'G' | 'T' | 'N') {
            return Err(FastaError::UnsupportedBase { base: character }.into());
        }
        sequence.push(base);
    }

    if sequence.is_empty() || sequence.len() > MAX_REFERENCE_LENGTH {
        return Err(FastaError::Length {
            length: sequence.len(),
            maximum: MAX_REFERENCE_LENGTH,
        }
        .into());
    }

    let sequence_sha256 = hex_sha256(sequence.as_bytes());
    Ok(Reference {
        name: name.to_owned(),
        sequence,
        topology,
        sequence_sha256,
    })
}

#[cfg(test)]
mod tests {
    use std::fs;
    use std::path::PathBuf;

    use tempfile::{TempDir, tempdir};

    use crate::error::FastaError;

    use super::*;

    fn write_reference(contents: &[u8]) -> Result<(TempDir, PathBuf)> {
        let directory = tempdir().map_err(|source| Error::Output {
            path: "temporary directory".into(),
            source,
        })?;
        let path = directory.path().join("ref.fa");
        fs::write(&path, contents).map_err(|source| Error::Output {
            path: path.clone(),
            source,
        })?;
        Ok((directory, path))
    }

    #[test]
    fn loads_one_record_and_preserves_dna_reference_semantics() -> Result<()> {
        let (_directory, path) =
            write_reference(b">rCRS Homo sapiens mitochondrial reference\nacgt\nn\n")?;

        let reference = load(&path, ReferenceTopology::Circular)?;

        assert_eq!(reference.name, "rCRS");
        assert_eq!(reference.sequence, "ACGTN");
        assert_eq!(reference.topology, ReferenceTopology::Circular);
        assert_eq!(reference.sequence_sha256, hex_sha256(b"ACGTN"));
        Ok(())
    }

    #[test]
    fn rejects_multiple_records() -> Result<()> {
        let (_directory, path) = write_reference(b">one\nACGT\n>two\nACGT\n")?;
        assert!(matches!(
            load(&path, ReferenceTopology::Linear),
            Err(Error::Fasta(FastaError::MultipleRecords))
        ));
        Ok(())
    }

    #[test]
    fn rejects_a_blank_header_line_as_unparseable() -> Result<()> {
        let (_directory, path) = write_reference(b">   \nACGT\n")?;
        assert!(matches!(
            load(&path, ReferenceTopology::Linear),
            Err(Error::Fasta(FastaError::Parse(_)))
        ));
        Ok(())
    }

    #[test]
    fn rejects_identifier_made_only_of_unicode_whitespace() -> Result<()> {
        let (_directory, path) = write_reference(">\u{a0}\nACGT\n".as_bytes())?;
        assert!(matches!(
            load(&path, ReferenceTopology::Linear),
            Err(Error::Fasta(FastaError::MissingIdentifier))
        ));
        Ok(())
    }

    #[test]
    fn rejects_empty_sequence() -> Result<()> {
        let (_directory, path) = write_reference(b">ref\n")?;
        assert!(matches!(
            load(&path, ReferenceTopology::Linear),
            Err(Error::Fasta(FastaError::Length { length: 0, .. }))
        ));
        Ok(())
    }

    #[test]
    fn rejects_unsupported_reference_base() -> Result<()> {
        let (_directory, path) = write_reference(b">ref\nACGR\n")?;
        assert!(matches!(
            load(&path, ReferenceTopology::Linear),
            Err(Error::Fasta(FastaError::UnsupportedBase { base: 'R' }))
        ));
        Ok(())
    }

    #[test]
    fn rejects_reference_above_length_cap() -> Result<()> {
        let sequence = "A".repeat(MAX_REFERENCE_LENGTH + 1);
        let fasta = format!(">ref\n{sequence}\n");
        let (_directory, path) = write_reference(fasta.as_bytes())?;

        assert!(matches!(
            load(&path, ReferenceTopology::Linear),
            Err(Error::Fasta(FastaError::Length { length, maximum }))
                if length == MAX_REFERENCE_LENGTH + 1 && maximum == MAX_REFERENCE_LENGTH
        ));
        Ok(())
    }

    #[test]
    fn rejects_non_utf8_reference() -> Result<()> {
        let (_directory, path) = write_reference(&[b'>', b'r', b'e', b'f', b'\n', 0xff, b'\n'])?;
        assert!(matches!(
            load(&path, ReferenceTopology::Linear),
            Err(Error::Fasta(FastaError::NotUtf8 {
                part: "reference",
                ..
            }))
        ));
        Ok(())
    }
}
