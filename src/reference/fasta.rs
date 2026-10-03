//! Strict one-record plain FASTA loading.

use std::fs;
use std::path::Path;

use crate::checksum::hex_sha256;
use crate::config::{MAX_REFERENCE_BYTES, MAX_REFERENCE_LENGTH};
use crate::error::{Error, Result};
use crate::model::reference::{Reference, ReferenceTopology};

/// Loads one normalized reference record.
pub(crate) fn load(path: &Path, topology: ReferenceTopology) -> Result<Reference> {
    let metadata = fs::metadata(path).map_err(|source| Error::Read {
        kind: "reference",
        path: path.to_path_buf(),
        source,
    })?;
    if metadata.len() == 0 || metadata.len() > MAX_REFERENCE_BYTES as u64 {
        return Err(Error::Fasta(format!(
            "reference file size {} is outside 1..={MAX_REFERENCE_BYTES} bytes",
            metadata.len()
        )));
    }
    let bytes = fs::read(path).map_err(|source| Error::Read {
        kind: "reference",
        path: path.to_path_buf(),
        source,
    })?;
    if bytes.is_empty() || bytes.len() > MAX_REFERENCE_BYTES {
        return Err(Error::Fasta(format!(
            "reference file size {} is outside 1..={MAX_REFERENCE_BYTES} bytes",
            bytes.len()
        )));
    }
    let text = std::str::from_utf8(&bytes)
        .map_err(|error| Error::Fasta(format!("reference must be UTF-8: {error}")))?;
    let mut lines = text.lines();
    let header = lines
        .next()
        .ok_or_else(|| Error::Fasta("reference is empty".into()))?;
    let name = header
        .strip_prefix('>')
        .and_then(|value| value.split_whitespace().next())
        .filter(|value| !value.is_empty())
        .ok_or_else(|| Error::Fasta("first line must contain a FASTA identifier".into()))?;
    let mut sequence = String::new();
    for line in lines {
        if line.starts_with('>') {
            return Err(Error::Fasta(
                "reference must contain exactly one record".into(),
            ));
        }
        for character in line.chars().filter(|character| !character.is_whitespace()) {
            let base = character.to_ascii_uppercase();
            if !matches!(base, 'A' | 'C' | 'G' | 'T' | 'N') {
                return Err(Error::Fasta(format!(
                    "unsupported reference base {character:?}"
                )));
            }
            sequence.push(base);
        }
    }
    if sequence.is_empty() || sequence.len() > MAX_REFERENCE_LENGTH {
        return Err(Error::Fasta(format!(
            "reference length {} is outside 1..={MAX_REFERENCE_LENGTH}",
            sequence.len()
        )));
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
        assert!(load(&path, ReferenceTopology::Linear).is_err());
        Ok(())
    }

    #[test]
    fn rejects_missing_identifier() -> Result<()> {
        let (_directory, path) = write_reference(b">   \nACGT\n")?;
        assert!(load(&path, ReferenceTopology::Linear).is_err());
        Ok(())
    }

    #[test]
    fn rejects_empty_sequence() -> Result<()> {
        let (_directory, path) = write_reference(b">ref\n")?;
        assert!(load(&path, ReferenceTopology::Linear).is_err());
        Ok(())
    }

    #[test]
    fn rejects_unsupported_reference_base() -> Result<()> {
        let (_directory, path) = write_reference(b">ref\nACGR\n")?;
        assert!(load(&path, ReferenceTopology::Linear).is_err());
        Ok(())
    }

    #[test]
    fn rejects_reference_above_length_cap() -> Result<()> {
        let sequence = "A".repeat(MAX_REFERENCE_LENGTH + 1);
        let fasta = format!(">ref\n{sequence}\n");
        let (_directory, path) = write_reference(fasta.as_bytes())?;

        assert!(load(&path, ReferenceTopology::Linear).is_err());
        Ok(())
    }

    #[test]
    fn rejects_non_utf8_reference() -> Result<()> {
        let (_directory, path) = write_reference(&[b'>', b'r', b'e', b'f', b'\n', 0xff, b'\n'])?;
        assert!(load(&path, ReferenceTopology::Linear).is_err());
        Ok(())
    }
}
