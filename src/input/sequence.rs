//! Sequence modality: reviewed consensus sequences in FASTA (ADR-0069).
//!
//! Each FASTA record is one read whose bases a reviewer already resolved. The
//! adapter vouches for every base up to the record's ends, masks nothing, and
//! raises no support veto. A canonical base carries a one-hot profile; an IUPAC
//! ambiguity code is an unresolved `N` whose profile shares its weight equally
//! among the bases the code admits, and `N` carries no profile.

use std::collections::BTreeSet;
use std::fs;
use std::path::{Path, PathBuf};

use noodles_fasta as fasta;

use crate::checksum::hex_sha256;
use crate::config::{Config, MAX_REFERENCE_LENGTH, MAX_SEQUENCE_BYTES, MAX_SEQUENCE_RECORDS};
use crate::error::{Error, Result, SequenceError};
use crate::model::reference::Reference;
use crate::profile::Profile;
use crate::read_evidence::{CallEvidence, EvidenceProfile, ReadEvidence, VetoSet};

use super::{load_config, load_profile, load_reference, require_regular_file};

/// One consensus sequence read from a FASTA record.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct SequenceRead {
    /// Record identifier, unique within the call.
    pub(crate) name: String,
    /// SHA-256 of the normalized upper-case sequence.
    pub(crate) sha256: String,
    /// Upper-case bases and IUPAC codes.
    pub(crate) sequence: String,
}

/// Scientific inputs for one core-only call over consensus sequences.
pub(crate) struct CallInputs {
    pub(crate) config: Config,
    pub(crate) profile: Profile,
    pub(crate) reads: Vec<SequenceRead>,
    pub(crate) reference: Reference,
}

/// Validates and loads consensus sequences, configuration, profile, and the
/// shared reference.
pub(crate) fn load_call(
    sequence_paths: &[PathBuf],
    reference_path: &Path,
    config_path: &Path,
) -> Result<CallInputs> {
    if sequence_paths.is_empty() {
        return Err(SequenceError::NoSequences.into());
    }
    for path in sequence_paths {
        require_regular_file(path, "sequence")?;
    }
    require_regular_file(reference_path, "reference")?;
    let config = load_config(config_path)?;
    let profile = load_profile(&config)?;
    let mut reads = Vec::new();
    for path in sequence_paths {
        reads.extend(load(path)?);
        if reads.len() > MAX_SEQUENCE_RECORDS {
            return Err(SequenceError::TooManyRecords {
                maximum: MAX_SEQUENCE_RECORDS,
            }
            .into());
        }
    }
    require_unique(&reads)?;
    let reference = load_reference(reference_path, &profile)?;
    Ok(CallInputs {
        config,
        profile,
        reads,
        reference,
    })
}

/// Loads every record of one FASTA file.
fn load(path: &Path) -> Result<Vec<SequenceRead>> {
    let metadata = fs::metadata(path).map_err(|source| Error::Read {
        kind: "sequence",
        path: path.to_path_buf(),
        source,
    })?;
    if metadata.len() > MAX_SEQUENCE_BYTES as u64 {
        return Err(SequenceError::FileSize {
            bytes: metadata.len(),
            maximum: MAX_SEQUENCE_BYTES,
        }
        .into());
    }
    let bytes = fs::read(path).map_err(|source| Error::Read {
        kind: "sequence",
        path: path.to_path_buf(),
        source,
    })?;
    parse(&bytes)
}

fn parse(bytes: &[u8]) -> Result<Vec<SequenceRead>> {
    if bytes.is_empty() || bytes.len() > MAX_SEQUENCE_BYTES {
        return Err(SequenceError::FileSize {
            bytes: u64::try_from(bytes.len()).unwrap_or(u64::MAX),
            maximum: MAX_SEQUENCE_BYTES,
        }
        .into());
    }
    std::str::from_utf8(bytes).map_err(|error| SequenceError::NotUtf8 {
        part: "sequence file",
        error,
    })?;
    let mut reader = fasta::io::Reader::new(bytes);
    let mut reads = Vec::new();
    for record in reader.records() {
        let record = record.map_err(SequenceError::Parse)?;
        if reads.len() == MAX_SEQUENCE_RECORDS {
            return Err(SequenceError::TooManyRecords {
                maximum: MAX_SEQUENCE_RECORDS,
            }
            .into());
        }
        let name = std::str::from_utf8(record.name())
            .map_err(|error| SequenceError::NotUtf8 {
                part: "sequence identifier",
                error,
            })?
            .trim()
            .to_owned();
        if name.is_empty() {
            return Err(SequenceError::MissingIdentifier.into());
        }
        let raw = std::str::from_utf8(record.sequence().as_ref()).map_err(|error| {
            SequenceError::NotUtf8 {
                part: "sequence",
                error,
            }
        })?;
        let mut sequence = String::with_capacity(raw.len());
        for symbol in raw.chars().filter(|symbol| !symbol.is_whitespace()) {
            let base = symbol.to_ascii_uppercase();
            if admitted(base).is_none() {
                return Err(SequenceError::UnsupportedSymbol { name, symbol }.into());
            }
            sequence.push(base);
        }
        if sequence.is_empty() || sequence.len() > MAX_REFERENCE_LENGTH {
            return Err(SequenceError::Length {
                name,
                length: sequence.len(),
                maximum: MAX_REFERENCE_LENGTH,
            }
            .into());
        }
        reads.push(SequenceRead {
            name,
            sha256: hex_sha256(sequence.as_bytes()),
            sequence,
        });
    }
    if reads.is_empty() {
        return Err(SequenceError::Empty.into());
    }
    Ok(reads)
}

fn require_unique(reads: &[SequenceRead]) -> Result<()> {
    let mut names = BTreeSet::new();
    let mut sequences = BTreeSet::new();
    for read in reads {
        if !names.insert(read.name.as_str()) {
            return Err(SequenceError::DuplicateName {
                name: read.name.clone(),
            }
            .into());
        }
        if !sequences.insert(read.sha256.as_str()) {
            return Err(SequenceError::DuplicateSequence {
                name: read.name.clone(),
            }
            .into());
        }
    }
    Ok(())
}

/// The A/C/G/T channels an upper-case base or IUPAC code admits; `N` admits
/// all four.
const fn admitted(symbol: char) -> Option<[bool; 4]> {
    Some(match symbol {
        'A' => [true, false, false, false],
        'C' => [false, true, false, false],
        'G' => [false, false, true, false],
        'T' => [false, false, false, true],
        'R' => [true, false, true, false],
        'Y' => [false, true, false, true],
        'S' => [false, true, true, false],
        'W' => [true, false, false, true],
        'K' => [false, false, true, true],
        'M' => [true, true, false, false],
        'B' => [false, true, true, true],
        'D' => [true, false, true, true],
        'H' => [true, true, false, true],
        'V' => [true, true, true, false],
        'N' => [true, true, true, true],
        _ => return None,
    })
}

/// Builds the core evidence of one consensus sequence.
pub(crate) fn read_evidence(read: &SequenceRead) -> Result<ReadEvidence> {
    let calls = read
        .sequence
        .chars()
        .map(|symbol| {
            let channels = admitted(symbol).unwrap_or([true; 4]);
            let count = channels.iter().filter(|&&channel| channel).count();
            let canonical = count == 1;
            let profile = (count < 4).then(|| {
                let share = 1.0 / count as f64;
                EvidenceProfile {
                    weights: channels.map(|channel| if channel { share } else { 0.0 }),
                }
            });
            CallEvidence {
                base: if canonical { symbol } else { 'N' },
                profile,
                mask: None,
                vetoes: VetoSet::default(),
            }
        })
        .collect::<Vec<_>>();
    let informative = 0..calls.len();
    Ok(ReadEvidence::new(calls, informative, Vec::new())?.with_vouched_ends())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_records_and_normalizes_case() -> Result<()> {
        let reads = parse(b">HV1 reviewed\nacgT\nNY\n>HV2|68-340\nGG\n")?;
        assert_eq!(reads.len(), 2);
        assert_eq!(reads[0].name, "HV1");
        assert_eq!(reads[0].sequence, "ACGTNY");
        assert_eq!(reads[0].sha256, hex_sha256(b"ACGTNY"));
        assert_eq!(reads[1].name, "HV2|68-340");
        Ok(())
    }

    #[test]
    fn rejects_invalid_records() {
        assert!(parse(b"").is_err());
        assert!(parse(b">HV1\nAC-GT\n").is_err());
        assert!(parse(b">HV1\n\n").is_err());
        assert!(parse(b">\nACGT\n").is_err());
        let reads = |names: [&str; 2], sequences: [&str; 2]| {
            names
                .iter()
                .zip(sequences)
                .map(|(name, sequence)| SequenceRead {
                    name: (*name).into(),
                    sha256: hex_sha256(sequence.as_bytes()),
                    sequence: sequence.into(),
                })
                .collect::<Vec<_>>()
        };
        assert!(require_unique(&reads(["a", "a"], ["AC", "GT"])).is_err());
        assert!(require_unique(&reads(["a", "b"], ["AC", "AC"])).is_err());
        assert!(require_unique(&reads(["a", "b"], ["AC", "GT"])).is_ok());
    }

    #[test]
    fn builds_vouched_evidence_with_iupac_profiles() -> Result<()> {
        let read = SequenceRead {
            name: "HV1".into(),
            sha256: String::new(),
            sequence: "AYN".into(),
        };
        let evidence = read_evidence(&read)?;
        assert_eq!(evidence.ends(), crate::read_evidence::ReadEnds::Vouched);
        assert_eq!(evidence.informative(), 0..3);
        let calls = evidence.calls();
        assert_eq!(calls[0].base, 'A');
        assert_eq!(
            calls[0].profile.map(|p| p.weights),
            Some([1.0, 0.0, 0.0, 0.0])
        );
        assert_eq!(calls[1].base, 'N');
        assert_eq!(
            calls[1].profile.map(|p| p.weights),
            Some([0.0, 0.5, 0.0, 0.5])
        );
        assert_eq!(calls[2].base, 'N');
        assert_eq!(calls[2].profile, None);
        assert!(calls.iter().all(|call| call.mask.is_none()));
        Ok(())
    }
}
