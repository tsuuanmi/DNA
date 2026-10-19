//! Run signals (`dna.signal/v1`) as input to `review` (PROP-0009).

use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::Path;

use serde::Deserialize;

use crate::error::{KernelError, Result, RunSignalError};
use dna_kernel::checksum::hex_sha256;
use dna_sanger::run_signal::{Histogram, Key, RunSignal, SELECTION};

use super::require_regular_file;

/// Largest accepted run signal document.
const MAX_SIGNAL_BYTES: u64 = 64 * 1024 * 1024;

/// The contract this reader accepts.
const SCHEMA_VERSION: &str = "dna.signal/v1";

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct RawSignal {
    schema_version: String,
    run_id: String,
    provenance: RawProvenance,
    traces: Vec<String>,
    reads: usize,
    #[serde(rename = "rejected_reads")]
    _rejected_reads: usize,
    calls: u64,
    keys: Vec<RawKey>,
}

#[derive(Deserialize)]
struct RawProvenance {
    selection: u32,
    flank: usize,
    resolution: u16,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct RawKey {
    context: String,
    base: char,
    second: char,
    shares: Vec<[u64; 2]>,
}

/// Reads and validates one run signal document, with its SHA-256.
pub(crate) fn load(path: &Path) -> Result<(RunSignal, String)> {
    require_regular_file(path, "run signal")?;
    let size = fs::metadata(path)
        .map_err(|source| KernelError::Read {
            kind: "run signal",
            path: path.to_path_buf(),
            source,
        })?
        .len();
    if size > MAX_SIGNAL_BYTES {
        return Err(KernelError::Path {
            kind: "run signal",
            path: path.to_path_buf(),
            reason: "file exceeds the run signal size limit",
        }
        .into());
    }
    let bytes = fs::read(path).map_err(|source| KernelError::Read {
        kind: "run signal",
        path: path.to_path_buf(),
        source,
    })?;
    let sha256 = hex_sha256(&bytes);
    Ok((parse(&bytes, path)?, sha256))
}

fn parse(bytes: &[u8], path: &Path) -> std::result::Result<RunSignal, RunSignalError> {
    let raw: RawSignal = serde_json::from_slice(bytes).map_err(|error| RunSignalError::Parse {
        path: path.to_path_buf(),
        source: Box::new(error),
    })?;
    if raw.schema_version != SCHEMA_VERSION {
        return Err(RunSignalError::UnsupportedSchema {
            found: raw.schema_version,
        });
    }
    if raw.provenance.selection != SELECTION {
        return Err(RunSignalError::Selection {
            found: raw.provenance.selection,
            expected: SELECTION,
        });
    }
    let RawProvenance {
        flank, resolution, ..
    } = raw.provenance;
    let canonical = |base: char| matches!(base, 'A' | 'C' | 'G' | 'T');
    let mut keys = BTreeMap::new();
    let mut total = 0_u64;
    for key in raw.keys {
        let context = key.context.chars().collect::<Vec<_>>();
        if context.len() != 2 * flank + 1
            || !context.iter().all(|&base| canonical(base))
            || context[flank] != key.base
            || !canonical(key.second)
            || key.second == key.base
        {
            return Err(RunSignalError::Key {
                context: key.context,
            });
        }
        let mut histogram = Histogram::new();
        for [share, count] in key.shares {
            let Ok(share) = u16::try_from(share) else {
                return Err(RunSignalError::Key {
                    context: key.context,
                });
            };
            if share > resolution || count == 0 {
                return Err(RunSignalError::Key {
                    context: key.context,
                });
            }
            total = total.saturating_add(count);
            *histogram.entry(share).or_default() += count;
        }
        keys.insert(
            Key {
                context: key.context,
                base: key.base,
                second: key.second,
            },
            histogram,
        );
    }
    if total != raw.calls.saturating_mul(3) {
        return Err(RunSignalError::Total {
            calls: raw.calls,
            total,
        });
    }
    Ok(RunSignal {
        run_id: raw.run_id,
        selection: raw.provenance.selection,
        flank,
        resolution,
        traces: raw.traces.into_iter().collect::<BTreeSet<_>>(),
        reads: raw.reads,
        calls: raw.calls,
        keys,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    const DOCUMENT: &str = r#"{"schema_version":"dna.signal/v1","run_id":"R","provenance":{"configuration_sha256":"0","plugins":[],"selection":2,"flank":2,"resolution":1000},"traces":["a"],"reads":1,"rejected_reads":0,"calls":1,"keys":[{"context":"ACGTA","base":"G","second":"A","shares":[[0,1]]},{"context":"ACGTA","base":"G","second":"C","shares":[[12,1]]},{"context":"ACGTA","base":"G","second":"T","shares":[[3,1]]}]}"#;

    // SRS-API-014 SRS-API-016
    #[test]
    fn reads_a_run_signal_and_refuses_another_selection_or_shape()
    -> std::result::Result<(), RunSignalError> {
        let signal = parse(DOCUMENT.as_bytes(), Path::new("r.json"))?;
        assert_eq!(signal.calls, 1);
        assert_eq!(signal.keys.len(), 3);
        let other = DOCUMENT.replace("\"selection\":2", "\"selection\":1");
        assert!(matches!(
            parse(other.as_bytes(), Path::new("r.json")),
            Err(RunSignalError::Selection { found: 1, .. })
        ));
        let short = DOCUMENT.replace("\"calls\":1", "\"calls\":2");
        assert!(matches!(
            parse(short.as_bytes(), Path::new("r.json")),
            Err(RunSignalError::Total { .. })
        ));
        // A context of another width than the flank's, or a share above the
        // resolution, is not of this run signal.
        let narrow = DOCUMENT.replace("\"flank\":2", "\"flank\":1");
        assert!(matches!(
            parse(narrow.as_bytes(), Path::new("r.json")),
            Err(RunSignalError::Key { .. })
        ));
        let coarse = DOCUMENT.replace("\"resolution\":1000", "\"resolution\":10");
        assert!(matches!(
            parse(coarse.as_bytes(), Path::new("r.json")),
            Err(RunSignalError::Key { .. })
        ));
        Ok(())
    }
}
