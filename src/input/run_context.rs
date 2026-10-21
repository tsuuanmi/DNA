//! Run contexts (`dna.run_context/v1`, PROP-0014) as input to `consensus` and
//! `review`.

use std::collections::BTreeSet;
use std::fs;
use std::path::Path;

use serde::Deserialize;

use crate::error::{KernelError, Result, RunContextError};
use crate::model::run_context::{LENGTH_METHOD, RunContext};
use dna_consensus::noise::{LengthStatistics, RunStatistics};
use dna_kernel::checksum::hex_sha256;

use super::require_regular_file;

/// Largest accepted run context document.
const MAX_CONTEXT_BYTES: u64 = 64 * 1024 * 1024;

/// The contract this reader accepts.
pub(crate) const SCHEMA_VERSION: &str = "dna.run_context/v1";

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct RawContext {
    schema_version: String,
    run_id: String,
    provenance: RawProvenance,
    reads: Vec<String>,
    rejected_reads: usize,
    #[serde(rename = "samples")]
    _samples: usize,
    length: RawLength,
}

/// The provenance fields a consumer checks; the rest is recorded only.
#[derive(Deserialize)]
struct RawProvenance {
    reference: RawReference,
}

#[derive(Deserialize)]
struct RawReference {
    sha256: String,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct RawLength {
    method: u32,
    windows: Vec<[usize; 2]>,
    noise_forward: [f64; 2],
    noise_reverse: [f64; 2],
    edges: Vec<f64>,
    runs: Vec<RawRun>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct RawRun {
    run: usize,
    pairs: usize,
    sum_z2: f64,
    sum_z4: f64,
    histogram: Vec<usize>,
}

/// Reads and validates the run context at `path`.
///
/// # Errors
///
/// When the file is unreadable or too large, not a `dna.run_context/v1`
/// document, of another length method, or internally inconsistent.
pub(crate) fn load(path: &Path) -> Result<RunContext> {
    require_regular_file(path, "run context")?;
    let read_error = |source| KernelError::Read {
        kind: "run context",
        path: path.to_path_buf(),
        source,
    };
    let length = fs::metadata(path).map_err(read_error)?.len();
    if length == 0 || length > MAX_CONTEXT_BYTES {
        return Err(RunContextError::FileSize {
            bytes: length,
            maximum: MAX_CONTEXT_BYTES,
        }
        .into());
    }
    let bytes = fs::read(path).map_err(read_error)?;
    let mut context = parse(&bytes, path)?;
    context.sha256 = Some(hex_sha256(&bytes));
    Ok(context)
}

fn parse(bytes: &[u8], path: &Path) -> Result<RunContext> {
    let raw = serde_json::from_slice::<RawContext>(bytes).map_err(|source| {
        KernelError::RunContextParse {
            path: path.to_path_buf(),
            source: Box::new(source),
        }
    })?;
    if raw.schema_version != SCHEMA_VERSION {
        return Err(RunContextError::UnsupportedSchema {
            found: raw.schema_version,
        }
        .into());
    }
    if raw.length.method != LENGTH_METHOD {
        return Err(RunContextError::Method {
            section: "length",
            found: raw.length.method,
            expected: LENGTH_METHOD,
        }
        .into());
    }
    let length = raw.length;
    let edges_valid = !length.edges.is_empty()
        && length.edges.iter().all(|edge| edge.is_finite())
        && length.edges.windows(2).all(|pair| pair[0] < pair[1]);
    let noise_valid = length
        .noise_forward
        .iter()
        .chain(&length.noise_reverse)
        .all(|noise| noise.is_finite() && (0.0..=0.5).contains(noise));
    if !edges_valid || !noise_valid {
        return Err(RunContextError::Inconsistent("length noise or edges are invalid").into());
    }
    let windows_valid = length
        .windows
        .iter()
        .all(|[first, last]| 1 <= *first && first <= last)
        && length
            .windows
            .windows(2)
            .all(|pair| pair[0][1] < pair[1][0]);
    if !windows_valid {
        return Err(RunContextError::Inconsistent("length windows are invalid").into());
    }
    let mut runs = length.runs;
    runs.sort_by_key(|run| run.run);
    let [first, second] = runs.as_slice() else {
        return Err(RunContextError::Inconsistent("the length section needs runs 1 and 2").into());
    };
    if first.run != 1 || second.run != 2 {
        return Err(RunContextError::Inconsistent("the length section needs runs 1 and 2").into());
    }
    let statistics = |run: &RawRun| -> Result<RunStatistics> {
        if run.histogram.len() != length.edges.len() + 1
            || run.histogram.iter().sum::<usize>() != run.pairs
            || !run.sum_z2.is_finite()
            || run.sum_z2 < 0.0
            || !run.sum_z4.is_finite()
            || run.sum_z4 < 0.0
        {
            return Err(RunContextError::Inconsistent(
                "a length run's histogram or sum is inconsistent",
            )
            .into());
        }
        Ok(RunStatistics {
            pairs: run.pairs,
            sum_z2: run.sum_z2,
            sum_z4: run.sum_z4,
            histogram: run.histogram.clone(),
        })
    };
    let reads = raw.reads.into_iter().collect::<BTreeSet<_>>();
    if reads
        .iter()
        .any(|read| read.len() != 64 || !read.bytes().all(|byte| byte.is_ascii_hexdigit()))
    {
        return Err(RunContextError::Inconsistent("a read identity is not a SHA-256").into());
    }
    Ok(RunContext {
        run_id: raw.run_id,
        reads,
        rejected_reads: raw.rejected_reads,
        reference_sha256: raw.provenance.reference.sha256,
        windows: length.windows.clone(),
        length: LengthStatistics {
            forward: length.noise_forward,
            reverse: length.noise_reverse,
            edges: length.edges.clone(),
            runs: [statistics(first)?, statistics(second)?],
        },
        sha256: None,
    })
}
