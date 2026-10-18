//! Consensus documents (`dna.consensus/v1`) as input to `call` (PROP-0006).
//!
//! The document's segments and sites are read back into the core's consensus
//! model, so `call` reports exactly the edits the consensus decided. The
//! document's own provenance stays with it; the call checks the sample and
//! the reference identity.

use std::fs;
use std::path::Path;

use serde::Deserialize;

use crate::config::Config;
use crate::error::{ConsensusDocumentError, KernelError, Result};
use crate::input::sequence::{self, SequenceRead};
use crate::report;
use dna_core::model::consensus::{
    Consensus, ConsensusRun, ConsensusSegment, ConsensusSite, ConsensusSummary, RunEvidence,
    SiteState,
};
use dna_kernel::model::reference::Reference;
use dna_kernel::profile::Profile;

use super::{load_config, load_profile, load_reference, require_regular_file};

/// Largest accepted consensus document.
const MAX_CONSENSUS_BYTES: usize = 16 * 1024 * 1024;

/// The contract this reader accepts.
const SCHEMA_VERSION: &str = "dna.consensus/v1";

/// Scientific inputs for calling one consensus document.
pub(crate) struct ConsensusInputs {
    pub(crate) config: Config,
    pub(crate) profile: Profile,
    pub(crate) reference: Reference,
    pub(crate) consensus: Consensus,
    /// The segments as the consensus FASTA writes them, in segment order.
    pub(crate) reads: Vec<SequenceRead>,
}

#[derive(Deserialize)]
struct RawDocument {
    schema_version: String,
    sample_id: String,
    provenance: RawProvenance,
    segments: Vec<RawSegment>,
    sites: Vec<RawSite>,
}

#[derive(Deserialize)]
struct RawProvenance {
    reference: RawReference,
}

#[derive(Deserialize)]
struct RawReference {
    sha256: String,
}

#[derive(Deserialize)]
struct RawSegment {
    reference: RawInterval,
    sequence: String,
}

#[derive(Deserialize)]
struct RawInterval {
    start: usize,
    end: usize,
}

#[derive(Deserialize)]
struct RawSite {
    start: usize,
    end: usize,
    reference: String,
    call: Option<String>,
    state: String,
    supporting_reads: Vec<String>,
    opposing_reads: Vec<String>,
    uninformative_reads: Vec<String>,
    runs: Vec<RawRun>,
}

#[derive(Deserialize)]
struct RawRun {
    base: char,
    length: usize,
    length_evidence: String,
}

/// The site state a document label names.
fn state(label: &str) -> Option<SiteState> {
    [
        SiteState::Called,
        SiteState::Contested,
        SiteState::Unresolved,
    ]
    .into_iter()
    .find(|state| state.label() == label)
}

/// The run evidence a document label names.
fn evidence(label: &str) -> Option<RunEvidence> {
    [
        RunEvidence::InPhase,
        RunEvidence::AnchoredEnd,
        RunEvidence::PhaseLoss,
        RunEvidence::ReferenceFrame,
    ]
    .into_iter()
    .find(|evidence| evidence.label() == label)
}

/// Validates paths, loads configuration, profile, and reference, and reads
/// the consensus document of `sample_id` against that reference.
pub(crate) fn load_call(
    sample_id: &str,
    document: &Path,
    reference_path: &Path,
    config_path: &Path,
) -> Result<ConsensusInputs> {
    require_regular_file(document, "consensus")?;
    require_regular_file(reference_path, "reference")?;
    let config = load_config(config_path)?;
    let profile = load_profile(&config)?;
    let reference = load_reference(reference_path, &profile)?;
    let read_error = |source| KernelError::Read {
        kind: "consensus",
        path: document.to_path_buf(),
        source,
    };
    let length = fs::metadata(document).map_err(read_error)?.len();
    if length == 0 || length > MAX_CONSENSUS_BYTES as u64 {
        return Err(ConsensusDocumentError::FileSize {
            bytes: length,
            maximum: MAX_CONSENSUS_BYTES,
        }
        .into());
    }
    let bytes = fs::read(document).map_err(read_error)?;
    let consensus = parse(&bytes, document, sample_id, &reference)?;
    let reads = sequence::from_fasta(&report::consensus_segments_fasta(sample_id, &consensus))?;
    Ok(ConsensusInputs {
        config,
        profile,
        reference,
        consensus,
        reads,
    })
}

fn parse(bytes: &[u8], path: &Path, sample_id: &str, reference: &Reference) -> Result<Consensus> {
    let raw: RawDocument =
        serde_json::from_slice(bytes).map_err(|error| KernelError::ConsensusParse {
            path: path.to_path_buf(),
            source: Box::new(error),
        })?;
    if raw.schema_version != SCHEMA_VERSION {
        return Err(ConsensusDocumentError::UnsupportedSchema {
            found: raw.schema_version,
        }
        .into());
    }
    if raw.sample_id != sample_id {
        return Err(ConsensusDocumentError::SampleMismatch {
            expected: sample_id.to_owned(),
            found: raw.sample_id,
        }
        .into());
    }
    if raw.provenance.reference.sha256 != reference.sequence_sha256 {
        return Err(ConsensusDocumentError::ReferenceMismatch.into());
    }
    let length = reference.len();
    let segments = raw
        .segments
        .into_iter()
        .map(|segment| {
            let RawInterval { start, end } = segment.reference;
            if start >= end || end > length || segment.sequence.is_empty() {
                return Err(ConsensusDocumentError::Inconsistent(
                    "a segment interval is empty, reversed, or outside the reference",
                ));
            }
            Ok(ConsensusSegment {
                start_0based: start,
                end_0based: end - 1,
                sequence: segment.sequence,
            })
        })
        .collect::<std::result::Result<Vec<_>, _>>()?;
    let sites = raw
        .sites
        .into_iter()
        .map(|site| {
            if site.start == 0 || site.start > site.end || site.end > length {
                return Err(ConsensusDocumentError::Inconsistent(
                    "a site interval is empty, reversed, or outside the reference",
                ));
            }
            let unknown = || {
                ConsensusDocumentError::Inconsistent("a site has an unknown state or run evidence")
            };
            Ok(ConsensusSite {
                first_0based: site.start - 1,
                last_0based: site.end - 1,
                reference: site.reference,
                call: site.call,
                state: state(&site.state).ok_or_else(unknown)?,
                supporting: site.supporting_reads,
                opposing: site.opposing_reads,
                uninformative: site.uninformative_reads,
                runs: site
                    .runs
                    .into_iter()
                    .map(|run| {
                        Ok(ConsensusRun {
                            base: run.base,
                            length: run.length,
                            evidence: evidence(&run.length_evidence).ok_or_else(unknown)?,
                        })
                    })
                    .collect::<std::result::Result<_, _>>()?,
            })
        })
        .collect::<std::result::Result<Vec<_>, _>>()?;
    Ok(Consensus {
        segments,
        sites,
        summary: ConsensusSummary::default(),
    })
}
