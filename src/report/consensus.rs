//! Projection of a sample consensus into `dna.consensus/v1` and FASTA.

use std::collections::BTreeMap;

use crate::model::consensus_result::{
    ConsensusReadResult, ConsensusRejectedReadResult, ConsensusResult, ConsensusSegmentResult,
    ConsensusSiteResult, ConsensusSummaryResult,
};
use crate::model::read_observation::RejectedRead;
use crate::model::result::{IntervalResult, ReferenceResult};
use crate::model::sample_result::SampleProvenanceResult;
use crate::report::json::{project_plugins, project_profile};
use crate::report::sample::{project_rejection, reviewer_read_names};
use dna_core::model::consensus::Consensus;
use dna_core::model::sample_evidence::SampleEvidence;
use dna_kernel::error::{ReportError, Result};
use dna_kernel::model::reference::Reference;
use dna_kernel::plugin::PluginDescriptor;
use dna_kernel::profile::ProfileIdentity;

/// FASTA sequence line width.
const FASTA_WIDTH: usize = 70;

/// Inputs consumed to build one immutable consensus document.
pub(crate) struct CompletedConsensus {
    pub(crate) sample_id: String,
    pub(crate) reference: Reference,
    pub(crate) profile: ProfileIdentity,
    /// The sample's read registry and rejected reads.
    pub(crate) evidence: SampleEvidence,
    /// Rejected reads with their causes, by content SHA-256.
    pub(crate) rejected: BTreeMap<String, RejectedRead>,
    pub(crate) consensus: Consensus,
    /// Plugins of the workflow, in execution order.
    pub(crate) plugins: &'static [&'static PluginDescriptor],
}

/// Builds `dna.consensus/v1` without filesystem side effects.
pub(crate) fn build(completed: CompletedConsensus) -> Result<ConsensusResult> {
    let CompletedConsensus {
        sample_id,
        reference,
        profile,
        evidence,
        rejected,
        consensus,
        plugins,
    } = completed;
    let (names, rejected_names) = reviewer_read_names(&evidence)?;
    let reviewer = evidence
        .reads
        .iter()
        .zip(&names)
        .map(|(read, name)| (read.input_name.clone(), name.clone()))
        .collect::<BTreeMap<_, _>>();
    let rename = |reads: &[String]| -> Result<Vec<String>> {
        reads
            .iter()
            .map(|read| {
                reviewer.get(read).cloned().ok_or_else(|| {
                    ReportError::Inconsistent("consensus site names an unregistered read").into()
                })
            })
            .collect()
    };
    let rejected_reads = evidence
        .rejected_reads
        .iter()
        .zip(rejected_names)
        .map(|(read, name)| {
            let cause = rejected
                .get(&read.input_sha256)
                .ok_or(ReportError::Inconsistent(
                    "rejected consensus read lacks its rejection cause",
                ))?;
            Ok(ConsensusRejectedReadResult {
                name,
                sha256: read.input_sha256.clone(),
                rejection: project_rejection(cause.cause),
            })
        })
        .collect::<Result<_>>()?;
    let sites = consensus
        .sites
        .iter()
        .map(|site| {
            Ok(ConsensusSiteResult {
                start: site.first_0based + 1,
                end: site.last_0based + 1,
                reference: site.reference.clone(),
                call: site.call.clone(),
                state: site.state.label(),
                supporting_reads: rename(&site.supporting)?,
                opposing_reads: rename(&site.opposing)?,
                uninformative_reads: rename(&site.uninformative)?,
            })
        })
        .collect::<Result<_>>()?;
    let segments = consensus
        .segments
        .iter()
        .map(|segment| ConsensusSegmentResult {
            name: format!(
                "{sample_id}_{}-{}",
                segment.start_0based + 1,
                segment.end_0based + 1
            ),
            reference: IntervalResult {
                start: segment.start_0based,
                end: segment.end_0based + 1,
            },
            sequence: segment.sequence.clone(),
        })
        .collect();
    Ok(ConsensusResult {
        schema_version: "dna.consensus/v1",
        sample_id,
        provenance: SampleProvenanceResult {
            reference: ReferenceResult {
                name: reference.name,
                topology: reference.topology,
                sha256: reference.sequence_sha256,
            },
            configuration_sha256: evidence.configuration_sha256,
            profile: project_profile(profile),
            plugins: project_plugins(plugins),
        },
        reads: evidence
            .reads
            .iter()
            .zip(names)
            .map(|(read, name)| ConsensusReadResult {
                name,
                sha256: read.input_sha256.clone(),
            })
            .collect(),
        rejected_reads,
        segments,
        sites,
        summary: ConsensusSummaryResult {
            called_positions: consensus.summary.called_positions,
            contested_sites: consensus.summary.contested_sites,
            unresolved_positions: consensus.summary.unresolved_positions,
        },
    })
}

/// The consensus segments as FASTA, one record per segment, usable as
/// `dna call` input.
pub(crate) fn fasta(result: &ConsensusResult) -> Vec<u8> {
    let mut text = String::new();
    for segment in &result.segments {
        text.push('>');
        text.push_str(&segment.name);
        text.push('\n');
        let bytes = segment.sequence.as_bytes();
        for line in bytes.chunks(FASTA_WIDTH) {
            text.push_str(&String::from_utf8_lossy(line));
            text.push('\n');
        }
    }
    text.into_bytes()
}
