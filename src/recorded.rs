//! A sample analysed from its reads' recorded evidence in memory
//! (PROP-0011): a trace's `dna.read_evidence/v1` document as `dna evidence`
//! writes it, and a sample analysed from such documents as
//! `dna sample --evidence` analyses it from files. `dna-eval`'s evidence
//! cache reads and writes through it.

use std::num::NonZeroUsize;
use std::path::Path;

use crate::error::Result;
use crate::input::evidence::load_sample_documents;
use crate::pipeline::plugins;
use crate::sample_analysis::reads::Reads;
use crate::sample_analysis::{CallerScope, Loaded, Outcome, SampleAnalysis, analyze};

/// The read evidence document (`dna.read_evidence/v1`) of one trace under
/// the configuration at `config`, exactly as `dna evidence` writes it, in
/// memory (PROP-0011). The read stage reads no profile and no reference, so
/// the document depends on the trace, the configuration, and the build alone.
///
/// # Errors
///
/// As `dna evidence`: an unreadable or invalid trace or configuration, or a
/// failure of the Sanger plugin other than setting the read aside, which the
/// document records.
pub fn read_evidence_document(trace: &Path, config: &Path) -> Result<Vec<u8>> {
    crate::pipeline::evidence_document(trace, config)
}

/// What one sample analysis from recorded evidence yields: the analysis, or,
/// when every read was set aside, how many reads were.
#[derive(Debug, Clone)]
pub enum EvidenceOutcome {
    /// At least one read was admitted.
    Analyzed(Box<SampleAnalysis>),
    /// Every read was set aside; the number of distinct reads.
    Rejected {
        /// Distinct reads set aside.
        reads: usize,
    },
}

/// Runs a sample analysis from its reads' read evidence `documents` in
/// memory, each named for its errors, as `dna sample --evidence` does from
/// files (SRS-API-019): the core places every recorded read exactly as the
/// read of its trace. Reads from recorded evidence have no Sanger products,
/// so the analysis's `reads` and `rejected` are empty; the sample evidence
/// lists the rejected reads.
///
/// # Errors
///
/// As [`analyze_sanger_outcome`](crate::sample_analysis::analyze_sanger_outcome), and when a document is malformed or was
/// not written by this build's plugins under this configuration.
///
/// # Panics
///
/// As [`analyze_sanger`](crate::sample_analysis::analyze_sanger).
pub fn analyze_evidence_outcome(
    documents: &[(String, Vec<u8>)],
    reference: &Path,
    config: &Path,
    workers: NonZeroUsize,
) -> Result<EvidenceOutcome> {
    let inputs = load_sample_documents(documents, reference, config, plugins::EVIDENCE)?;
    let loaded = Loaded {
        config: inputs.config,
        profile: inputs.profile,
        reference: inputs.reference,
    };
    Ok(
        match analyze(
            Reads::Recorded(&inputs.reads),
            loaded,
            &CallerScope::current(),
            workers,
        )? {
            Outcome::Analyzed(analysis) => EvidenceOutcome::Analyzed(analysis),
            Outcome::Rejected { rejections, .. } => EvidenceOutcome::Rejected {
                reads: rejections.len(),
            },
        },
    )
}
