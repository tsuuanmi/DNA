//! A sample analysed from its reads' recorded evidence in memory
//! (PROP-0011): a trace's `dna.read_evidence/v1` document as `dna evidence`
//! writes it, and a sample analysed from such documents as
//! `dna sample --evidence` analyses it from files. `dna-eval`'s evidence
//! cache reads and writes through it, and its placement cache through the
//! placed-read contract (`dna.placed_read/v1`), which `dna-align` alone
//! writes and reads back.

use std::num::NonZeroUsize;
use std::path::Path;

use dna_kernel::checksum::hex_sha256;

use crate::error::Result;
use crate::input::evidence::{Recorded, load_sample_documents};
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

/// The placed-read document (`dna.placed_read/v1`, PROP-0011) of one read
/// evidence `document`, named `name` for its errors: the read's placement as
/// the core searches it on the reference at `reference` under the
/// configuration at `config`, written by `dna-align`. It records no
/// placement for a read the Sanger plugin set aside or one the core cannot
/// place; analysing such a read searches again and rejects it as before.
///
/// # Errors
///
/// When the document, the reference, or the configuration is unreadable or
/// invalid, as [`analyze_evidence_outcome`].
pub fn placed_read_document(
    name: &str,
    document: &[u8],
    reference: &Path,
    config: &Path,
) -> Result<Vec<u8>> {
    let inputs = load_sample_documents(
        &[(name.to_owned(), document.to_vec())],
        reference,
        config,
        plugins::EVIDENCE,
    )?;
    let alignment = match inputs.reads.first() {
        Some(Recorded::Read(read)) => dna_align::align_best(
            &read.evidence,
            &inputs.reference,
            &inputs.config.core.alignment,
        )
        .ok(),
        _ => None,
    };
    Ok(dna_align::placed::write(
        alignment.as_ref(),
        dna_align::placed::Identity {
            evidence_sha256: &hex_sha256(document),
            reference_sha256: &inputs.reference.sequence_sha256,
            configuration_sha256: &inputs.config.source_sha256,
        },
    )?)
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
    let documents = documents
        .iter()
        .map(|(name, document)| ReadDocuments {
            name: name.clone(),
            evidence: document.clone(),
            placed: None,
        })
        .collect::<Vec<_>>();
    analyze_placed_outcome(&documents, reference, config, workers)
}

/// One read's recorded documents: its read evidence document, named for its
/// errors, and its placed-read document when one is given.
#[derive(Debug, Clone)]
pub struct ReadDocuments {
    /// The name the read's errors carry.
    pub name: String,
    /// Its `dna.read_evidence/v1` document.
    pub evidence: Vec<u8>,
    /// Its `dna.placed_read/v1` document, when one is given.
    pub placed: Option<Vec<u8>>,
}

/// As [`analyze_evidence_outcome`], each read with its placed-read document
/// when one is given (`dna.placed_read/v1`, [`placed_read_document`]): the
/// core starts from that placement instead of searching again, and the
/// analysis equals the one that searches.
///
/// # Errors
///
/// As [`analyze_evidence_outcome`], and when a placed-read document is
/// malformed, was written for another evidence document, reference, or
/// configuration, or does not hold for its read.
///
/// # Panics
///
/// As [`analyze_sanger`](crate::sample_analysis::analyze_sanger).
pub fn analyze_placed_outcome(
    documents: &[ReadDocuments],
    reference: &Path,
    config: &Path,
    workers: NonZeroUsize,
) -> Result<EvidenceOutcome> {
    let evidence = documents
        .iter()
        .map(|read| (read.name.clone(), read.evidence.clone()))
        .collect::<Vec<_>>();
    let mut inputs = load_sample_documents(&evidence, reference, config, plugins::EVIDENCE)?;
    for (recorded, documents) in inputs.reads.iter_mut().zip(documents) {
        if let (Recorded::Read(read), Some(placed)) = (recorded, &documents.placed) {
            read.placed = dna_align::placed::read(
                placed,
                &read.evidence,
                &inputs.reference,
                dna_align::placed::Identity {
                    evidence_sha256: &hex_sha256(&documents.evidence),
                    reference_sha256: &inputs.reference.sequence_sha256,
                    configuration_sha256: &inputs.config.source_sha256,
                },
            )?;
        }
    }
    let loaded = Loaded {
        config: inputs.config,
        profile: inputs.profile,
        reference: inputs.reference,
        run_context: None,
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
