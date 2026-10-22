//! Typed failures exposed by the DNA application boundary.
//!
//! [`Error`] names the failing stage; each stage owns a `#[non_exhaustive]`
//! failure vocabulary in a submodule. Failures that depend on observed values
//! carry them as structured data; fixed rules are distinct variants, some with
//! static rule text; arithmetic overflow and internal-consistency guards carry a
//! static description of the violated invariant. Stage failures render inline
//! (`"<stage prefix>: <failure>"`) and are reached by matching, not through
//! [`std::error::Error::source`].

mod alignment;
mod call_evidence;
mod config;
mod consensus_document;
mod evidence;
mod fasta;
mod profile;
mod report;
mod representation;
mod run_context;
mod sample;
mod sequence;
mod variant;
mod variants;

use std::path::PathBuf;

pub use alignment::AlignmentError;
pub use call_evidence::CallEvidenceError;
pub use config::ConfigError;
pub use consensus_document::ConsensusDocumentError;
pub use evidence::EvidenceError;
pub use fasta::FastaError;
pub use profile::ProfileError;
pub use report::ReportError;
pub use representation::{NomenclatureError, NormalizationError, RepresentationError};
pub use run_context::RunContextError;
pub use sample::SampleError;
pub use sequence::SequenceError;
pub use variant::VariantError;
pub use variants::VariantsError;

/// Result type returned by DNA operations.
pub type Result<T> = std::result::Result<T, Error>;

/// Type-erased failure from a third-party parser or serializer.
///
/// Erasure keeps dependency types out of the public API; the message and
/// source chain remain available through [`std::error::Error`].
pub type ForeignError = Box<dyn std::error::Error + Send + Sync>;

/// A failure in a validated analysis stage.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum Error {
    /// A required path is absent, invalid, or has the wrong filesystem type.
    #[error("invalid {kind} path {path}: {reason}")]
    Path {
        /// Role of the path, for example `"AB1"` or `"output"`.
        kind: &'static str,
        /// The rejected path.
        path: PathBuf,
        /// Why the path was rejected.
        reason: &'static str,
    },
    /// A file could not be read.
    #[error("failed to read {kind} file {path}: {source}")]
    Read {
        /// Role of the file, for example `"AB1"` or `"configuration"`.
        kind: &'static str,
        /// The unreadable file.
        path: PathBuf,
        /// Underlying I/O failure.
        #[source]
        source: std::io::Error,
    },
    /// Configuration bytes were not valid TOML.
    #[error("invalid configuration {path}: {source}")]
    ConfigParse {
        /// The configuration file.
        path: PathBuf,
        /// Underlying TOML parse failure.
        #[source]
        source: ForeignError,
    },
    /// A configuration value violates the scientific contract.
    #[error("invalid configuration value: {0}")]
    Config(ConfigError),
    /// Target-profile bytes were not valid TOML.
    #[error("invalid target profile {path}: {source}")]
    ProfileParse {
        /// The profile file.
        path: PathBuf,
        /// Underlying TOML parse failure.
        #[source]
        source: ForeignError,
    },
    /// A target profile is invalid or does not apply to the supplied reference.
    #[error("invalid target profile: {0}")]
    Profile(ProfileError),
    /// The reference FASTA is invalid.
    #[error("invalid reference FASTA: {0}")]
    Fasta(FastaError),
    /// A consensus-sequence FASTA is invalid.
    #[error("invalid sequence FASTA: {0}")]
    Sequence(SequenceError),
    /// Variants-document bytes were not valid JSON of the expected shape.
    #[error("invalid variants document {path}: {source}")]
    VariantsParse {
        /// The variants document.
        path: PathBuf,
        /// Underlying JSON parse failure.
        #[source]
        source: ForeignError,
    },
    /// A variants document cannot be used for the requested operation.
    #[error("invalid variants document: {0}")]
    Variants(VariantsError),
    /// Consensus-document bytes were not valid JSON of the expected shape.
    #[error("invalid consensus document {path}: {source}")]
    ConsensusParse {
        /// The consensus document.
        path: PathBuf,
        /// Underlying JSON parse failure.
        #[source]
        source: ForeignError,
    },
    /// A consensus document cannot be called.
    #[error("invalid consensus document: {0}")]
    ConsensusDocument(ConsensusDocumentError),
    /// Run-context bytes were not valid JSON of the expected shape.
    #[error("invalid run context {path}: {source}")]
    RunContextParse {
        /// The run context document.
        path: PathBuf,
        /// Underlying JSON parse failure.
        #[source]
        source: ForeignError,
    },
    /// A run context cannot be made or used.
    #[error("invalid run context: {0}")]
    RunContext(RunContextError),
    /// A modality's per-read evidence could not enter the core.
    #[error("read evidence failed: {0}")]
    Evidence(EvidenceError),
    /// Pairwise alignment failed or was not uniquely interpretable.
    #[error("alignment failed: {0}")]
    Alignment(AlignmentError),
    /// Variant extraction or configured eligibility failed.
    #[error("variant calling failed: {0}")]
    Variant(VariantError),
    /// Post-calling variant representation normalization failed.
    #[error("variant normalization failed: {0}")]
    VariantNormalization(NormalizationError),
    /// Target-specific variant nomenclature representation failed.
    #[error("variant nomenclature failed: {0}")]
    VariantNomenclature(NomenclatureError),
    /// Sample-level read evidence is inconsistent or invalid.
    #[error("sample evidence failed: {0}")]
    Sample(SampleError),
    /// A completed model could not be assembled consistently.
    #[error("failed to assemble analysis report: {0}")]
    Report(ReportError),
    /// A completed model could not be serialized.
    #[error("failed to serialize result JSON: {0}")]
    Serialize(#[source] ForeignError),
    /// A log directory or append-only log file operation failed.
    #[error("failed to access log path {path}: {source}")]
    Log {
        /// The log directory or file.
        path: PathBuf,
        /// Underlying I/O failure.
        #[source]
        source: std::io::Error,
    },
    /// An output file operation failed.
    #[error("failed to publish output {path}: {source}")]
    Output {
        /// The output or temporary file.
        path: PathBuf,
        /// Underlying I/O failure.
        #[source]
        source: std::io::Error,
    },
}

/// Wraps each stage failure vocabulary into its [`Error`] variant.
macro_rules! stage_failures {
    ($($variant:ident($failure:ty)),+ $(,)?) => {
        $(
            impl From<$failure> for Error {
                fn from(failure: $failure) -> Self {
                    Self::$variant(failure)
                }
            }
        )+
    };
}

stage_failures!(
    Config(ConfigError),
    Profile(ProfileError),
    Fasta(FastaError),
    Sequence(SequenceError),
    Variants(VariantsError),
    ConsensusDocument(ConsensusDocumentError),
    RunContext(RunContextError),
    Evidence(EvidenceError),
    Alignment(AlignmentError),
    Variant(VariantError),
    VariantNormalization(NormalizationError),
    VariantNomenclature(NomenclatureError),
    Sample(SampleError),
    Report(ReportError),
);

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn renders_stage_prefixed_messages_for_cli_and_logs() {
        let cases: Vec<(Error, &str)> = vec![
            (
                FastaError::UnsupportedBase { base: 'R' }.into(),
                "invalid reference FASTA: unsupported reference base 'R'",
            ),
            (
                ConfigError::NotFinitePositive {
                    key: "signal_processing.minimum_primary_snr",
                }
                .into(),
                "invalid configuration value: signal_processing.minimum_primary_snr must be finite and positive",
            ),
            (
                AlignmentError::LowIdentity {
                    identity: 0.5,
                    minimum: 0.8,
                }
                .into(),
                "alignment failed: alignment callable identity 0.5000 is below 0.8000",
            ),
            (
                VariantError::ReferenceAlleleMismatch { position: 9 }.into(),
                "variant calling failed: variant reference allele disagrees with the supplied reference at position 9",
            ),
            (
                SequenceError::UnsupportedSymbol {
                    name: "HV1".into(),
                    symbol: '-',
                }
                .into(),
                "invalid sequence FASTA: unsupported sequence symbol '-' in \"HV1\"",
            ),
            (
                EvidenceError::InvalidReason("read_end").into(),
                "read evidence failed: invalid modality reason label \"read_end\"",
            ),
            (
                SampleError::NoAdmittedReads { rejected: 2 }.into(),
                "sample evidence failed: all 2 reads were rejected and none could be analyzed",
            ),
            (
                SampleError::CallEvidence(CallEvidenceError::MissingCall { index: 7 }).into(),
                "sample evidence failed: variant references missing call index 7",
            ),
            (
                ReportError::DuplicateReadName {
                    name: "read".into(),
                }
                .into(),
                "failed to assemble analysis report: sample read name \"read\" is not unique",
            ),
            (
                NormalizationError::ReferenceIdentityMismatch.into(),
                "variant normalization failed: called variants do not match the supplied reference identity",
            ),
            (
                NomenclatureError::Representation(RepresentationError::ReferenceAlleleMismatch)
                    .into(),
                "variant nomenclature failed: called variant reference allele disagrees with the supplied reference",
            ),
        ];
        for (error, expected) in cases {
            assert_eq!(error.to_string(), expected);
        }
    }
}
