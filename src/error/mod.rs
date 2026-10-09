//! Typed failures exposed by the DNA application boundary.
//!
//! [`Error`] names the failing stage; each stage owns a `#[non_exhaustive]`
//! failure vocabulary in a submodule. Failures that depend on observed values
//! carry them as structured data; fixed rules are distinct variants, some with
//! static rule text; arithmetic overflow and internal-consistency guards carry a
//! static description of the violated invariant. Stage failures render inline
//! (`"<stage prefix>: <failure>"`) and are reached by matching, not through
//! [`std::error::Error::source`].

mod abif;
mod alignment;
mod basecalling;
mod call_evidence;
mod callability;
mod config;
mod fasta;
mod locus;
mod profile;
mod quality_control;
mod report;
mod representation;
mod sample;
mod signal;
mod variant;

use std::path::PathBuf;

pub use abif::{AbifError, Tag};
pub use alignment::AlignmentError;
pub use basecalling::BasecallingError;
pub use call_evidence::CallEvidenceError;
pub use callability::CallabilityError;
pub use config::ConfigError;
pub use fasta::FastaError;
pub use locus::LocusWindowError;
pub use profile::ProfileError;
pub use quality_control::QualityControlError;
pub use report::ReportError;
pub use representation::{NomenclatureError, NormalizationError, RepresentationError};
pub use sample::SampleError;
pub use signal::SignalError;
pub use variant::VariantError;

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
    /// The ABIF container or one of its required records is invalid.
    #[error("invalid ABIF input: {0}")]
    Abif(AbifError),
    /// The reference FASTA is invalid.
    #[error("invalid reference FASTA: {0}")]
    Fasta(FastaError),
    /// DNA-derived base re-calling failed.
    #[error("base re-calling failed: {0}")]
    Basecalling(BasecallingError),
    /// Observational signal-quality feature extraction failed.
    #[error("signal processing failed: {0}")]
    Signal(SignalError),
    /// Signal-derived read callability could not be derived.
    #[error("read callability failed: {0}")]
    Callability(CallabilityError),
    /// Quality scoring or end trimming failed.
    #[error("quality control failed: {0}")]
    QualityControl(QualityControlError),
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
    /// An operation failed and its terminal error record also could not be persisted.
    #[error("{operation}; additionally failed to persist the operation error log: {logging}")]
    OperationAndLog {
        /// The original operation failure.
        operation: Box<Error>,
        /// The failure to record it in the operation log.
        logging: Box<Error>,
    },
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
    Abif(AbifError),
    Fasta(FastaError),
    Basecalling(BasecallingError),
    Signal(SignalError),
    Callability(CallabilityError),
    QualityControl(QualityControlError),
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
                AbifError::MissingSignature.into(),
                "invalid ABIF input: missing ABIF signature",
            ),
            (
                AbifError::MissingTag {
                    tag: Tag(*b"PLOC"),
                    number: 2,
                }
                .into(),
                "invalid ABIF input: missing required tag PLOC.2",
            ),
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
                BasecallingError::LocusWindow(LocusWindowError::InvalidWindow {
                    start: 3,
                    end: 3,
                    position: 3,
                })
                .into(),
                "base re-calling failed: invalid locus window 3..3 for locus position 3",
            ),
            (
                SignalError::TooFewCalls {
                    calls: 4,
                    window: 5,
                }
                .into(),
                "signal processing failed: 4 calls are fewer than window_size_bases 5",
            ),
            (
                CallabilityError::TooFewPositions {
                    positions: 1,
                    minimum: 2,
                }
                .into(),
                "read callability failed: 1 call positions are fewer than the minimum 2",
            ),
            (
                CallabilityError::TooFewCallableCalls {
                    callable: 12,
                    minimum: 20,
                }
                .into(),
                "read callability failed: 12 callable calls are fewer than minimum_callable_calls 20",
            ),
            (
                QualityControlError::InvalidCallableSpan {
                    start: 4,
                    end: 2,
                    calls: 10,
                }
                .into(),
                "quality control failed: callable span 4..2 does not fit 10 calls",
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
                SampleError::NoAdmittedReads { rejected: 2 }.into(),
                "sample evidence failed: all 2 reads have too few callable calls to be analyzed",
            ),
            (
                SampleError::CallEvidence(CallEvidenceError::MissingPeakEvidence { index: 7 })
                    .into(),
                "sample evidence failed: variant call index 7 lacks primary-event peak evidence",
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
