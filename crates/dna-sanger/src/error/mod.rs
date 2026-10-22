//! Typed failures of the Sanger modality: decoding a trace, re-calling its
//! bases, processing its signal, its callability, its quality control, and
//! assembling its core evidence. Each renders with the stage prefix it has
//! always had, so command-line messages and logs do not depend on which
//! crate owns the failure.

mod abif;
mod basecalling;
mod callability;
mod locus;
mod quality_control;
mod read_evidence;
mod run_signal;
mod signal;

pub use abif::{AbifError, Tag};
pub use basecalling::BasecallingError;
pub use callability::CallabilityError;
pub use locus::LocusWindowError;
pub use quality_control::QualityControlError;
pub use read_evidence::ReadEvidenceError;
pub use run_signal::RunSignalError;
pub use signal::SignalError;

use dna_kernel::error::{ConfigError, EvidenceError};

/// Result type returned by the Sanger modality.
pub type Result<T> = std::result::Result<T, SangerError>;

/// A failure of the Sanger modality.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum SangerError {
    /// The trace is not a valid ABIF file.
    #[error("invalid ABIF input: {0}")]
    Abif(AbifError),
    /// Base re-calling failed.
    #[error("base re-calling failed: {0}")]
    Basecalling(BasecallingError),
    /// Signal processing failed.
    #[error("signal processing failed: {0}")]
    Signal(SignalError),
    /// Read callability failed.
    #[error("read callability failed: {0}")]
    Callability(CallabilityError),
    /// Quality control failed.
    #[error("quality control failed: {0}")]
    QualityControl(QualityControlError),
    /// A run signal cannot be used.
    #[error("invalid run signal: {0}")]
    RunSignal(RunSignalError),
    /// A recorded read evidence document cannot be read back.
    #[error("invalid read evidence: {0}")]
    ReadEvidence(ReadEvidenceError),
    /// The Sanger products of one read disagree while its core evidence is
    /// assembled.
    #[error("{0}")]
    Products(ProductsError),
    /// A failure of the shared kernel: configuration or evidence validation.
    #[error(transparent)]
    Kernel(#[from] dna_kernel::error::Error),
}

/// The Sanger products of one read disagree. Each renders as the stage the
/// core reported it under before the modality owned it.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum ProductsError {
    /// Signal loci and quality records disagree in count.
    #[error(
        "alignment failed: signal/quality call count mismatch: {loci} loci, {qualities} quality records"
    )]
    CallCountMismatch {
        /// Signal locus count.
        loci: usize,
        /// Quality record count.
        qualities: usize,
    },
    /// The trim interval does not fit the locus profiles.
    #[error("alignment failed: invalid trim interval {start}..{end} for {profiles} locus profiles")]
    InvalidTrim {
        /// Trim start (0-based).
        start: usize,
        /// Trim end (0-based, exclusive).
        end: usize,
        /// Locus profile count.
        profiles: usize,
    },
    /// The retained sequence and its profiles disagree in length.
    #[error(
        "alignment failed: retained sequence/profile length mismatch: {bases} bases, {profiles} profiles"
    )]
    RetainedLengthMismatch {
        /// Retained base count.
        bases: usize,
        /// Retained profile count.
        profiles: usize,
    },
    /// The callability mask and quality evidence disagree in call count.
    #[error("alignment failed: callability mask and quality evidence disagree in call count")]
    MaskCountMismatch,
    /// A call index has no base call.
    #[error("variant calling failed: variant filter references missing call index {index}")]
    MissingCall {
        /// Call index (0-based).
        index: usize,
    },
    /// A call index does not match its call and quality records.
    #[error(
        "variant calling failed: variant filter call index {index} does not match call/quality records"
    )]
    CallMismatch {
        /// Call index (0-based).
        index: usize,
    },
    /// A call has no channel peaks.
    #[error("variant calling failed: call index {index} has no channel peaks")]
    NoChannelPeaks {
        /// Call index (0-based).
        index: usize,
    },
    /// A reported call has no quality record.
    #[error("failed to assemble analysis report: variant references missing quality index {index}")]
    MissingQuality {
        /// Call index (0-based).
        index: usize,
    },
    /// A reported call's call and quality records disagree on their own
    /// index.
    #[error(
        "failed to assemble analysis report: variant call index {index} does not match call/quality records"
    )]
    IndexMismatch {
        /// Call index (0-based).
        index: usize,
    },
    /// A reported call has no primary-event peak evidence.
    #[error(
        "failed to assemble analysis report: variant call index {index} lacks primary-event peak evidence"
    )]
    MissingPeakEvidence {
        /// Call index (0-based).
        index: usize,
    },
}

impl From<AbifError> for SangerError {
    fn from(error: AbifError) -> Self {
        Self::Abif(error)
    }
}

impl From<BasecallingError> for SangerError {
    fn from(error: BasecallingError) -> Self {
        Self::Basecalling(error)
    }
}

impl From<SignalError> for SangerError {
    fn from(error: SignalError) -> Self {
        Self::Signal(error)
    }
}

impl From<CallabilityError> for SangerError {
    fn from(error: CallabilityError) -> Self {
        Self::Callability(error)
    }
}

impl From<QualityControlError> for SangerError {
    fn from(error: QualityControlError) -> Self {
        Self::QualityControl(error)
    }
}

impl From<RunSignalError> for SangerError {
    fn from(failure: RunSignalError) -> Self {
        Self::RunSignal(failure)
    }
}

impl From<ReadEvidenceError> for SangerError {
    fn from(failure: ReadEvidenceError) -> Self {
        Self::ReadEvidence(failure)
    }
}

impl From<ProductsError> for SangerError {
    fn from(error: ProductsError) -> Self {
        Self::Products(error)
    }
}

impl From<ConfigError> for SangerError {
    fn from(error: ConfigError) -> Self {
        Self::Kernel(error.into())
    }
}

impl From<EvidenceError> for SangerError {
    fn from(error: EvidenceError) -> Self {
        Self::Kernel(error.into())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn renders_stage_prefixed_messages_for_cli_and_logs() {
        let cases: Vec<(SangerError, &str)> = vec![
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
                ProductsError::InvalidTrim {
                    start: 5,
                    end: 4,
                    profiles: 10,
                }
                .into(),
                "alignment failed: invalid trim interval 5..4 for 10 locus profiles",
            ),
            (
                ProductsError::NoChannelPeaks { index: 3 }.into(),
                "variant calling failed: call index 3 has no channel peaks",
            ),
            (
                ConfigError::NotFinitePositive {
                    key: "signal_processing.minimum_primary_snr",
                }
                .into(),
                "invalid configuration value: signal_processing.minimum_primary_snr must be finite and positive",
            ),
        ];
        for (error, expected) in cases {
            assert_eq!(error.to_string(), expected);
        }
    }
}
