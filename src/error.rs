//! The application's typed failures.
//!
//! [`Error`] is a failure of the shared kernel or of the Sanger modality,
//! each rendering with the stage prefix it has always had, so command-line
//! messages and logs name the failing stage whichever crate owns it. The
//! kernel's and the modality's failure vocabularies are re-exported here.

/// The shared kernel's failure, which [`Error::Kernel`] wraps.
pub use dna_kernel::error::Error as KernelError;
pub use dna_kernel::error::*;
pub use dna_sanger::error::{
    AbifError, BasecallingError, CallabilityError, ProductsError, QualityControlError,
    RunSignalError, SangerError, SignalError, Tag,
};

/// Result type returned by DNA operations.
pub type Result<T> = std::result::Result<T, Error>;

/// A failure of an operation.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum Error {
    /// A failure of the shared kernel or the core.
    #[error(transparent)]
    Kernel(#[from] dna_kernel::error::Error),
    /// A failure of the Sanger modality's own stages.
    #[error(transparent)]
    Sanger(SangerError),
    /// The operation failed, and so did persisting its error log.
    #[error("{operation}; additionally failed to persist the operation error log: {logging}")]
    OperationAndLog {
        /// The operation failure.
        operation: Box<Error>,
        /// The logging failure.
        logging: Box<Error>,
    },
}

/// A kernel failure the modality met (configuration or evidence validation,
/// reading a trace) is the kernel's, whichever crate returned it.
impl From<SangerError> for Error {
    fn from(failure: SangerError) -> Self {
        match failure {
            SangerError::Kernel(failure) => Self::Kernel(failure),
            failure => Self::Sanger(failure),
        }
    }
}

impl From<RunSignalError> for Error {
    fn from(failure: RunSignalError) -> Self {
        Self::Sanger(SangerError::RunSignal(failure))
    }
}

/// Each kernel failure vocabulary converts through the kernel's error.
macro_rules! kernel_failures {
    ($($failure:ty),+ $(,)?) => {
        $(
            impl From<$failure> for Error {
                fn from(failure: $failure) -> Self {
                    Self::Kernel(failure.into())
                }
            }
        )+
    };
}

kernel_failures!(
    ConfigError,
    ProfileError,
    FastaError,
    SequenceError,
    VariantsError,
    ConsensusDocumentError,
    EvidenceError,
    AlignmentError,
    VariantError,
    NormalizationError,
    NomenclatureError,
    SampleError,
    ReportError,
);

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_kernel_failure_met_by_the_modality_is_the_kernels() {
        let error = Error::from(SangerError::from(ConfigError::Constraint("rule")));
        assert!(matches!(
            &error,
            Error::Kernel(KernelError::Config(ConfigError::Constraint("rule")))
        ));
        assert_eq!(error.to_string(), "invalid configuration value: rule");
    }

    #[test]
    fn a_modality_failure_keeps_its_stage_prefix() {
        let error = Error::from(SangerError::from(AbifError::UnequalChannels));
        assert!(matches!(&error, Error::Sanger(SangerError::Abif(_))));
        assert!(error.to_string().starts_with("invalid ABIF input: "));
    }
}
