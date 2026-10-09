//! Strict, reproducible scientific configuration.

mod defaults;
mod load;
mod types;

pub(crate) use defaults::{
    MAX_ABIF_BYTES, MAX_ALIGNMENT_CELLS, MAX_REFERENCE_BYTES, MAX_REFERENCE_LENGTH,
};
pub(crate) use load::{load_path, resolve_path};
pub(crate) use types::{
    AlignmentConfig, BasecallingConfig, Config, QualityControlConfig, SampleReconciliationConfig,
    SangerEvidenceConfig, SignalProcessingConfig, VariantCallingConfig,
};
