//! Strict, reproducible scientific configuration.

mod defaults;
mod load;
mod types;

pub(crate) use defaults::{
    MAX_AB1_BYTES, MAX_ALIGNMENT_CELLS, MAX_REFERENCE_BYTES, MAX_REFERENCE_LENGTH,
};
pub(crate) use load::load;
pub(crate) use types::{
    AlignmentConfig, BasecallingConfig, Config, DNAProcessingConfig, QualityControlConfig,
    SampleReconciliationConfig, VariantCallingConfig,
};
