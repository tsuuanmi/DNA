//! Bounded deterministic affine-gap alignment and strand selection.

mod config;

pub use config::{AlignmentConfig, MAX_ALIGNMENT_CELLS, RawAlignmentConfig};

mod canonical;
mod exact;
mod gotoh;
mod orient;
mod runs;
mod scoring;
mod traceback;

pub(crate) use orient::align_best;
