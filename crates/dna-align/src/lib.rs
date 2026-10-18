//! Bounded deterministic affine-gap alignment and strand selection.

pub mod model;

mod config;

pub use config::{AlignmentConfig, MAX_ALIGNMENT_CELLS, RawAlignmentConfig};

mod canonical;
mod exact;
mod gotoh;
mod orient;
mod runs;
mod scoring;
mod traceback;

pub use orient::align_best;
