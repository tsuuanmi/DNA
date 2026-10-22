//! DNA Sanger modality (ADR-0069): ABIF decoding, signal-derived basecalling,
//! signal processing, callability, quality control, and the adapter that turns
//! one read into `ReadEvidence`.
//!
//! It depends only on the kernel.

pub mod abif;
pub mod basecalling;
pub mod callability;
mod channel;
pub mod error;
pub mod locus;
mod median;
pub mod model;
pub mod quality_control;
pub mod read_processing;
pub mod run_signal;
pub mod second_cause;
pub mod signal_processing;
pub mod tract;
