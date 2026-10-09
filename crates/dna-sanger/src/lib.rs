//! DNA Sanger modality (ADR-0069): ABIF decoding, signal-derived basecalling,
//! signal processing, callability, quality control, and the adapter that turns
//! one read into `ReadEvidence`.
//!
//! It depends only on the kernel.

pub mod abif;
pub mod basecalling;
pub mod callability;
pub mod locus;
pub mod model;
pub mod quality_control;
pub mod read_processing;
pub mod signal_processing;
