//! DNA core caller (ADR-0069): evidence-profile alignment, per-read variant
//! calling, and sample aggregation over modality-neutral `ReadEvidence`.
//!
//! It depends only on the kernel; no modality is visible here.

pub mod alignment;
pub mod consensus;
pub mod model;
pub mod read_call;
pub mod sample;
pub mod variant_calling;
