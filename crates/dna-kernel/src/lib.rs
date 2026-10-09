//! DNA kernel: the contracts every plugin family shares (ADR-0069).
//!
//! It holds the modality → core evidence contract, the public called-variant
//! contracts, the one error type of every stage, target profiles, reference
//! loading, and the plugin descriptor types with their compile-time
//! validation. It depends on no other DNA crate.

pub mod bounds;
pub mod checksum;
pub mod error;
pub mod model;
pub mod plugin;
pub mod profile;
pub mod read_evidence;
pub mod reference;
pub mod variant;
