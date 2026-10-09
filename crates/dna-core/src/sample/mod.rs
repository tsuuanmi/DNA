//! Sample-level evidence aggregation in reference-coordinate and variant space.

mod config;

pub use config::{RawSampleReconciliationConfig, SampleReconciliationConfig};

mod aggregate;
mod call_evidence;
mod contribution;
mod coverage;
mod loci;
mod nucleotide_support;
mod overlap;
mod profile_geometry;
mod variants;

pub use aggregate::aggregate;
