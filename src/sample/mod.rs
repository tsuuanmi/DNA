//! Sample-level evidence aggregation in reference-coordinate and variant space.

mod aggregate;
mod call_evidence;
mod contribution;
mod coverage;
mod loci;
mod nucleotide_support;
mod overlap;
mod profile_geometry;
mod variants;

pub(crate) use aggregate::aggregate;
