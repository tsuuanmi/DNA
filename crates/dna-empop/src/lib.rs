//! DNA post-calling forensic format (ADR-0069, ADR-0077): reports
//! represented calls that break the profile's forensic (EMPOP-style) notation
//! conventions, without rewriting them, and exports an offline EMPOP query.
//!
//! It depends only on the kernel.

pub mod conformance;
pub mod export;
