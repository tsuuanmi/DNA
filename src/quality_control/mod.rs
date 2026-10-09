//! Uncalibrated relative quality scoring and the trim interval derived from
//! the read's callable span.

mod config;

pub(crate) use config::{QualityControlConfig, RawQualityControlConfig};

mod penalty;
mod quality;
mod trim;

pub(crate) use trim::analyze;
