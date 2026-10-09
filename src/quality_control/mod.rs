//! Uncalibrated relative quality scoring and the trim interval derived from
//! the read's callable span.

mod penalty;
mod quality;
mod trim;

pub(crate) use trim::analyze;
