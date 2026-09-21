//! Uncalibrated relative quality scoring and low-quality end trimming.

mod penalty;
mod quality;
mod trim;

pub(crate) use trim::analyze;
