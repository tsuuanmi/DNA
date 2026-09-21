//! Bounded deterministic affine-gap alignment and strand selection.

mod canonical;
mod gotoh;
mod orient;
mod scoring;
mod traceback;

pub(crate) use orient::align_best;
