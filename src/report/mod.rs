//! Typed deterministic JSON reporting and atomic publication.

mod atomic;
mod basecall;
mod call;
mod callability;
mod consensus;
mod json;
mod notation;
mod review;
mod sample;
mod sanger_call;
mod signal;
mod variant;

pub(crate) use atomic::publish;
pub(crate) use basecall::{CompletedBasecall, build as build_basecall};
pub(crate) use call::{CompletedCall, build as build_call};
pub(crate) use consensus::{
    CompletedConsensus, build as build_consensus, fasta as consensus_fasta,
    segments_fasta as consensus_segments_fasta,
};
pub(crate) use json::{CompletedAnalysis, build_analysis, serialize};
pub(crate) use notation::{
    CompletedNotation, NamedRead, ReadRepresentation, SampleNotation, build as build_notation,
    calls as notation_calls,
};
pub(crate) use review::{
    CompletedReview, PLUGIN as REVIEW, Reviewed, build as build_review, checks as review_checks,
};
pub(crate) use sample::{CompletedSampleEvidence, SangerSampleEvidence, build as build_sample};
