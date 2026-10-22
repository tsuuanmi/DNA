//! Typed deterministic JSON reporting and atomic publication.

mod aligned_calls;
mod analysis;
mod atomic;
mod basecall;
mod call;
mod consensus;
mod evidence;
mod json;
pub(crate) mod length;
mod notation;
mod review;
mod run_context;
mod run_signal;
mod sample;
mod sanger_call;
mod variant;

pub(crate) use analysis::{CompletedAnalysis, build_analysis};
pub(crate) use atomic::publish;
pub(crate) use basecall::{CompletedBasecall, build as build_basecall};
pub(crate) use call::{CompletedCall, build as build_call};
pub(crate) use consensus::{
    CompletedConsensus, build as build_consensus, fasta as consensus_fasta,
    segments_fasta as consensus_segments_fasta,
};
pub(crate) use evidence::{CompletedEvidence, Evidenced, build as build_evidence};
pub(crate) use json::serialize;
pub(crate) use notation::{
    CompletedNotation, NamedRead, ReadRepresentation, SampleNotation, build as build_notation,
    calls as notation_calls,
};
pub(crate) use review::{
    CompletedReview, PLUGIN as REVIEW, Reviewed, build as build_review, list as review_list,
};
pub(crate) use run_context::{CompletedRunContext, Skipped, build as build_run_context};
pub(crate) use run_signal::{CompletedRunSignal, build as build_run_signal};
pub(crate) use sample::{CompletedSampleEvidence, SangerSampleEvidence, build as build_sample};
