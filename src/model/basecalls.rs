//! DNA-derived call records and peak evidence.

use serde::Serialize;

use crate::model::nucleotide::Nucleotide;

/// How a channel value was selected inside a call window.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum PeakSource {
    /// A positive local maximum was found.
    LocalMaximum,
    /// No positive local maximum existed; the PLOC sample was used.
    PlocFallback,
}

/// Strongest evidence for one channel in one base window.
#[derive(Debug, Clone, Copy)]
pub struct ChannelPeak {
    pub(crate) base: Nucleotide,
    pub(crate) height: i32,
    pub(crate) position_0based: usize,
    pub(crate) source: PeakSource,
}

/// Raw A/C/G/T observations at the uniquely strongest primary event.
///
/// Channel order follows `Nucleotide::ALL`: A, C, G, T.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct PrimaryPeakEvidence {
    /// Selected sample position of the uniquely strongest primary channel.
    pub(crate) position_0based: usize,
    /// Raw analyzed A/C/G/T channel values at `position_0based`.
    pub(crate) channel_heights: [i32; 4],
}

/// One signal-derived base call.
#[derive(Debug, Clone)]
pub struct BaseCall {
    pub(crate) index_0based: usize,
    pub(crate) ploc_0based: usize,
    pub(crate) window_start_0based: usize,
    pub(crate) window_end_0based_exclusive: usize,
    pub(crate) peaks: [ChannelPeak; 4],
    pub(crate) primary_peak_evidence: Option<PrimaryPeakEvidence>,
    pub(crate) primary: char,
    pub(crate) ambiguity: char,
    pub(crate) qualifying_channels: Vec<Nucleotide>,
    pub(crate) vendor_agrees: Option<bool>,
}

/// Ordered calls and the primary sequence consumed by downstream stages.
#[derive(Debug, Clone)]
pub struct BaseCalls {
    pub(crate) calls: Vec<BaseCall>,
    pub(crate) primary_sequence: String,
}

impl BaseCalls {
    /// Number of call loci.
    pub fn len(&self) -> usize {
        self.calls.len()
    }

    /// Whether no calls were produced.
    pub fn is_empty(&self) -> bool {
        self.calls.is_empty()
    }
}
