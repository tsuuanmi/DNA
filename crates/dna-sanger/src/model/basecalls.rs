//! DNA-derived call records and peak evidence.

use serde::Serialize;

use dna_kernel::model::nucleotide::Nucleotide;

/// How a channel value was selected inside a call window.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum PeakSource {
    /// A positive local maximum was found.
    LocalMaximum,
    /// No positive local maximum existed; the canonical locus sample was used.
    LocusFallback,
}

/// Strongest evidence for one channel in one base window.
#[derive(Debug, Clone, Copy)]
pub struct ChannelPeak {
    /// Channel base.
    pub base: Nucleotide,
    /// Selected peak height.
    pub height: i32,
    /// Trace sample of the selected peak (0-based).
    pub position_0based: usize,
    /// How the peak was selected.
    pub source: PeakSource,
}

/// Raw A/C/G/T observations at the uniquely strongest primary event.
///
/// Channel order follows `Nucleotide::ALL`: A, C, G, T.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PrimaryPeakEvidence {
    /// Selected sample position of the uniquely strongest primary channel.
    pub position_0based: usize,
    /// Raw analyzed A/C/G/T channel values at `position_0based`.
    pub channel_heights: [i32; 4],
}

/// One signal-derived base call.
#[derive(Debug, Clone)]
pub struct BaseCall {
    /// Call index in trace order (0-based).
    pub index_0based: usize,
    /// Trace sample of the call's locus (0-based).
    pub locus_position_0based: usize,
    /// First trace sample of the call's locus window (0-based).
    pub window_start_0based: usize,
    /// Trace sample after the call's locus window (0-based, exclusive).
    pub window_end_0based_exclusive: usize,
    /// Selected peak of every channel, in A/C/G/T order.
    pub peaks: [ChannelPeak; 4],
    /// Channel heights at the primary peak, when the call has one.
    pub primary_peak_evidence: Option<PrimaryPeakEvidence>,
    /// Primary base, or `N` when unresolved.
    pub primary: char,
    /// IUPAC symbol of the qualifying channels.
    pub ambiguity: char,
    /// Channels whose signal qualifies at the primary peak.
    pub qualifying_channels: Vec<Nucleotide>,
    /// Whether the vendor base agrees, when the trace carries one.
    pub vendor_agrees: Option<bool>,
}

/// Ordered calls and the primary sequence consumed by downstream stages.
#[derive(Debug, Clone)]
pub struct BaseCalls {
    /// Calls in trace order.
    pub calls: Vec<BaseCall>,
    /// Primary bases of every call, in trace order.
    pub primary_sequence: String,
}

impl BaseCalls {
    /// Number of call loci.
    #[must_use]
    pub fn len(&self) -> usize {
        self.calls.len()
    }

    /// Whether no calls were produced.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.calls.is_empty()
    }
}
