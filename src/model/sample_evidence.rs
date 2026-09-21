//! Compact sample evidence aggregated across independently processed reads.

use serde::Serialize;

use crate::model::alignment::{Orientation, ReferenceSegment};
use crate::model::locus_evidence::EvidenceProfile;
use crate::model::signal::TraceIntegrity;
use crate::model::variant::{VariantCallRole, VariantExclusionReason, VariantKind};

/// Why a mapped read pair is not admitted as reliable overlap evidence.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub(crate) enum OverlapExclusionReason {
    ComparableBasesBelowMinimum,
    AgreementBelowMinimum,
}

/// Pairwise overlap evidence discovered after independent reference placement.
#[derive(Debug, Clone)]
pub(crate) struct ReadOverlapEvidence {
    pub(crate) left_read_index: usize,
    pub(crate) right_read_index: usize,
    pub(crate) shared_positions: usize,
    pub(crate) comparable_bases: usize,
    pub(crate) agreements: usize,
    pub(crate) conflicts: usize,
    pub(crate) agreement: Option<f64>,
    pub(crate) eligible: bool,
    pub(crate) exclusion_reasons: Vec<OverlapExclusionReason>,
}

/// One maximal reference interval with constant independently placed read depth.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct SampleCoverageEvidence {
    pub(crate) start_0based: usize,
    pub(crate) end_0based_exclusive: usize,
    pub(crate) read_depth: usize,
    pub(crate) forward_depth: usize,
    pub(crate) reverse_depth: usize,
}

/// How one read observes a locus retained because at least one read differs.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub(crate) enum LocusState {
    Reference,
    Alternate,
    Unresolved,
    Deletion,
}

/// Structural nucleotide-contribution eligibility for one retained locus observation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum NucleotideContribution {
    Eligible,
    MissingProfile,
    DeletionEvent,
}

/// Factorized topology of reads observing one retained differential locus.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct LocusSupportTopology {
    pub(crate) reads: usize,
    pub(crate) forward_reads: usize,
    pub(crate) reverse_reads: usize,
    pub(crate) reference_reads: usize,
    pub(crate) alternate_reads: usize,
    pub(crate) unresolved_reads: usize,
    pub(crate) deletion_reads: usize,
    pub(crate) profile_reads: usize,
    pub(crate) profile_forward_reads: usize,
    pub(crate) profile_reverse_reads: usize,
}

/// Concise evidence supporting one selected read placement.
#[derive(Debug, Clone)]
pub(crate) struct SampleReadAlignmentEvidence {
    pub(crate) orientation: Orientation,
    pub(crate) callable_bases: usize,
    pub(crate) identity: f64,
    pub(crate) gap_opens: usize,
    pub(crate) unresolved_bases: usize,
    pub(crate) reference_segments: Vec<ReferenceSegment>,
    pub(crate) wraps_origin: bool,
}

/// One read retained once at sample scope.
#[derive(Debug, Clone)]
pub(crate) struct SampleReadEvidence {
    pub(crate) input_name: String,
    pub(crate) input_sha256: String,
    pub(crate) integrity: TraceIntegrity,
    pub(crate) alignment: SampleReadAlignmentEvidence,
}

/// Reference-oriented signal evidence associated with one source call.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct CallDNAEvidence {
    pub(crate) corrected_amplitudes: [f64; 4],
    pub(crate) snrs: [f64; 4],
    pub(crate) profile: Option<EvidenceProfile>,
    pub(crate) in_noisy_region: bool,
}

/// One read observation at one retained sample reference locus.
#[derive(Debug, Clone)]
pub(crate) struct SampleLocusObservation {
    pub(crate) read_index: usize,
    pub(crate) state: LocusState,
    pub(crate) base: Option<char>,
    pub(crate) quality: Option<u8>,
    pub(crate) signal: Option<CallDNAEvidence>,
    pub(crate) nucleotide_contribution: NucleotideContribution,
}

/// Threshold-free heterogeneity geometry for one non-empty nucleotide-profile partition.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct ProfileHeterogeneity {
    pub(crate) within_profile_impurity: f64,
    pub(crate) between_profile_dispersion: f64,
    pub(crate) total_profile_heterogeneity: f64,
}

/// Unweighted eligible nucleotide-profile support retained at one differential locus.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct LocusNucleotideSupport {
    pub(crate) contributors: usize,
    pub(crate) forward_contributors: usize,
    pub(crate) reverse_contributors: usize,
    pub(crate) support: [f64; 4],
    pub(crate) forward_support: [f64; 4],
    pub(crate) reverse_support: [f64; 4],
    pub(crate) mean_profile: Option<EvidenceProfile>,
    pub(crate) forward_mean_profile: Option<EvidenceProfile>,
    pub(crate) reverse_mean_profile: Option<EvidenceProfile>,
    pub(crate) heterogeneity: Option<ProfileHeterogeneity>,
    pub(crate) forward_heterogeneity: Option<ProfileHeterogeneity>,
    pub(crate) reverse_heterogeneity: Option<ProfileHeterogeneity>,
    pub(crate) directional_profile_distance: Option<f64>,
}

/// All covering-read observations retained at one sample reference locus.
#[derive(Debug, Clone)]
pub(crate) struct SampleLocusEvidence {
    pub(crate) position_1based: usize,
    pub(crate) reference_base: char,
    pub(crate) support_topology: LocusSupportTopology,
    pub(crate) nucleotide_support: LocusNucleotideSupport,
    pub(crate) observations: Vec<SampleLocusObservation>,
}

/// One trace call directly supporting or flanking a normalized variant.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct VariantCallEvidence {
    pub(crate) role: VariantCallRole,
    pub(crate) base: char,
    pub(crate) peak_heights: [i32; 4],
    pub(crate) quality: u8,
    pub(crate) signal: CallDNAEvidence,
}

/// One read observing a normalized variant, with configured eligibility retained.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct VariantSupport {
    pub(crate) read_index: usize,
    pub(crate) eligible: bool,
    pub(crate) exclusion_reasons: Vec<VariantExclusionReason>,
    pub(crate) calls: Vec<VariantCallEvidence>,
}

/// Factorized topology of reads observing one normalized variant.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct VariantSupportTopology {
    pub(crate) reads: usize,
    pub(crate) eligible_reads: usize,
    pub(crate) forward_reads: usize,
    pub(crate) reverse_reads: usize,
    pub(crate) eligible_forward_reads: usize,
    pub(crate) eligible_reverse_reads: usize,
}

/// One normalized observed variant with factorized read support.
#[derive(Debug, Clone)]
pub(crate) struct VariantEvidence {
    pub(crate) position_1based: usize,
    pub(crate) reference: String,
    pub(crate) alternate: String,
    pub(crate) kind: VariantKind,
    pub(crate) support_topology: VariantSupportTopology,
    pub(crate) support: Vec<VariantSupport>,
}

/// Complete compact evidence for one sample.
#[derive(Debug, Clone)]
pub(crate) struct SampleEvidence {
    pub(crate) reference_sha256: String,
    pub(crate) configuration_sha256: String,
    pub(crate) reads: Vec<SampleReadEvidence>,
    pub(crate) coverage: Vec<SampleCoverageEvidence>,
    pub(crate) overlaps: Vec<ReadOverlapEvidence>,
    pub(crate) locus_differences: Vec<SampleLocusEvidence>,
    pub(crate) variants: Vec<VariantEvidence>,
}
