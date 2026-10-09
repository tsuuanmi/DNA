//! Compact sample evidence aggregated across independently processed reads.

use serde::Serialize;

use crate::model::alignment::{Orientation, ReferenceSegment};
use crate::model::variant::{VariantCallRole, VariantExclusionReason, VariantKind};
use dna_kernel::read_evidence::EvidenceProfile;

/// Why a mapped read pair is not admitted as reliable overlap evidence.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum OverlapExclusionReason {
    /// Fewer comparable canonical bases than `minimum_comparable_bases`.
    ComparableBasesBelowMinimum,
    /// Canonical-base agreement below `minimum_overlap_agreement`.
    AgreementBelowMinimum,
}

/// Pairwise overlap evidence discovered after independent reference placement.
#[derive(Debug, Clone)]
pub struct ReadOverlapEvidence {
    /// Registry index of the first read of the pair.
    pub left_read_index: usize,
    /// Registry index of the second read of the pair.
    pub right_read_index: usize,
    /// Reference positions both reads cover.
    pub shared_positions: usize,
    /// Shared positions where both reads call a canonical base.
    pub comparable_bases: usize,
    /// Comparable positions where both reads agree.
    pub agreements: usize,
    /// Comparable positions where the reads disagree.
    pub conflicts: usize,
    /// Fraction of comparable positions that agree, when any are comparable.
    pub agreement: Option<f64>,
    /// Whether the pair is admitted as overlap evidence.
    pub eligible: bool,
    /// Why the pair is not admitted; empty when eligible.
    pub exclusion_reasons: Vec<OverlapExclusionReason>,
}

/// One maximal reference interval with constant independently placed read depth.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SampleCoverageEvidence {
    /// First reference position of the segment (0-based, inclusive).
    pub start_0based: usize,
    /// Position after the segment (0-based, exclusive).
    pub end_0based_exclusive: usize,
    /// Reads covering every position of the segment.
    pub read_depth: usize,
    /// Forward-oriented reads covering the segment.
    pub forward_depth: usize,
    /// Reverse-oriented reads covering the segment.
    pub reverse_depth: usize,
}

/// How one read observes a locus retained because at least one read differs.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum LocusState {
    /// The read's canonical base matches the reference.
    Reference,
    /// The read's canonical base differs from the reference.
    Alternate,
    /// The read or reference base is not canonical.
    Unresolved,
    /// The read has a gap at the position.
    Deletion,
    /// Covered by a masked call: kept for review, never callable.
    Masked,
}

/// Structural nucleotide-contribution eligibility for one retained locus observation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NucleotideContribution {
    /// The observation carries an evidence profile and contributes nucleotide mass.
    Eligible,
    /// The observed call has no evidence profile.
    MissingProfile,
    /// A deletion has no called base to contribute.
    DeletionEvent,
    /// A masked call never contributes nucleotide mass.
    MaskedCall,
}

/// Factorized topology of reads observing one retained differential locus.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct LocusSupportTopology {
    /// Reads observing the locus.
    pub reads: usize,
    /// Forward-oriented observing reads.
    pub forward_reads: usize,
    /// Reverse-oriented observing reads.
    pub reverse_reads: usize,
    /// Reads observing the reference base.
    pub reference_reads: usize,
    /// Reads observing another canonical base.
    pub alternate_reads: usize,
    /// Reads observing an unresolved base.
    pub unresolved_reads: usize,
    /// Reads with a gap at the locus.
    pub deletion_reads: usize,
    /// Reads whose call at the locus is masked.
    pub masked_reads: usize,
    /// Reads whose call carries an evidence profile.
    pub profile_reads: usize,
    /// Forward-oriented reads with an evidence profile.
    pub profile_forward_reads: usize,
    /// Reverse-oriented reads with an evidence profile.
    pub profile_reverse_reads: usize,
}

/// Concise evidence supporting one selected read placement.
#[derive(Debug, Clone)]
pub struct SampleReadAlignmentEvidence {
    /// Selected read orientation.
    pub orientation: Orientation,
    /// Alignment columns pairing canonical read and reference bases.
    pub callable_bases: usize,
    /// Fraction of callable columns whose bases agree.
    pub identity: f64,
    /// Number of gap openings.
    pub gap_opens: usize,
    /// Unresolved read bases outside masked calls.
    pub unresolved_bases: usize,
    /// Alignment columns on masked calls.
    pub masked_bases: usize,
    /// Reference segments the read covers.
    pub reference_segments: Vec<ReferenceSegment>,
    /// Reference segments observed by unmasked calls.
    pub callable_reference_segments: Vec<ReferenceSegment>,
    /// Whether the read crosses the origin of a circular reference.
    pub wraps_origin: bool,
}

/// One read retained once at sample scope.
#[derive(Debug, Clone)]
pub struct SampleReadEvidence {
    /// Reviewer-facing name of the read's source.
    pub input_name: String,
    /// SHA-256 content identity of the read's source.
    pub input_sha256: String,
    /// Selected-alignment summary of the read.
    pub alignment: SampleReadAlignmentEvidence,
}

/// A read its modality rejected before the core, recorded at sample scope
/// without contributing alignment, coverage, overlaps, loci, or variants. The
/// modality's reason is joined by the report.
#[derive(Debug, Clone)]
pub struct RejectedSampleRead {
    /// Reviewer-facing name of the rejected read's source.
    pub input_name: String,
    /// SHA-256 content identity of the rejected read's source.
    pub input_sha256: String,
}

/// One read observation at one retained sample reference locus.
#[derive(Debug, Clone)]
pub struct SampleLocusObservation {
    /// Registry index of the observing read.
    pub read_index: usize,
    /// What the read observed at the locus.
    pub state: LocusState,
    /// Reference-strand base observed; absent for a deletion.
    pub base: Option<char>,
    /// Source call of the read; absent for a deletion.
    pub call_index_0based: Option<usize>,
    /// Reference-oriented evidence profile of the call.
    pub profile: Option<EvidenceProfile>,
    /// Whether the observation contributes nucleotide-profile mass.
    pub nucleotide_contribution: NucleotideContribution,
}

/// Threshold-free heterogeneity geometry for one non-empty nucleotide-profile partition.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ProfileHeterogeneity {
    /// Mean impurity of the individual profiles.
    pub within_profile_impurity: f64,
    /// Dispersion of the profiles around their mean.
    pub between_profile_dispersion: f64,
    /// Impurity of the mean profile: within plus between.
    pub total: f64,
}

/// Unweighted eligible nucleotide-profile support retained at one differential locus.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LocusNucleotideSupport {
    pub(crate) contributors: usize,
    pub(crate) forward_contributors: usize,
    pub(crate) reverse_contributors: usize,
    /// Summed A/C/G/T profile mass of all contributors.
    pub support: [f64; 4],
    pub(crate) forward_support: [f64; 4],
    pub(crate) reverse_support: [f64; 4],
    /// Mean profile of all contributors.
    pub mean_profile: Option<EvidenceProfile>,
    /// Mean profile of forward-oriented contributors.
    pub forward_mean_profile: Option<EvidenceProfile>,
    /// Mean profile of reverse-oriented contributors.
    pub reverse_mean_profile: Option<EvidenceProfile>,
    /// Profile heterogeneity of all contributors.
    pub heterogeneity: Option<ProfileHeterogeneity>,
    /// Profile heterogeneity of forward-oriented contributors.
    pub forward_heterogeneity: Option<ProfileHeterogeneity>,
    /// Profile heterogeneity of reverse-oriented contributors.
    pub reverse_heterogeneity: Option<ProfileHeterogeneity>,
    /// Total-variation distance between the forward and reverse mean profiles.
    pub directional_profile_distance: Option<f64>,
}

/// All covering-read observations retained at one sample reference locus.
#[derive(Debug, Clone)]
pub struct SampleLocusEvidence {
    /// 1-based reference position of the locus.
    pub position_1based: usize,
    /// Reference base at the locus.
    pub reference_base: char,
    /// Counts of observing reads by orientation and state.
    pub support_topology: LocusSupportTopology,
    /// Nucleotide-profile support of the locus.
    pub nucleotide_support: LocusNucleotideSupport,
    /// Every read's observation of the locus, in registry order.
    pub observations: Vec<SampleLocusObservation>,
}

/// One trace call directly supporting or flanking a normalized variant.
#[derive(Debug, Clone, PartialEq)]
pub struct VariantCallEvidence {
    /// Whether the call supports or flanks the variant.
    pub role: VariantCallRole,
    /// Reference-strand base.
    pub base: char,
    /// Source call of the read.
    pub call_index_0based: usize,
}

/// One read observing a normalized variant, with configured eligibility retained.
#[derive(Debug, Clone, PartialEq)]
pub struct VariantSupport {
    /// Registry index of the supporting read.
    pub read_index: usize,
    /// Whether the read's observation is eligible for reporting.
    pub eligible: bool,
    /// Why the observation is not eligible; empty when eligible.
    pub exclusion_reasons: Vec<VariantExclusionReason>,
    /// The read's public calls of the variant.
    pub calls: Vec<VariantCallEvidence>,
}

/// Factorized topology of reads observing one normalized variant.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct VariantSupportTopology {
    /// Reads observing the variant.
    pub reads: usize,
    /// Reads whose observation is eligible.
    pub eligible_reads: usize,
    /// Forward-oriented observing reads.
    pub forward_reads: usize,
    /// Reverse-oriented observing reads.
    pub reverse_reads: usize,
    /// Forward-oriented reads with an eligible observation.
    pub eligible_forward_reads: usize,
    /// Reverse-oriented reads with an eligible observation.
    pub eligible_reverse_reads: usize,
}

/// Admitted reads that callably observe a variant's reference span without
/// supporting the variant, in read-registry order.
///
/// Opposition is evidence, not a vote or a verdict: a forward and a reverse
/// read can share an assay artifact, and a read can oppose a true minor allele.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VariantOpposition {
    /// Registry indices of the opposing reads, in registry order.
    pub read_indices: Vec<usize>,
    /// Forward-oriented opposing reads.
    pub forward_reads: usize,
    /// Reverse-oriented opposing reads.
    pub reverse_reads: usize,
}

/// One normalized observed variant with factorized read support.
#[derive(Debug, Clone)]
pub struct VariantEvidence {
    /// 1-based reference position of the anchored variant.
    pub position_1based: usize,
    /// Anchored reference allele.
    pub reference: String,
    /// Anchored alternate allele.
    pub alternate: String,
    /// Variant type.
    pub kind: VariantKind,
    /// Counts of supporting reads by orientation and eligibility.
    pub support_topology: VariantSupportTopology,
    /// Every supporting read's observation, in registry order.
    pub support: Vec<VariantSupport>,
    /// Reads that callably observe the variant's span without supporting it.
    pub opposition: VariantOpposition,
}

/// Complete compact evidence for one sample.
#[derive(Debug, Clone)]
pub struct SampleEvidence {
    /// SHA-256 of the reference every read was placed on.
    pub reference_sha256: String,
    /// SHA-256 of the configuration every read was called with.
    pub configuration_sha256: String,
    /// Admitted reads in registry (content SHA-256) order.
    pub reads: Vec<SampleReadEvidence>,
    /// Reads the modality rejected, in content SHA-256 order.
    pub rejected_reads: Vec<RejectedSampleRead>,
    /// Run-length coverage segments of the admitted reads.
    pub coverage: Vec<SampleCoverageEvidence>,
    /// Pairwise overlap evidence of the admitted reads.
    pub overlaps: Vec<ReadOverlapEvidence>,
    /// Loci where at least one read differs from the reference.
    pub locus_differences: Vec<SampleLocusEvidence>,
    /// Normalized variants with their support and opposition.
    pub variants: Vec<VariantEvidence>,
}
