//! Canonical Sanger sequencing evidence after source-format decoding.

/// Optional calls and quality values stored by the ABI basecaller.
#[derive(Debug, Clone, Default)]
pub struct VendorEvidence {
    /// Vendor base calls (`PBAS.2`), when present.
    pub primary: Option<String>,
    /// Vendor qualities (`PCON.2`), when present.
    pub qualities: Option<Vec<u8>>,
}

/// Canonical analyzed Sanger chromatogram samples in A/C/G/T order.
#[derive(Debug, Clone)]
pub struct Chromatogram {
    /// Reviewer-facing name of the trace file.
    pub source_name: String,
    /// SHA-256 of the trace file bytes.
    pub source_sha256: String,
    pub(crate) channels: [Vec<i32>; 4],
    pub(crate) locus_positions: Vec<usize>,
    /// Optional vendor evidence, never authoritative.
    pub vendor: VendorEvidence,
}

impl Chromatogram {
    /// Number of samples in every channel.
    #[must_use]
    pub fn sample_count(&self) -> usize {
        self.channels[0].len()
    }

    /// Number of vendor-defined base loci.
    #[must_use]
    pub fn call_count(&self) -> usize {
        self.locus_positions.len()
    }
}
