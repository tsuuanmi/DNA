//! Canonical Sanger sequencing evidence after source-format decoding.

/// Optional calls and quality values stored by the ABI basecaller.
#[derive(Debug, Clone, Default)]
pub(crate) struct VendorEvidence {
    pub(crate) primary: Option<String>,
    pub(crate) qualities: Option<Vec<u8>>,
}

/// Canonical analyzed Sanger chromatogram samples in A/C/G/T order.
#[derive(Debug, Clone)]
pub(crate) struct Chromatogram {
    pub(crate) source_name: String,
    pub(crate) source_sha256: String,
    pub(crate) channels: [Vec<i32>; 4],
    pub(crate) locus_positions: Vec<usize>,
    pub(crate) vendor: VendorEvidence,
}

impl Chromatogram {
    /// Number of samples in every channel.
    pub(crate) fn sample_count(&self) -> usize {
        self.channels[0].len()
    }

    /// Number of vendor-defined base loci.
    pub(crate) fn call_count(&self) -> usize {
        self.locus_positions.len()
    }
}
