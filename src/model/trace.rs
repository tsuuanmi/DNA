//! Validated four-channel chromatogram and vendor evidence.

/// Optional calls and quality values stored by the ABI basecaller.
#[derive(Debug, Clone, Default)]
pub struct VendorEvidence {
    pub(crate) primary: Option<String>,
    pub(crate) qualities: Option<Vec<u8>>,
}

/// Decoded analyzed chromatogram samples in canonical A/C/G/T order.
#[derive(Debug, Clone)]
pub struct Chromatogram {
    pub(crate) source_name: String,
    pub(crate) source_sha256: String,
    pub(crate) channels: [Vec<i32>; 4],
    pub(crate) base_locations: Vec<usize>,
    pub(crate) vendor: VendorEvidence,
}

impl Chromatogram {
    /// Number of samples in every channel.
    pub fn sample_count(&self) -> usize {
        self.channels[0].len()
    }

    /// Number of vendor-defined base loci.
    pub fn call_count(&self) -> usize {
        self.base_locations.len()
    }
}
