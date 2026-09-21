//! Relative quality scores and auditable end-trim bounds.

/// Per-call quality-control evidence.
#[derive(Debug, Clone)]
pub struct CallQuality {
    pub(crate) index_0based: usize,
    pub(crate) penalty: i32,
    pub(crate) relative_quality_score: u8,
    pub(crate) vendor_quality_applies: bool,
}

/// Complete quality-control result.
#[derive(Debug, Clone)]
pub struct QualityControlResult {
    pub(crate) per_call: Vec<CallQuality>,
    pub(crate) trim_start_0based: usize,
    pub(crate) trim_end_0based_exclusive: usize,
    pub(crate) retained_sequence: String,
}
