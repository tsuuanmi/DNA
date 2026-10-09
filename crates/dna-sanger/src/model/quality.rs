//! Relative quality scores and auditable end-trim bounds.

/// Per-call quality-control evidence.
#[derive(Debug, Clone)]
pub struct CallQuality {
    /// Call index (0-based).
    pub index_0based: usize,
    /// Signal penalty of the call.
    pub penalty: i32,
    /// Uncalibrated relative quality score.
    pub relative_quality_score: u8,
    /// Whether a vendor quality was available for the call.
    pub vendor_quality_applies: bool,
}

/// Complete quality-control result.
#[derive(Debug, Clone)]
pub struct QualityControlResult {
    /// Quality of every call, in call order.
    pub per_call: Vec<CallQuality>,
    /// First retained call (0-based, inclusive).
    pub trim_start_0based: usize,
    /// Call after the retained interval (0-based, exclusive).
    pub trim_end_0based_exclusive: usize,
    /// Primary bases inside the trim interval.
    pub retained_sequence: String,
}
