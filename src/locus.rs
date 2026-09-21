//! Shared PLOC-defined locus-window geometry.

use crate::model::trace::Chromatogram;

/// Half-open sample window around one validated PLOC locus.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct LocusWindow {
    pub(crate) start: usize,
    pub(crate) end: usize,
}

/// Builds symmetric neighboring-midpoint windows around every PLOC locus.
///
/// This geometry is shared by basecalling and signal-evidence extraction. The
/// caller maps geometry failures into its own stage-specific error type.
pub(crate) fn windows(trace: &Chromatogram) -> std::result::Result<Vec<LocusWindow>, String> {
    let positions = &trace.base_locations;
    if positions.len() < 2 {
        return Err("at least two PLOC positions are required".into());
    }
    let sample_count = trace.sample_count();
    let mut output = Vec::with_capacity(positions.len());
    for index in 0..positions.len() {
        let start = if index == 0 {
            let spacing = positions[1] - positions[0];
            positions[0].saturating_sub(spacing / 2)
        } else {
            midpoint(positions[index - 1], positions[index])?
        };
        let end = if index + 1 == positions.len() {
            let spacing = positions[index] - positions[index - 1];
            positions[index]
                .checked_add(spacing.div_ceil(2))
                .unwrap_or(sample_count)
                .min(sample_count)
        } else {
            midpoint(positions[index], positions[index + 1])?
        };
        if start >= end || end > sample_count || positions[index] < start || positions[index] >= end
        {
            return Err(format!(
                "invalid locus window {start}..{end} for PLOC {}",
                positions[index]
            ));
        }
        output.push(LocusWindow { start, end });
    }
    Ok(output)
}

fn midpoint(left: usize, right: usize) -> std::result::Result<usize, String> {
    left.checked_add((right - left) / 2)
        .ok_or_else(|| "locus-window midpoint overflow".into())
}

#[cfg(test)]
mod tests {
    use crate::model::trace::{Chromatogram, VendorEvidence};

    use super::*;

    #[test]
    fn builds_neighbor_midpoint_windows() -> std::result::Result<(), String> {
        let trace = Chromatogram {
            source_name: "synthetic.ab1".into(),
            source_sha256: String::new(),
            channels: std::array::from_fn(|_| vec![0; 12]),
            base_locations: vec![2, 6, 10],
            vendor: VendorEvidence::default(),
        };

        assert_eq!(
            windows(&trace)?,
            vec![
                LocusWindow { start: 0, end: 4 },
                LocusWindow { start: 4, end: 8 },
                LocusWindow { start: 8, end: 12 },
            ]
        );
        Ok(())
    }
}
