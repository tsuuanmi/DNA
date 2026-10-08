//! Shared canonical Sanger locus-window geometry.

use crate::error::LocusWindowError;
use crate::model::sanger::Chromatogram;

/// Half-open sample window around one validated Sanger locus.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct LocusWindow {
    pub(crate) start: usize,
    pub(crate) end: usize,
}

/// Builds symmetric neighboring-midpoint windows around every Sanger locus.
///
/// This geometry is shared by basecalling and signal-evidence extraction. The
/// caller wraps geometry failures in its own stage-specific error.
pub(crate) fn windows(trace: &Chromatogram) -> Result<Vec<LocusWindow>, LocusWindowError> {
    let positions = &trace.locus_positions;
    if positions.len() < 2 {
        return Err(LocusWindowError::TooFewLoci);
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
            return Err(LocusWindowError::InvalidWindow {
                start,
                end,
                position: positions[index],
            });
        }
        output.push(LocusWindow { start, end });
    }
    Ok(output)
}

fn midpoint(left: usize, right: usize) -> Result<usize, LocusWindowError> {
    left.checked_add((right - left) / 2)
        .ok_or(LocusWindowError::MidpointOverflow)
}

#[cfg(test)]
mod tests {
    use crate::model::sanger::{Chromatogram, VendorEvidence};

    use super::*;

    #[test]
    fn rejects_a_trace_with_fewer_than_two_loci() {
        let trace = Chromatogram {
            source_name: "synthetic.ab1".into(),
            source_sha256: String::new(),
            channels: std::array::from_fn(|_| vec![0; 12]),
            locus_positions: vec![2],
            vendor: VendorEvidence::default(),
        };

        assert_eq!(windows(&trace), Err(LocusWindowError::TooFewLoci));
    }

    #[test]
    fn builds_neighbor_midpoint_windows() -> Result<(), LocusWindowError> {
        let trace = Chromatogram {
            source_name: "synthetic.ab1".into(),
            source_sha256: String::new(),
            channels: std::array::from_fn(|_| vec![0; 12]),
            locus_positions: vec![2, 6, 10],
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
