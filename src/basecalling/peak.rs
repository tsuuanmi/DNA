//! Channel-local peak selection inside shared PLOC locus windows.

use crate::error::{Error, Result};
use crate::locus::{self, LocusWindow};
use crate::model::basecalls::{ChannelPeak, PeakSource};
use crate::model::nucleotide::Nucleotide;
use crate::model::trace::Chromatogram;

/// Builds the shared PLOC-defined windows and maps geometry failures to basecalling.
pub(crate) fn windows(trace: &Chromatogram) -> Result<Vec<LocusWindow>> {
    locus::windows(trace).map_err(Error::Basecalling)
}

/// Finds one positive local peak per channel or samples PLOC explicitly.
pub(crate) fn peaks(trace: &Chromatogram, window: LocusWindow, ploc: usize) -> [ChannelPeak; 4] {
    std::array::from_fn(|channel_index| {
        let channel = &trace.channels[channel_index];
        let search_start = window.start.max(1);
        let search_end = window.end.min(channel.len().saturating_sub(1));
        let mut selected = None;
        for position in search_start..search_end {
            let value = channel[position];
            let local = (channel[position - 1] <= value && value > channel[position + 1])
                || (channel[position - 1] < value && value >= channel[position + 1]);
            if local && value > 0 && selected.is_none_or(|(_, best)| value > best) {
                selected = Some((position, value));
            }
        }
        let (position, height, source) = selected.map_or_else(
            || (ploc, channel[ploc], PeakSource::PlocFallback),
            |(position, height)| (position, height, PeakSource::LocalMaximum),
        );
        ChannelPeak {
            base: Nucleotide::ALL[channel_index],
            height,
            position_0based: position,
            source,
        }
    })
}

#[cfg(test)]
mod tests {
    use crate::model::trace::{Chromatogram, VendorEvidence};

    use super::*;

    #[test]
    fn selects_each_channel_independently_and_falls_back_to_ploc() {
        let trace = Chromatogram {
            source_name: "synthetic.ab1".into(),
            source_sha256: String::new(),
            channels: [
                vec![0, 1, 10, 1, 0, 0, 0],
                vec![0, 1, 2, 20, 1, 0, 0],
                vec![0, 0, 1, 2, 30, 1, 0],
                vec![0, 1, 2, 3, 4, 5, 6],
            ],
            base_locations: vec![3, 5],
            vendor: VendorEvidence::default(),
        };
        let selected = peaks(&trace, LocusWindow { start: 1, end: 6 }, 3);
        assert_eq!((selected[0].position_0based, selected[0].height), (2, 10));
        assert_eq!((selected[1].position_0based, selected[1].height), (3, 20));
        assert_eq!((selected[2].position_0based, selected[2].height), (4, 30));
        assert_eq!((selected[3].position_0based, selected[3].height), (3, 3));
        assert_eq!(selected[3].source, PeakSource::PlocFallback);
    }
}
