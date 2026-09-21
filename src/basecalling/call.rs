//! DNA-derived call orchestration at validated PLOC loci.

use crate::basecalling::iupac;
use crate::basecalling::peak;
use crate::config::BasecallingConfig;
use crate::error::{Error, Result};
use crate::model::basecalls::{BaseCall, BaseCalls, PrimaryPeakEvidence};
use crate::model::trace::Chromatogram;

/// Re-calls every vendor-defined locus from analyzed channel signals.
pub(crate) fn call(trace: &Chromatogram, config: &BasecallingConfig) -> Result<BaseCalls> {
    let windows = peak::windows(trace)?;
    let mut calls = Vec::with_capacity(trace.call_count());
    let mut primary_sequence = String::with_capacity(trace.call_count());

    for (index, (&ploc, window)) in trace.base_locations.iter().zip(windows).enumerate() {
        let peaks = peak::peaks(trace, window, ploc);
        if peaks
            .iter()
            .any(|peak| peak.position_0based < window.start || peak.position_0based >= window.end)
        {
            return Err(Error::Basecalling(format!(
                "selected peak escaped call window {}..{} at call {index}",
                window.start, window.end
            )));
        }
        let mut order = [0_usize, 1, 2, 3];
        order.sort_by(|left, right| {
            peaks[*right]
                .height
                .cmp(&peaks[*left].height)
                .then_with(|| left.cmp(right))
        });
        let top_height = peaks[order[0]].height;
        let tied_top = top_height > 0 && peaks[order[1]].height == top_height;
        let vendor_primary = trace
            .vendor
            .primary
            .as_ref()
            .and_then(|sequence| sequence.as_bytes().get(index))
            .map(|value| char::from(*value));

        let (primary, ambiguity, qualifying_channels, primary_peak_evidence) = if top_height <= 0
            || tied_top
        {
            ('N', 'N', Vec::new(), None)
        } else {
            let top_index = order[0];
            let primary_peak_position = peaks[top_index].position_0based;
            let primary_peak_evidence = PrimaryPeakEvidence {
                position_0based: primary_peak_position,
                channel_heights: std::array::from_fn(|channel| {
                    trace.channels[channel][primary_peak_position]
                }),
            };
            let qualifying_channels = order
                .iter()
                .filter(|channel| {
                    let selected_height = peaks[**channel].height;
                    let colocated_height = primary_peak_evidence.channel_heights[**channel];
                    reaches_ratio(selected_height, top_height, config.secondary_peak_ratio)
                        && reaches_ratio(colocated_height, top_height, config.secondary_peak_ratio)
                })
                .map(|channel| peaks[*channel].base)
                .collect::<Vec<_>>();
            let strongest = peaks[top_index].base.as_char();
            let (primary, ambiguity) = match qualifying_channels.len() {
                1 => (strongest, strongest),
                2 => (strongest, iupac::code(&qualifying_channels)),
                3 => (strongest, 'N'),
                _ => ('N', 'N'),
            };
            (
                primary,
                ambiguity,
                qualifying_channels,
                Some(primary_peak_evidence),
            )
        };

        primary_sequence.push(primary);
        calls.push(BaseCall {
            index_0based: index,
            ploc_0based: ploc,
            window_start_0based: window.start,
            window_end_0based_exclusive: window.end,
            peaks,
            primary_peak_evidence,
            primary,
            ambiguity,
            qualifying_channels,
            vendor_agrees: vendor_primary.map(|vendor| vendor == primary),
        });
    }

    Ok(BaseCalls {
        calls,
        primary_sequence,
    })
}

fn reaches_ratio(height: i32, top_height: i32, minimum_ratio: f64) -> bool {
    height > 0 && f64::from(height) / f64::from(top_height) >= minimum_ratio
}

#[cfg(test)]
mod tests {
    use crate::model::trace::{Chromatogram, VendorEvidence};

    use super::*;

    fn trace(channels: [Vec<i32>; 4]) -> Chromatogram {
        trace_at(channels, vec![2, 6])
    }

    fn trace_at(channels: [Vec<i32>; 4], base_locations: Vec<usize>) -> Chromatogram {
        Chromatogram {
            source_name: "synthetic.ab1".into(),
            source_sha256: String::new(),
            channels,
            base_locations,
            vendor: VendorEvidence::default(),
        }
    }

    #[test]
    fn calls_unambiguous_strongest_channel() -> Result<()> {
        let chromatogram = trace_at(
            [
                vec![0, 0, 0, 1, 100, 1, 0, 1, 100, 1],
                vec![0, 0, 0, 1, 20, 1, 0, 1, 20, 1],
                vec![0, 0, 0, 1, 3, 1, 0, 1, 3, 1],
                vec![0, 0, 0, 0, -7, 0, 0, 0, -7, 0],
            ],
            vec![4, 8],
        );
        let calls = call(
            &chromatogram,
            &BasecallingConfig {
                secondary_peak_ratio: 0.33,
            },
        )?;
        assert_eq!(calls.primary_sequence, "AA");
        assert_eq!(
            calls
                .calls
                .iter()
                .map(|call| call.ambiguity)
                .collect::<String>(),
            "AA"
        );
        assert_eq!(calls.calls[0].peaks[0].position_0based, 4);
        assert_eq!(
            calls.calls[0].primary_peak_evidence,
            Some(PrimaryPeakEvidence {
                position_0based: 4,
                channel_heights: [100, 20, 3, -7],
            })
        );
        Ok(())
    }

    #[test]
    fn exact_strongest_tie_is_unresolved() -> Result<()> {
        let chromatogram = trace([
            vec![0, 1, 20, 1, 0, 1, 20, 1],
            vec![0, 1, 20, 1, 0, 1, 20, 1],
            vec![0; 8],
            vec![0; 8],
        ]);
        let calls = call(
            &chromatogram,
            &BasecallingConfig {
                secondary_peak_ratio: 0.33,
            },
        )?;
        assert_eq!(calls.primary_sequence, "NN");
        assert!(
            calls
                .calls
                .iter()
                .all(|call| call.primary_peak_evidence.is_none())
        );
        Ok(())
    }

    #[test]
    fn non_positive_channels_have_no_primary_evidence() -> Result<()> {
        let chromatogram = trace([vec![0; 8], vec![0; 8], vec![0; 8], vec![0; 8]]);
        let calls = call(
            &chromatogram,
            &BasecallingConfig {
                secondary_peak_ratio: 0.33,
            },
        )?;
        assert_eq!(calls.primary_sequence, "NN");
        assert!(
            calls
                .calls
                .iter()
                .all(|call| call.primary_peak_evidence.is_none())
        );
        Ok(())
    }

    #[test]
    fn off_locus_secondary_peak_does_not_create_ambiguity() -> Result<()> {
        let chromatogram = trace([
            vec![0, 1, 100, 1, 0, 1, 100, 1],
            vec![0, 0, 1, 40, 0, 40, 1, 0],
            vec![0; 8],
            vec![0; 8],
        ]);
        let calls = call(
            &chromatogram,
            &BasecallingConfig {
                secondary_peak_ratio: 0.33,
            },
        )?;
        assert_eq!(calls.primary_sequence, "AA");
        assert_eq!(
            calls
                .calls
                .iter()
                .map(|call| call.ambiguity)
                .collect::<String>(),
            "AA"
        );
        assert!(
            calls
                .calls
                .iter()
                .all(|call| call.qualifying_channels.len() == 1
                    && call.qualifying_channels[0].as_char() == 'A')
        );
        assert_eq!(
            (
                calls.calls[0].peaks[1].height,
                calls.calls[0].peaks[1].position_0based
            ),
            (40, 3)
        );
        assert_eq!(
            calls.calls[0]
                .primary_peak_evidence
                .as_ref()
                .map(|evidence| evidence.channel_heights[1]),
            Some(1)
        );
        Ok(())
    }

    #[test]
    fn overlapping_offset_secondary_peak_still_qualifies() -> Result<()> {
        let chromatogram = trace([
            vec![0, 1, 100, 1, 0, 1, 100, 1],
            vec![0, 20, 40, 50, 0, 50, 40, 20],
            vec![0; 8],
            vec![0; 8],
        ]);
        let calls = call(
            &chromatogram,
            &BasecallingConfig {
                secondary_peak_ratio: 0.33,
            },
        )?;
        assert_eq!(calls.primary_sequence, "AA");
        assert_eq!(
            calls
                .calls
                .iter()
                .map(|call| call.ambiguity)
                .collect::<String>(),
            "MM"
        );
        assert_eq!(
            (
                calls.calls[0].peaks[1].height,
                calls.calls[0].peaks[1].position_0based
            ),
            (50, 3)
        );
        assert_eq!(
            calls.calls[0]
                .primary_peak_evidence
                .as_ref()
                .map(|evidence| evidence.channel_heights[1]),
            Some(40)
        );
        Ok(())
    }

    #[test]
    fn offset_secondary_signal_at_exact_ratio_threshold_qualifies() -> Result<()> {
        let chromatogram = trace([
            vec![0, 1, 100, 1, 0, 1, 100, 1],
            vec![0, 20, 50, 60, 0, 60, 50, 20],
            vec![0; 8],
            vec![0; 8],
        ]);
        let calls = call(
            &chromatogram,
            &BasecallingConfig {
                secondary_peak_ratio: 0.5,
            },
        )?;
        assert!(
            calls
                .calls
                .iter()
                .all(|call| call.ambiguity == 'M' && call.qualifying_channels.len() == 2)
        );
        Ok(())
    }

    #[test]
    fn ploc_fallback_primary_always_qualifies() -> Result<()> {
        let chromatogram = trace([(0..8).collect(), vec![0; 8], vec![0; 8], vec![0; 8]]);
        let calls = call(
            &chromatogram,
            &BasecallingConfig {
                secondary_peak_ratio: 1.0,
            },
        )?;
        assert_eq!(calls.primary_sequence, "AA");
        assert!(calls.calls.iter().all(|call| {
            call.peaks[0].source == crate::model::basecalls::PeakSource::PlocFallback
                && call.qualifying_channels.len() == 1
                && call.qualifying_channels[0].as_char() == 'A'
                && call.primary_peak_evidence.as_ref().is_some_and(|evidence| {
                    evidence.position_0based == call.ploc_0based
                        && evidence.channel_heights
                            == std::array::from_fn(|channel| {
                                chromatogram.channels[channel][call.ploc_0based]
                            })
                })
        }));
        Ok(())
    }

    #[test]
    fn ploc_fallback_secondary_must_be_colocated() -> Result<()> {
        let chromatogram = trace([
            vec![0, 100, 1, 0, 0, 100, 1, 0],
            vec![0, 1, 40, 41, 42, 43, 44, 45],
            vec![0; 8],
            vec![0; 8],
        ]);
        let calls = call(
            &chromatogram,
            &BasecallingConfig {
                secondary_peak_ratio: 0.33,
            },
        )?;
        assert_eq!(calls.calls[0].primary, 'A');
        assert_eq!(calls.calls[0].ambiguity, 'A');
        assert_eq!(
            calls.calls[0].peaks[1].source,
            crate::model::basecalls::PeakSource::PlocFallback
        );
        assert_eq!(calls.calls[0].qualifying_channels.len(), 1);
        Ok(())
    }

    #[test]
    fn three_qualifying_channels_keep_primary_but_not_ambiguity() -> Result<()> {
        let chromatogram = trace([
            vec![0, 1, 100, 1, 0, 1, 100, 1],
            vec![0, 1, 50, 1, 0, 1, 50, 1],
            vec![0, 1, 40, 1, 0, 1, 40, 1],
            vec![0; 8],
        ]);
        let calls = call(
            &chromatogram,
            &BasecallingConfig {
                secondary_peak_ratio: 0.33,
            },
        )?;
        assert_eq!(calls.primary_sequence, "AA");
        assert_eq!(
            calls
                .calls
                .iter()
                .map(|call| call.ambiguity)
                .collect::<String>(),
            "NN"
        );
        Ok(())
    }

    #[test]
    fn four_qualifying_channels_are_fully_unresolved() -> Result<()> {
        let chromatogram = trace([
            vec![0, 1, 100, 1, 0, 1, 100, 1],
            vec![0, 1, 80, 1, 0, 1, 80, 1],
            vec![0, 1, 60, 1, 0, 1, 60, 1],
            vec![0, 1, 40, 1, 0, 1, 40, 1],
        ]);
        let calls = call(
            &chromatogram,
            &BasecallingConfig {
                secondary_peak_ratio: 0.33,
            },
        )?;
        assert_eq!(calls.primary_sequence, "NN");
        assert_eq!(
            calls
                .calls
                .iter()
                .map(|call| call.ambiguity)
                .collect::<String>(),
            "NN"
        );
        assert!(
            calls
                .calls
                .iter()
                .all(|call| call.qualifying_channels.len() == 4)
        );
        Ok(())
    }
}
