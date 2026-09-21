//! Auditable low-quality end trimming.

use crate::config::QualityControlConfig;
use crate::error::{Error, Result};
use crate::model::basecalls::BaseCalls;
use crate::model::quality::{CallQuality, QualityControlResult};
use crate::model::trace::Chromatogram;
use crate::quality_control::penalty;
use crate::quality_control::quality;

/// Scores calls and selects one retained interval.
pub(crate) fn analyze(
    trace: &Chromatogram,
    calls: &BaseCalls,
    config: &QualityControlConfig,
) -> Result<QualityControlResult> {
    if calls.len() < config.minimum_retained_bases {
        return Err(Error::QualityControl(format!(
            "{} calls are fewer than minimum_retained_bases {}",
            calls.len(),
            config.minimum_retained_bases
        )));
    }
    let penalty = penalty::calculate(calls, config.trim_window_size, config.best_section_fraction)?;
    let scores = quality::relative_scores(&penalty.penalties, config.max_relative_quality_score);
    let threshold = config.trim_stringency * penalty.best_average * config.trim_window_size as f64;
    let mut trim_start = 0;
    for start in (0..penalty.best_start).rev() {
        let end = start
            .saturating_add(config.trim_window_size)
            .min(calls.len());
        let sum: i64 = penalty.penalties[start..end]
            .iter()
            .map(|value| i64::from(*value))
            .sum();
        if sum as f64 > threshold {
            trim_start = end.min(penalty.best_start);
            break;
        }
    }
    let mut trim_end = calls.len();
    for start in penalty.best_end..calls.len() {
        let end = start
            .saturating_add(config.trim_window_size)
            .min(calls.len());
        let sum: i64 = penalty.penalties[start..end]
            .iter()
            .map(|value| i64::from(*value))
            .sum();
        if sum as f64 > threshold {
            trim_end = start.max(penalty.best_end);
            break;
        }
    }
    if trim_end <= trim_start || trim_end - trim_start < config.minimum_retained_bases {
        return Err(Error::QualityControl(format!(
            "retained interval {trim_start}..{trim_end} is shorter than minimum {}",
            config.minimum_retained_bases
        )));
    }
    let retained_sequence = calls.primary_sequence[trim_start..trim_end].to_owned();
    let per_call = calls
        .calls
        .iter()
        .enumerate()
        .map(|(index, call)| {
            let vendor_quality = trace
                .vendor
                .qualities
                .as_ref()
                .and_then(|qualities| qualities.get(index))
                .copied();
            CallQuality {
                index_0based: index,
                penalty: penalty.penalties[index],
                relative_quality_score: scores[index],
                vendor_quality_applies: vendor_quality.is_some()
                    && call.vendor_agrees == Some(true),
            }
        })
        .collect();
    Ok(QualityControlResult {
        per_call,
        trim_start_0based: trim_start,
        trim_end_0based_exclusive: trim_end,
        retained_sequence,
    })
}

#[cfg(test)]
mod tests {
    use crate::model::basecalls::{BaseCall, ChannelPeak, PeakSource};
    use crate::model::nucleotide::Nucleotide;
    use crate::model::trace::VendorEvidence;

    use super::*;

    #[test]
    fn retains_a_uniform_read_and_applies_matching_vendor_quality() -> Result<()> {
        let locations = [2, 6, 10, 14];
        let calls = BaseCalls {
            calls: locations
                .iter()
                .enumerate()
                .map(|(index, &ploc)| BaseCall {
                    index_0based: index,
                    ploc_0based: ploc,
                    window_start_0based: ploc.saturating_sub(1),
                    window_end_0based_exclusive: ploc + 2,
                    peaks: std::array::from_fn(|channel| ChannelPeak {
                        base: Nucleotide::ALL[channel],
                        height: 1,
                        position_0based: ploc,
                        source: PeakSource::PlocFallback,
                    }),
                    primary_peak_evidence: None,
                    primary: 'A',
                    ambiguity: 'A',
                    qualifying_channels: vec![Nucleotide::A],
                    vendor_agrees: Some(true),
                })
                .collect(),
            primary_sequence: "AAAA".into(),
        };
        let trace = Chromatogram {
            source_name: "synthetic.ab1".into(),
            source_sha256: String::new(),
            channels: std::array::from_fn(|_| vec![0; 16]),
            base_locations: locations.to_vec(),
            vendor: VendorEvidence {
                primary: Some("AAAA".into()),
                qualities: Some(vec![40; 4]),
            },
        };
        let result = analyze(
            &trace,
            &calls,
            &QualityControlConfig {
                trim_window_size: 2,
                best_section_fraction: 0.5,
                max_relative_quality_score: 60,
                trim_stringency: 7.0,
                minimum_retained_bases: 4,
            },
        )?;
        assert_eq!(
            (result.trim_start_0based, result.trim_end_0based_exclusive),
            (0, 4)
        );
        assert_eq!(result.retained_sequence, "AAAA");
        assert!(result.per_call.iter().all(|quality| {
            quality.relative_quality_score == 60 && quality.vendor_quality_applies
        }));
        Ok(())
    }
}
