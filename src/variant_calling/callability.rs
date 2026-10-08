//! Read callability: positions of a read whose calls cannot support a variant.
//!
//! Sanger reads lose reliability at the ends of the aligned interval and in the
//! calls that follow a long homopolymer in the sequencing direction, where
//! polymerase slippage shifts phase. These rules use only the read's own calls,
//! in trace (sequencing) order, so they need no primer or filename knowledge.

use crate::config::VariantCallingConfig;
use crate::model::basecalls::BaseCalls;
use crate::model::nucleotide::is_canonical;
use crate::model::quality::QualityControlResult;
use crate::model::variant::{VariantCallMapping, VariantExclusionReason};

/// Untrusted call indices of one read.
pub(super) struct ReadCallability {
    /// First call index outside the leading end margin.
    trusted_start: usize,
    /// End (exclusive) of the calls outside the trailing end margin.
    trusted_end: usize,
    /// Per call index: inside a post-homopolymer window.
    post_homopolymer: Vec<bool>,
}

impl ReadCallability {
    /// Derives the untrusted calls of one read from its retained interval and
    /// primary calls.
    pub(super) fn new(
        calls: &BaseCalls,
        quality: &QualityControlResult,
        config: &VariantCallingConfig,
    ) -> Self {
        Self {
            trusted_start: quality
                .trim_start_0based
                .saturating_add(config.read_end_margin),
            trusted_end: quality
                .trim_end_0based_exclusive
                .saturating_sub(config.read_end_margin),
            post_homopolymer: post_homopolymer_windows(calls, config),
        }
    }

    /// Callability reasons against a variant's mapped calls, each reported once.
    pub(super) fn reasons(&self, mappings: &[VariantCallMapping]) -> Vec<VariantExclusionReason> {
        let indices = || mappings.iter().map(|mapping| mapping.call_index_0based);
        let mut reasons = Vec::new();
        if indices().any(|index| index < self.trusted_start || index >= self.trusted_end) {
            reasons.push(VariantExclusionReason::ReadEnd);
        }
        if indices().any(|index| self.post_homopolymer.get(index).copied().unwrap_or(false)) {
            reasons.push(VariantExclusionReason::PostHomopolymer);
        }
        reasons
    }
}

/// Marks, for every run of at least `homopolymer_min_length` identical
/// canonical primary calls, `post_homopolymer_window` calls starting with the
/// run's last call.
fn post_homopolymer_windows(calls: &BaseCalls, config: &VariantCallingConfig) -> Vec<bool> {
    let primary: Vec<char> = calls.calls.iter().map(|call| call.primary).collect();
    let mut flagged = vec![false; primary.len()];
    let mut start = 0;
    while start < primary.len() {
        let base = primary[start];
        let length = primary[start..]
            .iter()
            .take_while(|&&next| next == base)
            .count();
        let last = start + length - 1;
        if is_canonical(base) && length >= config.homopolymer_min_length {
            let end = last
                .saturating_add(config.post_homopolymer_window)
                .min(primary.len());
            flagged[last..end].fill(true);
        }
        start += length;
    }
    flagged
}

#[cfg(test)]
mod tests {
    use crate::model::basecalls::{BaseCall, ChannelPeak, PeakSource};
    use crate::model::nucleotide::Nucleotide;
    use crate::model::quality::{CallQuality, QualityControlResult};
    use crate::model::variant::{VariantCallMapping, VariantCallRole};

    use super::*;

    fn read(primary: &str, trim: (usize, usize)) -> (BaseCalls, QualityControlResult) {
        let calls = primary
            .chars()
            .enumerate()
            .map(|(index, base)| BaseCall {
                index_0based: index,
                locus_position_0based: index * 4,
                window_start_0based: index * 4,
                window_end_0based_exclusive: index * 4 + 4,
                peaks: Nucleotide::ALL.map(|channel| ChannelPeak {
                    base: channel,
                    height: 0,
                    position_0based: index * 4,
                    source: PeakSource::LocalMaximum,
                }),
                primary_peak_evidence: None,
                primary: base,
                ambiguity: base,
                qualifying_channels: Vec::new(),
                vendor_agrees: None,
            })
            .collect();
        (
            BaseCalls {
                calls,
                primary_sequence: primary.into(),
            },
            QualityControlResult {
                per_call: (0..primary.len())
                    .map(|index| CallQuality {
                        index_0based: index,
                        penalty: 0,
                        relative_quality_score: 40,
                        vendor_quality_applies: false,
                    })
                    .collect(),
                trim_start_0based: trim.0,
                trim_end_0based_exclusive: trim.1,
                retained_sequence: String::new(),
            },
        )
    }

    fn settings(read_end_margin: usize, minimum_run: usize, window: usize) -> VariantCallingConfig {
        VariantCallingConfig {
            max_indel_length: 50,
            minimum_peak_height: 150,
            relative_quality_threshold: 30,
            regions: vec![[1, 50_000]],
            read_end_margin,
            homopolymer_min_length: minimum_run,
            post_homopolymer_window: window,
        }
    }

    fn reasons_at(callability: &ReadCallability, index: usize) -> Vec<VariantExclusionReason> {
        callability.reasons(&[VariantCallMapping {
            role: VariantCallRole::Supporting,
            call_index_0based: index,
            reference_position_0based: Some(index),
        }])
    }

    #[test]
    fn flags_calls_within_the_margin_of_the_retained_interval_ends() {
        let (calls, quality) = read(&"ACGT".repeat(8), (5, 25));
        let callability = ReadCallability::new(&calls, &quality, &settings(3, 8, 0));

        for (index, read_end) in [(7, true), (8, false), (21, false), (22, true)] {
            assert_eq!(
                reasons_at(&callability, index) == [VariantExclusionReason::ReadEnd],
                read_end,
                "call {index}"
            );
        }
    }

    #[test]
    fn flags_the_end_of_a_long_homopolymer_and_the_window_after_it() {
        let primary = format!("ACGT{}TACGTACGTACGT", "C".repeat(8));
        let (calls, quality) = read(&primary, (0, primary.len()));
        let callability = ReadCallability::new(&calls, &quality, &settings(0, 8, 3));

        for (index, after_run) in [(10, false), (11, true), (12, true), (13, true), (14, false)] {
            assert_eq!(
                reasons_at(&callability, index) == [VariantExclusionReason::PostHomopolymer],
                after_run,
                "call {index}"
            );
        }
    }

    #[test]
    fn ignores_short_runs_and_runs_broken_by_unresolved_calls() {
        let primary = format!("A{}A{}A", "C".repeat(7), "CCCCNCCCC");
        let (calls, quality) = read(&primary, (0, primary.len()));
        let callability = ReadCallability::new(&calls, &quality, &settings(0, 8, 6));

        assert!((0..primary.len()).all(|index| reasons_at(&callability, index).is_empty()));
    }

    #[test]
    fn zero_settings_trust_every_retained_call() {
        let primary = format!("{}A", "C".repeat(12));
        let (calls, quality) = read(&primary, (0, primary.len()));
        let callability = ReadCallability::new(&calls, &quality, &settings(0, 8, 0));

        assert!((0..primary.len()).all(|index| reasons_at(&callability, index).is_empty()));
    }

    #[test]
    fn reports_each_reason_once_per_variant() {
        let primary = format!("{}AAAA", "C".repeat(8));
        let (calls, quality) = read(&primary, (0, primary.len()));
        let callability = ReadCallability::new(&calls, &quality, &settings(4, 8, 3));
        let mappings: Vec<_> = (7..12)
            .map(|index| VariantCallMapping {
                role: VariantCallRole::Flanking,
                call_index_0based: index,
                reference_position_0based: Some(index),
            })
            .collect();

        assert_eq!(
            callability.reasons(&mappings),
            [
                VariantExclusionReason::ReadEnd,
                VariantExclusionReason::PostHomopolymer
            ]
        );
    }
}
