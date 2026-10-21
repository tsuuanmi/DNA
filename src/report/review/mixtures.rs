//! Each admitted read's phase mixtures (ADR-0082) as the review reports
//! them: every component the modality asserts, by offset and fraction, an
//! explained mixture's placed components compared with the reference, and
//! the reference positions of the calls the components do not explain.

use std::collections::{BTreeMap, BTreeSet};

use super::round;
use crate::error::{Error, Result};
use crate::model::review_result::{ComparedResult, ComponentResult, MixtureResult};
use crate::sample_analysis::SampleAnalysis;
use dna_call::components::ComponentCall;

/// Each admitted read's phase mixtures (ADR-0082), by content SHA-256: every
/// component placed and called as a read of its own, its eligible edits
/// written in the profile notation over the reference interval it covers
/// (`notation`, which the pipeline's representation gives),
/// and the reference positions of the calls the components do not explain.
pub(crate) fn mixtures(
    analysis: &SampleAnalysis,
    notation: impl Fn(&ComponentCall) -> Result<Vec<String>>,
) -> Result<BTreeMap<String, Vec<MixtureResult>>> {
    let mut found = BTreeMap::new();
    for read in analysis.components() {
        let positions = analysis
            .called
            .iter()
            .find(|called| called.input_sha256 == read.sha256)
            .map(|called| {
                called
                    .alignment
                    .columns()
                    .iter()
                    .filter_map(|column| {
                        Some((
                            column.original_call_index_0based()?,
                            column.reference_index_0based()? + 1,
                        ))
                    })
                    .collect::<BTreeMap<_, _>>()
            })
            .unwrap_or_default();
        let mixtures = read
            .mixtures
            .iter()
            .map(|mixture| {
                // Every component the modality asserts, by offset and
                // fraction; an explained mixture's placed ones also compared
                // with the reference.
                let components = mixture
                    .mixture
                    .components
                    .iter()
                    .map(|component| {
                        let compared = mixture
                            .components
                            .iter()
                            .find(|placed| placed.offset == component.offset)
                            .filter(|_| mixture.mixture.explained)
                            .map(|placed| {
                                Ok::<_, Error>(ComparedResult {
                                    first: placed.span[0],
                                    last: placed.span[1],
                                    sequence: placed.sequence.clone(),
                                    calls: notation(placed)?,
                                })
                            })
                            .transpose()?;
                        Ok(ComponentResult {
                            offset: component.offset,
                            fraction: round(component.fraction),
                            compared,
                        })
                    })
                    .collect::<Result<Vec<_>>>()?;
                Ok(MixtureResult {
                    explained: mixture.mixture.explained,
                    components,
                    unexplained: mixture
                        .mixture
                        .unexplained
                        .iter()
                        .filter_map(|call| positions.get(call).copied())
                        .collect::<BTreeSet<_>>()
                        .into_iter()
                        .collect(),
                })
            })
            .collect::<Result<Vec<_>>>()?;
        found.insert(read.sha256, mixtures);
    }
    Ok(found)
}
