//! The placed-read step contract, `dna.placed_read/v1` (PROP-0011): one
//! read's placement as [`align_best`](crate::align_best) selected it, written
//! and read back by this crate alone, so a later step can start from a placed
//! read without searching again and `dna-align` stays the only producer of
//! alignments.
//!
//! The document names the read evidence document it places (its SHA-256),
//! the reference, and the configuration. Reading one back re-checks what the
//! search guarantees and what the read's evidence and the reference fix:
//! every column against the evidence's call and the reference's base, the
//! columns' order for the orientation, the segments against the columns, the
//! callable segments and masked bases recomputed from the evidence's masks,
//! and the metrics' own consistency. The score is recorded as selected; the
//! search that made it is identified by whoever keys the document (the
//! evaluation cache keys it by this crate's code).

use serde::{Deserialize, Serialize};

use dna_kernel::error::{AlignmentError, Result};
use dna_kernel::model::reference::Reference;
use dna_kernel::read_evidence::{Placement, ReadEvidence};

use crate::model::alignment::{
    Alignment, AlignmentColumn, AlignmentMetrics, Orientation, ReferenceSegment,
};
use crate::orient::callable_segments;

/// The contract's schema identifier.
pub const SCHEMA: &str = "dna.placed_read/v1";

#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct Document {
    schema_version: String,
    provenance: Provenance,
    /// The placement; `null` when the read's evidence places no read (a read
    /// the modality set aside).
    alignment: Option<Record>,
}

#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct Provenance {
    input: Input,
    reference_sha256: String,
    configuration_sha256: String,
}

/// The placed evidence document.
#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct Input {
    sha256: String,
}

#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct Record {
    orientation: Orientation,
    score: i64,
    wraps_origin: bool,
    reference_segments: Vec<[usize; 2]>,
    callable_segments: Vec<[usize; 2]>,
    metrics: Metrics,
    /// `[query, reference, call, position]` per column, in reference order.
    columns: Vec<(char, char, Option<usize>, Option<usize>)>,
}

#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct Metrics {
    exact_matches: usize,
    mismatches: usize,
    gap_opens: usize,
    callable_columns: usize,
    callable_identity: f64,
    unresolved_query_bases: usize,
    masked_query_bases: usize,
}

/// What a placed-read document is written for: the evidence document it
/// places, the reference, and the configuration, each by SHA-256.
#[derive(Debug, Clone, Copy)]
pub struct Identity<'a> {
    /// SHA-256 of the read evidence document.
    pub evidence_sha256: &'a str,
    /// SHA-256 of the reference sequence.
    pub reference_sha256: &'a str,
    /// SHA-256 of the configuration.
    pub configuration_sha256: &'a str,
}

fn segments(segments: &[ReferenceSegment]) -> Vec<[usize; 2]> {
    segments
        .iter()
        .map(|segment| [segment.start_0based, segment.end_0based_exclusive])
        .collect()
}

/// The `dna.placed_read/v1` document of `alignment`, or of no placement.
///
/// # Errors
///
/// Returns `AlignmentError::Inconsistent` when the document cannot be
/// serialized.
pub fn write(alignment: Option<&Alignment>, identity: Identity<'_>) -> Result<Vec<u8>> {
    let document = Document {
        schema_version: SCHEMA.to_owned(),
        provenance: Provenance {
            input: Input {
                sha256: identity.evidence_sha256.to_owned(),
            },
            reference_sha256: identity.reference_sha256.to_owned(),
            configuration_sha256: identity.configuration_sha256.to_owned(),
        },
        alignment: alignment.map(|alignment| Record {
            orientation: alignment.orientation,
            score: alignment.score,
            wraps_origin: alignment.wraps_origin,
            reference_segments: segments(&alignment.reference_segments),
            callable_segments: segments(&alignment.callable_segments),
            metrics: Metrics {
                exact_matches: alignment.metrics.exact_matches,
                mismatches: alignment.metrics.mismatches,
                gap_opens: alignment.metrics.gap_opens,
                callable_columns: alignment.metrics.callable_columns,
                callable_identity: alignment.metrics.callable_identity,
                unresolved_query_bases: alignment.metrics.unresolved_query_bases,
                masked_query_bases: alignment.metrics.masked_query_bases,
            },
            columns: alignment
                .columns
                .iter()
                .map(|column| {
                    (
                        column.query_base,
                        column.reference_base,
                        column.original_call_index_0based,
                        column.reference_index_0based,
                    )
                })
                .collect(),
        }),
    };
    serde_json::to_vec(&document)
        .map_err(|_| AlignmentError::Inconsistent("placed read does not serialize").into())
}

fn inconsistent<T>(reason: &'static str) -> Result<T> {
    Err(AlignmentError::Inconsistent(reason).into())
}

/// The placement a `dna.placed_read/v1` document records for `evidence` on
/// `reference`, `None` when it records none.
///
/// # Errors
///
/// Returns `AlignmentError::Inconsistent` when the document is malformed, was
/// written for another evidence document, reference, or configuration, or
/// does not hold for the evidence and the reference.
pub fn read(
    bytes: &[u8],
    evidence: &ReadEvidence,
    reference: &Reference,
    identity: Identity<'_>,
) -> Result<Option<Alignment>> {
    let Ok(document) = serde_json::from_slice::<Document>(bytes) else {
        return inconsistent("placed read is not a dna.placed_read/v1 document");
    };
    if document.schema_version != SCHEMA
        || document.provenance.input.sha256 != identity.evidence_sha256
        || document.provenance.reference_sha256 != identity.reference_sha256
        || document.provenance.configuration_sha256 != identity.configuration_sha256
    {
        return inconsistent("placed read was written for other inputs");
    }
    let Some(record) = document.alignment else {
        return Ok(None);
    };
    let calls = evidence.calls();
    let bases = reference.sequence.as_bytes();
    let columns = record
        .columns
        .iter()
        .map(|&(query, reference_base, call, position)| {
            AlignmentColumn::new(query, reference_base, call, position)
        })
        .collect::<Vec<_>>();
    // Each column against the read's call and the reference's base.
    let unresolved = |placement: Placement| matches!(placement, Placement::Unresolved(_));
    for column in &columns {
        let read_side = match column.original_call_index_0based {
            None => column.query_base == '-',
            Some(index) => calls.get(index).is_some_and(|call| {
                column.query_base
                    == if unresolved(call.placement) {
                        'N'
                    } else {
                        record.orientation.reference_base(call.base)
                    }
            }),
        };
        let reference_side = match column.reference_index_0based {
            None => column.reference_base == '-',
            Some(position) => {
                bases.get(position).map(|&base| char::from(base)) == Some(column.reference_base)
            }
        };
        if !read_side || !reference_side {
            return inconsistent("placed read column differs from its evidence or reference");
        }
    }
    // The columns' order: calls in the read's order for the orientation,
    // reference positions increasing, once across the origin when it wraps.
    let call_indices = columns
        .iter()
        .filter_map(AlignmentColumn::original_call_index_0based)
        .collect::<Vec<_>>();
    let calls_ordered = call_indices
        .windows(2)
        .all(|pair| match record.orientation {
            Orientation::Forward => pair[0] < pair[1],
            Orientation::Reverse => pair[0] > pair[1],
        });
    let positions = columns
        .iter()
        .filter_map(AlignmentColumn::reference_index_0based)
        .collect::<Vec<_>>();
    let wraps = positions
        .windows(2)
        .filter(|pair| pair[1] <= pair[0])
        .count();
    let positions_ordered = positions.windows(2).all(|pair| {
        pair[1] > pair[0] || (record.wraps_origin && pair[1] == 0 && pair[0] + 1 == bases.len())
    });
    if !calls_ordered || !positions_ordered || wraps != usize::from(record.wraps_origin) {
        return inconsistent("placed read columns are out of order");
    }
    // The segments against the columns, and the masks the evidence fixes.
    let to_segments = |pairs: &[[usize; 2]]| {
        pairs
            .iter()
            .map(|&[start, end]| ReferenceSegment {
                start_0based: start,
                end_0based_exclusive: end,
            })
            .collect::<Vec<_>>()
    };
    let reference_segments = to_segments(&record.reference_segments);
    let callable = to_segments(&record.callable_segments);
    let covered = |position: usize| {
        reference_segments.iter().any(|segment| {
            segment.start_0based <= position && position < segment.end_0based_exclusive
        })
    };
    let mask_of = |column: &AlignmentColumn| {
        column
            .original_call_index_0based
            .and_then(|index| calls.get(index))
            .map(|call| call.placement)
            .filter(|placement| *placement != Placement::Scored)
    };
    let recomputed = callable_segments(&columns, |column| {
        mask_of(column).is_none() && dna_kernel::model::nucleotide::is_canonical(column.query_base)
    });
    let masked = columns
        .iter()
        .filter(|column| mask_of(column).is_some())
        .count();
    let metrics = &record.metrics;
    if !positions.iter().all(|&position| covered(position))
        || recomputed != callable
        || masked != metrics.masked_query_bases
        || metrics.exact_matches + metrics.mismatches != metrics.callable_columns
    {
        return inconsistent("placed read segments or metrics do not hold for its evidence");
    }
    Ok(Some(Alignment {
        orientation: record.orientation,
        score: record.score,
        reference_segments,
        callable_segments: callable,
        wraps_origin: record.wraps_origin,
        metrics: AlignmentMetrics {
            exact_matches: metrics.exact_matches,
            mismatches: metrics.mismatches,
            gap_opens: metrics.gap_opens,
            callable_columns: metrics.callable_columns,
            callable_identity: metrics.callable_identity,
            unresolved_query_bases: metrics.unresolved_query_bases,
            masked_query_bases: metrics.masked_query_bases,
        },
        columns,
    }))
}

#[cfg(test)]
mod tests {
    use dna_kernel::model::reference::ReferenceTopology;
    use dna_kernel::read_evidence::EvidenceReason;

    use super::*;
    use crate::{AlignmentConfig, align_best};

    fn config() -> AlignmentConfig {
        AlignmentConfig {
            match_score: 3,
            mismatch_score: -5,
            ambiguous_score: 0,
            gap_open_score: -10,
            gap_extension_score: -4,
            minimum_callable_bases: 1,
            minimum_identity: 0.8,
            max_stretch: 64,
        }
    }

    fn reference(topology: ReferenceTopology) -> Reference {
        Reference {
            name: "ref".into(),
            sequence: "TTGACGTCAGTACGATCGTACCTGAGTACGAGGACTTGCAAGT".into(),
            topology,
            sequence_sha256: String::new(),
        }
    }

    const IDENTITY: Identity<'static> = Identity {
        evidence_sha256: "e",
        reference_sha256: "r",
        configuration_sha256: "c",
    };

    // PROP-0011
    #[test]
    fn reads_back_what_it_writes() -> Result<()> {
        let linear = reference(ReferenceTopology::Linear);
        // A forward read with a deletion and a masked call, a reverse read,
        // and a read across the origin of a circular reference.
        let mut masked = ReadEvidence::clean("TTGACGTCAGTACGATCGTCCTGAGTACGAGG");
        let call =
            masked.calls()[14].masked(Placement::Anchoring(EvidenceReason::new("dephased_signal")));
        masked = masked.with_call(14, call);
        let reverse = ReadEvidence::clean(&dna_kernel::model::nucleotide::reverse_complement(
            "GACGTCAGTACGATCGTACCTGAGTAC",
        ));
        let circular = reference(ReferenceTopology::Circular);
        let across = ReadEvidence::clean("GGACTTGCAAGTTTGACGTCAGTAC");
        for (evidence, reference) in [
            (&masked, &linear),
            (&reverse, &linear),
            (&across, &circular),
        ] {
            let alignment = align_best(evidence, reference, &config())?;
            let written = write(Some(&alignment), IDENTITY)?;
            let back = read(&written, evidence, reference, IDENTITY)?
                .ok_or(AlignmentError::Inconsistent("no placement read back"))?;
            assert_eq!(write(Some(&back), IDENTITY)?, written);
            assert_eq!(back.orientation(), alignment.orientation());
            assert_eq!(back.score(), alignment.score());
        }
        assert!(read(&write(None, IDENTITY)?, &masked, &linear, IDENTITY)?.is_none());
        Ok(())
    }

    // PROP-0011
    #[test]
    fn refuses_a_document_for_other_inputs_or_evidence() -> Result<()> {
        let reference = reference(ReferenceTopology::Linear);
        let evidence = ReadEvidence::clean("TTGACGTCAGTACGATCGTACCTGAG");
        let written = write(
            Some(&align_best(&evidence, &reference, &config())?),
            IDENTITY,
        )?;
        let other = Identity {
            evidence_sha256: "other",
            ..IDENTITY
        };
        assert!(read(&written, &evidence, &reference, other).is_err());
        // The same placement for a read with another base.
        let changed = ReadEvidence::clean("TTGACGTCAGTACGTTCGTACCTGAG");
        assert!(read(&written, &changed, &reference, IDENTITY).is_err());
        // A tampered column.
        let text = String::from_utf8(written.clone()).unwrap_or_default();
        let tampered = text.replacen("[\"T\",\"T\",0,0]", "[\"T\",\"T\",1,0]", 1);
        assert_ne!(tampered, text);
        assert!(read(tampered.as_bytes(), &evidence, &reference, IDENTITY).is_err());
        assert!(read(b"{}", &evidence, &reference, IDENTITY).is_err());
        Ok(())
    }
}
