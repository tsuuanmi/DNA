//! Each call of the selected alignment with its reference column and the
//! core's trust in it, so one analysis document shows the chain sanger →
//! align → call for its read (PROP-0011).

use dna_align::model::alignment::AlignmentColumn;
use dna_call::model::trust::{CallTrust, TrustScope, Withholder};

use crate::error::{ReportError, Result};
use crate::model::result::{
    AlignedCallResult, ColumnResult, ScopeResult, TrustResult, WithholderResult,
};
use crate::report::evidence::placement;

/// Projects every call the alignment `columns` place, in call order, with
/// its 1-based reference position and insertion ordinal and its `trust`.
///
/// An inserted call takes the last reference position before it in reference
/// order and its ordinal among the calls inserted there; one inserted before
/// any reference position has no column.
pub(super) fn project(
    columns: &[AlignmentColumn],
    trust: &[CallTrust],
) -> Result<Vec<AlignedCallResult>> {
    let mut previous: Option<usize> = None;
    let mut ordinal = 0;
    let mut calls = Vec::new();
    for column in columns {
        if let Some(reference) = column.reference_index_0based() {
            previous = Some(reference);
            ordinal = 0;
        } else {
            ordinal += 1;
        }
        let Some(index) = column.original_call_index_0based() else {
            continue;
        };
        let call = trust.get(index).ok_or(ReportError::Inconsistent(
            "analysis alignment places a call without trust",
        ))?;
        calls.push(AlignedCallResult {
            index,
            column: previous.map(|reference| ColumnResult {
                position: reference + 1,
                ordinal,
            }),
            trust: project_trust(call),
        });
    }
    calls.sort_by_key(|call| call.index);
    if calls.windows(2).any(|pair| pair[0].index == pair[1].index) {
        return Err(ReportError::Inconsistent("analysis alignment places a call twice").into());
    }
    Ok(calls)
}

fn project_trust(trust: &CallTrust) -> TrustResult {
    TrustResult {
        placement: placement(trust.placement),
        trusted: trust.trusted,
        exempt: trust.exempt,
        withheld_by: trust
            .withheld_by
            .iter()
            .map(|withholder| match withholder {
                Withholder::ReadEnd => WithholderResult::ReadEnd,
                Withholder::Untrimmed => WithholderResult::Untrimmed,
                Withholder::Masked(reason) => WithholderResult::Masked {
                    reason: reason.label(),
                },
            })
            .collect(),
        scope: match trust.scope {
            TrustScope::Full => ScopeResult::Full,
            TrustScope::Base => ScopeResult::Base,
        },
    }
}

#[cfg(test)]
mod tests {
    use dna_kernel::read_evidence::{EvidenceReason, Placement};

    use super::*;

    const DEPHASED: EvidenceReason = EvidenceReason::new("dephased_signal");

    fn trust(placement: Placement, withheld_by: Vec<Withholder>) -> CallTrust {
        CallTrust {
            placement,
            unresolved: false,
            vetoes: Vec::new(),
            trusted: withheld_by.is_empty(),
            exempt: false,
            withheld_by,
            scope: TrustScope::Full,
            second: None,
        }
    }

    fn column(call: Option<usize>, reference: Option<usize>) -> AlignmentColumn {
        AlignmentColumn::new('A', 'A', call, reference)
    }

    #[test]
    fn places_inserted_calls_after_the_last_reference_position() -> Result<()> {
        // A reverse read in reference order: a call inserted before any
        // reference position, a matched call, a deletion, two inserted calls,
        // and a matched call.
        let columns = [
            column(Some(5), None),
            column(Some(4), Some(9)),
            column(None, Some(10)),
            column(Some(3), None),
            column(Some(2), None),
            column(Some(1), Some(11)),
        ];
        let trust = vec![
            trust(Placement::Scored, Vec::new()),
            trust(Placement::Scored, vec![Withholder::ReadEnd]),
            trust(Placement::Scored, Vec::new()),
            trust(
                Placement::Unresolved(DEPHASED),
                vec![Withholder::Masked(DEPHASED)],
            ),
            trust(
                Placement::Scored,
                vec![Withholder::Masked(DEPHASED), Withholder::Untrimmed],
            ),
            trust(Placement::Scored, vec![Withholder::Untrimmed]),
        ];
        let calls = serde_json::to_value(project(&columns, &trust)?).unwrap_or_default();
        let scored = serde_json::json!({"kind": "scored"});
        assert_eq!(
            calls,
            serde_json::json!([
                {"index": 1, "column": {"position": 12, "ordinal": 0}, "trust": {
                    "placement": scored, "trusted": false, "exempt": false,
                    "withheld_by": [{"kind": "read_end"}], "scope": "full"}},
                {"index": 2, "column": {"position": 11, "ordinal": 2}, "trust": {
                    "placement": scored, "trusted": true, "exempt": false,
                    "withheld_by": [], "scope": "full"}},
                {"index": 3, "column": {"position": 11, "ordinal": 1}, "trust": {
                    "placement": {"kind": "unresolved", "reason": "dephased_signal"},
                    "trusted": false, "exempt": false,
                    "withheld_by": [{"kind": "masked", "reason": "dephased_signal"}],
                    "scope": "full"}},
                {"index": 4, "column": {"position": 10, "ordinal": 0}, "trust": {
                    "placement": scored, "trusted": false, "exempt": false,
                    "withheld_by": [
                        {"kind": "masked", "reason": "dephased_signal"},
                        {"kind": "untrimmed"}
                    ],
                    "scope": "full"}},
                {"index": 5, "trust": {
                    "placement": scored, "trusted": false, "exempt": false,
                    "withheld_by": [{"kind": "untrimmed"}], "scope": "full"}},
            ])
        );
        Ok(())
    }

    #[test]
    fn rejects_a_call_without_trust() {
        let columns = [column(Some(1), Some(0))];
        assert!(project(&columns, &[trust(Placement::Scored, Vec::new())]).is_err());
    }
}
