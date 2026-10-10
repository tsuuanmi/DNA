//! What one placed read observes at each reference position and junction.

use std::collections::BTreeMap;

use crate::model::called_read::CalledRead;
use crate::variant_calling::VariantCallingConfig;
use crate::variant_calling::eligibility::ReadEligibility;
use dna_kernel::read_evidence::{MaskedAlignment, VetoSet};

/// What one read contributes at one reference position or junction.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) enum SiteValue {
    /// A base at a reference position.
    Base(char),
    /// A deleted reference position.
    Deletion,
    /// The bases inserted at a junction; empty when the read spans the
    /// junction without an insertion.
    Insertion(String),
}

/// One observed base: a call aligned to a position or inserted at a junction.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) struct Char {
    pub(super) base: char,
    /// Unmasked, trusted (`read_end`), canonical, and no support veto.
    pub(super) clean: bool,
    /// Canonical and informative, unmasked or anchoring: it shows where a
    /// neighbouring run ends (ADR-0071).
    pub(super) resolving: bool,
}

/// One read's value at one site, and whether it may decide the site.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct Observation {
    pub(super) value: SiteValue,
    /// The calls behind the value are unmasked, trusted (`read_end`), and
    /// raise no support veto.
    pub(super) clean: bool,
    /// The observed bases: one for a base, none for a deletion, the inserted
    /// ones at a junction.
    pub(super) chars: Vec<Char>,
}

/// Every observation of one read, keyed by 0-based reference position; a
/// junction is keyed by the position before it.
#[derive(Debug, Clone)]
pub(super) struct ReadObservations {
    pub(super) name: String,
    pub(super) positions: BTreeMap<usize, Observation>,
    pub(super) junctions: BTreeMap<usize, Observation>,
}

/// Walks the read's alignment columns in reference order.
///
/// A base is clean when its call is clean and canonical. A deletion is clean
/// when the nearest calls on both sides are clean. A junction between two
/// adjacent covered positions is clean when every inserted call and the
/// nearest calls on both sides are clean; it carries the inserted bases, or
/// none.
pub(super) fn observe(
    read: &CalledRead,
    config: &VariantCallingConfig,
    reference_length: usize,
) -> ReadObservations {
    let calls = read.evidence.calls();
    let eligibility = ReadEligibility::new(&read.evidence, config);
    let clean_call = |index: usize| {
        calls.get(index).is_some_and(|call| {
            call.mask.is_none() && call.vetoes == VetoSet::default() && eligibility.trusted(index)
        })
    };
    let interval = read.evidence.informative();
    let observed_char = |base: char, call: usize| Char {
        base,
        clean: canonical(base) && clean_call(call),
        resolving: canonical(base)
            && interval.contains(&call)
            && calls.get(call).is_some_and(|evidence| {
                evidence
                    .mask
                    .is_none_or(|mask| mask.alignment == MaskedAlignment::Anchoring)
            }),
    };
    let columns = &read.alignment.columns;
    // Nearest call at or before, and at or after, every column.
    let mut before = Vec::with_capacity(columns.len());
    let mut last = None;
    for column in columns {
        last = column.original_call_index_0based.or(last);
        before.push(last);
    }
    let mut after = vec![None; columns.len()];
    let mut next = None;
    for (offset, column) in columns.iter().enumerate().rev() {
        next = column.original_call_index_0based.or(next);
        after[offset] = next;
    }
    let bounded = |left: Option<usize>, right: Option<usize>| {
        left.is_some_and(clean_call) && right.is_some_and(clean_call)
    };

    let mut positions = BTreeMap::new();
    let mut junctions = BTreeMap::new();
    let mut previous: Option<(usize, usize)> = None; // (reference position, column)
    let mut inserted: Vec<(char, usize)> = Vec::new();
    for (offset, column) in columns.iter().enumerate() {
        let Some(position) = column.reference_index_0based else {
            if let Some(call) = column.original_call_index_0based {
                inserted.push((column.query_base, call));
            }
            continue;
        };
        if let Some((previous_position, previous_column)) = previous
            && (previous_position + 1) % reference_length == position
        {
            let bases = inserted.iter().map(|(base, _)| *base).collect::<String>();
            let clean = inserted
                .iter()
                .all(|(base, call)| canonical(*base) && clean_call(*call))
                && bounded(before[previous_column], after[offset]);
            junctions.insert(
                previous_position,
                Observation {
                    value: SiteValue::Insertion(bases),
                    clean,
                    chars: inserted
                        .iter()
                        .map(|(base, call)| observed_char(*base, *call))
                        .collect(),
                },
            );
        }
        let observation = if column.query_base == '-' {
            Observation {
                value: SiteValue::Deletion,
                clean: bounded(
                    offset.checked_sub(1).and_then(|left| before[left]),
                    after.get(offset + 1).copied().flatten(),
                ),
                chars: Vec::new(),
            }
        } else {
            Observation {
                value: SiteValue::Base(column.query_base),
                clean: canonical(column.query_base)
                    && column.original_call_index_0based.is_some_and(clean_call),
                chars: column
                    .original_call_index_0based
                    .map(|call| observed_char(column.query_base, call))
                    .into_iter()
                    .collect(),
            }
        };
        positions.insert(position, observation);
        previous = Some((position, offset));
        inserted.clear();
    }
    // Indels beyond the caller's length cap are never clean evidence.
    let cap = config.max_indel_length;
    for observation in junctions.values_mut() {
        if observation.chars.len() > cap {
            observation.clean = false;
            for char in &mut observation.chars {
                char.clean = false;
            }
        }
    }
    let deleted = positions
        .iter()
        .filter(|(_, observation)| observation.value == SiteValue::Deletion)
        .map(|(&position, _)| position)
        .collect::<Vec<_>>();
    for run in deleted.chunk_by(|left, right| left + 1 == *right) {
        if run.len() > cap {
            for position in run {
                if let Some(observation) = positions.get_mut(position) {
                    observation.clean = false;
                }
            }
        }
    }
    ReadObservations {
        name: read.input_name.clone(),
        positions,
        junctions,
    }
}

impl ReadObservations {
    /// The read's observed bases over `first..=last`, insertions at both
    /// edges and inside included, framed by the bases at the two neighbouring
    /// positions (`None` in the frame for a deleted neighbour); `None` when
    /// the read does not cover the stretch and both its neighbours.
    pub(super) fn span(
        &self,
        first: usize,
        last: usize,
        length: usize,
    ) -> Option<(Option<Char>, Vec<Char>, Option<Char>)> {
        let before = (first + length - 1) % length;
        let after = (last + 1) % length;
        let frame = |position: usize| -> Option<Option<Char>> {
            Some(self.positions.get(&position)?.chars.first().copied())
        };
        let (left, right) = (frame(before)?, frame(after)?);
        let mut chars = self.junctions.get(&before)?.chars.clone();
        let mut position = first;
        loop {
            let observation = self.positions.get(&position)?;
            if observation.value == SiteValue::Deletion && !observation.clean {
                // An untrusted deletion hides the base like an unresolved call.
                chars.push(Char {
                    base: 'N',
                    clean: false,
                    resolving: false,
                });
            }
            chars.extend(observation.chars.iter().copied());
            chars.extend(self.junctions.get(&position)?.chars.iter().copied());
            if position == last {
                break;
            }
            position = (position + 1) % length;
        }
        Some((left, chars, right))
    }
}

const fn canonical(base: char) -> bool {
    matches!(base, 'A' | 'C' | 'G' | 'T')
}
