//! Helpers the review's unit tests share.

use crate::model::result::ColumnResult;
use crate::model::review_result::{Decision, ItemKind, ItemResult, LabelResult};
use dna_align::model::alignment::Orientation;

pub(super) const F: Orientation = Orientation::Forward;

pub(super) const R: Orientation = Orientation::Reverse;

pub(super) fn read(read: &str) -> String {
    read.to_owned()
}

pub(super) fn item(kind: ItemKind, labels: Vec<LabelResult>) -> ItemResult {
    ItemResult {
        first: ColumnResult {
            position: 1,
            ordinal: 0,
        },
        last: ColumnResult {
            position: 1,
            ordinal: 0,
        },
        kind,
        decision: Decision::Called,
        sequence: None,
        calls: Vec::new(),
        labels,
        evidence: Vec::new(),
    }
}
