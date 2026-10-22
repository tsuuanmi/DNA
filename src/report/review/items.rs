//! The review's items: the drafts over the consensus sites, the second bases,
//! and the gaps, and each draft finished into an item with its labels, kind,
//! and evidence, or dropped when it is DNA's own decision.

use super::Sample;
use super::labels::{key, priority};
use super::ranges::{inside, within};
use super::sites::{beside_not_analysed, called_near, one_strand};
use crate::error::Result;
use crate::model::result::ColumnResult;
use crate::model::review_result::{Decision, ItemKind, ItemResult, LabelResult};
use dna_consensus::model::consensus::ConsensusSite;
use std::collections::BTreeSet;

/// An item before its labels are complete.
pub(super) struct Draft<'a> {
    pub(super) first: usize,
    pub(super) last: usize,
    pub(super) site: Option<&'a ConsensusSite>,
    /// The decided sequence: a site's call, or the reference base of a
    /// position outside every site.
    pub(super) sequence: Option<String>,
    /// Whether the consensus decided the item without an `N`.
    pub(super) decided: bool,
    pub(super) labels: Vec<LabelResult>,
}

impl<'a> Sample<'a> {
    /// Every item, in priority order.
    pub(super) fn items(&self) -> Result<Vec<ItemResult>> {
        let in_site = |position: usize| {
            self.sites
                .iter()
                .any(|site| site.first_0based < position && position <= site.last_0based + 1)
        };
        let mut drafts = Vec::new();
        // A: every consensus site, clipped to the regions.
        for site in self.sites {
            if let Some([first, last]) =
                within(self.regions, site.first_0based + 1, site.last_0based + 1)
            {
                let mut labels = self.site_labels(site);
                let mut sequence = site.call.clone();
                if let Some(shown) = self.clipped_decision(site, first, last) {
                    self.clip_disagreement(&mut labels, first, last, &shown);
                    sequence = Some(shown);
                }
                drafts.push(Draft {
                    first,
                    last,
                    site: Some(site),
                    sequence,
                    decided: site.call.as_ref().is_some_and(|call| !call.contains('N')),
                    labels,
                });
            }
        }
        // B: callable positions outside every site where a second base shows;
        // C: maximal runs of positions outside every site and not callable.
        let mut run: Option<[usize; 2]> = None;
        for position in self.regions.iter().flat_map(|&[start, end]| start..=end) {
            let open = !in_site(position) && !inside(self.callable, position);
            match (&mut run, open) {
                (Some([_, last]), true) if *last + 1 == position => *last = position,
                (_, true) => {
                    if let Some([first, last]) = run.replace([position, position]) {
                        drafts.push(self.gap(first, last));
                    }
                }
                (_, false) => {
                    if let Some([first, last]) = run.take() {
                        drafts.push(self.gap(first, last));
                    }
                }
            }
            if !in_site(position)
                && inside(self.callable, position)
                && let Some(label) = self.heteroplasmy(position)
            {
                drafts.push(Draft {
                    first: position,
                    last: position,
                    site: None,
                    sequence: Some(self.base(position).to_string()),
                    decided: true,
                    labels: vec![label],
                });
            }
        }
        if let Some([first, last]) = run {
            drafts.push(self.gap(first, last));
        }
        // A call no other item holds, for example one the notation moved.
        let held = drafts
            .iter()
            .map(|draft| [draft.first, draft.last])
            .collect::<Vec<_>>();
        let lone = self
            .calls
            .iter()
            .map(|call| call.position)
            .filter(|&position| {
                inside(self.regions, position)
                    && !held
                        .iter()
                        .any(|&[first, last]| first <= position && position <= last)
            })
            .collect::<BTreeSet<_>>();
        for position in lone {
            drafts.push(Draft {
                first: position,
                last: position,
                site: None,
                sequence: None,
                decided: true,
                labels: Vec::new(),
            });
        }

        // D: every stretch the call does not analyse (SRS-CONS-012). It
        // leaves the region's reported range short, so an allele there would
        // be omitted from the report: a check, never a silent pass.
        for span in self.not_analysed {
            drafts.push(Draft {
                first: span.first,
                last: span.last,
                site: None,
                sequence: None,
                decided: false,
                labels: vec![LabelResult::NotAnalysed {
                    reasons: span.reasons.clone(),
                }],
            });
        }

        self.length_labels(&mut drafts);
        let mut items = drafts
            .into_iter()
            .map(|draft| self.finish(draft))
            .filter_map(Result::transpose)
            .collect::<Result<Vec<_>>>()?;
        items.sort_by_key(|item| (priority(&item.labels), item.first, item.last));
        Ok(items)
    }

    /// A run of positions outside every site that are not callable.
    fn gap(&self, first: usize, last: usize) -> Draft<'a> {
        let covered = self
            .covered
            .iter()
            .any(|&[start, end]| start <= last && first <= end);
        Draft {
            first,
            last,
            site: None,
            sequence: None,
            decided: false,
            labels: vec![if covered {
                LabelResult::Undecided
            } else {
                LabelResult::NotCovered
            }],
        }
    }

    /// Completes an item's calls, labels, and evidence; `None` for a site
    /// without a call or a label.
    fn finish(&self, draft: Draft<'_>) -> Result<Option<ItemResult>> {
        let Draft {
            first,
            last,
            site,
            sequence,
            decided,
            mut labels,
        } = draft;
        let calls = self
            .calls
            .iter()
            .filter(|call| first <= call.position && call.position <= last)
            .map(|call| call.text.clone())
            .collect::<Vec<_>>();
        // A difference: DNA's calls when it has them in the profile notation,
        // otherwise the consensus's own decision. A decided difference that
        // DNA does not call within a position of the item is `not_called`.
        let decided_difference = site
            .is_some_and(|site| decided && site.call.as_deref() != Some(site.reference.as_str()));
        let differs = if self.notation {
            !calls.is_empty()
        } else {
            decided_difference
        };
        if self.notation && decided_difference && !called_near(self.calls, first, last) {
            labels.push(LabelResult::NotCalled);
        }
        if let (Some(site), Some(sequence), true) = (site, sequence.as_deref(), decided) {
            labels.extend(self.disagreeing(site, first, last, sequence));
            let ends = self.end_disagreeing(first, last, site, sequence, &labels)?;
            labels.extend(ends);
        }
        if let (true, Some(site)) = (differs, site) {
            let indices = site
                .supporting
                .iter()
                .filter_map(|input| self.by_id.get(input.as_str()).copied())
                .collect::<Vec<_>>();
            let supporting = indices
                .iter()
                .map(|&index| &self.reads[index])
                .collect::<Vec<_>>();
            if let [read] = supporting.as_slice() {
                labels.push(LabelResult::OneRead {
                    read: read.name.to_owned(),
                });
            }
            if let Some(strand) = one_strand(
                &supporting
                    .iter()
                    .map(|read| read.strand)
                    .collect::<Vec<_>>(),
            ) {
                labels.push(LabelResult::OneStrand {
                    strand,
                    reads: supporting.iter().map(|read| read.name.to_owned()).collect(),
                });
                for &index in &indices {
                    labels.extend(self.reference_second(index, first, last)?);
                }
                if decided {
                    labels.extend(self.one_strand_insertion(site, first, last, strand, &indices));
                }
            }
        }
        if differs
            && beside_not_analysed(
                &self.reportable,
                &self
                    .not_analysed
                    .iter()
                    .map(|span| [span.first, span.last])
                    .collect::<Vec<_>>(),
                first,
                last,
                self.config.edge_zone,
            )
        {
            labels.push(LabelResult::BesideNotAnalysed);
        }
        labels.extend((first..=last).filter_map(|position| self.heteroplasmy(position)));
        let (read_labels, evidence) = self.reads_over(first, last)?;
        labels.extend(read_labels);
        // A label that leaves the item unassured makes it a check; the others
        // are context for the reviewer. An item with neither a call nor such a
        // label is DNA's own decision and is not listed.
        // A difference inside a read's poor signal is a check too: there the
        // reads' differences are often wrong (SRS-CALL-024).
        let poor = labels
            .iter()
            .any(|label| matches!(label, LabelResult::PoorSignal { .. }));
        let assured = !(labels.iter().any(LabelResult::unassured) || (differs && poor));
        let noted = labels.iter().any(|label| {
            matches!(
                label,
                LabelResult::PossibleHeteroplasmy {
                    explained: true,
                    ..
                } | LabelResult::TractLength { .. }
            )
        });
        let kind = match (assured, differs, noted) {
            (false, _, _) => ItemKind::Check,
            (true, true, _) => ItemKind::Variant,
            (true, false, true) => ItemKind::Note,
            (true, false, false) => return Ok(None),
        };
        labels.sort_by(|left, right| key(left).cmp(&key(right)));
        labels.dedup();
        let callable = (first..=last)
            .filter(|&position| inside(self.callable, position))
            .count();
        let decision = if callable == 0 {
            Decision::NotCallable
        } else if decided && callable == last - first + 1 {
            Decision::Called
        } else {
            Decision::PartlyCallable
        };
        Ok(Some(ItemResult {
            first: ColumnResult {
                position: first,
                ordinal: 0,
            },
            last: ColumnResult {
                position: last,
                ordinal: 0,
            },
            kind,
            decision,
            sequence,
            calls,
            labels,
            evidence,
        }))
    }
}
