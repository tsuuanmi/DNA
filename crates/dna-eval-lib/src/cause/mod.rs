//! The cause of a sample's discrepancy: whether DNA, the lab (its assay, its
//! runs, or its reviewers' truth), or the sequencing technology accounts for
//! it, read from the evidence of the reads, so each difference reaches the
//! one who can act on it.
//!
//! Every discordant position of the sample's class gives a [`Finding`] with
//! its cause, its [`Scope`], and the evidence for it; the sample's cause and
//! scope are those of its first defect, by `dna`, `lab`, `technology`,
//! `undetermined`, or of its first deferred item when it has no defect. The
//! truth is only read.

pub(crate) mod forms;
pub(crate) mod loci;

use std::collections::BTreeSet;

use crate::classes::Discrepancy;
use crate::compare::{self, Allele, Comparison, Failure};
use crate::geometry::Window;
use crate::lod;
use forms::{Forms, Majority};
use loci::{
    BESIDE_PHASE_LOSS, CONTESTED, COVERED, EDGE, INSIDE, Loci, MIXED_CALL, PHASE, RUN_UNREAD,
    SIGNAL, UNRESOLVED, WEAK_CALL, WEAK_REFERENCE,
};

/// Who accounts for a difference from the truth.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) enum Cause {
    /// DNA: its calling, callability, or notation.
    Dna,
    /// The lab: a missing or poor trace, or a truth the reads contradict.
    Lab,
    /// The technology: the lab did its part, but the read's signal shows
    /// dephasing, noise, or length heterogeneity the core cannot resolve.
    /// Always deferred, to signal processing.
    Technology,
    /// The evidence does not decide; adjudication does.
    Undetermined,
}

impl Cause {
    pub(crate) const fn label(self) -> &'static str {
        match self {
            Self::Dna => "dna",
            Self::Lab => "lab",
            Self::Technology => "technology",
            Self::Undetermined => "undetermined",
        }
    }
}

/// Whether a difference is a defect to fix or an item deferred by design.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) enum Scope {
    /// A defect to fix.
    Defect,
    /// Deferred by design: DNA writes no IUPAC code for a point
    /// heteroplasmy (SRS-VAR-009), the reviewers' convention does not record
    /// a call, or the signal needs processing the core does not do
    /// (every `technology` finding); deferrals the user made.
    Deferred,
}

impl Scope {
    pub(crate) const fn label(self) -> &'static str {
        match self {
            Self::Defect => "defect",
            Self::Deferred => "deferred",
        }
    }
}

/// One cause with its evidence.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) struct Finding {
    pub(crate) cause: Cause,
    pub(crate) scope: Scope,
    /// What shows it, such as `no_second_base`.
    pub(crate) evidence: &'static str,
    /// The position or window it concerns, when one does.
    pub(crate) at: Option<String>,
}

impl Finding {
    /// A finding of `cause`: deferred to signal processing for the
    /// technology, a defect otherwise.
    fn new(cause: Cause, evidence: &'static str, at: Option<String>) -> Self {
        Self {
            cause,
            scope: if cause == Cause::Technology {
                Scope::Deferred
            } else {
                Scope::Defect
            },
            evidence,
            at,
        }
    }

    fn deferred(self) -> Self {
        Self {
            scope: Scope::Deferred,
            ..self
        }
    }

    /// `cause:evidence@at`.
    pub(crate) fn text(&self) -> String {
        let at = self
            .at
            .as_deref()
            .map_or_else(String::new, |at| format!("@{at}"));
        format!("{}:{}{at}", self.cause.label(), self.evidence)
    }
}

/// What the cause of one sample is read from.
pub(crate) struct Sample<'a> {
    pub(crate) discrepancy: Discrepancy,
    pub(crate) failure: Option<&'a Failure>,
    /// The comparison under the class's measure.
    pub(crate) comparison: Option<&'a Comparison>,
    /// The primary comparison, whose uncovered alleles a coverage gap names.
    pub(crate) primary: Option<&'a Comparison>,
    /// The windows of the class's discordant positions.
    pub(crate) regions: &'a BTreeSet<Window>,
    pub(crate) dna: &'a [Allele],
    /// DNA's calls a profile convention labels `not_recorded`.
    pub(crate) not_recorded: Option<&'a BTreeSet<Allele>>,
    pub(crate) truth: &'a [Allele],
    pub(crate) loci: Option<&'a Loci>,
    pub(crate) signals: Option<&'a lod::Signals>,
}

/// The findings of one sample, sorted, and its cause and scope; `None` for
/// a sample without a discrepancy.
pub(crate) fn findings(
    sample: &Sample<'_>,
    forms: &Forms,
    reference: &str,
) -> Option<(Cause, Scope, Vec<Finding>)> {
    let mut found = match sample.discrepancy {
        Discrepancy::None => return None,
        Discrepancy::Failed => failed(sample.failure),
        Discrepancy::NotComparable => {
            vec![Finding::new(Cause::Undetermined, "not_comparable", None)]
        }
        Discrepancy::Coverage => coverage(sample),
        Discrepancy::Representation => representation(sample, forms, reference),
        // DNA reports where the reviewers did not look: only a review there
        // verifies it.
        Discrepancy::Unanalysed => sample
            .regions
            .iter()
            .map(|region| {
                Finding::new(
                    Cause::Undetermined,
                    "outside_truth_intervals",
                    Some(region.label().to_owned()),
                )
            })
            .collect(),
        Discrepancy::Heteroplasmy => heteroplasmy(sample),
        Discrepancy::Sequence => sequence(sample),
    };
    if found.is_empty() {
        found.push(Finding::new(Cause::Undetermined, "adjudicate", None));
    }
    found.sort();
    found.dedup();
    // A defect heads the sample before any deferred item, whatever its cause.
    let head = found
        .iter()
        .map(|finding| (finding.scope, finding.cause))
        .min()?;
    Some((head.1, head.0, found))
}

fn failed(failure: Option<&Failure>) -> Vec<Finding> {
    let (cause, evidence) = match failure {
        Some(Failure::TruthMissing) => (Cause::Lab, "truth_missing"),
        Some(Failure::Trace { .. }) => (Cause::Lab, "unreadable_trace"),
        Some(Failure::Rejected { .. }) => (Cause::Lab, "no_admitted_read"),
        Some(Failure::Dna(_) | Failure::Panic | Failure::Notation) => (Cause::Dna, "dna_failure"),
        Some(Failure::Reference { .. }) | None => (Cause::Undetermined, "failure"),
    };
    vec![Finding::new(cause, evidence, None)]
}

/// A truth allele where DNA has no callable base: no read there is the lab's
/// (a missing or short trace). A truth IUPAC code that a read shows is
/// judged as under heteroplasmy. Phase loss is the technology's, and so is a
/// run the consensus left undecided because no read shows both its ends;
/// weak signal is the run's; a tie of clean reads is for the reviewer. A
/// base read cleanly yet untrusted is the technology's beside a stretch
/// masked for phase loss or a read DNA ends at a repeat run; DNA's where a
/// weak-signal veto denies a call that reads the reference (a veto out of
/// its scope); undetermined where a call shows mixed channels; the lab's
/// (read quality) when every read has it within the margin of its own end
/// or with a weak-signal veto; and otherwise DNA's callability policy, as
/// is a base read as `N`.
fn coverage(sample: &Sample<'_>) -> Vec<Finding> {
    let Some(loci) = sample.loci else {
        return Vec::new();
    };
    sample
        .primary
        .map(|comparison| comparison.uncovered.as_slice())
        .unwrap_or_default()
        .iter()
        .filter_map(|allele| {
            let anchor = compare::anchor(&allele.position)?;
            let flags = loci.at(anchor);
            if flags & COVERED != 0
                && let Some(judged) = iupac(allele, anchor, loci, sample.signals)
            {
                return Some(judged);
            }
            let (cause, evidence) = if flags & COVERED == 0 {
                (Cause::Lab, "no_read")
            } else if flags & PHASE != 0 {
                (Cause::Technology, "phase_loss")
            } else if flags & SIGNAL != 0 {
                (Cause::Lab, "weak_signal")
            } else if flags & CONTESTED != 0 {
                (Cause::Undetermined, "contested_reads")
            } else if flags & RUN_UNREAD != 0 {
                (Cause::Technology, "run_unread")
            } else if flags & UNRESOLVED != 0 {
                (Cause::Dna, "unresolved_call")
            } else if flags & BESIDE_PHASE_LOSS != 0 {
                (Cause::Technology, "beside_phase_loss")
            } else if flags & WEAK_REFERENCE != 0 {
                (Cause::Dna, "veto_out_of_scope")
            } else if flags & MIXED_CALL != 0 {
                (Cause::Undetermined, "mixed_call")
            } else if flags & (EDGE | WEAK_CALL) != 0 && flags & INSIDE == 0 {
                if flags & WEAK_CALL != 0 {
                    (Cause::Lab, "weak_call")
                } else {
                    (Cause::Lab, "read_edge")
                }
            } else {
                (Cause::Dna, "untrusted_call")
            };
            Some(Finding::new(cause, evidence, Some(allele.position.clone())))
        })
        .collect()
}

/// A truth IUPAC code of two bases that DNA has no callable base for, as
/// under heteroplasmy: DNA's but deferred when a read shows the second base
/// at the limit of detection, the technology's when one shows it only
/// faintly; `None` otherwise, for the rest of the coverage rules.
fn iupac(
    allele: &Allele,
    position: i64,
    loci: &Loci,
    signals: Option<&lod::Signals>,
) -> Option<Finding> {
    let code = compare::bases(&allele.alternate);
    if !allele.is_iupac() || code.len() != 2 {
        return None;
    }
    let at = Some(allele.position.clone());
    if signals.is_some_and(|signals| signals.shown(position, code)) {
        Some(Finding::new(Cause::Dna, "iupac_not_reported", at).deferred())
    } else if code
        .iter()
        .any(|base| loci.faint.contains(&(position, *base)))
    {
        Some(Finding::new(Cause::Technology, "below_lod", at))
    } else {
        None
    }
}

/// The same sequence written differently: the truth's when DNA writes the
/// reviewers' majority form of the window sequence, DNA's when the truth
/// does. Outside the windows the majority is a proxy: the samples carrying
/// each side's calls.
fn representation(sample: &Sample<'_>, forms: &Forms, reference: &str) -> Vec<Finding> {
    let mut found = Vec::new();
    for region in sample.regions {
        let majority = if *region == Window::Other {
            let (dna, truth) = sample.comparison.map_or((&[][..], &[][..]), |comparison| {
                (
                    comparison.dna_unique.as_slice(),
                    comparison.truth_unique.as_slice(),
                )
            });
            forms.outside(dna, truth).map(|majority| (majority, true))
        } else {
            forms
                .window(region, sample.dna, sample.truth, reference)
                .map(|majority| (majority, false))
        };
        let Some((majority, proxy)) = majority else {
            continue;
        };
        let (cause, evidence) = match (majority, proxy) {
            (Majority::Dna, false) => (Cause::Lab, "truth_minority_form"),
            (Majority::Dna, true) => (Cause::Lab, "truth_minority_form_proxy"),
            (Majority::Truth, false) => (Cause::Dna, "dna_minority_form"),
            (Majority::Truth, true) => (Cause::Dna, "dna_minority_form_proxy"),
            (Majority::Neither, _) => (Cause::Undetermined, "no_majority_form"),
        };
        found.push(Finding::new(
            cause,
            evidence,
            Some(region.label().to_owned()),
        ));
    }
    found
}

/// A truth IUPAC code where DNA reports a base: DNA's when a read shows the
/// second base at the limit of detection, deferred since DNA writes no IUPAC
/// code by design; the
/// technology's when one shows it only faintly; the truth's, for a re-check,
/// when clean reads on both strands show no second base at all.
fn heteroplasmy(sample: &Sample<'_>) -> Vec<Finding> {
    let (Some(loci), Some(signals), Some(comparison)) =
        (sample.loci, sample.signals, sample.comparison)
    else {
        return Vec::new();
    };
    sample
        .truth
        .iter()
        .filter(|allele| allele.is_iupac() && comparison.positions.contains(&allele.position))
        .filter_map(|allele| {
            let position = allele.position.parse::<i64>().ok()?;
            let code = compare::bases(&allele.alternate);
            let flags = loci.at(position);
            let (cause, evidence) = if code.len() != 2 {
                (Cause::Undetermined, "iupac_code")
            } else if signals.shown(position, code) {
                return Some(
                    Finding::new(
                        Cause::Dna,
                        "iupac_not_reported",
                        Some(allele.position.clone()),
                    )
                    .deferred(),
                );
            } else if code
                .iter()
                .any(|base| loci.faint.contains(&(position, *base)))
            {
                (Cause::Technology, "below_lod")
            } else if flags & loci::FORWARD != 0 && flags & loci::REVERSE != 0 {
                (Cause::Lab, "no_second_base")
            } else {
                (Cause::Undetermined, "one_strand")
            };
            Some(Finding::new(cause, evidence, Some(allele.position.clone())))
        })
        .collect()
}

/// The sequences differ: when every DNA-only call is one a profile
/// convention labels `not_recorded` and the truth has no call of its own,
/// the reviewers' convention, deferred (SRS-NOM-017). Otherwise a DNA-only
/// call where clean reads on both strands show the reference is DNA's (a
/// phantom); a truth-only call there is the truth's, for a re-check.
/// Insertions are left to adjudication, since reads showing the reference at
/// the anchor do not exclude one.
fn sequence(sample: &Sample<'_>) -> Vec<Finding> {
    let (Some(loci), Some(comparison)) = (sample.loci, sample.comparison) else {
        return Vec::new();
    };
    if let Some(not_recorded) = sample.not_recorded
        && comparison.truth_unique.is_empty()
        && !comparison.dna_unique.is_empty()
        && comparison
            .dna_unique
            .iter()
            .all(|allele| not_recorded.contains(allele))
    {
        return sample
            .regions
            .iter()
            .map(|region| {
                Finding::new(
                    Cause::Lab,
                    "convention_not_recorded",
                    Some(region.label().to_owned()),
                )
                .deferred()
            })
            .collect();
    }
    let shown = |allele: &Allele| {
        !allele.position.contains('.')
            && comparison.positions.contains(&allele.position)
            && compare::anchor(&allele.position)
                .is_some_and(|anchor| loci.reference_on_both_strands(anchor))
    };
    let dna = comparison
        .dna_unique
        .iter()
        .filter(|allele| shown(allele))
        .map(|allele| Finding::new(Cause::Dna, "phantom", Some(allele.position.clone())));
    let truth = comparison
        .truth_unique
        .iter()
        .filter(|allele| shown(allele))
        .map(|allele| {
            Finding::new(
                Cause::Lab,
                "reads_show_reference",
                Some(allele.position.clone()),
            )
        });
    dna.chain(truth).collect()
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeSet;

    use super::*;
    use crate::compare::Concordance;
    use loci::{COVERED, FORWARD, OTHER, REVERSE};

    fn comparison(
        positions: &[&str],
        dna_unique: Vec<Allele>,
        truth_unique: Vec<Allele>,
        uncovered: Vec<Allele>,
    ) -> Comparison {
        Comparison {
            concordance: Concordance::Discordant,
            dna_alleles: 0,
            truth_alleles: 0,
            dna_unique,
            truth_unique,
            reasons: BTreeSet::new(),
            positions: positions
                .iter()
                .map(|position| (*position).to_owned())
                .collect(),
            iupac_only: false,
            uncovered,
            blocks: Vec::new(),
        }
    }

    fn texts(sample: &Sample<'_>) -> Vec<String> {
        findings(sample, &Forms::default(), "")
            .map(|(_, _, findings)| findings.iter().map(Finding::text).collect())
            .unwrap_or_default()
    }

    fn sample<'a>(
        discrepancy: Discrepancy,
        comparison: &'a Comparison,
        truth: &'a [Allele],
        loci: &'a Loci,
        signals: &'a lod::Signals,
        regions: &'a BTreeSet<Window>,
    ) -> Sample<'a> {
        Sample {
            discrepancy,
            failure: None,
            comparison: Some(comparison),
            primary: Some(comparison),
            regions,
            dna: &[],
            not_recorded: None,
            truth,
            loci: Some(loci),
            signals: Some(signals),
        }
    }

    #[test]
    fn a_coverage_gap_goes_to_whoever_left_the_position_unread() {
        // Positions 1-4: unread, phase loss, weak signal, read cleanly.
        // Positions 5-6: a tie of clean reads, a run no read shows whole;
        // 7-9: untrusted beside phase loss, at the read's end, inside; 10: a
        // mixed call at a truth IUPAC code no read shows the second base of;
        // 11: a weak call at a read end; 12: a weak call of the reference.
        let loci = Loci::of(
            1,
            vec![
                0,
                COVERED | PHASE,
                COVERED | SIGNAL,
                COVERED | FORWARD,
                COVERED | CONTESTED,
                COVERED | RUN_UNREAD,
                COVERED | FORWARD | BESIDE_PHASE_LOSS | EDGE,
                COVERED | FORWARD | EDGE,
                COVERED | FORWARD | REVERSE | EDGE | INSIDE,
                COVERED | FORWARD | MIXED_CALL,
                COVERED | FORWARD | WEAK_CALL | EDGE,
                COVERED | FORWARD | WEAK_REFERENCE | EDGE,
            ],
            &[],
        );
        let mut uncovered = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "11", "12"]
            .map(|position| Allele::new(position, "A", "G"))
            .to_vec();
        uncovered.push(Allele::new("10", "A", "R"));
        let comparison = comparison(&[], Vec::new(), Vec::new(), uncovered);
        let signals = lod::Signals::default();
        let regions = BTreeSet::new();
        let gap = sample(
            Discrepancy::Coverage,
            &comparison,
            &[],
            &loci,
            &signals,
            &regions,
        );
        assert_eq!(
            texts(&gap),
            [
                "dna:untrusted_call@4",
                "dna:untrusted_call@9",
                "dna:veto_out_of_scope@12",
                "lab:no_read@1",
                "lab:read_edge@8",
                "lab:weak_call@11",
                "lab:weak_signal@3",
                "technology:beside_phase_loss@7",
                "technology:phase_loss@2",
                "technology:run_unread@6",
                "undetermined:contested_reads@5",
                "undetermined:mixed_call@10"
            ]
        );
        assert_eq!(
            findings(&gap, &Forms::default(), "").map(|found| found.0),
            Some(Cause::Dna)
        );
    }

    #[test]
    fn a_heteroplasmy_is_the_truths_only_when_both_strands_show_no_second_base() {
        // 16189Y (C/T) where DNA calls T, at offsets 0-3 of a region from 16189.
        let both = COVERED | FORWARD | REVERSE;
        let loci = Loci::of(
            16189,
            vec![both, both, COVERED | FORWARD, both],
            &[(16190, 'C')],
        );
        let truth = ["16189", "16190", "16191", "16192"]
            .map(|position| Allele::new(position, "T", "Y"))
            .to_vec();
        let positions = ["16189", "16190", "16191", "16192"];
        let comparison = comparison(&positions, Vec::new(), Vec::new(), Vec::new());
        let signals = lod::Signals {
            calls: BTreeSet::new(),
            seconds: BTreeSet::from([(16192, 'C')]),
        };
        let regions = BTreeSet::new();
        let site = sample(
            Discrepancy::Heteroplasmy,
            &comparison,
            &truth,
            &loci,
            &signals,
            &regions,
        );
        assert_eq!(
            texts(&site),
            [
                "dna:iupac_not_reported@16192",
                "lab:no_second_base@16189",
                "technology:below_lod@16190",
                "undetermined:one_strand@16191"
            ]
        );
    }

    #[test]
    fn a_defect_heads_the_sample_before_a_deferred_item() -> Result<(), &'static str> {
        // 16189Y shown at the limit (DNA's, deferred) and 16190Y with no
        // second base on either clean strand (the truth's, a defect).
        let both = COVERED | FORWARD | REVERSE;
        let loci = Loci::of(16189, vec![both, both], &[]);
        let truth = [
            Allele::new("16189", "T", "Y"),
            Allele::new("16190", "T", "Y"),
        ];
        let comparison = comparison(&["16189", "16190"], Vec::new(), Vec::new(), Vec::new());
        let signals = lod::Signals {
            calls: BTreeSet::new(),
            seconds: BTreeSet::from([(16189, 'C')]),
        };
        let regions = BTreeSet::new();
        let site = sample(
            Discrepancy::Heteroplasmy,
            &comparison,
            &truth,
            &loci,
            &signals,
            &regions,
        );
        let (cause, scope, found) =
            findings(&site, &Forms::default(), "").ok_or("a discrepancy has findings")?;
        assert_eq!((cause, scope), (Cause::Lab, Scope::Defect));
        assert_eq!(found[0].text(), "dna:iupac_not_reported@16189");
        Ok(())
    }

    #[test]
    fn a_call_the_reads_contradict_on_both_strands_is_a_phantom_or_a_truth_to_recheck() {
        let both = COVERED | FORWARD | REVERSE;
        let loci = Loci::of(100, vec![both, both, both | OTHER], &[]);
        let comparison = comparison(
            &["100", "101", "102"],
            vec![Allele::new("100", "A", "G")],
            vec![Allele::new("101", "C", "T"), Allele::new("102", "G", "A")],
            Vec::new(),
        );
        let signals = lod::Signals::default();
        let regions = BTreeSet::new();
        let differs = sample(
            Discrepancy::Sequence,
            &comparison,
            &[],
            &loci,
            &signals,
            &regions,
        );
        // 102: a read shows something other than the reference.
        assert_eq!(
            texts(&differs),
            ["dna:phantom@100", "lab:reads_show_reference@101"]
        );
    }

    #[test]
    fn a_length_call_the_reviewers_do_not_record_is_their_convention() {
        let loci = Loci::of(16193, vec![COVERED | FORWARD | REVERSE], &[]);
        let length = Allele::new("16193.1", "-", "C");
        let not_recorded = BTreeSet::from([length.clone()]);
        let signals = lod::Signals::default();
        let regions = BTreeSet::from([Window::Named("hvs_i".to_owned())]);
        let texts_of = |dna_unique: Vec<Allele>, truth_unique: Vec<Allele>| {
            let comparison = comparison(&["16193.1"], dna_unique, truth_unique, Vec::new());
            let differs = Sample {
                not_recorded: Some(&not_recorded),
                ..sample(
                    Discrepancy::Sequence,
                    &comparison,
                    &[],
                    &loci,
                    &signals,
                    &regions,
                )
            };
            findings(&differs, &Forms::default(), "").map(|(cause, scope, findings)| {
                (
                    cause.label(),
                    scope.label(),
                    findings.iter().map(Finding::text).collect::<Vec<_>>(),
                )
            })
        };
        assert_eq!(
            texts_of(vec![length.clone()], Vec::new()),
            Some((
                "lab",
                "deferred",
                vec!["lab:convention_not_recorded@hvs_i".to_owned()]
            ))
        );
        // A truth-only call, or an unlabelled DNA-only call, is content.
        for (dna, truth) in [
            (vec![length.clone()], vec![Allele::new("16189", "T", "C")]),
            (vec![length, Allele::new("16194", "A", "C")], Vec::new()),
        ] {
            assert_ne!(
                texts_of(dna, truth).map(|(_, _, findings)| findings),
                Some(vec!["lab:convention_not_recorded@hvs_i".to_owned()])
            );
        }
    }

    #[test]
    fn a_technology_limit_is_deferred_to_signal_processing() {
        // 1: phase loss; 2: a clean call withheld by the read's own end.
        let loci = Loci::of(1, vec![COVERED | PHASE, COVERED | FORWARD | EDGE], &[]);
        let signals = lod::Signals::default();
        let regions = BTreeSet::new();
        let head = |positions: &[&str]| {
            let uncovered = positions
                .iter()
                .map(|position| Allele::new(position, "A", "G"))
                .collect();
            let comparison = comparison(&[], Vec::new(), Vec::new(), uncovered);
            findings(
                &sample(
                    Discrepancy::Coverage,
                    &comparison,
                    &[],
                    &loci,
                    &signals,
                    &regions,
                ),
                &Forms::default(),
                "",
            )
            .map(|(cause, scope, _)| (cause, scope))
        };
        assert_eq!(head(&["1"]), Some((Cause::Technology, Scope::Deferred)));
        // A lab defect heads a sample that also meets a technology limit.
        assert_eq!(head(&["1", "2"]), Some((Cause::Lab, Scope::Defect)));
    }

    #[test]
    fn a_failure_names_its_cause() {
        let failure = |failure: Failure| failed(Some(&failure))[0].text();
        assert_eq!(failure(Failure::TruthMissing), "lab:truth_missing");
        assert_eq!(
            failure(Failure::Rejected { reads: 2 }),
            "lab:no_admitted_read"
        );
        assert_eq!(failure(Failure::Panic), "dna:dna_failure");
    }
}
