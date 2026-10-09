//! Post-calling conformance: checks represented calls against the notation
//! conventions a target profile declares, such as the EMPOP/ISFG placement of
//! homopolymer length variants (ADR-0069).
//!
//! Findings are reported, never applied: the represented calls and the
//! reconstructed haplotype are unchanged. The checker knows no target; the
//! rules and the homopolymer length come from the profile, and runs are read
//! from the reference sequence.

use crate::model::reference::Reference;
use crate::profile::{Conformance, ConformanceRule};
use crate::variant::{Variant, VariantKind};

/// One represented variant that does not follow one rule.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct Finding {
    pub(crate) rule: ConformanceRule,
    pub(crate) variant: Variant,
}

/// Checks one read's represented variants, in variant order and then in the
/// profile's rule order.
pub(crate) fn check(
    reference: &Reference,
    conformance: &Conformance,
    variants: &[Variant],
) -> Vec<Finding> {
    let bases = reference.sequence.as_bytes();
    let mut findings = Vec::new();
    for variant in variants {
        let Some(insertion) = Insertion::of(variant) else {
            continue;
        };
        for &rule in &conformance.rules {
            let follows = match rule {
                ConformanceRule::InsertionAtRunEnd => {
                    insertion.at_run_end(bases, conformance.minimum_run_length)
                }
                ConformanceRule::InsertionMatchesRun => {
                    insertion.matches_run(bases, conformance.minimum_run_length)
                }
            };
            if !follows {
                findings.push(Finding {
                    rule,
                    variant: variant.clone(),
                });
            }
        }
    }
    findings
}

/// A left-anchored insertion: bases inserted after 1-based `anchor`.
struct Insertion<'a> {
    anchor: usize,
    inserted: &'a [u8],
}

impl<'a> Insertion<'a> {
    fn of(variant: &'a Variant) -> Option<Self> {
        if variant.kind != VariantKind::Ins || variant.reference.len() != 1 {
            return None;
        }
        let inserted = variant.alternate.as_bytes().get(1..)?;
        (variant.alternate.as_bytes().first() == variant.reference.as_bytes().first()
            && !inserted.is_empty())
        .then_some(Self {
            anchor: variant.position_1based,
            inserted,
        })
    }

    /// An insertion of one repeated base that touches a homopolymer of that
    /// base must follow the homopolymer's 3' base.
    fn at_run_end(&self, bases: &[u8], minimum: usize) -> bool {
        let base = self.inserted[0];
        if self.inserted.iter().any(|&inserted| inserted != base) {
            return true;
        }
        let touching = [self.anchor, self.anchor + 1]
            .into_iter()
            .find(|&position| reference_base(bases, position) == Some(base));
        let Some(position) = touching else {
            return true;
        };
        let (start, end) = run(bases, position);
        end - start + 1 < minimum || self.anchor == end
    }

    /// Bases inserted between two bases of one homopolymer must be that
    /// homopolymer's base.
    fn matches_run(&self, bases: &[u8], minimum: usize) -> bool {
        let (Some(left), Some(right)) = (
            reference_base(bases, self.anchor),
            reference_base(bases, self.anchor + 1),
        ) else {
            return true;
        };
        if left != right {
            return true;
        }
        let (start, end) = run(bases, self.anchor);
        end - start + 1 < minimum || self.inserted.iter().all(|&inserted| inserted == left)
    }
}

fn reference_base(bases: &[u8], position_1based: usize) -> Option<u8> {
    position_1based
        .checked_sub(1)
        .and_then(|index| bases.get(index))
        .copied()
}

/// The 1-based inclusive run of the base at `position` in the reference.
fn run(bases: &[u8], position: usize) -> (usize, usize) {
    let base = bases[position - 1];
    let mut start = position;
    while start > 1 && bases[start - 2] == base {
        start -= 1;
    }
    let mut end = position;
    while end < bases.len() && bases[end] == base {
        end += 1;
    }
    (start, end)
}

#[cfg(test)]
mod tests {
    use crate::checksum::hex_sha256;
    use crate::model::reference::ReferenceTopology;

    use super::*;

    /// rCRS 300-320: `AAACCCCCCCTCCCCCGCTTC` (an anchored C stretch at 303-315).
    fn reference() -> Reference {
        let sequence = format!("{}AAACCCCCCCTCCCCCGCTTC", "G".repeat(299));
        Reference {
            name: "rCRS".into(),
            sequence_sha256: hex_sha256(sequence.as_bytes()),
            sequence,
            topology: ReferenceTopology::Linear,
        }
    }

    fn conformance() -> Conformance {
        Conformance {
            rules: vec![
                ConformanceRule::InsertionAtRunEnd,
                ConformanceRule::InsertionMatchesRun,
            ],
            minimum_run_length: 3,
        }
    }

    fn insertion(anchor: usize, alternate: &str) -> Variant {
        let reference = reference();
        Variant {
            contig: "rCRS".into(),
            position_1based: anchor,
            reference: reference.sequence[anchor - 1..anchor].into(),
            alternate: alternate.into(),
            kind: VariantKind::Ins,
        }
    }

    fn rules(variants: &[Variant]) -> Vec<(usize, &'static str)> {
        check(&reference(), &conformance(), variants)
            .into_iter()
            .map(|finding| (finding.variant.position_1based, finding.rule.label()))
            .collect()
    }

    #[test]
    fn accepts_insertions_at_the_3_prime_end_of_a_run() {
        assert!(
            rules(&[
                insertion(309, "CC"),
                insertion(315, "CC"),
                insertion(315, "CCC")
            ])
            .is_empty()
        );
    }

    #[test]
    fn flags_a_run_insertion_away_from_its_3_prime_end() {
        assert_eq!(
            rules(&[insertion(305, "CC")]),
            [(305, "insertion_at_run_end")]
        );
        assert_eq!(
            rules(&[insertion(302, "AC")]),
            [(302, "insertion_at_run_end")]
        );
    }

    #[test]
    fn flags_another_base_inserted_inside_a_run() {
        assert_eq!(
            rules(&[insertion(313, "CG")]),
            [(313, "insertion_matches_run")]
        );
        assert!(rules(&[insertion(316, "CA")]).is_empty());
    }

    #[test]
    fn ignores_substitutions_deletions_and_short_runs() {
        let snv = Variant {
            kind: VariantKind::Snv,
            reference: "C".into(),
            alternate: "G".into(),
            ..insertion(314, "CG")
        };
        assert!(rules(&[snv]).is_empty());
        let short = Conformance {
            minimum_run_length: 8,
            ..conformance()
        };
        assert!(check(&reference(), &short, &[insertion(313, "CG")]).is_empty());
    }
}
