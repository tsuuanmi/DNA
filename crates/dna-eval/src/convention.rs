//! `dna-eval convention`: how far DNA's notation chain agrees with the way
//! the lab's reviewers write their own calls.
//!
//! Every truth sample's alleles (IUPAC-free samples only) go through DNA's
//! right-aligned normalization and the profile's nomenclature windows, as
//! DNA's own calls do. A sample whose written calls come back unchanged is a
//! fixed point of DNA's notation; a changed one shows a convention DNA writes
//! differently. Reviewer truth is only read, never rewritten or compared as
//! DNA output (ADR-0075).

use std::collections::BTreeMap;
use std::fmt::Write as _;

use dna_kernel::model::reference::Reference;
use dna_kernel::profile::Profile;
use dna_kernel::variant::{CalledVariantSet, Variant, VariantKind};
use dna_normalize::variant_nomenclature::{apply_with, from_normalization};
use dna_normalize::variant_normalization::{NormalizationPolicy, TradeBounds, normalize_read};

use crate::compare::Allele;
use crate::truth::TruthSample;

/// Largest number of change patterns the report lists.
const MAX_PATTERNS: usize = 200;

/// What DNA's notation does to the reviewers' written calls.
#[derive(Debug, Default)]
pub(crate) struct Convention {
    /// Samples whose calls went through the chain.
    pub(crate) checked: usize,
    /// Samples left out: with an IUPAC call, or calls DNA cannot apply.
    pub(crate) skipped: usize,
    /// Samples whose calls come back unchanged.
    pub(crate) unchanged: usize,
    /// Per change, the truth-only and DNA-only calls, and its samples.
    pub(crate) patterns: BTreeMap<(String, String), usize>,
}

/// Runs every sample's calls through the notation chain of `profile`.
pub(crate) fn check<'a>(
    samples: impl IntoIterator<Item = &'a TruthSample>,
    reference: &Reference,
    profile: &Profile,
) -> Convention {
    let mut convention = Convention::default();
    // As the sample workflow bounds a deletion trade, with every base read.
    let windows = profile
        .windows
        .iter()
        .map(|window| window.start..window.end())
        .collect::<Vec<_>>();
    let bounds = TradeBounds {
        windows: &windows,
        readable: None,
    };
    for sample in samples {
        let Some(variants) = variants(&sample.alleles, &reference.sequence, &reference.name) else {
            convention.skipped += 1;
            continue;
        };
        let called = CalledVariantSet {
            reference: reference.identity(),
            variants,
        };
        let represented = normalize_read(
            reference,
            &called,
            NormalizationPolicy::RightAligned,
            bounds,
        )
        .and_then(|normalized| apply_with(reference, profile, from_normalization(&normalized)));
        let Ok(represented) = represented else {
            convention.skipped += 1;
            continue;
        };
        convention.checked += 1;
        let written = texts(&sample.alleles);
        let rewritten = texts(&per_base(&represented.represented_variants));
        if written == rewritten {
            convention.unchanged += 1;
            continue;
        }
        let only = |left: &[String], right: &[String]| {
            left.iter()
                .filter(|call| !right.contains(call))
                .cloned()
                .collect::<Vec<_>>()
                .join(" ")
        };
        *convention
            .patterns
            .entry((only(&written, &rewritten), only(&rewritten, &written)))
            .or_default() += 1;
    }
    convention
}

impl Convention {
    /// `convention.tsv`: the change patterns, most samples first.
    pub(crate) fn table(&self) -> String {
        let mut text = String::from("samples\ttruth_calls\tdna_calls\n");
        for ((truth, dna), samples) in self.sorted().into_iter().take(MAX_PATTERNS) {
            let _ = writeln!(text, "{samples}\t{truth}\t{dna}");
        }
        text
    }

    /// One summary line.
    pub(crate) fn summary(&self) -> String {
        format!(
            "dna-eval convention: {} samples checked, {} unchanged ({:.4}), {} changed in {} patterns, {} skipped",
            self.checked,
            self.unchanged,
            if self.checked == 0 {
                0.0
            } else {
                self.unchanged as f64 / self.checked as f64
            },
            self.checked - self.unchanged,
            self.patterns.len(),
            self.skipped
        )
    }

    fn sorted(&self) -> Vec<(&(String, String), usize)> {
        let mut sorted = self
            .patterns
            .iter()
            .map(|(pattern, samples)| (pattern, *samples))
            .collect::<Vec<_>>();
        sorted.sort_by(|left, right| right.1.cmp(&left.1).then_with(|| left.0.cmp(right.0)));
        sorted
    }
}

/// The calls in the reviewer notation, sorted and without repeats.
fn texts(alleles: &[Allele]) -> Vec<String> {
    let mut calls = alleles
        .iter()
        .filter(|allele| allele.reference != allele.alternate)
        .map(Allele::text)
        .collect::<Vec<_>>();
    calls.sort();
    calls.dedup();
    calls
}

/// Anchored variants of per-base calls: runs of deleted bases and insertion
/// ordinals become one deletion or insertion after the preceding reference
/// base. `None` when a call is not a plain A/C/G/T call DNA could make.
fn variants(alleles: &[Allele], sequence: &str, contig: &str) -> Option<Vec<Variant>> {
    let base = |position: usize| sequence.get(position - 1..position);
    let mut substitutions = BTreeMap::new();
    let mut deletions = Vec::new();
    let mut insertions = BTreeMap::<usize, BTreeMap<usize, char>>::new();
    for allele in alleles
        .iter()
        .filter(|allele| allele.reference != allele.alternate)
    {
        let definite = |text: &str| text.len() == 1 && "ACGT".contains(text);
        let (anchor, ordinal) = match allele.position.split_once('.') {
            Some((anchor, ordinal)) => (
                anchor.parse::<usize>().ok()?,
                ordinal.parse::<usize>().ok()?,
            ),
            None => (allele.position.parse::<usize>().ok()?, 0),
        };
        if anchor == 0 || anchor > sequence.len() {
            return None;
        }
        match (ordinal, allele.alternate.as_str()) {
            (0, "-") => deletions.push(anchor),
            (0, alternate) if definite(alternate) => {
                substitutions.insert(anchor, alternate.to_owned());
            }
            (ordinal, alternate) if ordinal > 0 && definite(alternate) => {
                insertions
                    .entry(anchor)
                    .or_default()
                    .insert(ordinal, alternate.chars().next()?);
            }
            _ => return None,
        }
    }
    let mut variants = Vec::new();
    for (position, alternate) in substitutions {
        variants.push(Variant {
            contig: contig.to_owned(),
            position_1based: position,
            reference: base(position)?.to_owned(),
            alternate,
            kind: VariantKind::Snv,
        });
    }
    deletions.sort_unstable();
    deletions.dedup();
    let mut index = 0;
    while index < deletions.len() {
        let first = deletions[index];
        let mut last = first;
        while deletions.get(index + 1) == Some(&(last + 1)) {
            index += 1;
            last += 1;
        }
        index += 1;
        let anchor = first.checked_sub(1).filter(|anchor| *anchor > 0)?;
        variants.push(Variant {
            contig: contig.to_owned(),
            position_1based: anchor,
            reference: sequence.get(anchor - 1..last)?.to_owned(),
            alternate: base(anchor)?.to_owned(),
            kind: VariantKind::Del,
        });
    }
    for (anchor, ordinals) in insertions {
        // Ordinals must run 1, 2, ... without a gap.
        if ordinals.keys().copied().ne(1..=ordinals.len()) {
            return None;
        }
        let inserted = ordinals.values().collect::<String>();
        variants.push(Variant {
            contig: contig.to_owned(),
            position_1based: anchor,
            reference: base(anchor)?.to_owned(),
            alternate: format!("{}{inserted}", base(anchor)?),
            kind: VariantKind::Ins,
        });
    }
    variants.sort_by_key(|variant| variant.position_1based);
    Some(variants)
}

/// Per-base calls of anchored variants, as DNA's notation writes them.
fn per_base(variants: &[Variant]) -> Vec<Allele> {
    let mut alleles = Vec::new();
    for variant in variants {
        let position = variant.position_1based;
        match variant.kind {
            VariantKind::Ins => {
                for (ordinal, base) in variant.alternate.chars().skip(1).enumerate() {
                    alleles.push(Allele::new(
                        &format!("{position}.{}", ordinal + 1),
                        "-",
                        &base.to_string(),
                    ));
                }
            }
            VariantKind::Del => {
                for (offset, base) in variant.reference.chars().enumerate().skip(1) {
                    alleles.push(Allele::new(
                        &(position + offset).to_string(),
                        &base.to_string(),
                        "-",
                    ));
                }
            }
            _ => {
                for (offset, (reference, alternate)) in variant
                    .reference
                    .chars()
                    .zip(variant.alternate.chars())
                    .enumerate()
                {
                    if reference != alternate {
                        alleles.push(Allele::new(
                            &(position + offset).to_string(),
                            &reference.to_string(),
                            &alternate.to_string(),
                        ));
                    }
                }
            }
        }
    }
    alleles
}

#[cfg(test)]
mod tests {
    use std::path::Path;

    use dna_kernel::model::reference::ReferenceTopology;

    use super::*;

    fn rcrs() -> Result<(Reference, Profile), Box<dyn std::error::Error>> {
        let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
        let reference = dna_kernel::reference::load(
            &root.join("references/rCRS.fasta"),
            ReferenceTopology::Linear,
        )?;
        let profile = Profile::load(&root.join("config/profiles/human-mtdna-rcrs.toml"))?;
        Ok((reference, profile))
    }

    fn sample(calls: &[(&str, &str, &str)]) -> TruthSample {
        TruthSample {
            batch: "b".into(),
            status: None,
            alleles: calls
                .iter()
                .map(|(position, reference, alternate)| Allele::new(position, reference, alternate))
                .collect(),
            intervals: Vec::new(),
        }
    }

    #[test]
    fn counts_fixed_points_and_change_patterns() -> Result<(), Box<dyn std::error::Error>> {
        let (reference, profile) = rcrs()?;
        let lab = sample(&[("263", "A", "G"), ("309.1", "-", "C"), ("315.1", "-", "C")]);
        // A left-placed C-run insertion DNA writes at the run's 3' end.
        let left = sample(&[("16184.1", "-", "C")]);
        let iupac = sample(&[("16189", "T", "Y")]);
        let convention = check([&lab, &left, &iupac], &reference, &profile);
        assert_eq!(
            (convention.checked, convention.unchanged, convention.skipped),
            (2, 1, 1)
        );
        assert_eq!(
            convention.patterns,
            BTreeMap::from([(("16184.1C".to_owned(), "16188.1C".to_owned()), 1)])
        );
        assert!(convention.table().contains("1\t16184.1C\t16188.1C\n"));
        Ok(())
    }

    #[test]
    fn rebuilds_runs_of_deletions_and_insertion_ordinals() -> Result<(), Box<dyn std::error::Error>>
    {
        let (reference, _) = rcrs()?;
        let calls = [
            Allele::new("523", "A", "-"),
            Allele::new("524", "C", "-"),
            Allele::new("315.1", "-", "C"),
            Allele::new("315.2", "-", "C"),
        ];
        let built = variants(&calls, &reference.sequence, "rCRS").ok_or("variants")?;
        let shapes = built
            .iter()
            .map(|variant| {
                (
                    variant.position_1based,
                    variant.reference.as_str(),
                    variant.alternate.as_str(),
                )
            })
            .collect::<Vec<_>>();
        assert_eq!(shapes, [(315, "C", "CCC"), (522, "CAC", "C")]);
        assert_eq!(texts(&per_base(&built)), texts(&calls));
        // An ordinal gap is not a call DNA could make.
        assert!(
            variants(
                &[Allele::new("315.2", "-", "C")],
                &reference.sequence,
                "rCRS"
            )
            .is_none()
        );
        Ok(())
    }
}
