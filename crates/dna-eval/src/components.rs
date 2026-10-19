//! The components of every read's phase mixtures (ADR-0082), each reference
//! run cluster they cover measured against the truth: `components.tsv`.
//!
//! Reference homopolymers of at least [`MINIMUM_RUN`] bases one base apart
//! form a cluster (`303-309`, `310`, `311-315`), measured over a window with
//! one more base on each side. The component's haplotype there is the
//! reference with its edits applied, and the truth's the reference with its
//! alleles, so neither depends on how the aligner splits an edit. The run a
//! read enters in phase, the first after the window's first base for a
//! forward read and the last before its last base for a reverse one, is
//! written as its base and length for both (`entered`, `truth_entered`), and
//! `agrees` says whether the component's haplotype is the truth's over the
//! whole window. An edit across the window's edge, or an IUPAC truth symbol
//! in it, leaves it unmeasured (`.`).

use std::fmt::Write as _;

use dna::sample_analysis::SampleAnalysis;

use crate::compare::Allele;

/// Shortest reference homopolymer measured.
pub(crate) const MINIMUM_RUN: usize = 4;

/// One component of one phase mixture of one read.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Row {
    pub(crate) read: String,
    pub(crate) forward: bool,
    pub(crate) mixture: [usize; 2],
    pub(crate) explained: bool,
    pub(crate) components: usize,
    pub(crate) offset: i32,
    pub(crate) fraction: f64,
    pub(crate) residual: f64,
    pub(crate) unexplained: usize,
    /// 1-based reference interval, inclusive, crossing the origin when its
    /// first position exceeds its last; `None` for a mixture not decomposed.
    pub(crate) span: Option<[usize; 2]>,
    /// Its edits: 1-based anchor, reference and alternate alleles.
    pub(crate) edits: Vec<(usize, String, String)>,
}

// The fit gives finite fractions and residuals, so equality is reflexive.
impl Eq for Row {}

/// The component rows of an analysed sample: an explained mixture's placed
/// components with their calls, and the offsets and fractions alone of a
/// mixture none of whose components was compared (unexplained, or none
/// placed).
pub(crate) fn rows(analysis: &SampleAnalysis) -> Vec<Row> {
    let mut rows = Vec::new();
    for read in analysis.components() {
        for mixture in read.mixtures {
            let assertion = &mixture.mixture;
            let row = |offset: i32, fraction: f64, span, edits| Row {
                read: read.read.clone(),
                forward: read.forward,
                mixture: [assertion.call_start, assertion.call_end],
                explained: assertion.explained,
                components: assertion.components.len(),
                offset,
                fraction,
                residual: assertion.residual,
                unexplained: assertion.unexplained.len(),
                span,
                edits,
            };
            if mixture.components.is_empty() {
                rows.extend(
                    assertion.components.iter().map(|component| {
                        row(component.offset, component.fraction, None, Vec::new())
                    }),
                );
            } else {
                rows.extend(mixture.components.iter().map(|component| {
                    row(
                        component.offset,
                        component.fraction,
                        Some(component.span),
                        component
                            .variants
                            .iter()
                            .map(|observed| {
                                (
                                    observed.variant.position_1based,
                                    observed.variant.reference.clone(),
                                    observed.variant.alternate.clone(),
                                )
                            })
                            .collect(),
                    )
                }));
            }
        }
    }
    rows
}

/// The table's header.
pub(crate) const HEADER: &str = "sample_id\tbatch_date\tread\tstrand\tmixture\texplained\tcomponents\toffset\tfraction\tresidual\tunexplained\tspan\trun\tentered\ttruth_entered\tagrees\n";

/// Appends `sample`'s rows, each run cluster of `reference` they cover
/// measured against `truth`; a row covering none, or not decomposed, once
/// with its measures unset.
pub(crate) fn append(
    table: &mut String,
    sample: (&str, &str),
    rows: &[Row],
    truth: Option<&[Allele]>,
    reference: &[u8],
) {
    for row in rows {
        let head = format!(
            "{}\t{}\t{}\t{}\t{}-{}\t{}\t{}\t{}\t{:.3}\t{:.4}\t{}\t{}",
            sample.0,
            sample.1,
            row.read,
            if row.forward { "forward" } else { "reverse" },
            row.mixture[0],
            row.mixture[1],
            if row.explained { "yes" } else { "no" },
            row.components,
            row.offset,
            row.fraction,
            row.residual,
            row.unexplained,
            row.span
                .map_or_else(|| ".".to_owned(), |span| format!("{}-{}", span[0], span[1])),
        );
        let found = row.span.map_or_else(Vec::new, |span| {
            pieces(span, reference.len())
                .into_iter()
                .flat_map(|piece| clusters(reference, piece))
                .collect()
        });
        if found.is_empty() {
            let _ = writeln!(table, "{head}\t.\t.\t.\t.");
        }
        for cluster in found {
            let window = [cluster[0].0 - 1, cluster[cluster.len() - 1].1 + 1];
            let whole = haplotype(reference, window, &row.edits);
            let run = entered(window, row.forward, |window| {
                haplotype(reference, window, &row.edits)
            });
            let (known, truth_run) = truth.map_or((None, None), |truth| {
                (
                    truth_haplotype(reference, window, truth),
                    entered(window, row.forward, |window| {
                        truth_haplotype(reference, window, truth)
                    }),
                )
            });
            let agrees = match (&whole, &known) {
                (Some(whole), Some(known)) => {
                    if whole == known {
                        "yes"
                    } else {
                        "no"
                    }
                }
                _ => ".",
            };
            let text = |run: Option<(char, usize)>| {
                run.map_or_else(
                    || ".".to_owned(),
                    |(base, length)| format!("{base}{length}"),
                )
            };
            let _ = writeln!(
                table,
                "{head}\t{}-{}\t{}\t{}\t{agrees}",
                cluster[0].0,
                cluster[cluster.len() - 1].1,
                text(run),
                text(truth_run),
            );
        }
    }
}

/// The run a read enters in phase over the cluster `window`: the first block
/// after its first base for a `forward` read, the last before its last base
/// for a reverse one, of the haplotype `of` gives. A forward read's window
/// keeps its first base, where an insertion at the run's start is anchored,
/// and drops it from the haplotype; a reverse read's ends before its last.
fn entered(
    window: [usize; 2],
    forward: bool,
    of: impl Fn([usize; 2]) -> Option<String>,
) -> Option<(char, usize)> {
    if forward {
        let haplotype = of(window)?;
        blocks(haplotype.get(1..)?).first().copied()
    } else {
        blocks(&of([window[0], window[1] - 1])?).last().copied()
    }
}

/// The 1-based inclusive intervals of `span`: two when it crosses the origin
/// of a reference of `length` bases.
fn pieces(span: [usize; 2], length: usize) -> Vec<[usize; 2]> {
    if span[0] <= span[1] {
        vec![span]
    } else {
        vec![[span[0], length], [1, span[1]]]
    }
}

/// The reference homopolymers of at least [`MINIMUM_RUN`] bases, as 1-based
/// inclusive intervals with their base, grouped into clusters of runs one
/// base apart; only clusters whose runs lie in `span`.
fn clusters(reference: &[u8], span: [usize; 2]) -> Vec<Vec<(usize, usize, char)>> {
    let mut runs: Vec<(usize, usize, char)> = Vec::new();
    let mut position = span[0];
    // Step back to the start of a run the span opens inside.
    while position > 1 && reference.get(position - 2) == reference.get(position - 1) {
        position -= 1;
    }
    while position <= span[1] {
        let Some(&base) = reference.get(position - 1) else {
            break;
        };
        let mut last = position;
        while reference.get(last) == Some(&base) {
            last += 1;
        }
        if last + 1 - position >= MINIMUM_RUN {
            runs.push((position, last, char::from(base)));
        }
        position = last + 1;
    }
    let mut clusters: Vec<Vec<(usize, usize, char)>> = Vec::new();
    for run in runs {
        match clusters.last_mut() {
            Some(cluster) if cluster[cluster.len() - 1].1 + 2 == run.0 => cluster.push(run),
            _ => clusters.push(vec![run]),
        }
    }
    clusters.retain(|cluster| {
        cluster[0].0 > 1 && cluster[0].0 >= span[0] && cluster[cluster.len() - 1].1 <= span[1]
    });
    clusters
}

/// The reference over the 1-based inclusive `window` with `edits` (1-based
/// anchor, reference and alternate alleles) applied; `None` when an edit
/// crosses the window's edge or two overlap.
fn haplotype(
    reference: &[u8],
    window: [usize; 2],
    edits: &[(usize, String, String)],
) -> Option<String> {
    let mut inside = edits
        .iter()
        .filter(|(anchor, allele, _)| *anchor + allele.len() > window[0] && *anchor <= window[1])
        .collect::<Vec<_>>();
    inside.sort_by_key(|(anchor, _, _)| *anchor);
    let mut sequence = String::new();
    let mut position = window[0];
    for (anchor, allele, alternate) in inside {
        if *anchor < position || anchor + allele.len() - 1 > window[1] {
            return None;
        }
        sequence.push_str(std::str::from_utf8(reference.get(position - 1..anchor - 1)?).ok()?);
        sequence.push_str(alternate);
        position = anchor + allele.len();
    }
    sequence.push_str(std::str::from_utf8(reference.get(position - 1..window[1])?).ok()?);
    Some(sequence)
}

/// The truth's haplotype over `window`, as [`haplotype`]; `None` when an
/// IUPAC symbol or an unreadable position falls in it.
fn truth_haplotype(reference: &[u8], window: [usize; 2], truth: &[Allele]) -> Option<String> {
    let mut edits = Vec::new();
    for allele in truth {
        let (anchor, ordinal) = match allele.position.split_once('.') {
            Some((anchor, ordinal)) => (
                anchor.parse::<usize>().ok()?,
                ordinal.parse::<usize>().ok()?,
            ),
            None => (allele.position.parse::<usize>().ok()?, 0),
        };
        if !(window[0]..=window[1]).contains(&anchor) {
            continue;
        }
        if allele.is_iupac() {
            return None;
        }
        edits.push((anchor, ordinal, allele.alternate.replace('-', "")));
    }
    edits.sort();
    let mut merged: Vec<(usize, String, String)> = Vec::new();
    for (anchor, ordinal, bases) in edits {
        match merged.last_mut() {
            Some(last) if last.0 == anchor && ordinal > 0 => last.2.push_str(&bases),
            _ => {
                let own = char::from(*reference.get(anchor - 1)?).to_string();
                let first = if ordinal == 0 {
                    bases
                } else {
                    format!("{own}{bases}")
                };
                merged.push((anchor, own, first));
            }
        }
    }
    haplotype(reference, window, &merged)
}

/// Runs of equal symbols: symbol and length.
fn blocks(sequence: &str) -> Vec<(char, usize)> {
    let mut blocks: Vec<(char, usize)> = Vec::new();
    for symbol in sequence.chars() {
        match blocks.last_mut() {
            Some(block) if block.0 == symbol => block.1 += 1,
            _ => blocks.push((symbol, 1)),
        }
    }
    blocks
}

#[cfg(test)]
mod tests {
    use super::*;

    fn row(edits: &[(usize, &str, &str)]) -> Row {
        Row {
            read: "r".into(),
            forward: true,
            mixture: [0, 10],
            explained: true,
            components: 1,
            offset: 0,
            fraction: 1.0,
            residual: 0.0,
            unexplained: 0,
            span: Some([1, 20]),
            edits: edits
                .iter()
                .map(|&(anchor, reference, alternate)| (anchor, reference.into(), alternate.into()))
                .collect(),
        }
    }

    #[test]
    fn measures_a_cluster_whatever_the_edit_split() {
        // Flank A, C7 at 2-8, T at 9, C5 at 10-14, flank G.
        let reference = b"ACCCCCCCTCCCCCGA";
        assert_eq!(
            clusters(reference, [2, 16]),
            [vec![(2, 8, 'C'), (10, 14, 'C')]]
        );
        let truth = [
            Allele::new("8.1", "-", "C"),
            Allele::new("8.2", "-", "C"),
            Allele::new("14.1", "-", "C"),
        ];
        let window = [1, 15];
        let expected = "ACCCCCCCTCCCCCG";
        assert_eq!(haplotype(reference, window, &[]).as_deref(), Some(expected));
        let known = "ACCCCCCCCCTCCCCCCG";
        assert_eq!(
            truth_haplotype(reference, window, &truth).as_deref(),
            Some(known)
        );
        let truth_of = |window| truth_haplotype(reference, window, &truth);
        assert_eq!(entered(window, true, truth_of), Some(('C', 9)));
        assert_eq!(entered(window, false, truth_of), Some(('C', 6)));
        // The same sequence as an insertion of CCT and a substitution of T.
        let split = row(&[(8, "C", "CCCT"), (9, "T", "C")]);
        assert_eq!(
            haplotype(reference, window, &split.edits).as_deref(),
            Some(known)
        );
        // A lost T merges the runs.
        let lost = "ACCCCCCCCCCCCCG";
        assert_eq!(
            haplotype(reference, window, &row(&[(9, "T", "C")]).edits).as_deref(),
            Some(lost)
        );
        let lost_of = |window| haplotype(reference, window, &row(&[(9, "T", "C")]).edits);
        assert_eq!(entered(window, true, lost_of), Some(('C', 13)));
    }

    #[test]
    fn leaves_edge_crossing_edits_and_iupac_unmeasured() {
        let reference = b"ACCCCCCCTCCCCCGA";
        let window = [1, 15];
        assert_eq!(
            haplotype(reference, window, &row(&[(15, "GA", "G")]).edits),
            None
        );
        assert_eq!(
            haplotype(reference, window, &row(&[(15, "G", "GC")]).edits).as_deref(),
            Some("ACCCCCCCTCCCCCGC")
        );
        assert_eq!(
            truth_haplotype(reference, window, &[Allele::new("5", "C", "Y")]),
            None
        );
        assert_eq!(
            truth_haplotype(reference, window, &[Allele::new("5", "C", "-")]).as_deref(),
            Some("ACCCCCCTCCCCCG")
        );
    }
}
