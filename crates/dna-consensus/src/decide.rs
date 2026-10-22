//! Decision rule version 1: one site or stretch from its reads' observations.

use std::collections::BTreeMap;

use crate::model::consensus::{ConsensusRun, SiteState};

/// The decision at one site or stretch and the reads behind it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct Decision {
    pub(crate) value: Option<String>,
    pub(crate) state: SiteState,
    pub(crate) supporting: Vec<String>,
    pub(crate) opposing: Vec<String>,
    pub(crate) uninformative: Vec<String>,
    /// The decided runs, for a stretch decided run by run.
    pub(crate) runs: Vec<ConsensusRun>,
}

/// What a clean observation votes for: a sequence, or, for a read that
/// shows the site otherwise than every sequence voted for (`alone`), a
/// vote of its own that no other read shares and that never decides. Votes
/// order sequences first, then lone votes by read.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum Vote<'a> {
    Sequence(&'a str),
    Alone(&'a str),
}

/// Decides over clean observations only; each observation is a read's identity
/// and its sequence, `None` when it is not clean:
///
/// 1. all agree: that sequence;
/// 2. disagreement: a strict majority wins;
/// 3. one read against one read: the reference sequence, contested;
/// 4. any other tie: no sequence, contested;
/// 5. no clean observation: no sequence, unresolved.
pub(crate) fn decide(observations: &[(&str, Option<String>)], reference: &str) -> Decision {
    decide_with_lone(observations, &[], reference)
}

/// [`decide`], with the reads in `alone` each opposing every sequence with a
/// clean observation of its own that no other read shares. A lone vote never
/// decides: when one would lead, the site is contested.
pub(crate) fn decide_with_lone(
    observations: &[(&str, Option<String>)],
    alone: &[&str],
    reference: &str,
) -> Decision {
    let mut clean: BTreeMap<Vote<'_>, Vec<&str>> = BTreeMap::new();
    let mut uninformative = Vec::new();
    for (read, observation) in observations {
        match observation {
            Some(sequence) => clean
                .entry(Vote::Sequence(sequence.as_str()))
                .or_default()
                .push(read),
            None => uninformative.push((*read).to_owned()),
        }
    }
    for read in alone {
        clean.entry(Vote::Alone(read)).or_default().push(read);
    }
    let mut ranked: Vec<(Vote<'_>, usize)> = clean
        .iter()
        .map(|(vote, reads)| (*vote, reads.len()))
        .collect();
    ranked.sort_by_key(|(_, count)| std::cmp::Reverse(*count));
    let is_reference = |vote: &Vote<'_>| *vote == Vote::Sequence(reference);
    let (value, state) = match ranked.as_slice() {
        [] => (None, SiteState::Unresolved),
        [(Vote::Sequence(only), _)] => (Some((*only).to_owned()), SiteState::Called),
        [(Vote::Sequence(top), top_count), (_, second_count), ..] if top_count > second_count => {
            (Some((*top).to_owned()), SiteState::Called)
        }
        [(first, 1), (second, 1)] if is_reference(first) || is_reference(second) => {
            (Some(reference.to_owned()), SiteState::Contested)
        }
        _ => (None, SiteState::Contested),
    };
    let mut supporting = Vec::new();
    let mut opposing = Vec::new();
    for (vote, reads) in &clean {
        let target = if matches!(vote, Vote::Sequence(sequence) if value.as_deref() == Some(*sequence))
        {
            &mut supporting
        } else {
            &mut opposing
        };
        target.extend(reads.iter().map(|read| (*read).to_owned()));
    }
    Decision {
        value,
        state,
        supporting,
        opposing,
        uninformative,
        runs: Vec::new(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn decided(observations: &[(&str, Option<&str>)]) -> Decision {
        let owned = observations
            .iter()
            .map(|(read, sequence)| (*read, sequence.map(str::to_owned)))
            .collect::<Vec<_>>();
        decide(&owned, "C")
    }

    // SRS-CONS-005
    #[test]
    fn calls_agreement_and_strict_majority() {
        let agreed = decided(&[("a", Some("T")), ("b", Some("T"))]);
        assert_eq!(
            (agreed.value.as_deref(), agreed.state),
            (Some("T"), SiteState::Called)
        );
        let majority = decided(&[("a", Some("T")), ("b", Some("T")), ("c", Some("C"))]);
        assert_eq!(
            (majority.value.as_deref(), majority.state),
            (Some("T"), SiteState::Called)
        );
        assert_eq!(majority.opposing, ["c"]);
    }

    // SRS-CONS-005
    #[test]
    fn breaks_a_one_against_one_tie_toward_the_reference() {
        let tie = decided(&[("a", Some("T")), ("b", Some("C"))]);
        assert_eq!(
            (tie.value.as_deref(), tie.state),
            (Some("C"), SiteState::Contested)
        );
        assert_eq!(
            (tie.supporting, tie.opposing),
            (vec!["b".to_owned()], vec!["a".to_owned()])
        );
        let neither = decided(&[("a", Some("T")), ("b", Some("G"))]);
        assert_eq!((neither.value, neither.state), (None, SiteState::Contested));
        let wider = decided(&[
            ("a", Some("T")),
            ("b", Some("T")),
            ("c", Some("C")),
            ("d", Some("C")),
        ]);
        assert_eq!((wider.value, wider.state), (None, SiteState::Contested));
    }

    // SRS-CONS-005
    #[test]
    fn a_lone_vote_opposes_and_never_decides() {
        let owned = |observations: &[(&'static str, &str)]| {
            observations
                .iter()
                .map(|(read, sequence)| (*read, Some((*sequence).to_owned())))
                .collect::<Vec<_>>()
        };
        let majority = decide_with_lone(&owned(&[("a", "T"), ("b", "T")]), &["c"], "C");
        assert_eq!(
            (majority.value.as_deref(), majority.state),
            (Some("T"), SiteState::Called)
        );
        assert_eq!(majority.opposing, ["c"]);
        // One read against a lone read: the reference, contested.
        let tie = decide_with_lone(&owned(&[("a", "C")]), &["b"], "C");
        assert_eq!(
            (tie.value.as_deref(), tie.state),
            (Some("C"), SiteState::Contested)
        );
        let neither = decide_with_lone(&owned(&[("a", "T")]), &["b"], "C");
        assert_eq!((neither.value, neither.state), (None, SiteState::Contested));
        // A lone read alone decides nothing.
        let alone = decide_with_lone(&[], &["a"], "C");
        assert_eq!((alone.value, alone.state), (None, SiteState::Contested));
        assert_eq!(alone.opposing, ["a"]);
    }

    // SRS-CONS-003, SRS-CONS-005
    #[test]
    fn lets_only_clean_observations_decide() {
        let masked = decided(&[("a", Some("T")), ("b", None)]);
        assert_eq!(
            (masked.value.as_deref(), masked.state),
            (Some("T"), SiteState::Called)
        );
        assert_eq!(masked.uninformative, ["b"]);
        let none = decided(&[("a", None)]);
        assert_eq!((none.value, none.state), (None, SiteState::Unresolved));
    }
}
