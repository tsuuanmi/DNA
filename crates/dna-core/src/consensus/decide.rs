//! Decision rule version 1: one site or stretch from its reads' observations.

use std::collections::BTreeMap;

use crate::model::consensus::SiteState;

/// The decision at one site or stretch and the reads behind it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct Decision {
    pub(super) value: Option<String>,
    pub(super) state: SiteState,
    pub(super) supporting: Vec<String>,
    pub(super) opposing: Vec<String>,
    pub(super) uninformative: Vec<String>,
}

/// Decides over clean observations only; each observation is a read's name
/// and its sequence, `None` when it is not clean:
///
/// 1. all agree: that sequence;
/// 2. disagreement: a strict majority wins;
/// 3. one read against one read: the reference sequence, contested;
/// 4. any other tie: no sequence, contested;
/// 5. no clean observation: no sequence, unresolved.
pub(super) fn decide(observations: &[(&str, Option<String>)], reference: &str) -> Decision {
    let mut clean: BTreeMap<&str, Vec<&str>> = BTreeMap::new();
    let mut uninformative = Vec::new();
    for (read, observation) in observations {
        match observation {
            Some(sequence) => clean.entry(sequence.as_str()).or_default().push(read),
            None => uninformative.push((*read).to_owned()),
        }
    }
    let mut ranked: Vec<(&str, usize)> = clean
        .iter()
        .map(|(sequence, reads)| (*sequence, reads.len()))
        .collect();
    ranked.sort_by_key(|(_, count)| std::cmp::Reverse(*count));
    let (value, state) = match ranked.as_slice() {
        [] => (None, SiteState::Unresolved),
        [(only, _)] => (Some((*only).to_owned()), SiteState::Called),
        [(top, top_count), (_, second_count), ..] if top_count > second_count => {
            (Some((*top).to_owned()), SiteState::Called)
        }
        [(first, 1), (second, 1)] if *first == reference || *second == reference => {
            (Some(reference.to_owned()), SiteState::Contested)
        }
        _ => (None, SiteState::Contested),
    };
    let mut supporting = Vec::new();
    let mut opposing = Vec::new();
    for (sequence, reads) in &clean {
        let target = if value.as_deref() == Some(*sequence) {
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
