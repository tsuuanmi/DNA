//! Draft `region.minor_signal` rules: a position where the clean reads show a
//! second allele, which a single called base hides (point heteroplasmy).
//!
//! Every rule reads the minor channels of one position
//! ([`crate::evidence::minor`]) and counts a read when its second channel
//! reaches the rule's threshold (`>=`, in thousandths) and passes the rule's
//! shadow test:
//!
//! - `strands`: the same allele pair on both strands, shadows included;
//! - `unshadowed`: as `strands`, without reads whose second channel is the
//!   base called just before or after (a neighbour's shadow);
//! - `upstream`: as `strands`, without reads whose second channel is the
//!   base of the next reference position (a shadow in either strand's
//!   reading direction); a second channel equal to the previous reference
//!   base still counts;
//! - `noise`: as `strands`, with each read's second channel at least
//!   [`NOISE_RATIO`] times its local noise instead of the shadow test;
//! - `any`: one read without a shadow is enough;
//! - `lone`: `unshadowed`, or one read without a shadow where only one
//!   strand reads the position cleanly.

use crate::evidence::Position;
use crate::evidence::Strand;
use crate::evidence::minor::Minor;

/// `noise`: a read counts when its second channel is at least this many
/// times its local noise.
pub(crate) const NOISE_RATIO: u16 = 3;

/// How a minor-signal rule counts reads.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) enum Kind {
    Strands,
    Unshadowed,
    Upstream,
    Noise,
    Any,
    Lone,
}

impl Kind {
    pub(crate) const ALL: [Self; 6] = [
        Self::Strands,
        Self::Unshadowed,
        Self::Upstream,
        Self::Noise,
        Self::Any,
        Self::Lone,
    ];

    pub(crate) const fn label(self) -> &'static str {
        match self {
            Self::Strands => "strands",
            Self::Unshadowed => "unshadowed",
            Self::Upstream => "upstream",
            Self::Noise => "noise",
            Self::Any => "any",
            Self::Lone => "lone",
        }
    }
}

/// A minor-signal rule and its threshold.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) struct Rule {
    pub(crate) kind: Kind,
    /// The least second-channel weight a read needs, in thousandths.
    pub(crate) threshold: u16,
}

impl Rule {
    /// The rule the draft policy applies: on the development panel, the
    /// cheapest swept rule that removes more than half of the heteroplasmy
    /// leaks.
    pub(crate) const DEFAULT: Self = Self {
        kind: Kind::Upstream,
        threshold: 100,
    };

    /// `upstream@100`.
    pub(crate) fn label(self) -> String {
        format!("{}@{}", self.kind.label(), self.threshold)
    }

    /// Whether the rule marks a position with these minor channels.
    pub(crate) fn marks(self, minors: &[Minor], position: Position) -> bool {
        let counts = |minor: &&Minor| {
            minor.weight >= self.threshold
                && match self.kind {
                    Kind::Strands => true,
                    Kind::Unshadowed | Kind::Any | Kind::Lone => !minor.shadow(),
                    Kind::Upstream => !minor.downstream(),
                    Kind::Noise => {
                        u32::from(minor.weight) >= u32::from(NOISE_RATIO) * u32::from(minor.noise)
                    }
                }
        };
        let counted = || minors.iter().filter(counts);
        let both = || {
            counted().any(|forward| {
                forward.strand == Strand::Forward
                    && counted().any(|reverse| {
                        reverse.strand == Strand::Reverse && reverse.alleles() == forward.alleles()
                    })
            })
        };
        match self.kind {
            Kind::Strands | Kind::Unshadowed | Kind::Upstream | Kind::Noise => both(),
            Kind::Any => counted().next().is_some(),
            Kind::Lone => both() || (position.strands() == 1 && counted().next().is_some()),
        }
    }
}

#[cfg(test)]
mod tests {
    use dna_kernel::model::nucleotide::Nucleotide;

    use super::*;

    /// A read of the C/T pair (or `pair`) with a second channel of `weight`.
    fn read(strand: Strand, pair: [Nucleotide; 2], weight: u16) -> Minor {
        Minor {
            strand,
            pair,
            weight,
            before: false,
            after: false,
            noise: 20,
        }
    }

    const CT: [Nucleotide; 2] = [Nucleotide::C, Nucleotide::T];
    const TC: [Nucleotide; 2] = [Nucleotide::T, Nucleotide::C];

    fn two_strands() -> Position {
        Position {
            forward: 1,
            reverse: 1,
            ..Position::default()
        }
    }

    fn rule(kind: Kind, threshold: u16) -> Rule {
        Rule { kind, threshold }
    }

    #[test]
    fn cross_strand_rules_need_the_same_pair_on_both_strands_at_the_threshold() {
        let both = [
            read(Strand::Forward, CT, 150),
            read(Strand::Reverse, TC, 151),
        ];
        for kind in [
            Kind::Strands,
            Kind::Unshadowed,
            Kind::Upstream,
            Kind::Noise,
            Kind::Lone,
        ] {
            // At the threshold, and just below it.
            assert!(rule(kind, 150).marks(&both, two_strands()));
            assert!(!rule(kind, 151).marks(&both, two_strands()));
        }
        let one_strand = [
            read(Strand::Forward, CT, 300),
            read(Strand::Forward, CT, 300),
        ];
        assert!(!rule(Kind::Strands, 150).marks(&one_strand, two_strands()));
        assert!(rule(Kind::Any, 150).marks(&one_strand, two_strands()));
        let other_pair = [
            read(Strand::Forward, CT, 300),
            read(Strand::Reverse, [Nucleotide::C, Nucleotide::A], 300),
        ];
        assert!(!rule(Kind::Strands, 150).marks(&other_pair, two_strands()));
        assert!(!rule(Kind::Strands, 50).marks(&[], two_strands()));
    }

    #[test]
    fn shadows_count_only_for_strands_and_noise_needs_a_ratio() {
        let mut forward = read(Strand::Forward, CT, 200);
        forward.before = true;
        let reads = [forward, read(Strand::Reverse, CT, 200)];
        assert!(rule(Kind::Strands, 150).marks(&reads, two_strands()));
        assert!(!rule(Kind::Unshadowed, 150).marks(&reads, two_strands()));
        // Only the reverse read is no shadow.
        assert!(rule(Kind::Any, 150).marks(&reads, two_strands()));
        // The forward read's shadow is the previous reference base, which
        // `upstream` keeps; the next reference base (after a forward call,
        // before a reverse one) it drops.
        assert!(rule(Kind::Upstream, 150).marks(&reads, two_strands()));
        let mut downstream = reads;
        downstream[0].before = false;
        downstream[0].after = true;
        assert!(!rule(Kind::Upstream, 150).marks(&downstream, two_strands()));
        downstream[0].after = false;
        downstream[1].before = true;
        assert!(!rule(Kind::Upstream, 150).marks(&downstream, two_strands()));
        downstream[1].before = false;
        downstream[1].after = true;
        assert!(rule(Kind::Upstream, 150).marks(&downstream, two_strands()));
        // Noise 20: the ratio needs 60 or more; 70 is above, noise 24 is
        // exactly at it, noise 25 below.
        let mut noisy = [read(Strand::Forward, CT, 70), read(Strand::Reverse, CT, 72)];
        assert!(rule(Kind::Noise, 50).marks(&noisy, two_strands()));
        noisy[1].noise = 24;
        assert!(rule(Kind::Noise, 50).marks(&noisy, two_strands()));
        noisy[1].noise = 25;
        assert!(!rule(Kind::Noise, 50).marks(&noisy, two_strands()));
    }

    #[test]
    fn lone_takes_one_read_only_where_one_strand_reads_cleanly() {
        let reads = [read(Strand::Forward, CT, 200)];
        let one = Position {
            forward: 1,
            ..Position::default()
        };
        assert!(rule(Kind::Lone, 200).marks(&reads, one));
        assert!(!rule(Kind::Lone, 201).marks(&reads, one));
        assert!(!rule(Kind::Lone, 200).marks(&reads, two_strands()));
        assert!(!rule(Kind::Unshadowed, 200).marks(&reads, one));
    }

    #[test]
    fn labels_name_the_kind_and_threshold() {
        assert_eq!(Rule::DEFAULT.label(), "upstream@100");
        assert_eq!(
            Kind::ALL.map(Kind::label),
            ["strands", "unshadowed", "upstream", "noise", "any", "lone"]
        );
    }
}
