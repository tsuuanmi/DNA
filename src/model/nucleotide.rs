//! Canonical and ambiguous DNA symbols.

/// A canonical DNA base.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub(crate) enum Nucleotide {
    /// Adenine.
    A,
    /// Cytosine.
    C,
    /// Guanine.
    G,
    /// Thymine.
    T,
}

impl Nucleotide {
    /// Canonical channel order used throughout DNA.
    pub(crate) const ALL: [Self; 4] = [Self::A, Self::C, Self::G, Self::T];

    /// Parses an uppercase canonical base; every other symbol is `None`.
    pub(crate) const fn from_char(value: char) -> Option<Self> {
        match value {
            'A' => Some(Self::A),
            'C' => Some(Self::C),
            'G' => Some(Self::G),
            'T' => Some(Self::T),
            _ => None,
        }
    }

    /// Returns the uppercase nucleotide character.
    pub(crate) const fn as_char(self) -> char {
        match self {
            Self::A => 'A',
            Self::C => 'C',
            Self::G => 'G',
            Self::T => 'T',
        }
    }

    /// Returns the canonical A/C/G/T channel index.
    pub(crate) const fn channel_index(self) -> usize {
        match self {
            Self::A => 0,
            Self::C => 1,
            Self::G => 2,
            Self::T => 3,
        }
    }
}

/// Whether `value` is an uppercase canonical A/C/G/T base.
pub(crate) const fn is_canonical(value: char) -> bool {
    Nucleotide::from_char(value).is_some()
}

/// Complements an uppercase IUPAC DNA character.
pub(crate) const fn complement_iupac(value: char) -> char {
    match value {
        'A' => 'T',
        'C' => 'G',
        'G' => 'C',
        'T' => 'A',
        'R' => 'Y',
        'Y' => 'R',
        'S' => 'S',
        'W' => 'W',
        'K' => 'M',
        'M' => 'K',
        'B' => 'V',
        'D' => 'H',
        'H' => 'D',
        'V' => 'B',
        _ => 'N',
    }
}

/// Reverse-complements an uppercase DNA/IUPAC sequence.
pub(crate) fn reverse_complement(sequence: &str) -> String {
    sequence.chars().rev().map(complement_iupac).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_only_uppercase_canonical_bases_in_channel_order() {
        let parsed = ['A', 'C', 'G', 'T'].map(Nucleotide::from_char);
        assert_eq!(parsed, Nucleotide::ALL.map(Some));
        for (index, base) in Nucleotide::ALL.iter().enumerate() {
            assert_eq!(base.channel_index(), index);
        }
        for other in ['a', 'N', 'R', '-', 'Ł', '\0'] {
            assert_eq!(Nucleotide::from_char(other), None);
            assert!(!is_canonical(other));
        }
        assert!(['A', 'C', 'G', 'T'].into_iter().all(is_canonical));
    }
}
