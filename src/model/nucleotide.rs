//! Canonical and ambiguous DNA symbols.

/// A canonical DNA base.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Nucleotide {
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
    pub const ALL: [Self; 4] = [Self::A, Self::C, Self::G, Self::T];

    /// Returns the uppercase nucleotide character.
    pub const fn as_char(self) -> char {
        match self {
            Self::A => 'A',
            Self::C => 'C',
            Self::G => 'G',
            Self::T => 'T',
        }
    }

    /// Returns the canonical A/C/G/T channel index.
    pub const fn channel_index(self) -> usize {
        match self {
            Self::A => 0,
            Self::C => 1,
            Self::G => 2,
            Self::T => 3,
        }
    }
}

/// Complements an uppercase IUPAC DNA character.
pub const fn complement_iupac(value: char) -> char {
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
pub fn reverse_complement(sequence: &str) -> String {
    sequence.chars().rev().map(complement_iupac).collect()
}
