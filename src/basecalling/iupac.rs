//! Standard two-base IUPAC ambiguity mapping.

use crate::model::nucleotide::Nucleotide;

/// Returns the canonical or two-base IUPAC symbol for a sorted base set.
pub(crate) fn code(bases: &[Nucleotide]) -> char {
    let mut present = [false; 4];
    for base in bases {
        present[base.channel_index()] = true;
    }
    match present {
        [true, false, false, false] => 'A',
        [false, true, false, false] => 'C',
        [false, false, true, false] => 'G',
        [false, false, false, true] => 'T',
        [true, true, false, false] => 'M',
        [true, false, true, false] => 'R',
        [true, false, false, true] => 'W',
        [false, true, true, false] => 'S',
        [false, true, false, true] => 'Y',
        [false, false, true, true] => 'K',
        _ => 'N',
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn maps_two_base_codes() {
        assert_eq!(code(&[Nucleotide::A, Nucleotide::G]), 'R');
        assert_eq!(code(&[Nucleotide::C, Nucleotide::T]), 'Y');
        assert_eq!(code(&[Nucleotide::A, Nucleotide::C, Nucleotide::G]), 'N');
    }
}
