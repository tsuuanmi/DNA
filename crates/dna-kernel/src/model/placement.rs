//! The one primitive of 3'/right-most placement: a sequence-equivalent move
//! of an insertion or a deletion one base to the right (ADR-0047). Alignment
//! canonicalization in the core and variant normalization both move gaps
//! with it, so a read's indel and its called variant sit at the same place.

/// Whether an insertion or a deletion whose first inserted or deleted base is
/// `first` may move one base to the right across the reference base
/// `crossed` without changing the haplotype: the two bases are equal.
#[must_use]
pub const fn crosses(first: u8, crossed: u8) -> bool {
    first == crossed
}

/// The inserted bases of an insertion once it has moved one base to the
/// right across `crossed`: its first base leaves at the front and `crossed`
/// joins at the back.
#[must_use]
pub fn rotated(inserted: &[u8], crossed: u8) -> Vec<u8> {
    let mut rotated = inserted.get(1..).unwrap_or_default().to_vec();
    rotated.push(crossed);
    rotated
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn an_insertion_rotates_as_it_crosses_its_first_base() {
        // CA inserted before CACG moves across C, then A: AC, then CA.
        assert!(crosses(b'C', b'C'));
        assert_eq!(rotated(b"CA", b'C'), b"AC");
        assert_eq!(rotated(b"AC", b'A'), b"CA");
        assert!(!crosses(b'C', b'G'));
    }
}
