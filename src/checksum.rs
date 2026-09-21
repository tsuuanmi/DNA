//! Stable SHA-256 identities for validated input bytes.

use sha2::{Digest, Sha256};

/// Returns a lowercase hexadecimal SHA-256 digest.
pub(crate) fn hex_sha256(bytes: &[u8]) -> String {
    format!("{:x}", Sha256::digest(bytes))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hashes_bytes_deterministically() {
        assert_eq!(
            hex_sha256(b"dna"),
            "d0331ca087d6901424950a896dcc66d77de4758143396fa594bbdf343e09a843"
        );
    }
}
